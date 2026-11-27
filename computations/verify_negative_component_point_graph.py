"""Compute the two projected schemes of the generic negative Klein pairing."""
from itertools import combinations, product
import os
import subprocess

import sympy as s
from verify_klein_elliptic_completion import gram_twist


a,b,c,d,z = s.symbols('a b c d z')
x = s.symbols('x0:4')
K = range(4)
D = a*d-b*c
sample = {a:5,b:2,c:3,d:7,z:3}


def evaluation():
    X,Y,Z,T = x
    return s.Matrix([[a*X,c*Y,Z,-z*T],[b*X,d*Y,-Z,T],[-b*Y,-a*X,T,-Z],
                     [-d*Y,-c*X,-z*T,Z],[Z,z*T,X,Y],[T,Z,Y,X]])


def quadrics(right):
    X,Y,Z,T = x
    if right:
        H = D*(b*z-c)
        L,N = -(c+d)*(a*z+c),-(a+b)*(b*z+d)
        M,O = (a+c)*(c+d*z),(a+b*z)*(b+d)
    else:
        H = D*(a-d)
        L,N = -(b+d)*(c+d),-(a+b)*(a+c)
        M,O = (c+d*z)*(b*z+d),(a+b*z)*(a*z+c)
    assert all(v.subs(sample)!=0 for v in (L,M,N,O,L*O-M*N,H))
    return [H*X**2-L*Z**2-M*T**2,H*Y**2-N*Z**2-O*T**2]


def ideal(polys):
    return 'ideal('+','.join(str(p).replace('**','^') for p in polys)+')'


def verify(right):
    E = evaluation()
    co = [row.jacobian(x).T for row in (E.row(i) for i in range(6))]
    if right:
        E = s.Matrix.vstack(*[(m*s.Matrix(x)).T for m in co])
    # These six columns of the general 16 by 16 unfolding span the relations.
    A = s.Matrix([[a,b],[c,d]])
    W,_,_ = gram_twist(A,-A,z,s.S.NegativeOne)
    chosen = (0,5,1,4,2,3)
    assert W[:,list(chosen)]==s.Matrix.hstack(*[m.reshape(16,1) for m in co])
    X,Y,Z,T = x
    p0 = [Z,T,X*X-(z if right else 1)*Y*Y]
    p1 = [X,Y,Z*Z-(1 if right else z)*T*T]
    determinants = [E.extract(I,K).det(method='domain-ge') for I in combinations(range(6),4)]
    equations = quadrics(right)
    def rings(fiber):
        values = sample if fiber else {}
        def pol(fs):
            return ideal([f.subs(values) for f in fs])
        field = '0' if fiber else '(0,a,b,c,d,z)'
        commands = [f'ring r={field},(x0,x1,x2,x3),dp;',f'ideal I={pol(determinants)};',
                    f'ideal U=intersect(intersect({pol(equations)},{pol(p0)}),{pol(p1)});',
                    'print(size(reduce(I,std(U))));']
        if fiber:
            for j in K:
                chart = pol([x[j]-1])
                commands.append(f'print(size(reduce(U,std(I+{chart}))));')
            entries = [E.extract(I,J).subs(sample).det() for I in combinations(range(6),3)
                       for J in combinations(K,3)]
            commands.append(f'ideal minorThree={ideal(entries)};')
            for j in K:
                commands.append(f'print(std(minorThree+ideal(x{j}-1))[1]);')
        output = subprocess.run([os.environ.get('SINGULAR_BIN','Singular'),'-q'],
                                input='\n'.join(commands+['quit;']),text=True,capture_output=True,check=True)
        lines = [line for line in output.stdout.splitlines() if line and not line.startswith('//')]
        assert lines==(['0']*5+['1']*4 if fiber else ['0']),lines
    rings(False)
    rings(True)


if __name__=='__main__':
    for right in (False,True):
        verify(right)
    print('PASS: generic elliptic quartics and four distinct points on both projections, with rank-three graph')
