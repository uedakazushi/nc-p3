"""The full quotient, opposite ranks, and central bundles on the support-145 surface."""
from itertools import combinations, product
import os
import subprocess

from flint_rational_function import rational_function_field
from localized_unit_echelon import unit_echelon
from verify_klein_center_lift import central_equations,field_paths
from verify_parallel_excluded_shapes import path_dimension
from verify_parallel_invertible_rectangles import K,WORDS,X,completed,matrix,s,template
from verify_third_rectangle_frames import SINGULAR,frame_differences
from verify_third_rectangle_orbits import support_maps
from verify_transverse_rank_twelve_isomorphism import alternating_gradient

a,y = s.symbols('a y')
U = a*y*(a-1)*(y-1)


def window(a,y,t):
    return completed(112,t,a*y,1,t*a,1,y,t*(a+y-a*y))


def singular(commands,expected):
    source = SINGULAR.replace('(v,h,b)','(a,y)').replace('v*h*(h-1)*b*(b-1)','a*y*(a-1)*(y-1)')
    result = subprocess.run([os.environ.get('SINGULAR_BIN','Singular'),'-q'],
                            input=source+'\n'.join(commands)+'\nquit;\n',text=True,capture_output=True,check=True)
    out = [line.strip() for line in result.stdout.splitlines() if not line.startswith('//')]
    assert out==expected,(out,result.stderr)


def frames(t):
    w = window(a,y,t)
    commands = ['ideal G;']
    expected = []
    for pos,rows in frame_differences(w,(a,y),U,U):
        commands.append(f'G=pairIdeal(ideal({",".join(rows)}));attrib(G,"isSB",1);')
        if pos==(1,2,3):
            commands.extend(('print(reduce(a+1,G));','print(reduce(y+1,G));'))
            expected.extend(('0','0'))
            assert all(s.sympify(f.replace('^','**')).subs({a:-1,y:-1})==0 for f in rows)
        else:
            commands.append('print(G[1]);')
            expected.append('1')
    singular(commands,expected)
    w = {k:window_fixed for k,z in w.items() if (window_fixed:=z.subs({a:-1,y:-1}))}
    g = alternating_gradient({k:w.get(k,0) for k in WORDS})
    def spectra(p):
        matrices = [matrix(w,p,h)*matrix(g,p,h).T for h in K]
        return sorted(sum((m*m)[i,i] for m in matrices) for i in K)
    assert spectra(0)!=spectra(1)


def orbits(t):
    w = {k:f for k,f in window(a,y,t).items() if f}
    holes = [k for k,x in template(112).items() if x in X]
    thin = [k for k in w if k not in holes]
    from verify_klein_three_lines import character,smith
    for shift in K:
        z = {k[shift:]+k[:shift]:f for k,f in w.items()}
        perms = list(support_maps(z,w))
        assert len(perms)==(2 if shift%2==0 else 0)
        if not perms:
            continue
        target_parameters = ((a,y),(1/a,1/y)) if shift==0 else ((y,a),(1/y,1/a))
        for pa,py in target_parameters:
            target = window(pa,py,t)
            good = []
            for pi in perms:
                transformed = {tuple(p[i] for p,i in zip(pi,k)):f for k,f in z.items()}
                assert all(k in transformed and s.factor(transformed[k]**2-1)==0 for k in thin)
                keys = sorted(transformed)
                weights = s.Matrix([[int(k[p]==i) for p,i in product(K,repeat=2)] for k in keys])
                diag,left,_ = smith(weights)
                rows = left.tolist()[sum(bool(c) for c in diag.diagonal()):]
                ratios = [s.factor(transformed[k]/target[k]) for k in keys]
                if all(character(row,ratios)==1 for row in rows):
                    good.append(pi)
            assert len(good)==1
    # The inversion and half-turn are valid at the double boundary by the
    # independent support-142 field check in the curve verifier.


def ranks(t):
    w = window(a,y,t)
    b = a*y
    H = b**3-2*b*b*(a+y)+3*b*b+3*b-2*(a+y)+1
    targets = (0,(b-1)**2*(b+1)/U,0,0) if t==1 else (
        -4*(a+y)/((a-1)*(y-1)),-H/U,s.Integer(4),s.Integer(4))
    commands = []
    for h in K:
        m = s.Matrix([[w.get((i,j,i^h,j^h),0) for j in K] for i in K])
        assert s.factor(m.det()-targets[h])==0
        if t==1 and h!=1:
            assert m.rank()==2 and any(s.factor(m.extract(ii,jj).det()) in (-1,1)
                                        for ii,jj in product(combinations(K,2),repeat=2))
        check = ((2,1),(3,b-1)) if h==1 and t==1 else ((3,1),) if t==-1 and h<2 else ()
        for n,div in check:
            polys = [s.fraction(s.factor(m.extract(ii,jj).det()))[0] for ii,jj in product(combinations(K,n),repeat=2)]
            fs = [s.cancel(f/div) for f in polys]
            assert all(f.is_polynomial(a,y) for f in fs)
            commands.append('print(pairIdeal(ideal('+','.join(str(f).replace('**','^') for f in fs)+'))[1]);')
    singular(commands,['1']*len(commands))


def center_bundle():
    w = window(a,y,1)
    field = rational_function_field(('a','y'))
    factor = field.from_sympy(U*(a*y+1)).n
    def reduce(rows):
        return unit_echelon(rows,factor)
    rel = [[(h,{(j,j^h):field.from_sympy(s.factor(c[j])) for j in K if c[j]})
            for h in K for c in matrix(w,p,h).columnspace()] for p in K]
    basis,image = field_paths(rel,field.one,reduce)
    assert all([sum(len(basis[p,n,h]) for h in K) for n in range(6)]==[1,4,10,20,35,56] for p in K)
    for h in K:
        rows,width = central_equations(basis,image,h)
        assert width-len(reduce(rows))==(3,2,0,0)[h]
    bad = window(a,-1/a,1)
    assert [path_dimension(bad,0,4,h) for h in K]==[10,10,8,8]


if __name__=='__main__':
    for t in (1,-1):
        frames(t)
        orbits(t)
        ranks(t)
    center_bundle()
    print('PASS: full rectangular surface quotient, fixed point, ranks and globally free quartic center')
