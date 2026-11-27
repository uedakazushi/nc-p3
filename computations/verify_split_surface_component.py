"""Quartic centers, quantum boundary, and full tangent of the split Bell surface."""
from itertools import combinations, product
from subprocess import run

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from verify_isolated_split_component import PRIME, tangent_dimensions
from verify_klein_center_lift import central_equations, field_paths
from verify_split_surface_lifts import BINARY_PAIRS, U, V, WORDS, split_tensor


def verify_centers(w):
    field = sp.QQ.frac_field(U,V)
    tensor = {q:field.from_sympy(sp.sympify(value)) for q,value in w.items()}
    relations = []
    for i in range(4):
        def rotate(q):
            return q[-i:]+q[:-i] if i else q
        flat = DomainMatrix.from_list([[tensor.get(rotate(a+b),field.zero) for a in BINARY_PAIRS]
                                       for b in BINARY_PAIRS],field).to_sparse()
        matrix, pivots = flat.rref()
        assert len(pivots) == 6
        rows = []
        for row in matrix.to_list()[:6]:
            coeffs = {divmod(j,4):v for j,v in enumerate(row) if v}
            degrees = {a^b for a,b in coeffs}
            assert len(degrees) == 1
            rows.append((degrees.pop(),coeffs))
        relations.append(rows)
    basis, image = field_paths(relations,field.one)
    code = ['LIB "linalg.lib";ring r=(0,u,v),dummy,dp;']
    for color in range(4):
        rows,width = central_equations(basis,image,color)
        data = ','.join(str(row.get(j,0)).replace('**','^') for row in rows for j in range(width))
        code.append(f'matrix a{color}[{len(rows)}][{width}]={data};'
                    f'print("KERNEL {color} "+string({width}-rank(a{color})));')
    output = run(['Singular','-q'],input='\n'.join(code+['quit;']),text=True,
                 capture_output=True,check=True).stdout
    assert '?' not in output, output
    assert all(f'KERNEL {i} {v}' in output for i,v in enumerate((5,0,0,0))), output


def verify_quantum(w, sign):
    sample = {q:sp.simplify(sp.sympify(value).subs({U:sp.I,V:sign*sp.I})) for q,value in w.items()}
    sample = {q:value for q,value in sample.items() if value}
    flips = (0,0,1,1) if sign == 1 else (0,1,0,1)
    tensor = {tuple(c^flips[i] for i,c in enumerate(q)):v for q,v in sample.items()}
    assert len(tensor) == 24 and all(len(set(q))==4 for q in tensor)
    all_faces = []
    for i in range(4):
        def value(q):
            return tensor[tuple(q[-i:]+q[:-i] if i else q)]
        ratios = {}
        for a,b in combinations(range(4),2):
            c,d = [j for j in range(4) if j not in (a,b)]
            r = sp.cancel(-value([a,b,c,d])/value([b,a,c,d]))
            assert sp.simplify(r*value([b,a,d,c])+value([a,b,d,c])) == 0
            ratios[a,b] = r
        faces = tuple(sp.simplify(ratios[a,b]*ratios[b,c]/ratios[a,c])
                      for a,b,c in combinations(range(4),3))
        all_faces.append(faces)
    target = (sp.I,-sp.I,sp.I,-sp.I)
    assert len(set(all_faces)) == 1 and all_faces[0] in (target,tuple(-v for v in target))


def verify():
    x,z = sp.Matrix([[0,1],[1,0]]),sp.diag(1,-1)
    for h, sign in ((sp.eye(2),-1),(x,1),(z,-1),(x*z,1)):
        w = split_tensor(h)
        verify_quantum(w,sign)
        verify_centers(w)
        sample = {q:sp.Rational(sp.sympify(value).subs({U:2,V:3})) for q,value in w.items()}
        residue = {q:int(a.p)*pow(int(a.q),-1,PRIME)%PRIME for q,a in sample.items()}
        assert tangent_dimensions(residue) == [3,20,30,20,3]
        print("Five quartic centers, order-four quantum specialization and affine tangent dimension 76.")


if __name__ == '__main__':
    verify()
