"""Exact normalizations, Schur equations, point resultants and periods for the finite 0123 locus."""
from itertools import combinations, product

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

K = range(4)
A,B,C,D,H,F,X,Y,Z = sp.symbols('a b c d h f x y Z')


def rotate(w, k):
    return w[k:]+w[:k]


def monomial_chart():
    words = [w for w in product(K,repeat=4) if sum(w)==6
             and all(0<sum(rotate(w,i)[:2])<6 for i in K) and w not in ((1,2,1,2),(2,1,2,1))]
    ids = {w:i for i,w in enumerate(words)}
    identity = sp.eye(len(words))
    equations = []
    for cut in (0,1):
        for weight in (1,2,4,5):
            rows = [ij for ij in product(K,repeat=2) if sum(ij)==weight]
            cols = [ij for ij in product(K,repeat=2) if sum(ij)==6-weight]
            for p,q in combinations(rows,2):
                for r,s in combinations(cols,2):
                    def entry(u,v):
                        return identity[ids[rotate(u+v,(-cut)%4)],:]
                    equations.append(entry(p,r)+entry(q,s)-entry(p,s)-entry(q,r))
    lattice = sp.Matrix.vstack(*equations)
    colors = sp.Matrix([[int(w[i]==j) for i,j in product(K,K)] for w in words])
    pivots = colors.T.rref()[1]
    normalization = sp.Matrix.vstack(lattice,*(identity[j,:] for j in pivots))
    assert (len(words),lattice.rank(),len(pivots)) == (38,20,12)
    for mat in (lattice,colors[list(pivots),:],normalization):
        assert set(smith_normal_form(mat,domain=sp.ZZ).diagonal()) <= {-1,0,1}
    rows, fixed = normalization.rref()
    remaining = [i for i in range(len(words)) if i not in fixed]
    assert [words[i] for i in remaining] == [(3,0,3,0),(3,1,0,2),(3,1,1,1),(3,1,2,0),(3,2,0,1),(3,2,1,0)]
    values = dict(zip(remaining,(A,B,C,D,H,F)))
    for i,j in enumerate(fixed):
        values[j] = sp.prod(values[k]**(-rows[i,k]) for k in remaining)
    tensor = {w:values[i] for i,w in enumerate(words)}
    tensor.update({(1,2,1,2):X,(2,1,2,1):Y})
    exponents = [colors[list(pivots),:].T.gauss_jordan_solve(colors[j,:].T)[0] for j in remaining[:3]]
    return tensor,[words[i] for i in pivots],exponents


def schur(tensor):
    indices = ((0,3),(3,0),(1,2),(2,1))
    eq = []
    for cut in (0,1):
        m = sp.Matrix([[tensor.get(rotate(p+q,(-cut)%4),0) for q in indices] for p in indices])
        eq.extend(m[2:,2:]-m[2:,:2]*m[:2,:2].inv()*m[:2,2:])
    return list(map(sp.factor,eq))


def solve_chart(tensor):
    eq = schur(tensor)
    x = sp.solve(eq[0],X)[0]
    y = sp.solve(eq[3],Y)[0]
    assert sp.factor(x-H/C-1/A)==0 and sp.factor(y-B*D/C-D*H/(A*B))==0
    eq = [sp.fraction(sp.factor(v.subs({X:x,Y:y})))[0] for v in eq]
    d = sp.solve(eq[5],D)[0]
    h = sp.solve(eq[2].subs(D,d),H)[0]
    f = sp.solve(eq[7],F)[0]
    assert sp.factor(eq[1].subs({D:d,H:h}) - C*(A*B-C)*(A*(B-1)+C)) == 0
    branches = []
    for u,cc in ((1/B,C),(1/(1-B),sp.Integer(1))):
        sub = {A:C*u, D:sp.factor(d.subs(A,C*u)), H:sp.factor(h), F:sp.factor(f.subs({A:C*u,H:h}))}
        if u==1/(1-B):
            other = [sp.factor(eq[i].subs(sub,simultaneous=True)) for i in (4,6)]
            assert sp.factor(other[0]*(B-1)**3/(B*C*(2*B-1)**2*(C-1)))==B*C+B-C
            assert sp.factor(other[1].subs(C,B/(1-B))) != 0
            exceptional = sp.factor(other[1].subs(C,B/(1-B)))
            assert {p for p,_ in sp.factor_list(sp.fraction(exceptional)[0])[1]} <= {B,2*B-1}
        sub = {key:sp.factor(value.subs(C,cc)) for key,value in sub.items()}
        sub[C] = cc
        sub.update({X:sp.factor(x.subs(sub,simultaneous=True)),Y:sp.factor(y.subs(sub,simultaneous=True))})
        chart = {w:sp.factor(v.subs(sub,simultaneous=True)) for w,v in tensor.items()}
        assert schur(chart) == [0]*8
        branches.append(chart)
    return branches


def finite_resultant(w,cut,reverse):
    def val(i,j,k,l):
        word = (i,j,k,l)
        if reverse:
            word = tuple(3-v for v in word)
        return w.get(rotate(word,(-cut)%4),0)
    x = [1,1,Z,None]
    y = [1,None,None,None]
    y[1] = -val(1,0,2,3)/val(0,1,2,3)
    y[2] = sp.factor(-(val(1,1,1,3)*y[1]+val(2,0,1,3)*Z)/val(0,2,1,3))
    y[3] = sp.factor(-(val(1,2,0,3)*y[2]+val(2,1,0,3)*Z*y[1])/val(0,3,0,3))
    x[3] = sp.factor(-(val(1,2,3,0)*y[2]+val(2,1,3,0)*Z*y[1])/val(3,0,3,0))
    polys = []
    for total,(k,l) in ((4,(0,2)),(5,(0,1))):
        polys.append(sp.factor(sum(val(i,total-i,k,l)*x[i]*y[total-i] for i in K if 0<=total-i<4)))
    return sp.factor(sp.resultant(*polys,Z))


def verify():
    tensor,pivots,exps = monomial_chart()
    main,line = solve_chart(tensor)
    for chart,cc,disc,expected in (
        (main,C,4*B**2*C**2+4*B**2-4*B*C**2-4*B*C+C**2,
         ((1/B,1-B/C,(C-B)/(B*C)),(C/B,B,C))),
        (line,1,3*B**2-3*B+1,((1/B,1-B,1),(1/B,1-B,1)))):
        for cut,reverse in product(K,(False,True)):
            unit = sp.factor(finite_resultant(chart,cut,reverse)/disc)
            for part in sp.fraction(unit):
                assert {p for p,_ in sp.factor_list(part)[1]} <= {B,C,B-1,B-C}
        for trans,wanted in zip((lambda w:rotate(w,1),lambda w:tuple(3-v for v in w)),expected):
            def val(w):
                return chart[trans(w)]
            for word,exponent,target in zip(((3,0,3,0),(3,1,0,2),(3,1,1,1)),exps,wanted):
                actual = val(word)/sp.prod(val(p)**e for p,e in zip(pivots,exponent))
                assert sp.factor(actual-target)==0


if __name__ == '__main__':
    verify()
    print('Simple-weight chart: one surface and one curve, exact finite-point divisors, rotations and inversion.')
