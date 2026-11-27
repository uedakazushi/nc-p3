"""Six-parameter normal forms and the complete finite locus in the two quantum fibers."""
from functools import cache
from itertools import combinations, product

import sympy as sp

K = range(4)
P = [4*i+j for i,j in combinations(K, 2)]
W = (0, 3, 5, 9)
PARAMS = sp.symbols('a b c d e f')
FREE = ((-1,1,-1,1), (-1,-1,2,0), (0,-1,-1,2), (-1,-1,1,1),
        (-1,0,-1,2), (-1,-1,0,2))
LOW = ((-1,0,2,-1), (-1,2,-1,0), (1,-1,-1,1))
U, V, Z, U1, V1, Z1 = sp.symbols('u v z U V Z')
OFFSETS = (504*sp.I/78125, -7911*sp.I/12500000)


@cache
def words(h):
    return [q for q in product(K, repeat=4) if sum(W[i] for i in q) == sum(W)+h]


def cut(w, shift):
    def entry(i, j):
        q = (i//4, i%4, j//4, j%4)
        return w.get(q[-shift:]+q[:-shift] if shift else q, 0)
    return sp.Matrix(16, 16, entry)


def character(word):
    return tuple(word.count(i)-1 for i in K)


def normalized_tensor(zeta):
    crossing = {(1,2): zeta, (1,3): 1/zeta, (2,3): zeta}
    volume = {q: sp.LeviCivita(*q)*sp.prod(crossing.get((b,a),1)
              for i,a in enumerate(q) for b in q[i+1:] if a>b) for q in words(0)}
    tensor = dict(volume)
    pieces = []
    for shift in (0,1):
        m = cut(volume,shift)
        left = sp.Matrix.hstack(*m.T.nullspace()).T
        right = sp.Matrix.hstack(*m.nullspace())
        pieces.append((shift,left,right,{0:m[P,P].inv()},{0:m}))
    free = dict(zip(FREE,PARAMS))
    if zeta == sp.I:
        for j,offset in zip((3,4),OFFSETS):
            free[FREE[j]] += offset
    seen = set()
    for h in range(1,20):
        source = words(h)
        rhs = []
        for shift,left,right,inv,mats in pieces:
            contribution = sum((left*mats[i][:,P]*inv[h-i-j]*mats[j][P,:]*right
                                for i in range(1,h) for j in range(1,h-i+1)),sp.zeros(10))
            rhs.extend(contribution.applyfunc(sp.expand))
        rhs = sp.Matrix(rhs)
        term = {}
        if source:
            def column(q):
                return sp.Matrix.vstack(*[(l*cut({q:1},shift)*r).reshape(100,1)
                                         for shift,l,r,_,_ in pieces])
            linear = sp.Matrix.hstack(*(column(q) for q in source))
            piv = list(linear.rref()[1])
            rows = list(linear[:,piv].T.rref()[1])
            solution = sp.zeros(len(source),1)
            for j,value in zip(piv,linear[rows,piv].inv()*rhs[rows,:]):
                solution[j] = sp.expand(value)
            assert (linear*solution-rhs).applyfunc(sp.expand).is_zero_matrix
            kernel = linear.nullspace()
            if kernel:
                k = sp.Matrix.hstack(*kernel)
                solution = (solution-k*(k.H*k).inv()*k.H*solution).applyfunc(sp.expand)
            for vec in kernel:
                multi, = {character(q) for q,v in zip(source,vec) if v}
                if multi.count(-1)==2:
                    assert multi not in seen
                    seen.add(multi)
                    solution += vec*(1 if multi in LOW else free[multi])
                else:
                    assert sorted(multi)==[-1,0,0,1]
            term = dict(zip(source,solution))
            tensor.update(term)
        else:
            assert rhs.is_zero_matrix
        for shift,_,_,inv,mats in pieces:
            mats[h] = cut(term,shift)
            inv[h] = (-inv[0]*sum((mats[j][P,P]*inv[h-j] for j in range(1,h+1)),sp.zeros(6))).applyfunc(sp.expand)
    assert seen == set(FREE+LOW)
    return {q:c for q,c in tensor.items() if c}


def affine_cubics(tensor,shift):
    m = cut(tensor,shift)
    edge = (m[:,P]*m[P,P].inv()).applyfunc(sp.expand)
    x, y = (Z,U,V,1), (Z1,U1,V1,1)
    eq = [sp.expand(sum(x[i]*y[j]*edge[4*i+j,n] for i,j in product(K,repeat=2))) for n in range(6)]
    replacements = {}
    for j,letter in [(5,V1),(4,U1),(3,Z1),(2,Z)]:
        linear = sp.Poly(eq[j].subs(replacements),letter)
        assert linear.degree()==1 and not linear.nth(1).free_symbols
        replacements[letter] = sp.expand(-linear.nth(0)/linear.nth(1))
    return [sp.expand(q.subs(replacements)) for q in eq[:2]]


def finite_test(tensor,shift,zeta):
    t, intercept, h = sp.symbols('t intercept h')
    f,g = affine_cubics(tensor,shift)
    ff,gg = (sp.Poly(q,U,V).homogenize(h).as_expr() for q in (f,g))
    # A common curve cannot reach the transverse infinity point [1,0,0].
    jac = sp.Matrix([[q.diff(v) for v in (V,h)] for q in (ff,gg)])
    transverse = sp.expand(jac.det().subs({U:1,V:0,h:0}))
    assert transverse != 0 and not transverse.free_symbols
    f3,g3 = (sp.Poly(q.subs({U:t,V:1,h:0}),t) for q in (ff,gg))
    assert f3.degree()==2 and g3.degree()==1
    root = -g3.nth(0)/g3.nth(1)
    if zeta == -1:
        assert sp.expand(f3.eval(root)) == sp.Rational(1,2)
        return
    assert sp.expand(f3.eval(root))==0
    second = sp.Poly(g.subs(U,root*V+intercept),V).nth(2)
    line = root*V+sp.solve(second,intercept)[0]
    equations = [v for q in (f,g) for v in sp.Poly(q.subs(U,line),V).all_coeffs()]
    basis = sp.groebner(equations,*PARAMS,order='lex',extension=True)
    assert [p.as_expr().expand() for p in basis.polys] == list(PARAMS[1:])
    if shift == 0:
        a,b,c,d,e,f0 = PARAMS
        resultant = sp.Poly(sp.resultant(ff.subs({U:t,V:1}),gg.subs({U:t,V:1}),t),h)
        assert resultant.degree()==8
        expected = ([b], [b,d], [b,d,f0], [b,d,f0], [b,c*e,d,f0],
                    [a*c**2-e**2,b,c*e,d,e**3,f0],
                    [a*c**2-e**2,b,c**3,c*e,d,e**3,f0],
                    [a*c**2-e**2,b,c**3,c*e,d,e**3,f0])
        for j,gens in enumerate(expected,1):
            ideal = sp.groebner([resultant.nth(k) for k in range(1,j+1)],
                                *PARAMS,order='lex',domain=sp.QQ_I)
            assert [p.as_expr() for p in ideal.polys]==gens


def check_shift(tensor,zeta):
    # The rotation sends the background to its negative. The three unit
    # coefficients fix the compensating common diagonal transformation.
    def rotated(q):
        return q[1:]+q[:1]
    ratios = []
    for multi in LOW:
        q = next(q for q in words(1) if character(q)==multi and tensor.get(q,0))
        ratios.append(sp.simplify(-tensor[rotated(q)]/tensor[q]))
    basis = sp.Matrix([m[:3] for m in LOW]).T
    assert basis.det()==1
    exps = (2,3,3,0,0,1)
    for sign in (1,-1):
        substitute = {a:zeta**(sign*j)*a for a,j in zip(PARAMS,exps)}
        for q,value in tensor.items():
            coords = basis.inv()*sp.Matrix(character(q)[:3])
            change = sp.prod(r**(-c) for r,c in zip(ratios,coords))
            if sp.expand(-change*tensor.get(rotated(q),0)-value.subs(substitute))!=0:
                break
        else:
            return
    raise AssertionError('The rotation failed both orientation conventions.')


if __name__ == '__main__':
    for zeta in (-sp.Integer(1),sp.I):
        tensor = normalized_tensor(zeta)
        check_shift(tensor,zeta)
        for shift in K:
            finite_test(tensor,shift,zeta)
    print('Six free coordinates, four cuts, one exceptional line and all rotation eigenspaces verified.')
