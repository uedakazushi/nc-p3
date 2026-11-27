"""All Jordan completions with two distinct squared Bell eigenvalues."""
from collections import Counter
from functools import cache
from itertools import permutations, product

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from verify_split_inverse_windows import BITS, PAIRS, bell, reciprocal, transport
from verify_split_spectral_line import A, B, F, I, face_values, multiply, reconstruct, torus


def jordan(s, tip):
    p = bell((sp.Integer(1), sp.Integer(1), s, s))
    p[0,0,1,1] = 1
    if tip:
        p[tip] = 1
    return p


PLUS = (0,1,1,0)
MINUS = (1,0,0,1)
UNITS = (1, sp.I, -1, -sp.I)


def labeled(w):
    return {tuple(x^int(i>=2) for i,x in enumerate(q)):v for q,v in w.items()}


def decoration(w, edge):
    u,v = edge
    b,d = [j for j in range(4) if j not in edge]
    return sp.factor(w[u,v,u,v]*w[v,u,b,d]*w[b,d,v,u]/
                     (w[v,u,v,u]*w[u,v,b,d]*w[b,d,u,v]))


def diagonal_invariants(s, a, tip):
    w = labeled(reconstruct(jordan(s,tip), multiply(torus(a,1/(s**2*a)),F)))
    ordinary = {q:v for q,v in w.items() if len(set(q))==4}
    faces = face_values(labeled(ordinary))
    assert all(sp.cancel(x-y)==0 for x,y in zip(faces,(-a*s,-a/s,-a*s,-a*s**3)))
    additions = [q for q in w if q not in ordinary]
    assert all(q[:2]==q[2:] for q in additions)
    edges = {frozenset(q) for q in additions}
    expected = {frozenset((0,2))}
    if tip:
        expected.add(frozenset((0,3) if tip==PLUS else (1,2)))
    assert edges == expected
    for j in range(4):
        rotated = {q[j:]+q[:j]:v for q,v in w.items()}
        for pair in edges:
            for e in permutations(pair):
                assert decoration(rotated,e) == -1


def equations(p, flip, a, b):
    maps = multiply(torus(a,b),flip)
    q0 = transport(reciprocal(p,0),(I,)*4,PAIRS[0])
    q1 = transport(reciprocal(p,1),maps,PAIRS[1])
    return [sp.factor(v*q1[BITS[0]]-q1[q]*q0[BITS[0]]) for q,v in q0.items()]


@cache
def character_kernel(words):
    m = sp.Matrix([[int(q[i]==a) for q in words] for i,a in product(range(4),repeat=2)])
    _, _, right = smith_normal_decomp(DomainMatrix.from_Matrix(m).convert_to(sp.ZZ))
    return right.to_Matrix()[:,m.rank():]


def monomial_iso(w, v):
    def valences(tensor):
        return [Counter(q[i] for q in tensor) for i in range(4)]
    x,y = valences(w),valences(v)
    perms = [[p for p in permutations(range(4)) if all(x[i][j]==y[i][p[j]] for j in range(4))]
             for i in range(4)]
    words = tuple(sorted(w))
    for choices in product(*perms):
        target = [tuple(choices[i][j] for i,j in enumerate(q)) for q in words]
        if set(target)!=set(v):
            continue
        exponents = sp.Matrix([[UNITS.index(sp.cancel(sp.sympify(v[t])/w[q])) for t,q in zip(target,words)]])
        if all(j%4==0 for j in exponents*character_kernel(words)):
            return True
    return False


def exceptional():
    a,b = sp.symbols('a b', nonzero=True)
    representatives = []
    for tip in (PLUS,MINUS):
        p = jordan(sp.I,tip)
        for flip in ((I,)*4, A, B, F):
            eq = equations(p,flip,a,b)
            polys = [f.as_numer_denom()[0] for f in eq if f!=0]
            # Saturate the torus equations before comparing exact finite solutions.
            z = sp.Symbol('z')
            basis = sp.groebner(polys+[z*a*b-1],z,a,b,extension=sp.I)
            solutions = [] if list(basis)==[1] else sp.solve_poly_system(list(basis),z,a,b)
            points = {(sp.simplify(i),sp.simplify(j)) for _,i,j in solutions}
            expected = set()
            for e in (sp.Integer(1),sp.Integer(-1)):
                if flip==F:
                    expected.add((e,-e))
                elif tip==PLUS and flip==A:
                    expected.add((e,sp.I*e))
                elif tip==MINUS and flip==B:
                    expected.add((sp.I*e,e))
            assert points==expected
            if flip!=F:
                for i,j in points:
                    representatives.append(labeled(reconstruct(p,multiply(torus(i,j),flip))))
    w = representatives[0]
    assert all(monomial_iso(w,v) for v in representatives)
    assert monomial_iso(w, {q[1:]+q[:1]:v for q,v in w.items()})


if __name__ == '__main__':
    s,a = sp.symbols('s a',nonzero=True)
    diagonal_invariants(s,a,None)
    for e in (sp.Integer(1),sp.Integer(-1)):
        diagonal_invariants(s,e,PLUS)
        diagonal_invariants(s,e/s**2,MINUS)
    exceptional()
    print('Quantum faces, decoration ratios, eight exceptional torus systems and the mixed one-shift verified.')
