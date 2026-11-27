"""Coordinate-plane point orbits and the generic unmarked quotient of the period-four sheet."""
from itertools import combinations, product

import sympy as s

from verify_exchange_singular_corners import index_row
from verify_full_thin_exchange_cover import tensor
from verify_negative_component_quotient import K, relabel, thin, words
from verify_parallel_invertible_rectangles import smith

a,b,c,d = parameters = s.symbols("a b c d")
D = a*d-b*c
y = s.symbols("y0:4")


def relations():
    rows = [(0,[D,0,d-c,-c-d]),(0,[0,D,a-b,a+b]),
            (1,[D,0,-a-b,-a+b]),(1,[0,D,-c-d,-c+d]),
            (2,[1,1,1,-1]),(3,[1,1,1,1])]
    return [s.Matrix(4,4,lambda i,j:v[i] if j==i^h else 0) for h,v in rows]


def point_projection(F):
    rel = relations()
    gram = s.Matrix(16,16,lambda i,j:F.get((i//4,i%4,j//4,j%4),0))
    assert s.Matrix.hstack(gram,*(r.reshape(16,1) for r in rel)).rank()==gram.rank()==6
    E = s.Matrix.vstack(*[(r*s.Matrix(y)).T for r in rel])
    coeffs = [s.factor(s.Poly(E.extract(I,K).det(method='domain-ge'),*y).coeff_monomial(s.prod(y)))
              for I in combinations(range(6),4)]
    assert all(s.factor(v/(2*(b-c)*D)) in (0,1,-1) for v in coeffs)
    assert any(coeffs)


def actions(F):
    E = s.Matrix([index_row(q) for q in thin])
    H,L,R = smith(E)
    rank = E.rank()
    found = {0:[],2:[]}
    for shift,shear,i in product((0,2),range(2),K):
        for j,k in product((i&2,(i&2)^1),repeat=2):
            g = shift,shear,(i,j,k,i^j^k)
            G = {q:F[relabel(q,*g)] for q in words}
            rhs = L*s.Matrix([2*int(G[q]!=F[q]) for q in thin])
            if any(v%4 for v in rhs[rank:]):
                continue
            phase = R*s.Matrix([rhs[t]/H[t,t] if t<rank else 0 for t in range(16)])
            for q in words:
                G[q] *= s.I**(int(sum(phase[4*p+m] for p,m in enumerate(q)))%4)
            assert all(G[q]==F[q] for q in thin)
            A = [s.factor(G[m,m,n,n]) for m,n in product(range(2),repeat=2)]
            B = [s.factor(G[n,n^1,m,m^1]) for m,n in product(range(2),repeat=2)]
            assert B==[-A[3],A[1],A[2],-A[0]]
            p,q,r,w = s.symbols('p q r w')
            sub = {a:(p+q)/2,d:(p-q)/2,b:(r+w)/2,c:(r-w)/2}
            vals = [s.factor((u/(A[0]+A[3])).subs(sub)*p) for u in
                    (A[0]-A[3],A[1]+A[2],A[1]-A[2])]
            found[shift].append(tuple(vals))
    for shift,group in found.items():
        expected={(i*x,j*y,(-1 if shift else 1)*i*j*w)
                  for x,y in ((q,r),(r,q)) for i,j in product((1,-1),repeat=2)}
        assert len(group)==16 and set(group)==expected


if __name__=='__main__':
    F = tensor(s.Matrix([[a,b],[c,d]]),s.Matrix([[-d,b],[c,-a]]),-s.S.One,s.S.One)
    point_projection(F)
    actions(F)
    print("PASS: odd projection product coefficient, six relation rows and both unmarked action tables")
