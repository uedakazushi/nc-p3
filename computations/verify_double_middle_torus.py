"""Four mixed Schur complements and binary cubic resultants on the finite 0112 chart."""
import sympy as sp

K = range(4)
X,Y = sp.symbols('X Y')
ALPHA,BETA,P,Q,R,T = sp.symbols('a b p q r t', nonzero=True)
SWAP = sp.Matrix([[0,1],[1,0]])
GENERAL = sp.Matrix([[P,Q],[R,T]])


def factored(mat):
    return mat.applyfunc(sp.factor)


def forms(kind, alpha, beta):
    if kind == 'diagonal':
        pairing = sp.eye(2)
        def form(u,v):
            return sp.diag(u,v)
        bases = [sp.Matrix([[1,1],[1,h]]) for h in (alpha,beta)]
    else:
        pairing = SWAP
        def form(u,v):
            return sp.Matrix([[u,v],[v,0]])
        bases = [sp.Matrix([[0,-h],[1,1]]) for h in (1,beta)]
    bases += [pairing*m.inv().T for m in bases]
    return [form(*m[:,0]) for m in bases],[form(*m[:,1]) for m in bases]


def recurrence(A, B, first):
    mixed = [first]
    for i in K:
        mixed.append(factored(B[(i+3)%4].T*mixed[i].inv().T*A[(i+1)%4]))
    return mixed


def cubics(A,B,mixed,i):
    u = sp.Matrix([X,Y])
    v = -mixed[i].inv().T*A[(i+1)%4]*u
    z = sp.factor((u.T*B[(i+2)%4]*v)[0])
    zp = sp.factor((u.T*A[(i+2)%4]*v)[0])
    return factored(z*mixed[(i+2)%4]*v+zp*B[(i+1)%4]*u)


def resultant(polys):
    coeffs = [[sp.expand(p).coeff(X,3-i).coeff(Y,i) for i in K] for p in polys]
    rows = [[0]*k+r+[0]*(2-k) for r in coeffs for k in range(3)]
    return sp.factor(sp.Matrix(rows).det(method='domain-ge'))


def verify():
    for kind in ('diagonal','tangent'):
        A,B = forms(kind,ALPHA,BETA)
        final = recurrence(A,B,GENERAL)[-1]
        if kind=='diagonal':
            l = (BETA-1)/(ALPHA-1)
            m = ALPHA*l/BETA
            wanted = sp.diag(l,m)*GENERAL*sp.diag(m,l)
        else:
            delta = BETA-1
            left = sp.Matrix([[1,delta],[0,1]])
            right = sp.Matrix([[1,0],[-delta,1]])
            wanted = BETA**2*left*GENERAL*right
        assert factored(final-wanted).is_zero_matrix
    for kind,alpha,beta,first in [('diagonal',ALPHA,ALPHA,GENERAL),
                                  ('tangent',1,1,GENERAL),('tangent',1,-1,SWAP)]:
        A,B = forms(kind,alpha,beta)
        mixed = recurrence(A,B,first)
        assert factored(mixed[-1]-first).is_zero_matrix
        polys = cubics(A,B,mixed,0)
        assert sp.Poly(sp.gcd(*polys),X,Y).total_degree() >= 1
    A,B = forms('diagonal',ALPHA,1/ALPHA)
    mixed = recurrence(A,B,sp.eye(2))
    assert factored(mixed[-1]-sp.eye(2)).is_zero_matrix
    for low,high,opposite in ((A,B,mixed),(B,A,[mixed[(i+2)%4].T for i in K])):
        for i in K:
            res = resultant(cubics(low,high,opposite,i))
            unit = sp.factor(res/(ALPHA+1)**2)
            assert unit != 0
            for part in sp.fraction(unit):
                assert {p for p,_ in sp.factor_list(part)[1]} <= {ALPHA,ALPHA-1}


if __name__ == '__main__':
    verify()
    print('Double-middle chart: tangent and equal-ratio completions have moving points; the reciprocal family is finite.')
