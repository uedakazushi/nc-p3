"""Exact companion completions, smooth collision and a split twenty-point fiber."""
from itertools import combinations, combinations_with_replacement, product
import sympy as sp

INDICES = range(4)


def cross_blocks(u, v):
    x = sp.Matrix(2, 2, lambda a, c: v[c, a]*u[a, c ^ 1])
    y = sp.Matrix(2, 2, lambda c, a: u[a, c]*v[c, a ^ 1])
    assert y.det() == -x.det()
    return x, y


def other_completion(a, x, y):
    b = y*a.inv()*x
    i, j, t = sp.symbols('i j t')
    candidate = sp.Matrix([[i, a[0,1]], [a[1,0], j]])
    error = y*candidate.adjugate()*x - t*b
    coeff, rhs = sp.linear_eq_to_matrix([error[0,1], error[1,0]], [i, j])
    i1,j1 = coeff.inv()*rhs
    pencil = candidate.subs({i:i1, j:j1})
    quotient,remainder = sp.div(sp.Poly(pencil.det()-t,t), sp.Poly(t-a.det(),t))
    assert remainder.is_zero and quotient.degree() == 1
    return pencil.subs(t,-quotient.nth(0)/quotient.nth(1))


def gram(u,v,a,second):
    x,y = cross_blocks(u,v)
    def thick(block):
        return block.row_join(x).col_join(y.row_join(y*block.inv()*x))
    mats = [thick(a),thick(second).T]
    def thin(i,j):
        return u[i,j-2] if i < 2 else v[i-2,j]
    def value(i,j,k,l):
        if i ^ j ^ k ^ l:
            return 0
        h = i ^ j
        return mats[h][i,k] if h < 2 else thin(i,j)*thin(k,l ^ 1)
    w = sp.Matrix(16,16,lambda i,j:value(i//4,i%4,j//4,j%4))
    assert w.rank() == 6
    assert all(value(i,j,k,l) == value(k,j,i,l) for i,j,k,l in product(INDICES,repeat=4))
    return w


def opposite(w, symmetric):
    even = list(combinations_with_replacement(INDICES,2))
    odd = even if symmetric else list(product(INDICES,repeat=2))
    return sp.Matrix([[w[4*i+j,4*k+(l^1)] for j,l in odd] for i,k in even])


def check_collision():
    u = sp.Matrix([[1,2],[2,1]])
    v = sp.Matrix([[3,5],[5,3]])
    a = sp.Matrix([[2,3],[sp.Rational(44,87),7]])
    x,y = cross_blocks(u,v)
    assert other_completion(a,x,y) == a
    m = opposite(gram(u,v,a,a),True)
    assert m.rank() == 8
    pairs = list(combinations_with_replacement(INDICES,2))
    diag = [pairs.index((i,i)) for i in INDICES]
    for matrix in [m,m.T]:
        pencil = sp.Matrix.hstack(*matrix.nullspace()).T
        assert all(pencil[:,i].is_zero_matrix for i in range(10) if i not in diag)
        assert all(pencil[:,[i,j]].det() for i,j in combinations(diag,2))


def modular_points(w, prime):
    coords = sp.symbols('x:4')
    ev = sp.Matrix.vstack(*[(sp.Matrix(4,4,col).T*sp.Matrix(coords)).T
                            for col in w.columnspace()])
    minors = [sp.Poly(ev[list(rows),:].det(method='domain-ge'),*coords)
              for rows in combinations(range(6),4)]
    integral = []
    for poly in minors:
        denom, numerator = poly.clear_denoms(convert=True)
        assert denom % prime
        integral.append(numerator.as_expr())
    def triangular(equations, variables):
        if not variables:
            return [()] if all(int(f)%prime == 0 for f in equations) else []
        g = sp.groebner(equations,*variables,modulus=prime,order='lex')
        if list(g) == [1]:
            return []
        assert g.is_zero_dimensional
        last = variables[-1]
        eliminant = sp.Poly(g.polys[-1].as_expr(),last,modulus=prime)
        coefficients = [int(c) for c in eliminant.all_coeffs()]
        answers = []
        for root in range(prime):
            value = 0
            for c in coefficients:
                value = (value*root+c)%prime
            if value:
                continue
            factor = sp.Poly(last-root,last,modulus=prime)
            while eliminant.degree() and eliminant.eval(root)==0:
                eliminant = eliminant.exquo(factor)
            rest = [eq.subs(last,root) for eq in g]
            answers.extend(pt+(root,) for pt in triangular(rest,variables[:-1]))
        assert eliminant.degree()==0  # All geometric roots were found over the prime field.
        return answers
    points = []
    for i in INDICES:
        prefix = (0,)*i+(1,)
        eq = [f.subs(dict(zip(coords,prefix))) for f in integral]
        points.extend(prefix+pt for pt in triangular(eq,coords[i+1:]))
    return points


def trisecants(points,p):
    def minor(a,b,c,i,j,k):
        return (a[i]*(b[j]*c[k]-b[k]*c[j])-a[j]*(b[i]*c[k]-b[k]*c[i])
                +a[k]*(b[i]*c[j]-b[j]*c[i]))%p
    return [triple for triple in combinations(points,3)
            if all(minor(*triple,*ijk)==0 for ijk in combinations(INDICES,3))]


if __name__ == '__main__':
    check_collision()
    u = sp.Matrix([[1,2],[3,5]])
    v = sp.Matrix([[2,3],[5,7]])
    a = sp.Matrix([[2,1],[5,3]])
    second = other_completion(a,*cross_blocks(u,v))
    assert second == sp.Matrix([[sp.Rational(123,196),1],[5,sp.Rational(1335,176)]])
    w = gram(u,v,a,second)
    assert opposite(w,False).rank() == 10
    prime = 12619
    pts = modular_points(w,prime)
    assert len(set(pts)) == len(pts) == 20
    assert not trisecants(pts,prime)
    print('Smooth rank-eight collision; other sheet rank ten. Twenty split points, collinear triples:',
          len(trisecants(pts,prime)))
