"""A smooth ordinary Pauli point on the unrestricted fifteen-dimensional component."""

from itertools import product

from flint import nmod, nmod_mat

from verify_klein_elliptic_completion import coeffs
from verify_rank_fifteen_root_exchange import coefficients


P = 53
K = range(4)
WORDS = list(product(K, repeat=4))


def vector(values):
    first = next(x for x in values if x % P)
    return tuple(int(x) * pow(int(first), -1, P) % P for x in values)


def matrices():
    matrix = [[4, 5, 16], [23, -3, 19], [-21, 4, -1]]
    values = {(i, i ^ b, i ^ b ^ e, i ^ e): (-1)**((i >> 1)*((b ^ e) & 1))*matrix[b][e]
              for i, b, e in product(K, range(3), range(3))}
    def cut(index):
        return nmod_mat([[values.get(index(i,j,k,l),0) for k,l in product(K,repeat=2)]
                         for i,j in product(K,repeat=2)], P)
    return cut(lambda i,j,k,l: (i,j,k,l)), cut(lambda k,j,i,l: (i,j,k,l))


def point_certificate(first, second):
    basis, n = first.transpose().rref()
    assert n == 6
    relations = [nmod_mat(4,4,[int(basis[j,i]) for i in range(16)],P) for j in range(6)]
    points = {tuple(int(i == t) for i in K) for t in K}
    seed = (1,3,10,2)
    points.update(vector([(-1)**((c & i).bit_count()) * seed[i ^ x] for i in K])
                  for x,c in product(K,repeat=2))
    assert len(points) == 20
    conormal = []
    for a in sorted(points):
        av = nmod_mat([list(a)],P)
        left = nmod_mat([list((av*r).entries()) for r in relations], P)
        ker, d = left.nullspace()
        assert d == 1
        b = vector([int(ker[i,0]) for i in K])
        bv = nmod_mat([list(b)],P).transpose()
        jacobian = nmod_mat([list((r*bv).entries()) + list((av*r).entries()) for r in relations], P)
        assert jacobian.rank() == 6
        jb = nmod_mat([[int(k==r)*b[j] for r in K] for k,j in product(K,repeat=2)], P)
        ja = nmod_mat([[int(l==r)*a[i] for r in K] for i,l in product(K,repeat=2)], P)
        y, x = jb.transpose()*second, second*ja
        c, nc = y.transpose().nullspace()
        d, nd = x.nullspace()
        assert nc == nd == 1 and (y*ja).rank() == 0
        conormal.append([int(c[k,0])*b[j]*a[i]*int(d[l,0]) for k,j,i,l in product(K,repeat=4)])
    assert nmod_mat(conormal,P).rank() == 20


def tangent_characters(first, second):
    equations = []
    for side, matrix in enumerate((first,second)):
        left,n = matrix.transpose().nullspace()
        right,m = matrix.nullspace()
        assert n == m == 10
        for a,b in product(range(n),range(m)):
            def value(i,j,k,l):
                if side:
                    i,k = k,i
                return int(left[4*i+j,a]*right[4*k+l,b])
            equations.append({w:value(*w) for w in WORDS})
    result = []
    for total, character in product(K,repeat=2):
        representatives = [(0,j,k,total ^ j ^ k) for j,k in product(K,repeat=2)]
        rows = [[sum((-1)**((character & x).bit_count() + (x >> 1)*(w[2] & 1))
                     * row[tuple(a ^ x for a in w)] for x in K) for w in representatives]
                for row in equations]
        result.append(16 - nmod_mat(rows,P).rank())
    assert result == [4]*4 + [5]*12
    assert sum(result) == 76


def elliptic_chart():
    s,t,n = [nmod(i,P) for i in (25,35,28)]
    z, lam = s**4,(t/s)**2
    k,m = coefficients(z,lam)
    a,b,e,f = 1+n,k+m/n,1-n,k-m/n
    _, px, qy = coeffs(a,b,z*b,a,z,lam)
    x,y = a*a-z*f*f,e*e-z*b*b
    assert all(v != 0 for v in (a*a-z*b*b,e*e-z*f*f,z*z-1,lam,z,px[0],qy[0]))
    assert px[0]*x*x+px[1]*x+px[2] == qy[0]*y*y+qy[1]*y+qy[2] == 0
    assert 2*px[0]*x+px[1] == -(2*qy[0]*y+qy[1]) == 7


def verify_period_branches():
    lattice = [(h, c) for h in (1, 2, 3) for c in K]
    flips = {(h,c): (h,c ^ (((h & 1) << 1) + (h >> 1))) for h,c in lattice}
    exponents = []
    for x,z in product(K,repeat=2):
        supported = []
        for h,c in lattice:
            _, other = flips[h,c]
            if c >= other:
                continue
            first = ((x & c).bit_count() + (z & h).bit_count()) % 2
            second = ((x & other).bit_count() + (z & h).bit_count()) % 2
            if first == second:
                supported.append((h,c,first))
        exponents.append(tuple(supported))
        if x:
            assert {h for h,c,sign in supported} == {x}
            assert {sign for h,c,sign in supported} == {0,1}
    distinct = set(exponents)
    assert len([branch for branch in distinct if len(branch) == 6]) == 4
    assert len([branch for branch in distinct if len(branch) == 2]) == 6


if __name__ == "__main__":
    first, second = matrices()
    elliptic_chart()
    point_certificate(first, second)
    tangent_characters(first, second)
    verify_period_branches()
    print("PASS: six swapped pairs, one nine-dimensional and three five-dimensional quotient germs")
