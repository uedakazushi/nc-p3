"""Check the linear ideal filtrations in k^4 + M_2 for two opposite four-arrow spaces."""

from itertools import product
from random import Random

from flint import nmod, nmod_mat

P = 53
K = range(4)
PAIRS = ((0,1),(2,3))


def rowspace(vectors):
    mat = nmod_mat(vectors,P) if vectors else nmod_mat(0,8,P)
    reduced,rank = mat.rref()
    return reduced.tolist()[:rank]


def product_row(a,b):
    x,y = nmod_mat(2,2,a[4:],P),nmod_mat(2,2,b[4:],P)
    return [a[i]*b[i] for i in K]+(x*y).entries()


def point_rows(tangent):
    rng = Random(227)
    while True:
        basis = [[rng.randrange(P) for _ in K] for _ in K]
        if tangent:
            basis[:2] = [[1,0,0,1],[1,1,0,1]]
        x = nmod_mat(basis,P)
        if not x.det():
            continue
        y = 2*x.transpose().inv()
        sides = [[[nmod(int(i==j),P) for j in K]+list(t) for i,t in enumerate(rows)]
                 for rows in (x.tolist(),[[r[0],r[2],r[1],r[3]] for r in y.tolist()])]
        if any(not nmod_mat(2,2,r[4:],P).det() for side in sides for r in side):
            continue
        roots = [[[] for _ in PAIRS] for _ in sides]
        for p,side in enumerate(sides):
            for a,(i,j) in enumerate(PAIRS):
                for c in range(1,P):
                    r = [s+c*t for s,t in zip(side[i],side[j],strict=True)]
                    if not nmod_mat(2,2,r[4:],P).det():
                        roots[p][a].append(r)
        if all(secant for side in roots for secant in side) and (not tangent or len(roots[0][0])==1):
            return sides,roots


def verify(tangent):
    sides,roots = point_rows(tangent)
    full = [[nmod(int(i==j),P) for j in range(8)] for i in range(8)]
    def span_product(a,b):
        return rowspace([product_row(s,t) for s,t in product(a,b)])
    def same(a,b):
        return rowspace(a)==rowspace(b)
    halves = [span_product(sides[p],sides[1-p]) for p in (0,1)]
    assert len(halves[0]) == 7 and same(*halves)
    hyper = halves[0]
    assert all(2*sum(r[:4]) == r[4]+r[7] for r in hyper)
    for p in (0,1):
        start,nxt = sides[p],sides[1-p]
        assert same(span_product(hyper,start),full)
        for a,secant in enumerate(roots[p]):
            x = secant[0]
            x1 = start[PAIRS[a][0]]
            xm,xm1 = (nmod_mat(2,2,row[4:],P) for row in (x,x1))
            assert any(not any(product_row(x,r)) for r in roots[1-p][1-a])
            def correct_image(row):
                prod = xm1*nmod_mat(2,2,row[4:],P)
                return nmod_mat([list(aa)+list(bb) for aa,bb in zip(xm.tolist(),prod.tolist(),strict=True)],P).rank()==1
            y = next(r for r in roots[1-p][1-a] if correct_image(r))
            ym = nmod_mat(2,2,y[4:],P)
            # A second generator of the same matrix column ideal, with full scalar support.
            bs = nmod_mat([r[4:] for r in nxt],P).inv()
            y1 = None
            for t in range(P):
                mat = ym*nmod_mat([[1,t],[2*t,1]],P)
                co = (nmod_mat([mat.entries()],P)*bs).entries()
                if all(co) and len(rowspace([y,co+mat.entries()]))==2:
                    y1 = co+mat.entries()
                    break
            assert y1 is not None
            for first,second,side in ((x,x1,nxt),(y,y1,start)):
                assert [len(span_product([first],b)) for b in (side,hyper,full)] == [3,4,4]
                assert [len(span_product([first,second],b)) for b in (side,hyper,full)] == [5,6,6]
        for j in K:
            sub = [r for i,r in enumerate(start) if i!=j]
            assert [len(span_product(sub,b)) for b in (nxt,hyper,full)]==[6,7,7]


if __name__ == "__main__":
    verify(False)
    verify(True)
    print("PASS: both linear filtrations, complementary secants and tangents, principal annihilators and all tail ranks")
