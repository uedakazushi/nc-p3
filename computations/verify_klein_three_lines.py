"""Classify two (3,1,1,1) cuts, including sparse and exceptional fibers."""
from functools import lru_cache
from itertools import combinations, permutations, product
from pathlib import Path
from subprocess import check_output, run
from tempfile import TemporaryDirectory

import sympy as s

from verify_klein_mixed_supports import character, smith
from verify_klein_single_missing_supports import admissible

K = range(4)
WORDS = [(i,j,k,i^j^k) for i,j,k in product(K,repeat=3)]


def entry(p,b,i,k):
    w = (i,i^b,k,k^b)
    return w[-p:]+w[:-p] if p else w


def sparse():
    source = Path(__file__).with_name('klein_three_lines_supports.cpp')
    with TemporaryDirectory() as temp:
        binary = Path(temp)/'supports'
        run(['c++','-std=c++17','-O3',str(source),'-o',str(binary)],check=True)
        masks = set(map(int,check_output([str(binary)],text=True).splitlines()))
    assert len(masks)==5144
    def bit(w):
        i,j,k,_ = w
        return 1<<(16*(i^j)+4*i+k)
    @lru_cache(None)
    def matching(mask,indices):
        return max(sum(mask>>(4*i+pi[i])&1 for i in indices) for pi in permutations(K))
    def bound(mask,r,indices):
        return min(r,matching(mask,tuple(indices)))
    n = 0
    while masks:
        mask = min(masks)
        words = [w for w in WORDS if mask&bit(w)]
        orbit = set()
        for seq in (words,[w[::-1] for w in words]):
            for _ in K:
                seq = [w[1:]+w[:1] for w in seq]
                for pi in permutations(K):
                    orbit.add(sum(bit(tuple(pi[i] for i in w)) for w in seq))
        assert orbit<=masks
        masks -= orbit
        assert not admissible(words,lambda p,b: 1 if b else 3,bound)
        n += 1
    assert n==52


def thin_system(colors):
    thin = [w for w in WORDS if any((w[p]^w[p+1]) in colors[p] for p in range(2))]
    rows = []
    for p,group in enumerate(colors):
        for b in group:
            for i,k in product(range(1,4),repeat=2):
                r = [0]*len(thin)
                for a,j,sign in [(0,0,1),(i,k,1),(0,k,-1),(i,0,-1)]:
                    r[thin.index(entry(p,b,a,j))] += sign
                rows.append(r)
    matrix = s.Matrix(rows)
    diag,_,right = smith(matrix)
    torsion = [abs(d) for d in diag.diagonal() if d]
    gauge = s.Matrix([[int(w[p]==i) for p,i in product(K,repeat=2)] for w in thin])
    assert matrix*gauge==s.zeros(matrix.rows,16)
    coords = right.inv()*gauge
    assert coords[:len(torsion),:]==s.zeros(len(torsion),16)
    return thin,torsion,right,coords[len(torsion):,:]


def twist(w):
    return s.Integer(-1)**sum((a>>1)*(b&1) for a,b in combinations(w,2))


def dense():
    thin,torsion,right,gauge = thin_system([(1,2,3)]*2)
    assert len(thin)==60 and torsion==[1]*46+[2] and gauge.rank()==13
    signs = [twist(w) for w in thin]
    coords = [character(row,signs) for row in right.inv().tolist()[:47]]
    assert coords==[1]*46+[-1]
    assert all(twist(w)==twist(w[1:]+w[:1]) for w in WORDS)
    p = s.symbols('p:4')
    blocks = [s.Matrix([[twist((i,j,i^b,j^b)) for j in K] for i in K]) for b in K]
    assert all(abs(m.det())==16 for m in blocks)
    H = blocks[0]
    assert H**2==4*s.eye(4) and (H-2*s.eye(4)).rank()==1
    B = s.eye(4)-s.diag(*p)*H
    elementary3 = sum(s.prod(t) for t in combinations(p,3))
    assert s.expand(B.det()-(1-sum(p)+4*elementary3-16*s.prod(p)))==0
    for i in K:
        I = [j for j in K if j!=i]
        assert s.expand(B[I,I].det()-(1-sum(p[j] for j in I)+4*s.prod(p[j] for j in I)))==0
        special = {x:s.Rational(1 if i!=j else -1,2) for j,x in enumerate(p)}
        assert B.subs(special).rank()==2


if __name__=='__main__':
    sparse()
    dense()
    print('PASS: all 39461472 support candidates, 52 sparse orbits, two dense orbits, opposite rank strata')
