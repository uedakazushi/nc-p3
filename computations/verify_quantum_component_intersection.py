"""Quantum windows of canonical corner PI degree four with five quartic central sections."""

from itertools import combinations, permutations, product
from math import gcd

K = range(4)
TRIPLES = tuple(combinations(range(1,4),2))


def commutators(triangle):
    q = [[0]*4 for _ in K]
    for (i,j),v in zip(TRIPLES,triangle,strict=True):
        q[i][j],q[j][i] = v,-v
    return q


def triangles(q, p=tuple(K)):
    return tuple((q[p[0]][p[i]]+q[p[i]][p[j]]+q[p[j]][p[0]])%4 for i,j in TRIPLES)


def center(triangle):
    q = commutators(triangle)
    result = []
    for e in product(range(5),repeat=4):
        if sum(e)!=4:
            continue
        # Pair the degree-zero exponent with each e_j-e_0.
        characters = {sum((e[i]-1)*q[i][j] for i in K)%4 for j in K}
        if len(characters)==1:
            result.append(e)
    return result


def verify():
    remaining = {t for t in product(range(4),repeat=3) if gcd(4,*t)==1}
    table = []
    while remaining:
        t = min(remaining)
        q = commutators(t)
        orbit = {triangles(q,p) for p in permutations(K)}
        assert orbit <= remaining
        counts = {len(center(a)) for a in orbit}
        assert len(counts)==1
        table.append((t,len(orbit),counts.pop()))
        remaining -= orbit
    assert table == [((0,0,1),12,3),((0,1,2),24,2),((1,1,1),6,3),
                     ((1,1,2),12,1),((1,3,1),2,5)]
    quantum = [[0]*4 for _ in K]
    for (i,j),e in zip(combinations(K,2),(0,0,3,3,0,2),strict=True):
        quantum[i][j],quantum[j][i] = e,-e
    assert (1,3,1) in {triangles(quantum,p) for p in permutations(K)}
    assert len(center((1,3,1))) == 5


if __name__=="__main__":
    verify()
    print("PASS: 56 primitive degree-four triples, five permutation orbits, one surviving indexed class")
