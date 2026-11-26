"""Enumerate rectangular supports of the two squarefree cuts and their growth obstruction."""
from itertools import combinations, permutations, product

import sympy as sp


def blocks():
    words = list(permutations(range(4)))
    pairs = list(combinations(range(4),2))
    result = []
    for p in range(4):
        result.append([])
        for pair in pairs:
            other = sorted(set(range(4))-set(pair))
            block = []
            for ab,cd in product(permutations(pair),permutations(other)):
                local = ab+cd
                key = local[-p:]+local[:-p] if p else local
                block.append(1<<words.index(key))
            result[-1].append(block)
    return result,pairs


def main():
    groups,pairs = blocks()
    rectangles = (1,2,4,8,3,12,5,10,15)
    def masks(group):
        return [sum(group[i] for i in range(4) if rect>>i&1) for rect in rectangles]
    choices = [masks(group) for group in groups[0]]
    tests = [(sum(group),set(masks(group))) for group in groups[1]]
    directed = [(group[2*s]+group[2*s+1],pair[s],pair[1-s])
                for phase in groups for pair,group in zip(pairs,phase) for s in (0,1)]
    valid = 0
    survivors = []
    for bits in product(*choices):
        w = sum(bits)
        if any(w&mask not in allowed for mask,allowed in tests):
            continue
        valid += 1
        rows = [0]*4
        for mask,a,b in directed:
            if not w&mask:
                rows[a] |= 1<<b
        for via in range(4):
            for a in range(4):
                if rows[a]>>via&1:
                    rows[a] |= rows[via]
        if all(not row>>a&1 for a,row in enumerate(rows)):
            survivors.append(w)
    assert valid==22785 and survivors==[(1<<24)-1]

    # With full permutation support, a three-minor forces reciprocal diagonal entries to vanish together.
    a,b,c,d,e,f = sp.symbols('a b c d e f')
    minor = sp.Matrix([[a,0,c],[0,b,d],[e,f,0]]).det()
    assert sp.expand(minor+a*d*f+b*c*e)==0
    graphs = set()
    for bits in product((0,1),repeat=6):
        edges = {p for p,v in zip(pairs,bits) if v}
        if any(sum(i in p for p in edges)==3 for i in range(4)):
            continue
        if not any(set(p).isdisjoint(q) for p,q in combinations(edges,2)):
            continue
        graphs.add(tuple(sorted(sum(i in p for p in edges) for i in range(4))))
    assert graphs=={(1,1,1,1),(1,1,2,2),(2,2,2,2)}
    print('PASS: 531441 supports, 22785 compatible cuts, only the full permutation orbit survives')


if __name__=='__main__':
    main()
