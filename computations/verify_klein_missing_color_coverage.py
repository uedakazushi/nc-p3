"""Exact rank-two support enumeration after ordering the free coordinate transitions."""
from itertools import product


def bit(i,e):
    return 1 << (4*(e-1)+i)


def slices(b):
    result = {}
    for mask in range(4096):
        rows = [0]*4
        cols = [0]*4
        other = [0]*3
        for i,e in product(range(4),range(1,4)):
            if mask & bit(i,e):
                k = i^b^e
                rows[i] |= 1<<k
                cols[k] |= 1<<i
                other[e-1] |= bit(i^b,b)
        if any(i>(i^b) and (not rows[i] or not cols[i]) for i in range(4)):
            continue
        parent = list(range(8))
        def root(i):
            while parent[i]!=i:
                i = parent[i]
            return i
        for i,k in product(range(4),repeat=2):
            if rows[i] and cols[k] and not rows[i]>>k&1:
                parent[root(i)] = root(k+4)
        if any(rows[i]>>k&1 and root(i)==root(k+4) for i,k in product(range(4),repeat=2)):
            continue
        if any(len({root(i+d) for i in range(4) if side[i]})<2 for d,side in ((0,rows),(4,cols))):
            continue
        result[mask] = (other,not all(rows+cols))
    return result


def main():
    tables = [slices(b) for b in (1,2,3)]
    assert [len(t) for t in tables]==[171]*3
    count = 0
    for a,(u,free_a) in tables[0].items():
        for b,(v,free_b) in tables[1].items():
            ab = tuple(x|y for x,y in zip(u,v))
            for c,(w,free_c) in tables[2].items():
                key = (ab[0]|w[0],ab[1]|w[1],ab[2]|w[2])
                if any(key[i] not in tables[i] for i in range(3)):
                    continue
                count += 1
                assert not (free_a or free_b or free_c or any(tables[i][key[i]][1] for i in range(3)))
    assert count==108
    print('PASS: 171 supports per slice, 108 simultaneous supports, no acyclic nonempty set of free edges')


if __name__=='__main__':
    main()
