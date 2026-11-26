"""Support partitions and the descending mixed-pair pivots for the unmarked Klein chart."""
from itertools import product

import sympy as sp


def main():
    labels = range(4)
    partitions = {tuple(t.index(x) for x in t) for t in product(labels,repeat=4)}
    assert len(partitions)==15
    for b in (1,2,3):
        for t in partitions:
            if len(set(t))==1:
                continue
            def present(i,k):
                return t[i]!=t[k^b]
            assert all(present(i,k)==present(k^b,i^b) for i,k in product(labels,repeat=2))
            for e in set((1,2,3))-{b}:
                if all(not present(i,i^b^e) for i in labels):
                    assert all(present(i,i) for i in labels)
    for i,b,e in product(labels,range(1,4),range(1,4)):
        word = (i,i^b,i^b^e,i^e)
        reflected = {word,word[::-1],(word[1],word[0],word[3],word[2]),word[2:]+word[:2]}
        if b!=e:
            assert reflected == {tuple(x^j for x in word) for j in labels}
        else:
            assert len(reflected)==2
    for c in (1,2,3):
        small = {0,c}
        j,k = sorted(set(labels)-small)
        assert j^k==c
        for b in set((1,2,3))-{c}:
            # Mixed BA rows in the cut; use the two columns (j+b,k+b).
            colors = [[a^b^z for z in (j^b,k^b)] for a in (j,k)]
            assert colors==[[0,c],[c,0]]
    x,y = sp.symbols('x y')
    assert sp.Matrix([[0,x],[y,0]]).det()==-x*y
    print('PASS: all binary line partitions, both support reflections, and all complementary PBW pivots')


if __name__=='__main__':
    main()
