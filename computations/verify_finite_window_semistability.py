"""Verify the sharp weight bounds, mixed support obstruction and uniform PBW pivots."""
from itertools import product


def main():
    for d,g,s,t in product(range(6),repeat=4):
        if not d or s>t:
            continue
        mean = 3*d+s+t+g
        upper = (d, d+min(d,s), d+s+min(d,g), min(d+t+g,2*d+s+t),
                 2*d+t+g,2*d+s+t+g)
        compatible = all(upper[j]+upper[-j-1]>=mean for j in range(6))
        assert compatible == (s<=d<=g<=d+s)
        if not compatible:
            continue
        assert all(upper[j]+upper[-j-1]==mean for j in range(6))
        a,b = (0,d,d+s,d+s+g),(0,d,d+t,d+t+g)
        if s<t:
            sub = [(i,j) for i,j in product(range(4),repeat=2) if j<3 and a[i]+b[j]>=upper[4]]
            if s<g:
                assert sub==[(3,2)]
            else:
                assert a[3]+b[1] in upper if t==2*d else a[3]+b[1] not in upper
            continue
        if s==0:
            assert g==d
            continue
        for level,pivot in zip(upper,( (1,0),(2,0),(2,1),(3,0),(3,1),(3,2) )):
            block = [(i,j) for i,j in product(range(4),repeat=2) if a[i]+b[j]==level]
            assert pivot in block
            if g!=d or level!=upper[2]:
                assert max(block)==pivot
            else:
                assert block==[(0,3),(1,2),(2,1),(3,0)]
    # Every orientation of 02 and 13 occurs in the adjacent cut of a word
    # with an extreme 01 pair and the complementary 23 pair.
    pairs = set()
    for ab,cd in product(((0,1),(1,0)),((2,3),(3,2))):
        pairs.update(((ab[1],cd[0]),(cd[1],ab[0])))
    assert {(0,2),(2,0),(1,3),(3,1),(0,3),(3,0),(1,2),(2,1)}==pairs
    print('Extremal equality bounds, mixed-profile contradiction and six PBW leading pairs verified.')


if __name__ == '__main__':
    main()
