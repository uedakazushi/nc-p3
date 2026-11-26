"""Check extremal shadows and the split five-letter block ranks over a prime field."""
from itertools import combinations, combinations_with_replacement, product
from random import Random

from flint import nmod_mat

p = 101
for degree in (1,2,3):
    mons = [tuple(a.count(i) for i in range(4)) for a in combinations_with_replacement(range(4),degree)]
    for four in combinations(mons,4):
        shadow = {tuple(m[i]+(i==j) for i in range(4)) for j in range(4) for m in four}
        assert len(shadow) >= 10
        if len(shadow)==10:
            assert sum(min(m[i] for m in four) for i in range(4)) == degree-1

rng = Random(28)
K = [[[rng.randrange(p) for _ in range(4)] for _ in range(2)] for _ in range(3)]
A = [[[rng.randrange(p) for _ in range(4)] for _ in range(2)] for _ in range(3)]
left_rows = [k[0]+k[1] for k in K]
right_rows = []
for a in A:
    for shift in range(3):
        right_rows.append([a[j-shift][u] if j-shift in (0,1) else 0
                           for j,u in product(range(4),repeat=2)])
def quotient(rows):
    matrix = nmod_mat(rows,p)
    assert matrix.rank() == len(rows)
    null,dim = matrix.nullspace()
    return [[int(null[i,j]) for j in range(dim)] for i in range(matrix.ncols())]
left,right = quotient(left_rows),quotient(right_rows)
assert len(left[0]) == 5 and len(right[0]) == 7
def insertion(matrix,c,j):
    l = nmod_mat(left[4*c:4*c+4],p)
    r = nmod_mat(right[4*j:4*j+4],p)
    return (l.transpose()*nmod_mat(matrix,p)*r).entries()

restricted = [[[sum(K[h][b][z]*A[h][i][u] for h in range(3)) % p
                for u in range(4)] for z in range(4)] for b,i in product(range(2),repeat=2)]
rows = [insertion(mat,c,j) for mat,c,j in product(restricted,range(2),range(4))]
assert nmod_mat(rows,p).rank() <= 15
for _ in range(2):
    mat = [[rng.randrange(p) for u in range(4)] for z in range(4)]
    rows.extend(insertion(mat,c,j) for c,j in product(range(2),range(4)))
assert 35-nmod_mat(rows,p).rank() >= 4
print("PASS: four-monomial shadows and the projected conic block with two extra middle relations")
