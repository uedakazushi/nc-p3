"""Coefficient-independent PBW flags and exclusions on the confined parallel supports."""
from itertools import combinations, permutations, product

import sympy as s
from sympy.polys.matrices import DomainMatrix

from verify_klein_mixed_positive import bit, matrix
from verify_klein_single_missing_supports import lines
from verify_klein_sparse_transverse import primes, torus_forms
from verify_klein_three_lines import K, WORDS
from verify_klein_two_doublet_supports import candidates

PBW = (50,51,52,54,55,56,57,58,59,62,63,64,65,66,76,77,78,112,113,114,115,116,117)
EMPTY = (77,113,114,116)
# For each torus component use a starting position whose degree-four quotient is excessive.
# A negative entry means that the localized minors generate the unit ideal.
WITNESS = {5:(-1,0),13:(-1,2),23:(1,0),27:(-1,1),35:(3,0),36:(-1,2),
           38:(1,),40:(2,),41:(2,),43:(1,),44:(1,),45:(2,),46:(2,),47:(-1,2),
           49:(1,0),60:(2,),61:(3,0)}


def flag(mats):
    edges = []
    for p in K:
        edge = {}
        for pi in permutations(K):
            forbidden = [set() for _ in K]
            ambiguous = False
            for b,m in enumerate(mats[p]):
                rows,classes = lines(m)
                pivots = []
                for i in reversed(pi):
                    if not rows[i]:
                        continue
                    if pivots:
                        prev = pivots[0]
                        if classes[i]==classes[prev]:
                            continue
                        if rows[i]==rows[prev]:
                            ambiguous = True
                            break
                    forbidden[i].add(i^b)
                    pivots.append(i)
                    if len(pivots)==(2 if b<2 else 1):
                        break
                if ambiguous:
                    break
            sequence = [forbidden[i] for i in pi]+[set(K)]
            if not ambiguous and all(len(a)==j and a<=sequence[j+1] for j,a in enumerate(sequence[:-1])):
                edge[pi]=tuple(next(iter(b-a)) for a,b in zip(sequence,sequence[1:]))
        edges.append(edge)
    for pi in permutations(K):
        nxt = pi
        for edge in edges:
            nxt = edge.get(nxt)
        if nxt==pi:
            return True
    return False


def selected_columns(m,r):
    if r==1:
        return [next(m[:,j] for j in K if m[:,j]!=s.zeros(4,1))]
    for j,k in combinations(K,2):
        for i,l in combinations(K,2):
            minor = s.factor(m[i,j]*m[l,k]-m[i,k]*m[l,j])
            if minor!=0 and not any(z.is_Add for z in s.preorder_traversal(minor)):
                return [m[:,j],m[:,k]]
    raise ValueError('No pair with a Laurent unit minor')


def relaxed_fourth(w,p):
    cols = {(i,j,k,i^j^k):n for n,(i,j,k) in enumerate(product(K,repeat=3))}
    rows = []
    for l,b in product(range(3),K):
        for v in selected_columns(matrix(w,(p+l)%4,b),2 if b<2 else 1):
            for j in K:
                out = (j,j^b)
                row = [0]*64
                for i in K:
                    row[cols[out[:l]+(i,i^b)+out[l:]]]=v[i]
                rows.append(row)
    return 64-DomainMatrix.from_Matrix(s.Matrix(rows)).to_field().rank()


def main():
    flagged,confined,live = [],[],[]
    for n,mask in enumerate(candidates(1)):
        words = [t for t in WORDS if mask&bit(t)]
        mats = [[sum(1<<(4*t[p]+t[(p+2)%4]) for t in words if t[p]^t[(p+1)%4]==b) for b in K] for p in K]
        if any(len({i//2 for i in K if m&(15<<(4*i))})==1 for cut in mats for m in cut[2:]):
            confined.append(n)
        if flag(mats):
            flagged.append(n)
            for w,y,eq in torus_forms(mask,1):
                if eq==[1]:
                    continue
                if not primes(eq,y):
                    continue
                live.append(n)
                break
        elif n<67:
            forms = [(w,y,eq) for w,y,eq in torus_forms(mask,1) if eq!=[1]]
            witnesses = WITNESS.get(n,(0,))
            assert len(forms)==len(witnesses)
            for (w,y,eq),start in zip(forms,witnesses):
                if start<0:
                    assert not primes(eq,y)
                else:
                    assert relaxed_fourth(w,start)>9
    assert confined==list(range(67)) and tuple(flagged)==PBW
    assert live==[n for n in PBW if n not in EMPTY]
    print('PASS: 67 confined supports completely decided, 23 PBW flags and 19 realizable supports overall')


if __name__=='__main__':
    main()
