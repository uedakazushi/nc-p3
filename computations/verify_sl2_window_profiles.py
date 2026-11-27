"""All SL2 arrow and relation characters, and the singular repeated-standard exclusions."""
from collections import Counter
from itertools import permutations, product, combinations_with_replacement
from math import prod

import sympy as sp

from flint_rational_function import rational_function_field
from verify_klein_center_lift import echelon

K = range(4)
PERMS = [(p,(-1)**sum(p[i]>p[j] for i in K for j in range(i+1,4))) for p in permutations(K)]
EPS = sp.Matrix([[0,1],[-1,0]])


def partitions(n,low=1):
    if not n:
        yield ()
    for i in range(low,n+1):
        for tail in partitions(n-i,i):
            yield (i,)+tail


def character(rep,t,prime):
    return sum(pow(t,j,prime) for dim in rep for j in range(1-dim,dim,2))%prime


def tensor_product(left,right):
    return Counter(k for a in left for b in right for k in range(abs(a-b)+1,a+b,2))


def euler(arrows,relations,t,y,prime):
    a = [character(r,t,prime) for r in arrows]
    matrix = []
    for i in K:
        terms = (1+y,-a[i],character(relations[i%2],t,prime),-a[(i+3)%4])
        matrix.append([terms[(j-i)%4]*(y if j<i else 1) for j in K])
    return sum(sign*prod(matrix[i][p[i]] for i in K) for p,sign in PERMS)%prime


def check_characters():
    slots, rels = list(partitions(4)), list(partitions(6))
    prime,t,y = 1000003,7,11
    possible = {prod(1-y*pow(t,n,prime) for n in (i,-i,j,-j))%prime:(i,j)
                for i,j in combinations_with_replacement(range(13),2)}
    found = []
    for k in product(range(5),repeat=4):
        if k != min(k[i:]+k[:i] for i in K):
            continue
        arrows = [slots[i] for i in k]
        cuts = [tensor_product(arrows[i],arrows[(i+1)%4]) for i in K]
        choices = [[r for r in rels if Counter(r)<=cuts[i] and Counter(r)<=cuts[i+2]] for i in (0,1)]
        for pair in product(*choices):
            powers = possible.get(euler(arrows,pair,t,y,prime))
            if powers is not None:
                found.append((k,pair,powers))
    expected = [((1,1,1,1,1,1),(0,0)),((1,1,2,2),(0,4)),((3,3),(0,8)),
                ((1,1,1,3),(4,4)),((1,5),(4,12))]
    assert found == [((i,)*4,(r,r),p) for i,(r,p) in enumerate(expected)]


def repeated_tensor(maps,a,b):
    r = [EPS*m.T for m in maps]
    def value(i,j,k,l):
        u,v,w,x = (h//2 for h in (i,j,k,l))
        p,q,s,t = (h%2 for h in (i,j,k,l))
        return (a*r[1][v,w]*r[3][x,u]*EPS[p,q]*EPS[s,t]
                +b*r[0][u,v]*r[2][w,x]*EPS[q,s]*EPS[t,p])
    return {q:value(*q) for q in product(K,repeat=4)}


def flatten(w,start):
    def value(i,j):
        q = (i//4,i%4,j//4,j%4)
        return w[q[-start:]+q[:-start] if start else q]
    return sp.Matrix(16,16,value)


def fourth_dimension(w,start,t):
    rel = []
    for i in range(3):
        mat = flatten(w,(start+i)%4)
        rel.append([{divmod(j,4):v[j] for j in range(16) if v[j]}
                    for v in mat.columnspace()])
    field = rational_function_field((str(t),))
    dimension = 0
    for height in range(5):
        words = [q for q in product(K,repeat=4) if sum(i%2 for i in q)==height]
        ids = {q:i for i,q in enumerate(words)}
        rows = []
        for i,edge in enumerate(rel):
            for q in product(K,repeat=2):
                for r in edge:
                    if sum(j%2 for j in q+next(iter(r)))==height:
                        rows.append({ids[q[:i]+v+q[i:]]:field.from_sympy(c) for v,c in r.items()})
        dimension += len(words)-len(echelon(rows))
    return dimension


def check_repeated():
    unit,line = sp.eye(2),sp.diag(1,0)
    a,b,c,d,A,B,t = sp.symbols('a b c d A B t')
    m = sp.Matrix([[a,b],[c,d]])
    r = [EPS*g.T for g in (unit,unit,unit,m)]
    for start in (0,1):
        def value(i,j):
            q=(i//2,i%2,j//2,j%2)
            x,y,z,w = q[-start:]+q[:-start] if start else q
            return ((2-start)*A*r[1][y,z]*r[3][w,x]+(1+start)*B*r[0][x,y]*r[2][z,w])
        observed = sp.Matrix(4,4,value).det()
        expected = (-8*A**3*m.det()*(2*A*m.det()+B*m.trace()) if start==0
                    else -8*B**3*(A*m.trace()+2*B))
        assert sp.expand(observed-expected)==0
    w = repeated_tensor((unit,unit,unit,line),1,-sp.Rational(1,2))
    assert fourth_dimension(w,1,t) == 39
    for g,rank in ((line,6),(sp.diag(0,1),5),(sp.Matrix([[0,1],[0,0]]),5),
                   (sp.Matrix([[0,0],[1,0]]),6)):
        w = repeated_tensor((unit,unit,line,g),1,t)
        assert flatten(w,0).rank()==flatten(w,1).rank()==rank
        if rank==6:
            assert fourth_dimension(w,1,t) == 52
    for m,rank in ((unit,6),(sp.Matrix([[1,1],[0,1]]),12)):
        w = repeated_tensor((unit,unit,unit,m),1,-1)
        assert sp.Matrix(16,16,lambda i,j:w[i//4,j//4,i%4,j%4]).rank()==rank


def adjoint_tensor_rank():
    def value(q):
        omitted = [i for i,v in enumerate(q) if v == 0]
        if len(omitted)==1:
            i = omitted[0]
            return sp.LeviCivita(*(q[(i+j)%4] for j in (1,2,3)))
        if omitted==[0,2]:
            return int(q[1]==q[3])
        if omitted==[1,3]:
            return -int(q[0]==q[2])
        return 0
    opposite = sp.Matrix(16,16,lambda i,j:value((i//4,j//4,i%4,j%4)))
    assert opposite.rank()==8


if __name__ == '__main__':
    check_characters()
    check_repeated()
    adjoint_tensor_rank()
    print('Five uniform SL2 profiles. Singular multiplicity matrices fail ranks or the quartic Hilbert function.')
