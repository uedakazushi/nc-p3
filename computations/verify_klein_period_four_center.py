"""Rational four-periodic path spaces and the quartic center at the doublet sample."""
from functools import cache, reduce
from itertools import product
from math import gcd
from operator import xor
import sympy as sp
from verify_klein_two_doublets import cross_blocks, gram, other_completion

K = range(4)


def sample_relations():
    u,v = sp.Matrix([[1,2],[3,5]]),sp.Matrix([[2,3],[5,7]])
    a = sp.Matrix([[2,1],[5,3]])
    w = gram(u,v,a,other_completion(a,*cross_blocks(u,v)))
    return relations_from_gram(w)


def relations_from_gram(w):
    flip = [4*j+i for i,j in product(K,repeat=2)]
    relations = []
    for matrix in (w,w.T):
        rel = sp.Matrix.hstack(*matrix.columnspace()).T
        homogeneous = [(i^j, row) for row in rel.tolist()
                       for i,j in [divmod(next(k for k,x in enumerate(row) if x),4)]]
        relations += [homogeneous,[(h,[row[j] for j in flip]) for h,row in homogeneous]]
    return relations


def paths(relations,maximum):
    basis = {(p,0,h): [()] if h==0 else [] for p,h in product(K,repeat=2)}
    prepend = {}
    @cache
    def word_image(p,word):
        if not word:
            return sp.ones(1,1)
        h = reduce(xor,word,0)
        return prepend[p,len(word),h,word[0]]*word_image((p+1)%4,word[1:])
    for n in range(1,maximum+1):
        for p,h in product(K,repeat=2):
            cols = [(i,)+w for i in K for w in basis[(p+1)%4,n-1,h^i]]
            offsets = [sum(w[0]<i for w in cols) for i in range(5)]
            rows = []
            if n>=2:
                for g,row in relations[p]:
                    for tail in basis[(p+2)%4,n-2,g^h]:
                        vector = sp.zeros(1,len(cols))
                        for index,value in enumerate(row):
                            if value:
                                i,j = divmod(index,4)
                                vector[:,offsets[i]:offsets[i+1]] += value*word_image((p+1)%4,(j,)+tail).T
                        rows.append(vector)
            mat = sp.Matrix.vstack(*rows) if rows else sp.zeros(0,len(cols))
            reduced,piv = mat.to_DM().rref()
            free = [i for i in range(len(cols)) if i not in piv]
            q = sp.zeros(len(free),len(cols))
            for a,i in enumerate(free):
                q[a,i]=1
            if piv and free:
                rref = reduced.to_Matrix()
                for j,col in enumerate(piv):
                    q[:,col] = -rref[j,free].T
            basis[p,n,h]=[cols[i] for i in free]
            for i in K:
                prepend[p,n,h,i]=q[:,offsets[i]:offsets[i+1]]
    return basis,word_image


def commuting_space(basis,image,h,char,period=1):
    dims = [len(basis[p,4,h]) for p in K]
    starts = [sum(dims[:p]) for p in range(5)]
    rows = []
    for p,i in product(K,repeat=2):
        z = (p+1)%4
        left = sp.Matrix.hstack(*(image(p,w+(i,)) for w in basis[p,4,h]))
        right = sp.Matrix.hstack(*(image(p,(i,)+w) for w in basis[z,4,h]))
        block = sp.zeros(left.rows,starts[-1])
        block[:,starts[p]:starts[p+1]] = left
        block[:,starts[z]:starts[z+1]] = -((-1)**((char&i).bit_count()))*(period if z==0 else 1)*right
        rows.append(block)
    return sp.Matrix.vstack(*rows).to_DM().nullspace().to_Matrix().T


def artinian_quotient(relations,basis,image):
    prime = 1009
    def residue(x):
        return int(x.p)*pow(int(x.q),-1,prime)%prime
    rules = [[(h,{divmod(i,4):residue(x) for i,x in enumerate(row) if x})
              for h,row in rs] for rs in relations]
    for h,choices in [(0,[[1,3,5],[2,7,11]]),(1,[[1,0],[0,1]])]:
        center = commuting_space(basis,image,h,0)
        for coeff in choices:
            q = center*sp.Matrix(coeff)
            q *= sp.ilcm(*(x.q for x in q))
            q /= reduce(gcd,map(int,q))
            start = 0
            for p in K:
                words = basis[p,4,h]
                rule = {w:residue(q[start+j]) for j,w in enumerate(words) if q[start+j]}
                rules[p].append((h,rule))
                start += len(words)

    bases = {(p,0,h):[()] if h==0 else [] for p,h in product(K,repeat=2)}
    multiply = {}
    @cache
    def evaluate(p,word):
        if not word:
            return [1]
        mat = multiply[p,len(word),reduce(xor,word,0),word[0]]
        tail = evaluate((p+1)%4,word[1:])
        return [sum(x*y for x,y in zip(row,tail))%prime for row in mat]
    expected = [1,4,10,20,31,40,44,40,31,20,10,4,1,0]
    for n in range(1,len(expected)):
        for p,h in product(K,repeat=2):
            words = [(i,)+w for i in K for w in bases[(p+1)%4,n-1,h^i]]
            off = [sum(w[0]<i for w in words) for i in range(5)]
            rows = []
            for g,rule in rules[p]:
                length = len(next(iter(rule)))
                if n<length:
                    continue
                for tail in bases[(p+length)%4,n-length,g^h]:
                    row = [0]*len(words)
                    for term,c in rule.items():
                        i = term[0]
                        for j,v in enumerate(evaluate((p+1)%4,term[1:]+tail),off[i]):
                            row[j] += c*v
                    rows.append([x%prime for x in row])
            pivots = []
            for j in range(len(words)):
                r = len(pivots)
                pivot = next((k for k in range(r,len(rows)) if rows[k][j]),None)
                if pivot is None:
                    continue
                rows[r],rows[pivot] = rows[pivot],rows[r]
                inv = pow(rows[r][j],-1,prime)
                rows[r] = [x*inv%prime for x in rows[r]]
                for k,row in enumerate(rows):
                    if k!=r and row[j]:
                        rows[k]=[(x-row[j]*y)%prime for x,y in zip(row,rows[r])]
                pivots.append(j)
            free = [i for i in range(len(words)) if i not in pivots]
            bases[p,n,h] = [words[i] for i in free]
            projection = []
            for i in free:
                row = [0]*len(words)
                row[i]=1
                for k,j in enumerate(pivots):
                    row[j]=-rows[k][i]%prime
                projection.append(row)
            for i in K:
                multiply[p,n,h,i]=[row[off[i]:off[i+1]] for row in projection]
        assert all(sum(len(bases[p,n,h]) for h in K)==expected[n] for p in K)
    print('Central quotient through degree 13 at all four starting positions:',expected)


if __name__=='__main__':
    relations = sample_relations()
    basis,image = paths(relations,5)
    assert all([sum(len(basis[p,n,h]) for h in K) for n in range(6)]==[1,4,10,20,35,56] for p in K)
    for char,expected in [(0,[3,2,0,0]),(1,[0]*4),(2,[0,2,0,0]),(3,[0]*4)]:
        ds = [commuting_space(basis,image,h,char).cols for h in K]
        assert ds==expected,(char,ds)
        print('Quartic character',char,'dimensions',ds,flush=True)
    artinian_quotient(relations,basis,image)
