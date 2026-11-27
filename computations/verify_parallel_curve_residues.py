"""Certify irreducible characteristic-zero parameter factors at split residue points."""
from functools import cache, reduce
from operator import xor

from flint import fmpz_poly, nmod_mat
import numpy as np

from verify_parallel_third_rectangle import K, matrix, product, s, third_rectangle


def rules(value,prime):
    w = third_rectangle(s.Integer(value),s.Rational(2,value),2)
    def residue(v):
        a,b = s.fraction(v)
        return int(a)*pow(int(b),-1,prime)%prime
    return [[(h,{(j,j^h):residue(c[j]) for j in K if c[j]}) for h in K
             for c in matrix(w,p,h).columnspace()] for p in K]


def paths(relations,prime,limit):
    basis = {(p,0,h):[()] if h==0 else [] for p,h in product(K,repeat=2)}
    prepend = {}
    @cache
    def image(p,word):
        if not word:
            return np.ones(1,dtype=np.int64)
        return image((p+1)%4,word[1:])@prepend[p,len(word),reduce(xor,word),word[0]]%prime
    for n in range(1,limit+1):
        for p,h in product(K,repeat=2):
            words = [(i,)+w for i in K for w in basis[(p+1)%4,n-1,h^i]]
            offsets = np.cumsum([0]+[len(basis[(p+1)%4,n-1,h^i]) for i in K])
            size = len(words)
            rows = []
            for g,rule in relations[p]:
                length = len(next(iter(rule)))
                if n<length:
                    continue
                for tail in basis[(p+length)%4,n-length,h^g]:
                    row = np.zeros(size,dtype=np.int64)
                    for term,c in rule.items():
                        i = term[0]
                        row[offsets[i]:offsets[i+1]]+=c*image((p+1)%4,term[1:]+tail)
                    rows.append((row%prime).tolist())
            M = nmod_mat(rows,prime) if rows else nmod_mat(0,size,prime)
            reduced,rank = M.rref()
            projection = np.zeros((size,size-rank),dtype=np.int64)
            if rank:
                z = np.array(reduced.tolist()[:rank],dtype=np.int64)
                pivots = z.astype(bool).argmax(axis=1)
                free = np.delete(np.arange(size),pivots)
                projection[pivots,:]=(-z[:,free])%prime
            else:
                free = list(range(size))
            projection[free,:]=np.eye(size-rank,dtype=np.int64)
            for i in K:
                prepend[p,n,h,i]=projection[offsets[i]:offsets[i+1],:]
            basis[p,n,h]=[words[i] for i in free]
    return basis,image


def residue_test(value,prime,make_rules=rules):
    rel = make_rules(value,prime)
    basis,image = paths(rel,prime,5)
    assert all(sum(len(basis[p,n,h]) for h in K)==d for p in K for n,d in enumerate((1,4,10,20,35,56)))
    cut = [list(r) for r in rel]
    for h in K:
        offsets = np.cumsum([0]+[len(basis[p,4,h]) for p in K])
        rows = []
        for p,i in product(K,repeat=2):
            j = (p+1)%4
            block = np.zeros((len(basis[p,5,h^i]),offsets[-1]),dtype=np.int64)
            for k,w in enumerate(basis[p,4,h]):
                block[:,offsets[p]+k]=image(p,w+(i,))
            for k,w in enumerate(basis[j,4,h]):
                block[:,offsets[j]+k]=-image(p,(i,)+w)
            rows.extend((block%prime).tolist())
        kernel,dim = nmod_mat(rows,prime).nullspace()
        assert dim==(3,2,0,0)[h]
        if not dim:
            continue
        vectors = np.array(kernel.tolist(),dtype=np.int64)[:,:dim]
        weights = ((1,3,5),(2,7,11)) if h==0 else ((1,0),(0,1))
        for coeff in weights:
            vector = vectors@np.array(coeff,dtype=np.int64)%prime
            for p in K:
                section = {w:int(vector[offsets[p]+j]) for j,w in enumerate(basis[p,4,h]) if vector[offsets[p]+j]}
                if section:
                    cut[p].append((h,section))
    quotient,_ = paths(cut,prime,13)
    profile = (1,4,10,20,31,40,44,40,31,20,10,4,1,0)
    assert all(sum(len(quotient[p,n,h]) for h in K)==d for p in K for n,d in enumerate(profile))


def certify(exceptional,make_rules=rules):
    for coefficients in exceptional:
        f = fmpz_poly(list(coefficients))
        for prime in (1009,1013,1019,1021,1031,1033,1039,1049,1051,1061):
            root = next((j for j in range(1,prime) if int(f(j))%prime==0),None)
            if root is not None and int(f.leading_coefficient())%prime:
                residue_test(root,prime,make_rules)
                break
        else:
            raise AssertionError(('missing split residue',coefficients))
