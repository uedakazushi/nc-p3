"""The torus, tetrahedral action, rotation, and maximal supports of the quantum attraction fiber."""

from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from math import gcd

from flint import nmod_mat
from sympy import Matrix, ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from verify_quantum_ore_incidence import K, WORDS, tangent_equations

WEIGHTS = sorted(set(permutations((2,0,-1,-1))) | set(permutations((1,1,-1,-1))))
ROTATION = (3,2,1,1,3,2,2,3,1,0,0,0,0,0,0,0,0,0)


def verify_representation():
    rows,volume = tangent_equations()
    blocks = defaultdict(list)
    for i,word in enumerate(WORDS):
        blocks[tuple(word.count(j)-1 for j in K)].append(i)
    for weight,idx in blocks.items():
        mat = nmod_mat([[r[i] for i in idx] for r in rows],97)
        line,dim = mat.nullspace()
        if weight not in WEIGHTS:
            expected = {(-1,-1,-1,3):0,(-1,0,0,1):4,(0,0,0,0):13}
            assert dim == expected[tuple(sorted(weight))]
            continue
        assert dim==1
        co = {WORDS[i]:int(line[j,0]) for j,i in enumerate(idx)}
        scale = pow(22,ROTATION[WEIGHTS.index(weight)],97)
        assert all((-v-scale*co[w[1:]+w[:1]])%97==0 for w,v in co.items())
        if weight==(1,1,-1,-1):
            perm = (1,0,3,2)
            def move(w):
                return tuple(perm[i] for i in w)
            def phase(w):
                return pow(22,sum(j*(0,0,3,1)[i] for j,i in enumerate(w))%4,97)
            assert all((v*phase(w)+volume[move(w)])%97==0 for w,v in volume.items())
            assert all((v*phase(w)+co[move(w)])%97==0 for w,v in co.items())


def dot(a,b):
    return sum(x*y for x,y in zip(a,b,strict=True))


def cross(a,b):
    return tuple(a[(i+1)%3]*b[(i+2)%3]-a[(i+2)%3]*b[(i+1)%3] for i in range(3))


def verify_chambers():
    # Work modulo the diagonal cocharacter by setting its last coordinate to zero.
    normals = [w[:3] for w in WEIGHTS]
    chambers = set()
    for a,b in permutations(normals,2):
        ray = cross(a,b)
        if not any(ray):
            continue
        edge = cross(a,ray)
        for i,j,k in product((-1,1),repeat=3):
            def positive(w):
                return next(t for t in (i*dot(w,ray),j*dot(w,edge),k*dot(w,a)) if t)>0
            chambers.add(frozenset(n for n,w in enumerate(normals) if positive(w)))
    assert len(chambers)==168
    maximal = {s for s in chambers if not any(s<t for t in chambers)}
    assert Counter(map(len,maximal))=={9:36,11:12}
    rotations = [p for p in permutations(K) if sum(p[i]>p[j] for i,j in combinations(K,2))%2==0]
    def orbit(s):
        return {tuple(sorted(WEIGHTS.index(tuple(WEIGHTS[i][j] for j in p)) for i in s)) for p in rotations}
    reps={min(orbit(s)) for s in maximal}
    expected = {tuple(range(11)),(0,1,2,3,4,5,7,9,12),(0,1,2,3,4,5,9,10,12),(0,1,2,3,4,6,8,10,13)}
    assert reps==expected
    for s in reps:
        assert len(orbit(s))==12
        m=Matrix([[WEIGHTS[i][j] for i in s] for j in range(3)])
        diagonal,_,right=smith_normal_decomp(DomainMatrix.from_Matrix(m).convert_to(ZZ))
        assert list(diagonal.to_Matrix().diagonal())==[1,1,1]
        kernel=right.to_Matrix()[:,3:]
        assert m*kernel==Matrix.zeros(3,len(s)-3)
        phases=Matrix([[ROTATION[i] for i in s]])*kernel
        assert gcd(4,*map(int,phases))==1


if __name__=="__main__":
    verify_representation()
    verify_chambers()
    print("PASS: 18 weight lines, 168 chambers, four maximal support orbits and their exact periods")
