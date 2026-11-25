"""Weight and orbit checks for all distinct and single repeated quantum faces."""

from itertools import combinations, permutations, product

K = tuple(range(4))
MOD = 12


def triangles(q,i,j,k):
    return (q[i][j]+q[j][k]+q[k][i])%MOD


def verify():
    checked = [0,0]
    for a,b,c in product(range(MOD),repeat=3):
        q = [[0]*4 for _ in K]
        for (i,j),v in zip(combinations(K[1:],2),(a,b,c),strict=True):
            q[i][j],q[j][i] = v,-v
        face = tuple(((-1)**m*triangles(q,*(i for i in K if i!=m)))%MOD for m in K)
        assert sum(face)%MOD==0
        if len(set(face)) not in (3,4):
            continue
        repeated = [x for x in set(face) if face.count(x)==2]
        if repeated == [0]:
            continue
        checked[bool(repeated)] += 1
        for k,l in combinations(K,2):
            i,j = (v for v in K if v not in (k,l))
            for s in range(3):
                alpha = [-1]*4
                alpha[i],alpha[j]=s,2-s
                chi = [sum(q[t][u]*alpha[t] for t in K)%MOD for u in K]
                assert (chi[i]==chi[j]) == (face[k]==face[l])
        if repeated:
            pair = {i for i in K if face.count(face[i])==2}
            perms = []
            for p in permutations(K):
                sign = (-1)**sum(p[i]>p[j] for i,j in combinations(K,2))
                if all((face[p[i]]-sign*face[i])%MOD==0 for i in K):
                    perms.append(p)
            assert len(perms)==(2 if 2*repeated[0]%MOD==0 else 1)
            assert all(p[i]==i for p in perms for i in pair)
    assert all(checked)


if __name__=="__main__":
    verify()
    print("PASS: face-pair resonance, oriented permutation stabilizers, both multiplicity strata")
