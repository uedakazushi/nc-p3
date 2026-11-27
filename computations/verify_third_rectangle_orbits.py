"""Unmarked quotient, rotations and opposite ranks on all three support strata."""
from itertools import combinations, permutations, product

from verify_parallel_third_rectangle import K, s, smith, third_rectangle
from verify_third_rectangle_frames import check as check_frames


def support_maps(source,target):
    for p,q,r in product(permutations(K),repeat=3):
        last = [None]*4
        for i,j,k,l in source:
            image = p[i]^q[j]^r[k]
            if last[l] not in (None,image):
                break
            last[l] = image
        else:
            if len(set(last))==4 and {(p[i],q[j],r[k],last[l]) for i,j,k,l in source}==set(target):
                yield p,q,r,tuple(last)


def check():
    v,h,b = s.symbols('v h b',nonzero=True)
    general = third_rectangle(v,h/v,b)
    holes = [k for k in general if k[0]^k[1]<2 and k[1]^k[2]<2]
    thin = [k for k in general if k not in holes]
    def weight(k):
        return [int(k[p]==i) for p,i in product(K,repeat=2)]
    diag,_,right = smith(s.Matrix(list(map(weight,thin))))
    assert [n for n in diag.diagonal() if n]==[1]*11+[2]
    residual = s.Matrix(list(map(weight,holes)))*right
    assert all(a%2==0 for a in residual[:,11]) and any(residual[0,12:])
    for j in range(16):
        assert residual[j,12:]==(1 if j<8 else -1)*residual[0,12:]

    parameters = {'both':(s.Integer(2),s.Integer(2)),
                  'E':(b,b),'C':(b/(b-1),b),'open':(h,b)}
    tensors = {}
    for tag,(h0,b0) in parameters.items():
        raw = third_rectangle(v,h0/v,b0)
        tensors[tag] = {k:f for k,f in raw.items() if f}
    for tag,w in tensors.items():
        h0,b0 = parameters[tag]
        inverse_tag = {'E':'C','C':'E'}.get(tag,tag)
        for dest,target in tensors.items():
            if len(w)!=len(target):
                continue
            for shift in K:
                source = {k[shift:]+k[:shift]:f for k,f in w.items()}
                actual = []
                for perms in support_maps(source,target):
                    transported = {tuple(pi[i] for pi,i in zip(perms,k)):f for k,f in source.items()}
                    assert all(transported[k]==1 for k in thin)
                    z = transported[holes[2]]
                    normalized = [s.factor(transported.get(k,0)*(1/z if j<8 else z))
                                  for j,k in enumerate(holes)]
                    a,bb,d = normalized[0],normalized[1],normalized[7]
                    new = third_rectangle(a,d,bb)
                    assert all(s.factor(f-new[k])==0 for k,f in zip(holes,normalized))
                    actual.append(tuple(map(s.factor,(a,a*d,bb))))
                expected = []
                if shift%2==0:
                    a = v if shift==0 else h0/v
                    if tag==dest:
                        expected.append((a,h0,b0))
                    if inverse_tag==dest:
                        expected.append((-a,h0,b0/(b0-1)))
                assert sorted(actual,key=str)==sorted([tuple(map(s.factor,row)) for row in expected],key=str)

    rank_conditions = [b*(v+h/v-2)-h+1,b*(v+h/v+h+1)-h+1]
    dets = [-rank_conditions[0]/((b-1)*(h-1)),rank_conditions[1]/(h-1),-1,-1]
    for q in K:
        m = s.Matrix([[general.get((i,j,i^q,j^q),0) for j in K] for i in K])
        assert s.factor(m.det()-dets[q])==0
        assert any(s.factor(m.extract(ii,jj).det()) in (-1,1)
                   for ii,jj in product(combinations(K,3),repeat=2))
    print('PASS: all unmarked rectangle orbits, period strata and opposite-rank loci')


if __name__=='__main__':
    check_frames()
    check()
