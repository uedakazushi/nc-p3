"""One-return two-sided Koszul complex for the four-layer quadric truncation.

The caller supplies a cubic frame, a six-vector B2 dual frame, and four
commutative quadrics. The complex computes the internal-degree-five,
vertex-(position1,position2) summand of Tor^{A^e}(S^e, DA[degree4]).
"""
from itertools import combinations


def relation_frame(quadrics,K):
    result=[]
    for a,b in combinations(range(4),2):
        result.append({(a,b):K.one,(b,a):-K.one})
    result.extend(quadrics)
    assert len(result)==10
    return result


def one_return_complex(w,b2_dual,quadrics,K,relation_basis=None,coordinate_pairs=None,start_position=1,end_position=2):
    """Return sparse differentials d_h:C_h -> C_(h-1), h=2,3,4."""
    zero,one=K.zero,K.one
    q2=relation_frame(quadrics,K) if relation_basis is None else relation_basis
    assert len(q2)==10
    qdims=(1,4,10);bdims=(1,4,6,4)
    free2=coordinate_pairs or [(0,1),(0,2),(1,2),(1,3),(2,3),(3,3)]
    pmax,qmax=3-start_position,end_position
    assert pmax<=2 and qmax<=2
    hmin=max(0,end_position-start_position);hmax=pmax+qmax
    def bdegree(h):return h+start_position-end_position
    # C_i is dual to multiplication by e_i in the commutative ring B.
    contractions={}
    for r in (1,2,3):
        for i in range(4):
            for j in range(bdims[r]):
                if r==1:column={0:one} if i==j else {}
                elif r==2:column={a:b2_dual[j].get((i,a),zero) for a in range(4)}
                else:column={a:w.get(64*i+16*b+4*j+c,zero) for a,(b,c) in enumerate(free2)}
                contractions[r,i,j]={a:x for a,x in column.items() if x}
    bases={}
    for h in range(hmin,hmax+1):
        bases[h]=[(p,h-p,a,j,b) for p in range(pmax+1) if 0<=h-p<=qmax
                  for a in range(qdims[p]) for j in range(bdims[bdegree(h)])
                  for b in range(qdims[h-p])]
    indexes={h:{v:i for i,v in enumerate(base)} for h,base in bases.items()}
    differentials={}
    for h in range(hmin+1,hmax+1):
        columns=[];r=bdegree(h)
        for p,q,a,j,b in bases[h]:
            column={}
            def add(key,value):
                row=indexes[h-1][key];total=column.get(row,zero)+value
                if total:column[row]=total
                else:column.pop(row,None)
            if p:
                pieces={(0,a):one} if p==1 else q2[a]
                for (prefix,arrow),coefficient in pieces.items():
                    for m,value in contractions[r,arrow,j].items():
                        add((p-1,q,prefix,m,b),coefficient*value)
            if q:
                pieces={(b,0):one} if q==1 else q2[b]
                sign=one if p%2==0 else -one
                for (arrow,suffix),coefficient in pieces.items():
                    for m,value in contractions[r,arrow,j].items():
                        add((p,q-1,a,m,suffix),sign*coefficient*value)
            columns.append(column)
        differentials[h]=columns
    if (start_position,end_position)==(1,2):assert [len(bases[h]) for h in range(1,5)]==[8,144,480,400]
    if (start_position,end_position)==(1,1):assert [len(bases[h]) for h in range(4)]==[1,32,156,160]
    return differentials,bases


def point_quotient(differentials,bases):
    """Quotient by the point-evaluation subcomplex at p=(0,0,0,1).

    Requires B1 dual coordinate3=p, B2 dual coordinate5=p^2 and B3
    dual coordinate3=p^3. The known point complex has d3 rank29 and
    d4 rank51; H3(point)=0 gives the full-rank reconstruction.
    """
    remove={1:set(range(len(bases[1]))),
            2:{i for i,(_,_,a,j,b) in enumerate(bases[2]) if j==3},
            3:{i for i,(_,_,a,j,b) in enumerate(bases[3]) if j==5},
            4:{i for i,(_,_,a,j,b) in enumerate(bases[4]) if j==3}}
    keep={h:[i for i in range(len(bases[h])) if i not in remove[h]] for h in range(1,5)}
    indexes={h:{old:new for new,old in enumerate(keep[h])} for h in range(1,5)}
    result={}
    for h in (2,3,4):
        # Verify preservation of the point subcomplex before quotienting.
        assert all(not any(i not in remove[h-1] for i in differentials[h][j]) for j in remove[h])
        result[h]=[{indexes[h-1][i]:x for i,x in differentials[h][j].items() if i in indexes[h-1]}
                   for j in keep[h]]
    assert [len(keep[h]) for h in range(1,5)]==[0,108,400,300]
    return result
