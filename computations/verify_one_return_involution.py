#!/usr/bin/env python3
"""Check signed reversal on the point-quotient Tor complex over exact QQ."""
import argparse,json
from pathlib import Path
import sympy as s
from flint_rational_function import rational_function_field
from product_net_boundary import NAMES,tensor,terms,relation_windows
from trivial_extension_bitor import one_return_complex,point_quotient
from quintic_field_analysis import combine
import generic_boundary_field as field


def signed_reversal(bases,quadrics,K):
    from itertools import combinations
    pairs=list(combinations(range(4),2))
    tau=[{i:-K.one} for i in range(6)]
    for j,q in enumerate(quadrics):
        col={6+j:K.one}
        for i,pair in enumerate(pairs):
            value=q.get(pair,K.zero)
            if value:col[i]=-value
        tau.append(col)
    result={}
    for h,base in bases.items():
        lookup={v:i for i,v in enumerate(base)};columns=[]
        for p,q,a,j,b in base:
            left=tau[b] if q==2 else {b:K.one}
            right=tau[a] if p==2 else {a:K.one}
            sign=-K.one if p*q%2 else K.one
            columns.append({lookup[q,p,c,j,d]:sign*x*y for c,x in left.items() for d,y in right.items()})
        result[h]=columns
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path)
    parser.add_argument('--fixture',choices=('boolean','generic'),default='boolean');args=parser.parse_args()
    parameters=(0,0,0,1,1,0,0,0) if args.fixture=='boolean' else (2,5,11,3,7,17,23,31)
    K=rational_function_field(('z',));field.K,field.zero,field.one=K,K.zero,K.one
    class Point:
        zero,one=K.zero,K.one
        @staticmethod
        def from_sympy(e):return K.from_sympy(s.sympify(e).subs(dict(zip(s.symbols(NAMES),parameters))))
    point=Point();w=tensor(point);quadrics=terms(point)
    d,bases=one_return_complex(w,relation_windows(w,point)[0],quadrics,K)
    sigma=signed_reversal(bases,quadrics,K)
    # Same coordinate removal for d and sigma.
    remove={1:set(range(8)),2:{i for i,v in enumerate(bases[2]) if v[3]==3},
            3:{i for i,v in enumerate(bases[3]) if v[3]==5},4:{i for i,v in enumerate(bases[4]) if v[3]==3}}
    keep={h:[i for i in range(len(base)) if i not in remove[h]] for h,base in bases.items()}
    indexes={h:{old:new for new,old in enumerate(keep[h])} for h in bases}
    quotient_sigma={h:[{indexes[h][i]:x for i,x in sigma[h][j].items() if i in indexes[h]} for j in keep[h]] for h in bases}
    for h in bases:
        for j,column in enumerate(sigma[h]):assert combine(sigma[h],column,K.zero)=={j:K.one}
    for h in (2,3,4):
        for j in range(len(d[h])):
            assert combine(d[h],sigma[h][j],K.zero)==combine(sigma[h-1],d[h][j],K.zero)
    qd=point_quotient(d,bases);blocks={}
    for sign in (1,-1):
        eigen={}
        for h in (2,3,4):
            columns=[]
            for j,col in enumerate(quotient_sigma[h]):
                column={i:sign*x for i,x in col.items()}
                total=column.get(j,K.zero)+K.one
                if total:column[j]=total
                else:column.pop(j,None)
                columns.append(column)
            eigen[h]=list(field.echelon(columns).values())
        ranks={h:field.sparse_full_pivot_rank([combine(qd[h],v,K.zero) for v in eigen[h]]) for h in (3,4)}
        dims=[len(eigen[h]) for h in (2,3,4)]
        homology=[dims[0]-ranks[3],dims[1]-ranks[3]-ranks[4],dims[2]-ranks[4]]
        blocks[str(sign)]={'dimensions_C2_C3_C4':dims,'ranks_d3_d4':[ranks[3],ranks[4]],'homology_H2_H3_H4':homology}
    assert blocks['-1']['dimensions_C2_C3_C4']==[60,200,135]
    assert blocks['1']['dimensions_C2_C3_C4']==[48,200,165]
    if args.fixture=='boolean':assert sum(blocks[str(sign)]['homology_H2_H3_H4'][1] for sign in (1,-1))==7
    assert blocks['-1']['homology_H2_H3_H4'][1]>=5
    report={'result':'PASS','fixture':args.fixture,'parameters':list(parameters),'arithmetic':'exact QQ','involution_squared_identity':True,'differentials_commute_with_involution':True,
            'blocks':blocks,'universal_negative_euler_characteristic':-5,
            'universal_H3_lower_bound':5,
            'limits':['The universal bound follows from the proved signed-reversal dimensions, not this pointwise rank calculation.','The exact eight-variable H3 dimension is not determined.']}
    if args.output:args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
