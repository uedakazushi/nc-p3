"""Shared exact quintic Koszul matrices for generic four-window relations."""
import time
import generic_boundary_field as field


def combine(columns,vector,zero):
    result={}
    for j,coefficient in vector.items():
        for i,value in columns[j].items():
            total=result.get(i,zero)+coefficient*value
            if total:result[i]=total
            else:result.pop(i,None)
    return result


def analyze_quintic(relations,K,square_index=3,rank_strategy="columns",d2_only=False):
    field.K,field.zero,field.one=K,K.zero,K.one
    zero,one=K.zero,K.one
    start=time.monotonic()
    q3,free3=field.quotient(relations,2,3)
    q4,free4=field.quotient(relations,1,4)
    assert len(q3)==20 and len(q4)==35
    print('PASS quotient dimensions 20,35',flush=True)
    d2=[]
    for relation in relations[0]:
        for word in free3:
            column={}
            for (a,b),coefficient in relation.items():
                for j,row in enumerate(q4):
                    value=coefficient*row.get(64*b+word,zero)
                    if value:column[35*a+j]=column.get(35*a+j,zero)+value
            d2.append({i:x for i,x in column.items() if x})
    rank2=(field.sparse_full_pivot_rank(d2,lambda r:print("d2 pivot progress",r,flush=True))
           if rank_strategy=="full" else len(field.echelon(d2)))
    print('Exact quintic d2 rank:',rank2,flush=True)
    if d2_only:
        return {'d2_rank':rank2,'domain_dimension':120,'codomain_dimension':140,
                'C3_dimension':20,'C4_dimension':35,
                'elapsed_seconds':round(time.monotonic()-start,3)}
    dual=field.dual_windows(relations)
    w3=field.nullspace(field.inserted(dual,0,3),64)
    assert len(w3)==4
    ordered=sorted((min(4*a+b for a,b in relation),i) for i,relation in enumerate(relations[0]))
    def coordinates(vector):
        remainder=dict(vector);result={}
        for pivot,j in ordered:
            coefficient=remainder.get(pivot,zero)
            if not coefficient:continue
            result[j]=coefficient
            for (a,b),value in relations[0][j].items():
                key=4*a+b;total=remainder.get(key,zero)-coefficient*value
                if total:remainder[key]=total
                else:remainder.pop(key,None)
        assert not remainder
        return result
    syzygies=[]
    for z in w3:
        result={}
        for b in range(4):
            prefix={i//4:x for i,x in z.items() if i%4==b}
            for j,x in coordinates(prefix).items():result[4*j+b]=x
        syzygies.append(result)
    q2,free2=field.quotient(relations,3,2)
    assert len(q2)==10
    d3=[]
    for z in syzygies:
        for word in free2:
            column={}
            for index,coefficient in z.items():
                relation_index,b=divmod(index,4)
                for j,row in enumerate(q3):
                    value=coefficient*row.get(16*b+word,zero)
                    if value:
                        key=20*relation_index+j
                        column[key]=column.get(key,zero)+value
            d3.append({i:x for i,x in column.items() if x})
    rank3=len(field.echelon(d3))
    assert all(not combine(d2,col,zero) for col in d3)
    multiplication=[{j:row[64*square_index+word] for j,row in enumerate(q4) if row.get(64*square_index+word,zero)}
                    for word in free3]
    multiplication_rank=len(field.echelon(multiplication))
    multiplication_rows=[{j:col[i] for j,col in enumerate(multiplication) if i in col} for i in range(35)]
    annihilator=field.nullspace(multiplication_rows,20)
    assert len(annihilator)==20-multiplication_rank
    square=coordinates({5*square_index:one})
    cycles=[{20*i+j:x*y for i,x in square.items() for j,y in col.items()} for col in annihilator]
    assert all(not combine(d2,col,zero) for col in cycles)
    combined_rank=len(field.echelon(d3+cycles))
    print('Exact d3, square kernel, combined ranks:',rank3,len(annihilator),combined_rank,flush=True)
    return {'d2_rank':rank2,'d3_rank':rank3,'H2_dimension':120-rank2-rank3,
            'multiplication_rank':multiplication_rank,'square_kernel_dimension':len(annihilator),
            'combined_rank':combined_rank,'square_H2_dimension':combined_rank-rank3,
            'remaining_H2_quotient_dimension':120-rank2-combined_rank,
            'elapsed_seconds':round(time.monotonic()-start,3)}
