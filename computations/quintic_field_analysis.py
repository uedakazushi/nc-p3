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


def analyze_quintic(relations,K,square_index=3,rank_strategy="columns",d2_only=False,quotient_strategy="columns"):
    field.K,field.zero,field.one=K,K.zero,K.one
    zero,one=K.zero,K.one
    start=time.monotonic()
    q3,free3=field.quotient(relations,2,3,method=quotient_strategy)
    q4,free4=field.quotient(relations,1,4,method=quotient_strategy)
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
    q2,free2=field.quotient(relations,3,2,method=quotient_strategy)
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


def analyze_d2_normal_forms(relations,K):
    """Project only the 120 required columns; do not build all dual C4 vectors."""
    field.K,field.zero,field.one=K,K.zero,K.one
    zero=K.zero;start=time.monotonic()
    equations3,free3=field.full_pivot_equations(field.inserted(relations,2,3),64)
    assert len(equations3)==44 and len(free3)==20
    print('PASS C3 monomial section dimension20',flush=True)
    equations4,free4=field.full_pivot_equations(field.inserted(relations,1,4),256,
        lambda r:print('degree4 relation pivot progress',r,flush=True))
    assert len(equations4)==221 and len(free4)==35
    print('PASS C4 monomial quotient dimension35',flush=True)
    free_index={word:i for i,word in enumerate(free4)}
    def normal_form(vector):
        result=dict(vector)
        for pivot,row in equations4:
            coefficient=result.pop(pivot,zero)
            if coefficient:
                for column,value in row.items():
                    if column==pivot:continue
                    total=result.get(column,zero)-coefficient*value
                    if total:result[column]=total
                    else:result.pop(column,None)
        assert all(word in free_index for word in result)
        return result
    normal_cache={}
    def monomial_normal_form(word):
        if word not in normal_cache:
            normal_cache[word]=normal_form({word:K.one})
            if len(normal_cache)%20==0:print('normal monomial images',len(normal_cache),flush=True)
        return normal_cache[word]
    columns=[]
    for relation in relations[0]:
        for tail in free3:
            column={}
            for (a,b),coefficient in relation.items():
                for word,value in monomial_normal_form(64*b+tail).items():
                    row=35*a+free_index[word]
                    total=column.get(row,zero)+coefficient*value
                    if total:column[row]=total
                    else:column.pop(row,None)
            columns.append(column)
    print('PASS direct d2 quotient projection',flush=True)
    rank=field.sparse_full_pivot_rank(columns,lambda r:print('d2 pivot progress',r,flush=True))
    assert rank<=84
    print('Exact quintic d2 rank:',rank,flush=True)
    return {'d2_rank':rank,'domain_dimension':120,'codomain_dimension':140,
            'C3_dimension':20,'C4_dimension':35,'projection':'direct relation normal forms',
            'elapsed_seconds':round(time.monotonic()-start,3)}


def _reference_pivots(vectors,K):
    """Choose a pivot trace over an exact rational reference field."""
    remaining=[(j,{i:x for i,x in vector.items() if x}) for j,vector in enumerate(vectors)]
    remaining=[item for item in remaining if item[1]]
    trace=[];normalized_rows=[]
    while remaining:
        frequency={}
        for _,row in remaining:
            for i in row:frequency[i]=frequency.get(i,0)+1
        best=None;selected=None
        for j,(_,row) in enumerate(remaining):
            for i in row:
                score=((len(row)-1)*(frequency[i]-1),len(row))
                if best is None or score<best:best,selected=score,(j,i)
        j,i=selected;tag,row=remaining.pop(j);pivot=row[i]
        row={column:x/pivot for column,x in row.items()}
        updated=[]
        for other_tag,other in remaining:
            factor=other.pop(i,K.zero)
            if factor:
                for column,value in row.items():
                    if column==i:continue
                    total=other.get(column,K.zero)-factor*value
                    if total:other[column]=total
                    else:other.pop(column,None)
            if other:updated.append((other_tag,other))
        remaining=updated;trace.append((tag,i));normalized_rows.append((i,row))
    return trace,normalized_rows


def _reduce(vector,equations,zero):
    result=dict(vector)
    for pivot,row in equations:
        factor=result.pop(pivot,zero)
        if factor:
            for i,value in row.items():
                if i==pivot:continue
                total=result.get(i,zero)-factor*value
                if total:result[i]=total
                else:result.pop(i,None)
    return result


def _guided_equations(vectors,trace,K):
    result=[]
    for j,pivot in trace:
        row=_reduce(vectors[j],result,K.zero)
        coefficient=row[pivot]
        assert coefficient
        row={i:x/coefficient for i,x in row.items()}
        result.append((pivot,row))
        if len(result)%64==0:print('guided pivot progress',len(result),flush=True)
    return result


def _projected_d2(relations,free3,equations4,free4,K):
    index={word:i for i,word in enumerate(free4)};cache={}
    def normal(word):
        if word not in cache:
            cache[word]=_reduce({word:K.one},equations4,K.zero)
            assert all(i in index for i in cache[word])
        return cache[word]
    result=[]
    for relation in relations[0]:
        for tail in free3:
            column={}
            for (a,b),coefficient in relation.items():
                for word,value in normal(64*b+tail).items():
                    row=35*a+index[word]
                    total=column.get(row,K.zero)+coefficient*value
                    if total:column[row]=total
                    else:column.pop(row,None)
            result.append(column)
    return result


def analyze_d2_specialized_basis(relations,reference_relations,K,reference_K):
    """Exact generic rank, with bases chosen by nonzero rational reference minors.

    Generic C3=20 and C4=35 must already be proved by the nonempty
    Frobenius/quartic-exact open. Only basis choices use the reference point;
    every residual generic d2 column is reduced over the full function field.
    """
    start=time.monotonic()
    reference3=list(field.inserted(reference_relations,2,3))
    trace3,_=_reference_pivots(reference3,reference_K)
    assert len(trace3)==44
    free3=sorted(set(range(64))-{i for j,i in trace3})
    reference4=list(field.inserted(reference_relations,1,4))
    trace4,old_equations4=_reference_pivots(reference4,reference_K)
    assert len(trace4)==221
    free4=sorted(set(range(256))-{i for j,i in trace4})
    old_d2=_projected_d2(reference_relations,free3,old_equations4,free4,reference_K)
    trace_d2,_=_reference_pivots(old_d2,reference_K)
    assert len(trace_d2)==77
    print('PASS exact reference minors of sizes44,221,77',flush=True)
    generic4=list(field.inserted(relations,1,4))
    equations4=_guided_equations(generic4,trace4,K)
    print('PASS generic degree4 relation basis from the nonzero reference minor',flush=True)
    d2=_projected_d2(relations,free3,equations4,free4,K)
    print('PASS generic d2 projection',flush=True)
    pivots=_guided_equations(d2,trace_d2,K)
    selected={j for j,i in trace_d2}
    residual=[_reduce(column,pivots,K.zero) for j,column in enumerate(d2) if j not in selected]
    field.K,field.zero,field.one=K,K.zero,K.one
    extra=field.sparse_full_pivot_rank([column for column in residual if column])
    rank=77+extra
    assert rank<=84
    print('Exact generic d2 rank:',rank,'nonzero residual columns:',sum(bool(x) for x in residual),flush=True)
    return {'d2_rank':rank,'domain_dimension':120,'codomain_dimension':140,
            'C3_dimension':20,'C4_dimension':35,'reference_minor_sizes':[44,221,77],
            'generic_residual_zero':not any(residual),'extra_rank':extra,
            'projection':'exact specialization-guided relation normal forms',
            'generic_dimension_reason':'proved nonempty Frobenius/quartic-exact open, not inferred from numerical samples',
            'elapsed_seconds':round(time.monotonic()-start,3)}
