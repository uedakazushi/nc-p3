#!/usr/bin/env python3
"""Exact QQ(u,v) quintic ranks; uniform exclusion uses the manuscript lemma."""
import argparse,json,time
from pathlib import Path
import sympy as s
import generic_boundary_field as field
from flint_rational_function import K,zero,one


def backend_crosscheck():
    u,v=s.symbols('u v')
    for i in range(1,9):
        a=(i*u+v+1)/(u*v+i)
        b=(i*v-u+1)/(u*u+v+i)
        x,y=K.from_sympy(a),K.from_sympy(b)
        for actual,expected in [(x+y,a+b),(x-y,a-b),(x*y,a*b),(x/y,a/b)]:
            assert s.cancel(actual.as_expr()-expected)==0


def combine(columns,vector):
    result={}
    for j,coefficient in vector.items():
        for i,value in columns[j].items():
            total=result.get(i,zero)+coefficient*value
            if total:result[i]=total
            else:result.pop(i,None)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    backend_crosscheck()
    field.K,field.zero,field.one=K,zero,one
    w=field.tensor();relations=field.windows(w)
    print('PASS exact cubic frame and four generic pair ranks',flush=True)
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
    assert len(field.echelon(d2))==77
    print('PASS generic quintic d2 rank 77',flush=True)
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
    assert len(field.echelon(d3))==36
    assert all(not combine(d2,col) for col in d3)
    multiplication=[{j:row[192+word] for j,row in enumerate(q4) if row.get(192+word,zero)}
                    for word in free3]
    assert len(field.echelon(multiplication))==14
    multiplication_rows=[{j:col[i] for j,col in enumerate(multiplication) if i in col} for i in range(35)]
    annihilator=field.nullspace(multiplication_rows,20)
    assert len(annihilator)==6
    square=coordinates({15:one})
    cycles=[{20*i+j:x*y for i,x in square.items() for j,y in col.items()} for col in annihilator]
    assert all(not combine(d2,col) for col in cycles)
    assert len(field.echelon(d3+cycles))==42
    print('PASS d3 rank36, square annihilator6, combined rank42',flush=True)
    certificate={'result':'PASS','arithmetic':'exact rational functions QQ(u,v), FLINT rational multivariate polynomial gcd',
        'backend_crosscheck':'32 exact operations compared with SymPy',
        'quadric_coefficients':['1','u','v'],'cubic_normalization_determinant':'-2*u*v',
        'pair_ranks':[6]*4,'C_2_5_dimension':20,'C_1_5_dimension':35,
        'quintic_d2_rank':77,'quintic_d3_rank':36,'d2_d3_zero':True,
        'generic_quintic_H2_dimension':7,
        'square_multiplication_rank':14,'square_annihilator_dimension':6,
        'd3_plus_square_rank':42,'remaining_generic_H2_quotient_dimension':1,
        'uniform_conclusion':'Every parameter u*v!=0 in this cubic-frame family fails the standard Koszul AS conditions: outside the Frobenius/quartic-exact open it already fails necessary filters; on that open H2>=7 by rank specialization.',
        'limits':['This is a two-dimensional coefficient subfamily, not every four-quadric system with a basepoint.',
                  'Finite-Frobenius and quartic-exact exceptional strata are not completely classified.',
                  'The generic mixed seventh class has not been put into a compact closed formula.'],
        'elapsed_seconds':round(time.monotonic()-start,3)}
    if args.output:args.output.write_text(json.dumps(certificate,indent=2)+'\n')
    print(json.dumps(certificate,indent=2))


if __name__=='__main__':main()
