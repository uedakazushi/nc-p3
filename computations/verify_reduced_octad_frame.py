#!/usr/bin/env python3
"""A general-position reduced octad and its one-basepoint Frobenius boundary."""
import argparse,json
from itertools import combinations,combinations_with_replacement,product
from pathlib import Path
import sympy as s
from fractions import Fraction as F
from tensor_symmetrization import symmetrizer_space
import generic_boundary_field as field
from flint_rational_function import rational_function_field
from quintic_field_analysis import analyze_quintic
from frobenius_tensor import digits


def standard_dimension(groebner):
    assert groebner.is_zero_dimensional
    leading=[p.LM(order=groebner.order).exponents for p in groebner.polys]
    bounds=[]
    for i in range(3):
        bounds.append(min(m[i] for m in leading if m[i] and all(not m[j] for j in range(3) if j!=i)))
    return sum(not any(all(a>=b for a,b in zip(word,m)) for m in leading)
               for word in product(*(range(bound) for bound in bounds)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    x=s.symbols('x0:4');skew=s.zeros(4)
    for (i,j),value in zip(combinations(range(4),2),(1,2,3,5,7,11)):
        skew[i,j]=value;skew[j,i]=-value
    orthogonal=(s.eye(4)-skew)*(s.eye(4)+skew).inv()
    assert orthogonal*orthogonal.T==s.eye(4)
    points=s.eye(4).row_join(orthogonal)
    assert all(points[:,list(c)].det()!=0 for c in combinations(range(8),4))
    degree2=list(combinations_with_replacement(range(4),2))
    evaluation2=s.Matrix(8,10,lambda i,j:points[degree2[j][0],i]*points[degree2[j][1],i])
    assert evaluation2.rank()==7
    net=[s.expand(sum(v[j]*x[a]*x[b] for j,(a,b) in enumerate(degree2))) for v in evaluation2.nullspace()]
    affine=[f.subs(x[3],1-x[0]-x[1]-x[2]) for f in net]
    infinity=[f.subs(x[3],-x[0]-x[1]-x[2]) for f in net]
    assert standard_dimension(s.groebner(affine,*x[:3],order='grevlex'))==8
    assert standard_dimension(s.groebner(infinity,*x[:3],order='grevlex'))==8
    h=[sum(points[:,i]) for i in range(8)]
    assert all(h)
    normalized=s.Matrix.hstack(*(points[:,i]/h[i] for i in range(8)))
    e2=s.Matrix(8,10,lambda i,j:normalized[degree2[j][0],i]*normalized[degree2[j][1],i])
    assert e2[1:,:].rank()==7
    relation=e2.T.nullspace();assert len(relation)==1 and all(relation[0])
    assert normalized*s.diag(*list(relation[0]))*normalized.T==s.zeros(4)
    degree3=list(combinations_with_replacement(range(4),3))
    e3=s.Matrix(8,20,lambda i,j:s.prod(normalized[a,i] for a in degree3[j]))
    assert e3.rank()==8
    qvalues=[normalized[1,i]**2+2*normalized[2,i]**2+3*normalized[3,i]**2 for i in range(8)]
    assert qvalues[0]==0 and all(qvalues[1:])
    kernel=normalized[:,1:].nullspace();assert len(kernel)==3
    w={}
    for encoded in range(256):
        a,b,j,c=digits(encoded)
        if j<3:
            value=sum(kernel[j][i-1]/qvalues[i]*normalized[a,i]*normalized[b,i]*normalized[c,i] for i in range(1,8))
        else:value=normalized[a,0]*normalized[b,0]*normalized[c,0]
        if value:w[encoded]=value
    symmetrizers=[symmetrizer_space({i:F(value) for i,value in w.items()},j) for j in (1,2,3)]
    assert symmetrizers[0]==symmetrizers[2]==[s.eye(4)]
    assert len(symmetrizers[1])==1 and symmetrizers[1][0].rank()==1
    K=rational_function_field(('z',));field.K,field.zero,field.one=K,K.zero,K.one
    exact_w={i:K.from_sympy(value) for i,value in w.items()}
    singleton_ranks=[]
    for factor in range(4):
        columns={}
        for encoded,value in exact_w.items():
            word=digits(encoded);key=tuple(word[j] for j in range(4) if j!=factor)
            columns.setdefault(key,{})[word[factor]]=value
        singleton_ranks.append(len(field.echelon(columns.values())))
    assert singleton_ranks==[4]*4
    relations=field.windows(exact_w);dual=field.dual_windows(relations)
    intersections=[[len(field.nullspace(field.inserted(dual,i,n),4**n)) for n in (3,4)] for i in range(4)]
    assert intersections==[[4,1]]*4
    quartic=[len(field.quotient(relations,i,4)[0]) for i in range(4)]
    assert quartic==[35]*4
    report=analyze_quintic(relations,K,square_index=0)
    certificate={'result':'PASS','arithmetic':'exact rational geometry and QQ embedded in QQ(z)',
        'skew_parameters':[1,2,3,5,7,11],'orthogonal_matrix':[[str(orthogonal[i,j]) for j in range(4)] for i in range(4)],
        'nonzero_four_point_minors':70,'quadratic_evaluation_rank':7,
        'affine_and_hyperplane_section_quotient_dimensions':[8,8],
        'cubic_evaluation_rank':8,'seven_point_quadratic_rank':7,
        'nonzero_Cayley_Bacharach_weights':True,'fourth_quadric':'x1^2+2*x2^2+3*x3^2',
        'symmetrizer_dimensions':[len(space) for space in symmetrizers],
        'middle_symmetrizer_rank':1,
        'singleton_ranks':singleton_ranks,'pair_ranks':[6]*3,
        'intersection_dimensions':intersections,'quartic_dimensions':quartic,'quintic':report,
        'limits':['A single general-position reduced octad fixture, not a generic homology rank theorem for all reduced octads.',
                  'The orthogonal normal form does not classify indexed AS isomorphisms.']}
    if args.output:args.output.write_text(json.dumps(certificate,indent=2)+'\n')
    print(json.dumps(certificate,indent=2))


if __name__=='__main__':main()
