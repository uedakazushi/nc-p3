#!/usr/bin/env python3
"""A nonzero Grassmannian Jacobian minor proves dominance, not AS regularity."""
import argparse,json
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from product_net_boundary import terms,NAMES


def dominance_certificate():
    x=s.symbols('x0:4');parameters=s.symbols(NAMES)
    base=dict(zip(parameters,(2,3,5,2,3,5,7,11)))
    forms=[sum(value*x[i]*x[j] for (i,j),value in q.items()) for q in terms()]
    monomials=list(combinations_with_replacement(range(4),2));monomials.remove((3,3))
    def coefficients(f):return s.Matrix([s.Poly(s.expand(f),*x).coeff_monomial(x[i]*x[j]) for i,j in monomials])
    space=s.Matrix.hstack(*(coefficients(f.subs(base)) for f in forms));assert space.rank()==4
    quotient=s.Matrix.hstack(*space.T.nullspace()).T
    columns=[]
    for i in range(4):
        for j in range(4):
            if i<3 and j==3:continue
            projected=quotient*s.Matrix.hstack(*(coefficients((x[j]*s.diff(f,x[i])).subs(base)) for f in forms))
            columns.append(s.Matrix([projected[k,l] for l in range(4) for k in range(5)]))
    for parameter in parameters:
        projected=quotient*s.Matrix.hstack(*(coefficients(s.diff(f,parameter).subs(base)) for f in forms))
        columns.append(s.Matrix([projected[k,l] for l in range(4) for k in range(5)]))
    jacobian=s.Matrix.hstack(*columns)
    pivots=jacobian.rref()[1]
    assert jacobian.shape==(20,21) and len(pivots)==20
    determinant=jacobian[:,list(pivots)].det()
    assert determinant==-s.Rational(696,78125)
    al,be,ga=map(s.Rational,(2,3,5));denominator=1-al*be*ga
    points=[(0,0,0,1),(1,0,0,1),(0,1,0,1),(0,0,1,1),(1+al,1,0,1),(1,0,1+ga,1),(0,1+be,1,1),
            ((1+al+al*be)/denominator,(1+be+be*ga)/denominator,(1+ga+ga*al)/denominator,1)]
    assert len(set(points))==8
    for point in points:
        substitution=dict(zip(x,point))
        assert all(s.expand(f.subs(base).subs(substitution))==0 for f in forms[:3])
    values=[s.factor(forms[3].subs(base).subs(dict(zip(x,point)))) for point in points]
    assert values[0]==0 and all(values[1:])
    return {'result':'PASS','arithmetic':'exact rational tangent and evaluation matrices',
        'target':'Gr(4,9), quadrics through p=(0,0,0,1)','target_dimension':20,
        'stabilizer_GL_dimension':13,'coefficient_parameters':8,'tangent_shape':[20,21],
        'tangent_rank':20,'minor_columns':list(pivots),'nonzero_minor':str(determinant),
        'checked_parameters':[2,3,5,2,3,5,7,11],
        'other_seven_quadric_values':[str(z) for z in values[1:]],
        'conclusion':'The stabilizer action on the eight-parameter product-net chart is dominant in Gr(4,9).',
        'limits':['Generic coverage of quadric spaces, not indexed-AS moduli dimensions.',
                  'No regularity conclusion follows from the Jacobian test.']}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path);args=parser.parse_args()
    certificate=dominance_certificate()
    if args.output:args.output.write_text(json.dumps(certificate,indent=2)+'\n')
    print(json.dumps(certificate,indent=2))


if __name__=='__main__':main()
