#!/usr/bin/env python3
"""Exact adjacent-window overlap conditions and periodic Ore normalization."""
import argparse
import json
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement
from sympy import Matrix,Poly,expand,symbols,zeros,eye,kronecker_product
from verify_polynomial_ore import primal_window
from verify_quartic_obstruction import insertion_columns,exact_rank,relation_matrix

P=list(combinations_with_replacement(range(3),2))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    x=symbols('x0:3');d=symbols('d0:18')
    delta=[sum(d[6*i+j]*x[a]*x[b] for j,(a,b) in enumerate(P)) for i in range(3)]
    cubic=list(combinations_with_replacement(range(3),3))
    rows=[]
    for a,b in combinations(range(3),2):
        expression=Poly(expand(x[a]*delta[b]-x[b]*delta[a]),*x)
        for word in cubic:
            monomial=x[word[0]]*x[word[1]]*x[word[2]]
            value=expression.coeff_monomial(monomial)
            rows.append([value.coeff(variable) for variable in d])
    condition=Matrix(rows)
    assert condition.shape==(30,18) and condition.rank()==15
    radial=zeros(18,3)
    for i in range(3):
        for j in range(3):
            radial[6*i+P.index(tuple(sorted((i,j)))),j]=1
    assert condition*radial==zeros(30,3)
    cycle=zeros(120,72)
    for start in range(4):
        for i in range(30):
            for j in range(18):
                cycle[30*start+i,18*start+j]+=condition[i,j]
                cycle[30*start+i,18*((start+1)%4)+j]-=condition[i,j]
    assert cycle.rank()==45
    # Four periodic windows, with nonzero average radial correction.
    base=Matrix([1]+[0]*17)
    ell=[Matrix(v) for v in ([0,0,0],[1,0,0],[0,1,0],[0,0,1])]
    fields=[base+radial*v for v in ell]
    mean=sum(fields,zeros(18,1))/4
    meanell=sum(ell,zeros(3,1))/4
    residual=[v-meanell for v in ell]
    shifts=[zeros(3,1)]
    for v in residual:
        shifts.append(shifts[-1]-v)
    assert shifts[4]==shifts[0]
    windows=[primal_window(list(field)) for field in fields]
    target=relation_matrix(primal_window(list(mean)))
    annihilator=Matrix.hstack(*target.T.nullspace()).T
    maps=[]
    for shift in shifts[:4]:
        transform=eye(4)
        for i in range(3):transform[i,3]=shift[i]
        maps.append(transform)
    for n in range(4):
        image=kronecker_product(maps[n],maps[(n+1)%4])*relation_matrix(windows[n])
        assert annihilator*image==zeros(10,6)
        assert 64-exact_rank(insertion_columns(windows,n,3))==20
    # A genuine nonradial adjacent jump destroys the standard cubic dimension.
    bad=[primal_window([0]*18) for _ in range(4)]
    invalid=[0]*18;invalid[P.index((1,1))]=1
    bad[1]=primal_window(invalid)
    dimensions=[64-exact_rank(insertion_columns(bad,n,3)) for n in range(4)]
    assert any(value<20 for value in dimensions)
    certificate={'result':'PASS','arithmetic':'exact rational polynomial/linear algebra',
        'adjacent_condition_shape':[30,18],'adjacent_condition_rank':15,
        'four_window_condition_shape':[120,72],'four_window_condition_rank':45,
        'regular_coefficient_dimension':27,'periodic_shift_witness':'PASS',
        'nonregular_control_cubic_dimensions':dimensions,
        'limits':['27 is a coefficient-space dimension, not indexed moduli dimension.',
                  'The theorem covers the specified commuting Ore presentation only.']}
    if args.output:
        args.output.write_text(json.dumps(certificate,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(certificate,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
