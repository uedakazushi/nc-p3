#!/usr/bin/env python3
"""Exact classification of degree-one derivations of the elliptic Clifford B.

The tests use polynomial arithmetic in t,h, without specialization.
They verify a basis for ordinary derivations and for parity-twisted
sigma-derivations, together with a nonconstant Hesse j-invariant.
"""
from functools import lru_cache
from itertools import combinations_with_replacement
from pathlib import Path
import argparse
import json
import sympy as sp

T, H = sp.symbols('t h')
U, V, Z2 = sp.symbols('U V T')
CENTRAL = (U,V,Z2)
BASIS2 = tuple(combinations_with_replacement(range(3), 2))
# The cubic Clifford basis is (W,Y,Z) multiplied by the three central squares,
# followed by WYZ. All cubic products have unique coordinates in this basis.
CUBIC = tuple((i, z) for i in range(3) for z in CENTRAL) + ((3, None),)
RELS = (
    {(1,0): 1, (0,1): 1, (2,2): -T},
    {(2,0): 1, (0,2): 1, (1,1): -H},
    {(2,1): 1, (1,2): 1, (0,0): -1},
)
CORRECTION = {(1,0): T*Z2, (2,0): H*V, (2,1): U}

@lru_cache(None)
def normal(word):
    """Clifford normal form in the central-square algebra."""
    for i in range(len(word)-1):
        a, b = word[i:i+2]
        if a == b:
            return {m: sp.expand(CENTRAL[a]*c)
                    for m,c in normal(word[:i]+word[i+2:]).items()}
        if a > b:
            out = {}
            for m,c in normal(word[:i]+(b,a)+word[i+2:]).items():
                out[m] = out.get(m,0)-c
            for m,c in normal(word[:i]+word[i+2:]).items():
                out[m] = out.get(m,0)+CORRECTION[(a,b)]*c
            return {m: sp.expand(c) for m,c in out.items() if sp.expand(c)!=0}
    return {word:sp.Integer(1)}


def cubic(poly):
    output = {}
    for w,c in poly.items():
        for m,cc in normal(w).items():
            output[m] = sp.expand(output.get(m,0)+c*cc)
    entries = []
    for i,z in CUBIC:
        e = sp.expand(output.get((i,),0)).coeff(z) if z is not None else output.get((0,1,2),0)
        entries.append(sp.expand(e))
    return sp.Matrix(entries)


def relation_matrix(twisted=False):
    ind = sp.symbols('x0:18')
    images=[{BASIS2[j]:ind[6*i+j] for j in range(6)} for i in range(3)]
    rows=[]
    for rel in RELS:
        poly={}
        for word,c in rel.items():
            for position,letter in enumerate(word):
                for mon,value in images[letter].items():
                    name=word[:position]+mon+word[position+1:]
                    poly[name]=poly.get(name,0)+c*value*((-1)**position if twisted else 1)
        rows.extend(cubic(poly))
    return sp.Matrix(rows).jacobian(ind)


def coordinates(images):
    return sp.Matrix([images[i].get(m,0) for i in range(3) for m in BASIS2])


def ordinary_basis():
    half=sp.Rational(1,2)
    adW=[{}, {(0,1):2,(2,2):-T}, {(0,2):2,(1,1):-H}]
    adY=[{(0,1):-2,(2,2):T},{},{(1,2):2,(0,0):-1}]
    adZ=[{(0,2):-2,(1,1):H},{(1,2):-2,(0,0):1},{}]
    da=[{(1,2):1,(0,0):-half},{(0,1):1,(2,2):-T/2},{}]
    db=[{(0,1):-H,(2,2):H*T/2},{(0,2):1,(1,1):-H/2},{}]
    dc=[{(0,2):-T,(1,1):H*T/2},{},{(0,1):1,(2,2):-T/2}]
    return sp.Matrix.hstack(*[coordinates(images) for images in (adW,adY,adZ,da,db,dc)])


def sign_basis():
    cols=[]
    for i in range(3):
        for mon in ((0,0),(1,1),(2,2)):
            entries=[{} for _ in range(3)]
            entries[i][mon]=1
            cols.append(coordinates(entries))
    return sp.Matrix.hstack(*cols)


def verify():
    ordinary = relation_matrix()
    sign = relation_matrix(twisted=True)
    assert ordinary.shape == sign.shape == (30,18)
    D = ordinary_basis()
    P = sign_basis()
    assert ordinary*D == sp.zeros(30,6)
    assert sign*P == sp.zeros(30,9)

    # An everywhere nonzero minor of the ordinary constraint matrix,
    # and a direct-summand minor of its proposed solution basis.
    ord_rows=(0,1,2,3,4,5,9,10,11,12,19,20)
    ord_cols=(0,1,2,3,4,5,6,9,11,12,15,17)
    ord_det=sp.factor(ordinary.extract(ord_rows,ord_cols).det(method='domain-ge'))
    assert ord_det == -2048
    deriv_rows=(0,1,2,7,8,13)
    deriv_det=sp.factor(D.extract(deriv_rows,range(6)).det(method='domain-ge'))
    assert deriv_det == 4

    # Two minors certify that sigma-derivations have rank exactly 9
    # everywhere; their zero loci are disjoint, even before inverting 1+th.
    sign_cols=(1,2,4,7,8,10,13,14,16)
    first=(0,1,3,6,7,13,16,23,25)
    second=(0,1,2,3,4,5,6,7,12)
    det1=sp.factor(sign.extract(first,sign_cols).det(method='domain-ge'))
    det2=sp.factor(sign.extract(second,sign_cols).det(method='domain-ge'))
    assert sp.expand(det1-4*(H*T-2)*(H*T+4)**2)==0
    assert sp.expand(det2-3*H*T**3*(H*T-8)**2)==0
    assert P.extract((0,3,5,6,9,11,12,15,17),range(9)).det() in (1,-1)

    # An x=th coordinate is recovered up to finite ambiguity from the
    # intrinsic elliptic point-scheme cubic's j-invariant.
    x = sp.symbols('x')
    k3 = -(2-x)**3/(27*x)
    j = sp.factor(27*k3*(k3+8)**3/(k3-1)**3)
    expected=(x-2)**3*(x**3-6*x**2+228*x-8)**3/(x*(x-8)**3*(x+1)**6)
    assert sp.cancel(j-expected)==0
    a,b=sp.fraction(j)
    assert sp.degree(a,x)==12 and sp.degree(b,x)==10 and sp.gcd(a,b)==1

    return {
       'field':'char 0; on B(t,h) require 1+t*h != 0',
       'degree_one_derivations':{
          'ordinary':{'dimension':6,'constraint_rank':12,
                      'constant_constraint_minor':str(ord_det),
                      'basis':['ad_W','ad_Y','ad_Z','d_a','d_b','d_c'],
                      'constant_basis_minor':str(deriv_det)},
          'parity_sigma':{'dimension':9,'constraint_rank':9,
                          'solutions':'arbitrary central quadratic forms p_W,p_Y,p_Z',
                          'covering_minors':[str(det1),str(det2)]}
       },
       'ore_slice':{'raw_p_coefficients':9,'raw_d_coefficients':6,
                    'even_X_shift_gauge':3,'odd_X_shift_gauge':3,
                    'fixed_B_effective':9,'ternary_parameters':2,
                    'diagonal_scaling_gauge':2,'normal_slice_parameters':9},
       'elliptic_j':{'x':'t*h','formula':str(j),'map_degree':12,
           'usage':'necessary generic indexed isomorphism condition, not sufficient'},
       'limitations':['Not a global AS-regular Z-algebra classification',
           'Completeness only for the chosen Clifford B and specified Ore form',
           'No assertion of nine-family irreducible component maximality']
    }

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    data=verify()
    out=json.dumps(data,indent=2,ensure_ascii=False)+'\n'
    if args.output:
       args.output.parent.mkdir(parents=True,exist_ok=True)
       args.output.write_text(out,encoding='utf-8')
    print(out)
