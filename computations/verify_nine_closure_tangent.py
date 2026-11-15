#!/usr/bin/env python3
"""Exact rank certificate for the nine-parameter Ore orbit closure at Clifford point.

No finite-order truncation can prove an irreducible-component claim.  Here two
independent tangent vectors are produced by explicit algebraic curves in
seven-family orbit closure and the nine-family tangent rank comparison is exact.
"""
from pathlib import Path
import argparse
import json
import sympy as sp
from pbw_complex import PAIRS, MON2, clifford0, matrices, rank
from verify_elliptic_family import PARAMETERS, BASE_POINT, rules as nine_rules
from verify_algebraic_family import PARAMETERS as SEVEN_PARAMETERS, BASE_POINT as SEVEN_BASE, rules as seven_rules


def full_matrix(rules, symbols):
    return sp.Matrix([rules.get((b,a),{}).get(m,0) for a,b in PAIRS for m in MON2]).jacobian(symbols)


def verify():
    N, M = matrices(clifford0())
    orbit = N[0].row_join(N[1])
    gamma = full_matrix(nine_rules(), PARAMETERS)
    J = gamma.subs(dict(zip(PARAMETERS, BASE_POINT)))
    seven = full_matrix(seven_rules(), SEVEN_PARAMETERS).subs(dict(zip(SEVEN_PARAMETERS,SEVEN_BASE)))
    # Boundary algebra C_d: YW+WY = (-1)^i A ZX + d(Z^2+X^2/V^2).
    # At the Clifford specialization K=V=B=S=1 the d and A tangent vectors are:
    eta = sp.Matrix([int((a,b)==(0,1) and m in ((2,2),(3,3)))
                     for a,b in PAIRS for m in MON2])
    xi = sp.Matrix([int((a,b)==(0,1) and m==(2,3))
                    for a,b in PAIRS for m in MON2])
    assert N[0].shape == (60,16) and J.shape == (60,9)
    assert (M[0]+M[1])*eta == sp.zeros(80,1)
    assert (M[0]-M[1])*xi == sp.zeros(80,1)
    Q9=orbit.row_join(J)
    Q7=orbit.row_join(seven)
    values={
        'clifford_orbit':rank(orbit),
        'seven_orbit_plus_parameters':rank(Q7),
        'nine_orbit_plus_parameters':rank(Q9),
        'nine_plus_eta':rank(Q9.row_join(eta)),
        'nine_plus_xi':rank(Q9.row_join(xi)),
        'nine_plus_eta_plus_xi':rank(Q9.row_join(eta).row_join(xi)),
    }
    assert list(values.values()) == [30,37,39,40,40,41], values
    even=J[:,[0,1,2,3,7,8]]
    odd=J[:,[4,5,6]]
    parity={
      'even_orbit':rank(N[0]+N[1]),
      'odd_orbit':rank(N[0]-N[1]),
      'even_plus_nine':rank((N[0]+N[1]).row_join(even)),
      'odd_plus_nine':rank((N[0]-N[1]).row_join(odd)),
      'even_plus_eta':rank((N[0]+N[1]).row_join(even).row_join(eta)),
      'odd_plus_xi':rank((N[0]-N[1]).row_join(odd).row_join(xi))
    }
    assert list(parity.values()) == [15,15,21,18,22,19],parity
    lam,aa,bb,cc=sp.symbols('lambda a b c')
    oo=(lam,-bb,-aa,-bb,sp.Integer(0),cc)
    Q=sp.expand(oo[0]*oo[5]+oo[1]*oo[2]-oo[2]*oo[3]+oo[3]*oo[4])
    assert Q==lam*cc
    return {
      'result':'PASS',
      'base_point':dict(zip(map(str,PARAMETERS),map(str,BASE_POINT))),
      'arithmetic':'exact QQ rational Gaussian elimination, no random/numeric fitting',
      'ranks':values,
      'fourier_even_odd_ranks':parity,
      'quadratic_obstruction_on_algebraic_tangent_span':str(Q),
      'geometric_conclusion':{
        'G':'PGL4 x PGL4, effective dimension 30',
        'Z9':'closure of G-orbit saturation of nine-parameter elliptic Ore relation locus',
        'dim_Z9':39,
        'tangent_Z9_at_Clifford_lower_bound':41,
        'singular_at_Clifford':True,
        'integrable_single_curve_classes_even':7,
        'integrable_single_curve_classes_odd':4,
        'integrable_single_curve_total':11,
      },
      'dependencies':['Boundary eta and xi genuinely tangent to Z7, since v14 explicit algebraic flat specialization and singular generic conjugation',
        'Z7 subset Z9 as t=h=0 subfamily',
        'The elliptic nine-family is a smooth open nine-parameter family, with generically independent relation-class parameters'],
      'not_claimed':['All eleven tangent classes integrate simultaneously in an 11-parameter smooth family',
        'The entire AS-regular moduli has dimension 39 or is singular',
        'The nine-parameter family is an irreducible component',
        'A new irreducible component outside Z9 was found']
    }

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    d=verify()
    msg=json.dumps(d,ensure_ascii=False,indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(msg,encoding='utf-8')
    print(msg)
