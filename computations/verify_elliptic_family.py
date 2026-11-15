#!/usr/bin/env python3
"""Exact certificates for the nine-parameter elliptic Clifford--Ore family.

All checks use rational/polynomial arithmetic.  The all-fiber AS argument,
elliptic-divisor nonperiodicity proof, and generic-finiteness argument are in
nc-p3.tex.  No truncated power series or numerical root approximations occur.
"""
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path
import argparse
import json
import sympy as sp
from pbw_complex import PAIRS, MON2, add, clifford0, matrices, rank

PARAMETERS = sp.symbols('u v r s a b c t h')
u, v, r, s, a, b, c, t, h = PARAMETERS
BASE_POINT = (0, 1, 0, 0, 0, 0, 0, 0, 0)
CENTERS = sp.symbols('U V T')


def rules(epsilon=1):
    """Ordered degree-two coordinates; NOT a claimed quadratic Groebner basis."""
    e = sp.sympify(epsilon)
    out = clifford0()
    out[(1, 0)] = {(0, 1): -1, (2, 2): t}
    out[(2, 0)] = {(0, 2): -1, (1, 1): h}
    out[(3, 0)] = {(0, 3): -1, (1, 1): u+e*c*t*h/2,
                   (2, 2): v+e*b*h*t/2, (1, 2): e*a,
                   (0, 0): -e*a/2, (0, 1): -e*b*h, (0, 2): -e*c*t}
    out[(3, 1)] = {(1, 3): -1, (0, 0): 1, (2, 2): 1-e*a*t/2,
                   (1, 1): r-e*b*h/2, (0, 1): e*a, (0, 2): e*b}
    out[(3, 2)] = {(2, 3): -1, (1, 1): s, (0, 1): e*c,
                   (2, 2): -e*c*t/2}
    return {pair: {m: sp.sympify(x) for m, x in rhs.items() if x != 0}
            for pair, rhs in out.items()}


def clifford_checks():
    """Independent finite Clifford normal forms, central squares not guessed."""
    U, V, T = CENTERS
    correction = {(1, 0): t*T, (2, 0): h*V, (2, 1): U}

    @lru_cache(None)
    def normal(word):
        for j in range(len(word)-1):
            left, right = word[j:j+2]
            if left == right:
                return {m: sp.expand(CENTERS[left]*co)
                        for m, co in normal(word[:j]+word[j+2:]).items()}
            if left > right:
                result = {}
                for m, co in normal(word[:j]+(right, left)+word[j+2:]).items():
                    add(result, m, -co)
                for m, co in normal(word[:j]+word[j+2:]).items():
                    add(result, m, correction[(left, right)]*co)
                return {m: sp.expand(co) for m, co in result.items() if sp.expand(co) != 0}
        return {word: sp.Integer(1)}

    def reduce_polynomial(polynomial):
        result = {}
        for word, coefficient in polynomial.items():
            for mon, value in normal(word).items():
                add(result, mon, sp.expand(coefficient*value))
        return {mon: sp.expand(co) for mon, co in result.items() if sp.expand(co) != 0}

    R = [ {(1, 0): 1, (0, 1): 1, (2, 2): -t},
          {(2, 0): 1, (0, 2): 1, (1, 1): -h},
          {(2, 1): 1, (1, 2): 1, (0, 0): -1} ]
    d = [{(1, 2): a, (0, 0): -a/2, (0, 1): -b*h,
          (2, 2): b*h*t/2, (0, 2): -c*t, (1, 1): c*t*h/2},
         {(0, 1): a, (2, 2): -a*t/2, (0, 2): b, (1, 1): -b*h/2},
         {(0, 1): c, (2, 2): -c*t/2}]

    def derivative(polynomial, images, twisted=False):
        result = {}
        for word, coefficient in polynomial.items():
            for j, letter in enumerate(word):
                for mon, value in images[letter].items():
                    sign = (-1)**j if twisted else 1
                    add(result, word[:j]+mon+word[j+1:], coefficient*value*sign)
        return reduce_polynomial(result)

    assert all(derivative(relation, d) == {} for relation in R)
    omega = reduce_polynomial({(0, 1, 2): 1, (2, 2, 2): -t/2,
                               (1, 1, 1): h/2, (0, 0, 0): -sp.Rational(1, 2)})
    for j, scale in enumerate((2*a, -2*b, 2*c)):
        actual = derivative({(j, j): 1}, d)
        expected = {m: sp.expand(scale*co) for m, co in omega.items()}
        assert actual == expected
    for j, k in product(range(3), repeat=2):
        if j != k:
            assert reduce_polynomial({(j, j, k): 1, (k, j, j): -1}) == {}
    # A completely arbitrary central quadratic assignment is a sign derivation.
    p = sp.symbols('p0:9')
    images = [{(j, j): p[3*i+j] for j in range(3)} for i in range(3)]
    assert all(derivative(relation, images, twisted=True) == {} for relation in R)
    # Cliff(R,Q) itself has a confluent square-free normal form over R[U,V,T].
    reductions = {(j, j): {(): CENTERS[j]} for j in range(3)}
    reductions.update({pair: {tuple(reversed(pair)): -1, (): co}
                       for pair, co in correction.items()})
    overlaps = 0
    for word in product(range(3), repeat=3):
        if word[:2] not in reductions or word[1:] not in reductions:
            continue
        answers = []
        for j in (0, 1):
            pol = {}
            for mon, co in reductions[word[j:j+2]].items():
                add(pol, word[:j]+mon+word[j+2:], co)
            answers.append(reduce_polynomial(pol))
        assert answers[0] == answers[1]
        overlaps += 1
    return {'clifford_cubic_overlaps': overlaps,
            'ordinary_derivation_relations': 3,
            'arbitrary_central_sign_derivation_relations': 3,
            'center_derivatives': ['2*a*Omega', '-2*b*Omega', '2*c*Omega']}


def parameter_checks():
    N, M = matrices(clifford0())
    coordinates = sp.Matrix([rules().get((j, i), {}).get(mon, 0)
                             for i, j in PAIRS for mon in MON2])
    base = dict(zip(PARAMETERS, BASE_POINT))
    J = coordinates.jacobian(PARAMETERS).subs(base)
    even = J[:, [0, 1, 2, 3, 7, 8]]
    odd = J[:, [4, 5, 6]]
    assert rank(N[0].row_join(N[1])) == 30
    assert rank(N[0].row_join(N[1]).row_join(J)) == 39
    assert rank((N[0]+N[1]).row_join(even)) == 21
    assert rank((N[0]-N[1]).row_join(odd)) == 18
    assert (M[0]+M[1])*even == sp.zeros(80, 6)
    assert (M[0]-M[1])*odd == sp.zeros(80, 3)
    # The previous family is retained literally, not only up to isomorphism.
    from verify_algebraic_family import rules as old_rules
    for e in (1, -1):
        restricted = {pair: {mon: sp.expand(co.subs({t: 0, h: 0}))
                             for mon, co in rhs.items() if sp.expand(co.subs({t: 0, h: 0})) != 0}
                      for pair, rhs in rules(e).items()}
        assert restricted == old_rules(epsilon=e)
    return {'single_relation_space_orbit_rank': 30,
            'orbit_and_nine_parameter_rank': 39,
            'even_parameter_classes': 6, 'alternating_parameter_classes': 3,
            'seven_family_recovered_at': 't=h=0'}


def curve_checks():
    w, y, z = sp.symbols('w y z')
    p = sp.Matrix([w, y, z])
    K = sp.Matrix([[y, w, -t*z], [z, -h*y, w], [-w, z, y]])
    F = w**3+h*y**3+t*z**3+(2-t*h)*w*y*z
    assert sp.expand(K.det()+F) == 0
    assert sp.factor((2-t*h)**3+27*t*h) == -(t*h+1)**2*(t*h-8)
    # A generic rational point test of the graph/involution, not a period-one test.
    K1 = K.subs({t: 1, h: 1})
    p1 = sp.Matrix([1, -1, 0]); q1 = sp.Matrix([1, 1, -1])
    assert K1.subs(dict(zip((w, y, z), p1)))*q1 == sp.zeros(3, 1)
    assert K1.subs(dict(zip((w, y, z), q1)))*p1 == sp.zeros(3, 1)
    # Derivation section is the displayed nonzero invariant quadratic in p cross q.
    q = sp.Matrix(sp.symbols('qw qy qz'))
    cross = p.cross(q); A = cross[0]/2; B = cross[2]/2; C = -cross[1]/2
    ds = sp.Matrix([a*A-b*h*B-c*t*C, a*B+b*C, c*B])
    formula = (a*(cross[0]**2+cross[1]*cross[2])
               -b*(h*cross[0]*cross[2]+cross[1]**2)
               +c*(t*cross[0]*cross[1]+cross[2]**2))/2
    assert sp.expand(cross.dot(ds)-formula) == 0
    return {'point_cubic': str(F),
            'smooth_elliptic_open': 't*h*(1+t*h)*(t*h-8) != 0',
            'fixed_point_free_involution': 'proved in manuscript from K(p)*p=0',
            'invariant_derivation_section_verified': True}


def boundary_checks():
    """A necessary exclusion: the original ternary algebra jumps in degree 3."""
    words = list(product(range(3), repeat=3))
    def degree_three_rank(tv, hv):
        rels = [{(1, 0): 1, (0, 1): 1, (2, 2): -tv},
                {(2, 0): 1, (0, 2): 1, (1, 1): -hv},
                {(2, 1): 1, (1, 2): 1, (0, 0): -1}]
        cols = []
        for relation in rels:
            for letter in range(3):
                for left in (True, False):
                    pol = {}
                    for mon, coefficient in relation.items():
                        word = (letter,)+mon if left else mon+(letter,)
                        add(pol, word, coefficient)
                    cols.append(sp.Matrix([pol.get(word, 0) for word in words]))
        return rank(sp.Matrix.hstack(*cols))
    regular = degree_three_rank(1, 1)
    excluded = degree_three_rank(1, -1)
    assert (regular, excluded) == (17, 15)
    return {'regular_sample_t_h': [1, 1], 'regular_relation_rank': regular,
            'excluded_sample_t_h': [1, -1], 'excluded_relation_rank': excluded,
            'degree_three_dimensions': [27-regular, 27-excluded]}


def verify():
    return {'date': '2026-11-15', 'version': 10, 'result': 'PASS',
            'arithmetic': 'exact QQ and polynomial identities',
            'parameter_base': 'Spec k[u,v,r,s,a,b,c,t,h,(1+t*h)^(-1)]',
            'clifford_and_derivations': clifford_checks(),
            'parameter_ranks': parameter_checks(),
            'elliptic_geometry': curve_checks(),
            'excluded_boundary': boundary_checks(),
            'scope_limits': [
                'No claim of a full moduli irreducible component or dimension 14.',
                'Generic period two uses the manuscript elliptic-divisor proof, not a sampled numerical test.',
                'No claim of a quadratic Groebner basis when t and h are both nonzero.',
                'No independent peer review, formalization, or confirmed novelty.']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    certificate = verify()
    text = json.dumps(certificate, ensure_ascii=False, indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding='utf-8')
    print(text)
