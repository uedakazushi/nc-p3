#!/usr/bin/env python3
"""Exact checks for the seven-parameter algebraic AS-regular Z-family.

This verifies identities over QQ[u,v,r,s,a,b,c], not truncated power series.
The global AS argument and the open-set nonperiodicity proof are in papers/nc-p3.tex.
"""
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import argparse
import json
import sympy as sp
from pbw_complex import PAIRS, MON2, TRIPLES, add, clifford0, matrices, rank
from verify_moduli_jets import quotient_basis

PARAMETERS = sp.symbols('u v r s a b c')
BASE_POINT = (0, 1, 0, 0, 0, 0, 0)
SAMPLE = tuple(map(sp.Rational, ('13/3', '-13/6', '-1/4', '1/2',
                               '13/2', '13/5', '-13/5')))


def rules(parameters=PARAMETERS, epsilon=1):
    u, v, r, s, a, b, c = map(sp.sympify, parameters)
    out = clifford0()
    out[(3, 0)] = {(0, 3): -1, (1, 1): u, (2, 2): v,
                   (1, 2): epsilon*a, (0, 0): -epsilon*a/2}
    out[(3, 1)] = {(1, 3): -1, (0, 0): 1, (2, 2): 1,
                   (1, 1): r, (0, 1): epsilon*a, (0, 2): epsilon*b}
    out[(3, 2)] = {(2, 3): -1, (1, 1): s, (0, 1): epsilon*c}
    return {p: {m: sp.sympify(x) for m, x in rhs.items() if x != 0}
            for p, rhs in out.items()}


def check_polynomial_pbw():
    rr = [rules(epsilon=e) for e in (1, -1)]
    weights = (1, 2, 2, 4)
    key = lambda word: (sum(weights[j] for j in word), word)
    for rs in rr:
        for pair, rhs in rs.items():
            assert all(key(mon) < key(pair) for mon in rhs)
    @lru_cache(None)
    def normal(word, start):
        for j in range(len(word)-1):
            if word[j] > word[j+1]:
                out = {}
                for mon, coefficient in rr[(start+j) % 2][word[j:j+2]].items():
                    for new, value in normal(word[:j]+mon+word[j+2:], start).items():
                        add(out, new, sp.expand(coefficient*value))
                return {m: sp.expand(x) for m, x in out.items() if sp.expand(x) != 0}
        return {word: sp.Integer(1)}
    count = 0
    for start in (0, 1):
        for triple in TRIPLES:
            word = tuple(reversed(triple)); difference = {}
            for j, sign in ((0, 1), (1, -1)):
                for mon, coefficient in rr[(start+j) % 2][word[j:j+2]].items():
                    for new, value in normal(word[:j]+mon+word[j+2:], start).items():
                        add(difference, new, sp.expand(sign*coefficient*value))
            assert all(sp.expand(x) == 0 for x in difference.values())
            count += 1
    assert rules(tuple(map(sp.Integer, BASE_POINT))) == clifford0()
    return count


def check_parameter_rank():
    N, M = matrices(clifford0())
    alpha = sp.Matrix([rules().get((b, a), {}).get(mon, 0)
                       for a, b in PAIRS for mon in MON2])
    J = alpha.jacobian(PARAMETERS).subs(dict(zip(PARAMETERS, BASE_POINT)))
    assert rank(N[0].row_join(N[1])) == 30
    assert rank(N[0].row_join(N[1]).row_join(J)) == 37
    assert rank((N[0]+N[1]).row_join(J[:, :4])) == 19
    assert rank((N[0]-N[1]).row_join(J[:, 4:])) == 18
    assert (M[0]+M[1])*J[:, :4] == sp.zeros(80, 4)
    assert (M[0]-M[1])*J[:, 4:] == sp.zeros(80, 3)
    H = quotient_basis(N, M, -1)
    G = sp.Matrix.hstack(*(N[0]-N[1]).columnspace())
    coordinates = G.row_join(H).gauss_jordan_solve(J[:, 4:])[0][-6:, :]
    expected = sp.Matrix([[0, 0, 0], [0, -1, 0], [-1, 0, 0],
                          [0, -1, 0], [0, 0, 0], [0, 0, 1]])
    assert coordinates == expected
    a, b, c = PARAMETERS[4:]
    o = expected*sp.Matrix([a, b, c])
    assert sp.expand(o[0]*o[5]+o[1]*o[2]-o[2]*o[3]+o[3]*o[4]) == 0
    return {'pair_orbit_rank': 30, 'pair_orbit_and_parameter_rank': 37,
            'parameter_classes': 7, 'even_classes': 4, 'alternating_classes': 3,
            'odd_coordinates': ['0', '-b', '-a', '-b', '0', 'c']}


def point_matrix(point, parameters, epsilon):
    w, y, z, x = point
    u, v, r, s, a, b, c = parameters
    e = epsilon
    return sp.Matrix([[y, w, 0, 0], [z, 0, w, 0], [-w, z, y, 0],
                      [x+e*a*w/2, -u*y, -v*z-e*a*y, w],
                      [-w, x-r*y-e*a*w, -z-e*b*w, y],
                      [0, -s*y-e*c*w, x, z]])


def sextic(t, parameters, epsilon=1):
    _, _, r, s, a, b, c = parameters
    return 8*s*t**6+8*epsilon*c*t**5+4*r*t**4+4*epsilon*a*t**3-4*t**2-2*epsilon*b*t+1


def check_point_geometry():
    u, v, r, s, a, b, c = PARAMETERS
    t = sp.symbols('T', nonzero=True)
    z = -1/(2*t)
    for e in (1, -1):
        obstruction = z*(1-z**2-r*t**2)+s*t**3+e*(a/2-b*z**2+c*t**2)
        assert sp.expand(8*t**3*obstruction-sextic(t, PARAMETERS, e)) == 0
    w, y, zz = sp.symbols('w y z')
    mat = point_matrix((w, y, zz, sp.Symbol('x')), PARAMETERS, 1)
    assert sp.expand(mat[:3, :3].det()+w*(w**2+2*y*zz)) == 0
    roots = [sp.Rational(1, 5), sp.Rational(1, 2), sp.Integer(1),
             sp.Integer(5), sp.Integer(-1), sp.Rational(-1, 2)]
    assert sp.expand(sextic(t, SAMPLE)-4*sp.prod(t-q for q in roots)) == 0
    U, V, R, S, A, B, C = SAMPLE
    assert R**3+8*S == sp.Rational(255, 64)
    assert A**2+4*U*V == sp.Rational(169, 36)
    mon2 = list(combinations_with_replacement(range(3), 2))
    counts = []
    for e in (1, -1):
        projected = []
        for root in roots:
            yy = e*root; z = -1/(2*yy)
            total = -U*yy**2-V*z**2
            difference = (1-z**2-R*yy**2-e*(A*yy+B*z))/yy
            x = (total-difference)/2; xp = (total+difference)/2
            p = (1, yy, z, x); q = sp.Matrix([1, -yy, -z, xp])
            mat = point_matrix(p, SAMPLE, e)
            assert mat*q == sp.zeros(6, 1) and rank(mat) == 3
            projected.append(p[:3])
        for ell in (sp.Integer(e), sp.Integer(2*e)):
            assert V*ell**2+e*A*ell-U == 0
            total = R-ell**2; difference = S/ell
            x = (total-difference)/2; xp = (total+difference)/2
            p = (0, 1, ell, x); q = sp.Matrix([0, 1, -ell, xp])
            mat = point_matrix(p, SAMPLE, e)
            assert mat*q == sp.zeros(6, 1) and rank(mat) == 3
            projected.append(p[:3])
        evaluation = sp.Matrix([[p[i]*p[j] for i, j in mon2] for p in projected])
        conics = [ids for ids in combinations(range(8), 6)
                  if evaluation[list(ids), :].det() == 0]
        assert conics == [tuple(range(6))]
        assert rank(evaluation[:6, :]) == 5
        counts.append({'epsilon': e, 'simple_points_away_from_X': 8,
                       'six_point_conic_subsets': [list(ids) for ids in conics]})
    # Inversion T -> rho/T would force r^3 + 8s = 0; checked symbolically.
    rho2 = -r/(8*s)
    assert sp.factor(4*r*rho2**2+1/(2*s)) == (r**3+8*s)/(16*s**2)
    return {'sample': dict(zip(map(str, PARAMETERS), map(str, SAMPLE))),
            'conic_roots_plus': list(map(str, roots)),
            'line_roots_plus': ['1', '2'], 'line_roots_minus': ['-1', '-2'],
            'inversion_obstruction': '255/64', 'conic_checks': counts}


def verify():
    return {'date': '2026-11-15', 'version': 9, 'result': 'PASS',
            'arithmetic': 'exact rational polynomial identities',
            'base_point': dict(zip(map(str, PARAMETERS), BASE_POINT)),
            'polynomial_overlaps': check_polynomial_pbw(),
            'weights_WYZX': [1, 2, 2, 4],
            'parameter_rank': check_parameter_rank(),
            'point_geometry': check_point_geometry(),
            'proved_with_manuscript': [
                'Free PBW family over affine seven-space; all fibers are standard AS-regular.',
                'A nonempty Zariski-open subset has minimal period two.',
                'The family has seven-dimensional image in indexed isomorphism classes.'],
            'not_claimed': [
                'The family is an irreducible component of the full moduli space.',
                'All six alternating directions have been algebraically integrated.',
                'The full Clifford formal local ring or all period-one structures are determined.',
                'Novelty in the literature or independent peer review.']}


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
