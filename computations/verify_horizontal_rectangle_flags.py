"""Exact evaluation and pure-relation tests on the three horizontal rectangles."""
from functools import lru_cache
from itertools import combinations, product
import os
import subprocess

import sympy as s

from verify_parallel_invertible_rectangles import K, completed, matrix
from verify_parallel_horizontal_moduli import window as second_window
from verify_parallel_third_rectangle import third_rectangle


u, v, z, w = s.symbols('u v z w')
a = s.symbols('a:d')


def family(shape):
    if shape == 3:
        return third_rectangle(u, v, z), (u, v, z), (u, v, z, z-1, u*v-1)
    t = (u+v)/(1+u*v)
    x = u*(v+z)/(1+u*v)
    common = (u, v, u+v, 1+u*v)
    def factors(q):
        return (v+q, q-u, u*q-1, v*q+1)
    if shape == 1:
        W = completed(76, t, u*t, v*t, x, 1, 1, v*(1-u*z)/(u+v))
        return W, (u, v, z), (*common, z, *factors(z))
    y = u*(w+v)/(1+u*v)
    return second_window(t, u*t, v*t, x, y), (u, v, z, w), (
        *common, z*w-1, *factors(z), *factors(w))


def fmt(seq):
    return ','.join(str(f).replace('**', '^') for f in seq)


def singular_unit(equations, variables, units):
    # A vanishing unit can often be detected before saturation.
    code = (f'LIB "elim.lib";ring r=0,({fmt(variables)}),dp;short=0;'
            f'ideal H=std(ideal({fmt(equations)}));ideal N={fmt(units)};int j;'
            'for(j=1;j<=size(N);j++){if(reduce(N[j],H)==0){print("PASS");quit;};};'
            f'H=std(sat(H,ideal({fmt([s.prod(units)])}))[1]);'
            'if(H[1]==1){print("PASS");};quit;')
    result = subprocess.run([os.getenv('SINGULAR_BIN', 'Singular'), '-q'],
                            input=code, text=True, capture_output=True, check=True)
    assert 'PASS' in result.stdout.splitlines() and '?' not in result.stdout, result.stdout


def relation_frames(W, pos, unit):
    frames = []
    for h in K:
        m = matrix(W, pos, h)
        d = 2 if h < 2 else 1
        J, I = next((J, I) for J, I in product(combinations(K, d), repeat=2)
                    if unit(m.extract(I, J).det()))
        normalized = (m[:, J]*m.extract(I, J).inv()).applyfunc(s.factor)
        assert (normalized*m[I, :] - m).applyfunc(s.factor) == s.zeros(4)
        frames.append(normalized)
    return frames


def eval_support(frames, reverse, nz, params, units):
    x = [1 if j == nz[0] else a[j] if j in nz else 0 for j in K]
    extra = tuple(a[j] for j in nz[1:])
    variables = (*extra, *params)
    units = (*units, *extra)
    def monic(f):
        return s.Poly(f, *variables).monic().as_expr()
    inverted = set(map(monic, units))
    def factors(value):
        num = s.fraction(s.factor(value))[0]
        return frozenset(monic(g) for g, _ in s.factor_list(num)[1]) - inverted
    M = s.Matrix([[col[j if reverse else j^h]*x[j^h] for j in K]
                  for h, frame in enumerate(frames) for col in frame.columnspace()])
    # Once an invertible minor is found, its bordering minors suffice.
    pivot = next(((I, J) for I, J in product(combinations(range(6), 2), combinations(K, 2))
                  if (det := M.extract(I, J).det()) != 0 and s.factor(det) != 0
                  and not factors(det)), None)
    constraints = set()
    for I, J in product(combinations(range(6), 3), combinations(K, 3)):
        if pivot and not (set(pivot[0]) <= set(I) and set(pivot[1]) <= set(J)):
            continue
        value = s.factor(M.extract(I, J).det())
        if value:
            residual = factors(value)
            if not residual:
                return
            constraints.add(residual)
    assert constraints
    # The common zero set is the union over minimal choices of vanishing factors.
    cases = {frozenset()}
    for residual in sorted(constraints, key=len):
        cases = {old | {f} for old in cases for f in residual}
        cases = {case for case in cases if not any(other < case for other in cases)}
    for case in cases:
        singular_unit(sorted(case, key=str), variables, units)


def pure_third(frames, params, units, unit):
    x, y = s.symbols('x:4'), s.symbols('y:4')
    equations = []
    for h, frame in enumerate(frames):
        I = next(i for i in combinations(K, frame.cols) if unit(frame[i, :].det()))
        color = s.Matrix([x[j]*y[j^h] for j in K])
        residual = color-frame*frame[I, :].inv()*color[list(I), :]
        for val in residual:
            num, denom = s.fraction(s.factor(val))
            assert unit(denom)
            if num:
                equations.append(num)
    for i, j in product(K, repeat=2):
        chart = dict.fromkeys((*x[:i], *y[:j]), 0) | {x[i]: 1, y[j]: 1}
        singular_unit([f.xreplace(chart) for f in equations],
                      (*x[i+1:], *y[j+1:], *params), units)


def check(shape):
    W, params, units = family(shape)
    def monic(f):
        return s.Poly(f, *params).monic().as_expr()
    inverted = set(map(monic, units))
    @lru_cache(None)
    def unit(f):
        f = s.factor(f)
        return f != 0 and all(monic(g) in inverted for term in s.fraction(f)
                               for g, _ in s.factor_list(term)[1])
    for pos in (K if shape == 3 else (0,)):
        frames = relation_frames(W, pos, unit)
        if shape == 3:
            pure_third(frames, params, units, unit)
        else:
            c, d = frames[2][:, 0], frames[3][:, 0]
            # The two cross-block determinants are linear in the squared coefficients.
            assert unit(c[0]*c[1]*d[2]*d[3]-c[2]*c[3]*d[0]*d[1])
            for frame in frames[:2]:
                assert unit(frame[:2, :].det()) and unit(frame[2:, :].det())
        for reverse in (False, True):
            for n in range(1, 5):
                for nz in combinations(K, n):
                    eval_support(frames, reverse, nz, params, units)
    print(f'PASS: horizontal rectangle {shape}', flush=True)


if __name__ == '__main__':
    for shape in (1, 2, 3):
        check(shape)
