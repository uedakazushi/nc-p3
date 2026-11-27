"""Test cubic interpolation at the thirty nodes without choosing or splitting them."""
from itertools import product

import numpy as np
import sympy as sp
from flint import nmod_mpoly_ctx

from verify_thirty_central_nodes import DIGITS, MOD, central_form, lower_rank


def exponents(degree):
    return [e for e in product(range(degree + 1), repeat=5) if sum(e) == degree]


def certify(form):
    variables = sp.symbols("a b c d z")
    context = nmod_mpoly_ctx.get(tuple(map(str, variables)), MOD)
    quartic = context.from_dict({e: int(c) for e, c in sp.Poly(form, *variables).terms()})
    linear = sum((i + 1) * context.gen(i) for i in range(5))
    generators = [(quartic.derivative(i), 3) for i in range(5)] + [(linear**4, 4)]
    target = {exp: i for i, exp in enumerate(exponents(7))}
    rows = []
    for generator, degree in generators:
        terms = list(generator.to_dict().items())
        for multiplier in exponents(7 - degree):
            row = np.zeros(len(target), dtype=np.int64)
            for term, coeff in terms:
                row[target[tuple(i + j for i, j in zip(term, multiplier))]] = int(coeff)
            rows.append(row)
    counts = lower_rank(np.array(rows), DIGITS)
    assert counts == [298, 0, 32]
    assert sum(counts) == len(target)
    print("Cubic interpolation at all thirty nodes, degree-seven pivot counts:", counts)


if __name__ == "__main__":
    certify(central_form())
