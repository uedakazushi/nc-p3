"""Finite AS certification over four one-sided, quantum-commuting parameters."""

from fractions import Fraction
from itertools import combinations, product

from two_window_ext import PathRow, add_scaled, quotient_images
from verify_central_window_cutoff import apply_word, koszul_exact


def product_section(arrows, first, second):
    value = {}
    for (u, a), (v, b) in product(first.items(), second.items()):
        add_scaled(value, apply_word(arrows, 0, {0: 1}, u+v), a*b)
    return value


def parameter_relations(windows, sections, degrees, exchanges):
    length = max(degrees[i]+degrees[j] for i, j in combinations(range(4), 2))
    for start in range(4):
        arrows = PathRow(windows, start, length).arrows
        for i, j in combinations(range(4), 2):
            q = exchanges[i, j]
            xy = product_section(arrows, sections[i][start], sections[j][(start+degrees[i]) % 4])
            yx = product_section(arrows, sections[j][start], sections[i][(start+degrees[j]) % 4])
            add_scaled(yx, xy, -q)
            if not q or yx:
                return False
    return True


def right_quotient(windows, sections, degrees, start, last):
    # Only a parameter at the first position imposes a new right-module relation.
    dims, arrows = [1], []
    for n in range(1, last+1):
        rows = []
        def append_rule(rule, at, index):
            row = {}
            for word, c in rule.items():
                for j, v in apply_word(arrows, at, {index: 1}, word[:-1]).items():
                    k = 4*j + word[-1]
                    row[k] = row.get(k, 0) + c*v
            rows.append({k: v for k, v in row.items() if v})
        for section, degree in zip(sections, degrees):
            if degree == n:
                append_rule(section[start], 0, 0)
        if n >= 2:
            for index, rule in product(range(dims[n-2]), windows[(start+n-2) % 4]):
                append_rule(rule, n-2, index)
        size, images = quotient_images(4*dims[-1], rows, Fraction(1))
        dims.append(size)
        arrows.append(images)
    return dims


def certify(windows, sections, degrees, exchanges):
    """The standard finite Frobenius check on the quadratic dual is assumed."""
    if len(sections) != 4 or len(degrees) != 4 or min(degrees) < 1:
        raise ValueError("four positive parameter degrees are required")
    for section, degree in zip(sections, degrees):
        if len(section) != 4 or any(len(w) != degree for part in section for w in part):
            raise ValueError("each section must have the stated degree at all four starts")
    if not parameter_relations(windows, sections, degrees, exchanges):
        return "parameter relations failed"
    expected = [1]
    for degree in degrees:
        coeffs = [0]*(len(expected)+degree-1)
        for i, v in enumerate(expected):
            for j in range(degree):
                coeffs[i+j] += v
        expected = coeffs
    dims = [right_quotient(windows, sections, degrees, i, len(expected)) for i in range(4)]
    if any(d[-1] for d in dims):
        return "finite parameter cut unverified"
    if any(d != expected+[0] for d in dims):
        return "non-AS"
    seen = []
    for i in range(4):
        rotation = windows[i:] + windows[:i]
        if rotation not in seen:
            if not koszul_exact(rotation, 0, sum(degrees)):
                return "non-AS"
            seen.append(rotation)
    return "AS"


def verify():
    # A tensor product of the anticommuting plane and an infinite-order quantum plane.
    pairs = list(combinations(range(4), 2))
    q = {p: Fraction(1) for p in pairs}
    q[0, 1], q[2, 3] = -1, 2
    relations = [{(j, i): Fraction(1), (i, j): -q[i, j]} for i, j in pairs]
    windows = [relations]*4
    sections = [[part]*4 for part in ({(0,): 1, (1,): 1}, {(1, 1): 1}, {(2,): 1}, {(3,): 1})]
    exchanges = {pair: Fraction(1) for pair in pairs}
    exchanges[2, 3] = 2
    degrees = (1, 2, 1, 1)
    assert certify(windows, sections, degrees, exchanges) == 'AS'
    rejected = dict(exchanges)
    rejected[2, 3] = 3
    assert certify(windows, sections, degrees, rejected) == 'parameter relations failed'
    repeated = [*sections[:2], sections[2], sections[2]]
    assert certify(windows, repeated, degrees, {p: 1 for p in pairs}) == 'finite parameter cut unverified'
    assert right_quotient(windows, sections, degrees, 0, 2) == [1, 1, 0]


if __name__ == '__main__':
    verify()
    print('PASS: a nonnormal parameter in an infinite-order quantum example, degree-five AS and rejection controls')
