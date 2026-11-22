"""Identify the 32 deepest Pauli boundary points and their quantum polynomial window."""

from itertools import combinations, product

from sympy.polys.domains import QQ_I

from verify_rank_fifteen_parameter_orbits import dot, pattern, template_group

ONE = QQ_I.one
I = QQ_I(0, 1)
UNITS = (ONE, I, -ONE, -I)
PAIRS = tuple(combinations(range(4), 2))
GAMMA = pattern(0, 0, ONE, -I, I, ONE, I, ONE)


def verify_orbit():
    def normalize(row):
        return tuple(c / row[4] for c in row)

    stratum = {
        normalize(pattern(0, 0, ONE, m, a, -m/a, c, -m/c))
        for m, a, c in product((I, -I), UNITS, UNITS)
    }
    orbit = {
        normalize([I**e * GAMMA[j] for j, e in zip(perm, exp)])
        for perm, exp in template_group()
    }
    assert len(stratum) == 32 and orbit == stratum
    assert sum(bool(c) for c in GAMMA) == 12


def transformed_tensor():
    def entry(a, b, c, d):
        if a ^ c != b ^ d:
            return QQ_I.zero
        return sum(
            GAMMA[(a ^ c) + 4*v] * (-1)**(dot(v, c) + dot((0, 3, 2, 1)[v], d))
            for v in range(4)
        )

    tensor = {}
    for word in product(range(4), repeat=4):
        total = QQ_I.zero
        for bits in product(range(2), repeat=4):
            letters = [(a & 2) + b for a, b in zip(word, bits)]
            sign = (-1)**sum((a & 1)*b for a, b in zip(word, bits))
            total += sign * entry(*letters)
        if total:
            value = total / 8
            assert value.x**2 + value.y**2 == 2
            tensor[word] = value
    assert len(tensor) == 24
    return tensor


def quantum_rows(tensor):
    parameters = []
    for p in range(4):
        rows = {}
        for word, value in tensor.items():
            a, b, c, d = word[p:] + word[:p]
            a ^= 3*(p % 2)
            b ^= 3*(1 - p % 2)
            rows.setdefault((c, d), {})[a, b] = value
        side = {}
        for row in rows.values():
            (a, b), reverse = sorted(row)
            assert a < b and reverse == (b, a)
            q = -row[b, a] / row[a, b]
            assert q in UNITS and side.get((a, b), q) == q
            side[a, b] = q
        assert set(side) == set(PAIRS)
        parameters.append(side)
    return parameters


def verify_gauge(parameters):
    q = parameters[0]
    assert [q[pair] for pair in PAIRS] == [ONE, ONE, -I, -I, ONE, -ONE]
    period = [ONE]*4
    for side in parameters:
        ratios = [ONE] + [q[0, j] / side[0, j] for j in (1, 2, 3)]
        assert all(side[a, b]*ratios[b]/ratios[a] == q[a, b] for a, b in PAIRS)
        period = [a*b for a, b in zip(period, ratios)]
    assert period == [ONE]*4
    full = dict(q)
    full.update({(b, a): 1/v for (a, b), v in q.items()})
    full.update({(a, a): ONE for a in range(4)})
    central = set()
    for powers in product(range(5), repeat=4):
        if sum(powers) != 4:
            continue
        characters = [
            full[0, j]**powers[0] * full[1, j]**powers[1]
            * full[2, j]**powers[2] * full[3, j]**powers[3]
            for j in range(4)
        ]
        if len(set(characters)) == 1:
            central.add(powers)
    assert central == {(1, 1, 1, 1)} | {
        tuple(4*int(k == j) for k in range(4)) for j in range(4)
    }


if __name__ == "__main__":
    verify_orbit()
    verify_gauge(quantum_rows(transformed_tensor()))
    print("32 boundary points in one orbit; four-periodic quantum gauge and five quartic sections.")
