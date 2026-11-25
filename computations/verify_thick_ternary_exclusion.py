"""Mixed Ore restrictions on the supporting line of each exceptional thick cubic."""

import sympy as sp
from ternary_ore_parameters import intersection, kernel, presentation_pencil


def verify(gamma, a, q, r, expected, vanishing):
    e = sp.eye(9)

    def word(i, j):
        return e[:, 3 * i + j]

    x, y, z = range(3)
    relations = sp.Matrix.hstack(
        word(x, z) - gamma * word(z, x) - a * word(z, z),
        word(y, z) - gamma * word(z, y),
        word(x, y) - word(y, x) - q * word(z, y) - r * word(z, z),
    )
    left, right, _, lift = presentation_pencil(relations)
    core = sp.eye(18)
    while True:
        pairs = kernel((left * core).row_join(right * core))
        new = intersection(core * pairs[: core.cols, :], core * pairs[core.cols :, :])
        if new.cols == core.cols:
            break
        core = new
    assert core.cols == expected
    u, v = sp.symbols("u v")
    evaluate = sp.Matrix([[u**2, u * v, 0, u * v, v**2, 0, 0, 0, 0]]) * lift
    g = evaluate * core[12:, :]
    h = u * evaluate * core[6:12, :] - v * evaluate * core[:6, :]
    assert (g if vanishing == "g" else h).applyfunc(sp.expand).is_zero_matrix
    print(
        "PASS",
        gamma,
        a,
        q,
        r,
        "core",
        expected,
        vanishing,
        "vanishes on the complete line",
    )


if __name__ == "__main__":
    verify(-1, 0, 0, 1, 9, "h")
    verify(1, 0, 0, 1, 12, "g")
    verify(1, 1, 1, 1, 9, "g")
    verify(1, 1, 2, 0, 12, "g")
