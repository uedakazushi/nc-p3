"""Check every possible two-coordinate evaluation defect on an exceptional exchange sheet."""

from itertools import combinations, product

import sympy as sp


def data():
    a, b, c, d, z = sp.symbols("a b c d z")
    lines = [
        ({a: b, c: b}, b),
        ({a: z * b, c: z * z * b}, b),
        ({c: b, d: b}, b),
        ({d: z * b, c: z * z * b}, b),
        ({a: b, d: z * b}, b),
        ({a: z * b, d: b}, b),
        ({c: z * a, d: z * a}, z * a),
        ({a: z * d, c: z * d}, z * d),
    ]
    return (a, b, c, d, z), lines


def check():
    (a, b, c, d, z), lines = data()
    coords = (a, b, c, d)
    first = sp.Matrix([[a, b], [c, d]])
    D = first.det()
    L = -first.inv().T * sp.Matrix([[1, z], [1, 1]])
    K = -first.inv() * sp.Matrix([[1, 1], [1, z]])
    quadric = (a + d) * (z * b + c) - (z + 1) * (a * d + b * c)
    r = sp.Symbol("r", nonzero=True)
    free = sp.symbols("y:4")
    for num, (right, i, j) in enumerate(product(range(2), repeat=3)):
        axis = sp.zeros(4, 1)
        axis[i] = 1
        axis[2 + j] = r
        (u, v, h, k), (p, q, g, n) = (free, axis) if right else (axis, free)
        thick = [
            *(sp.Matrix([u * p, v * q]) - L * sp.Matrix([h * g, k * n])),
            *(sp.Matrix([u * q, v * p]) - K * sp.Matrix([h * n, k * g])),
        ]
        thin = [u * g + v * n + h * p + z * k * q, u * n + v * g + h * q + k * p]
        low = sp.solve(thin, free[:2], dict=True)
        assert len(low) == 1 and len(low[0]) == 2
        eq = sp.Matrix([t.subs(low[0]) for t in thick]).jacobian(free[2:])
        eq = [sp.factor(e * r * D) for e in eq if e != 0]
        restrictions, value = lines[num]
        constants = [
            sp.factor(f / r**2) for f in eq if sp.cancel(f / r**2).has(r) == False
        ]
        linear = sp.linear_eq_to_matrix(constants, coords)[0]
        assert linear.rank() == 2 and all(
            sp.expand(f.xreplace(restrictions)) == 0 for f in constants
        )
        assert sp.expand(quadric.xreplace(restrictions)) == 0
        after = [sp.factor(f.xreplace(restrictions)) for f in eq]
        assert all(sp.cancel(f.subs(r**2, value)) == 0 for f in after)
        det = sp.factor(D.xreplace(restrictions))
        assert det != 0 and value != 0
        assert any(
            sp.cancel(f * value / (det * (r**2 - value)))
            in (1, -1, z, -z, 1 / z, -1 / z)
            for f in after
            if f != 0
        )
    # Each line lies in two of the cross-difference planes and in the double quadric.
    planes = (a - b, a - c, a - z * b, z * a - c, b - d, c - d, z * b - d, c - z * d)
    incidences = [
        tuple(j for j, f in enumerate(planes) if sp.expand(f.xreplace(sub)) == 0)
        for sub, _ in lines
    ]
    assert len(set(incidences)) == 8 and all(len(t) == 2 for t in incidences)
    print(
        "All eight evaluation defects are exactly the specified lines and every one lies on the ramification quadric."
    )


if __name__ == "__main__":
    check()
