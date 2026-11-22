"""Every ordinary lift with rational x=K+H residue on the 9/16 level at 101."""

from itertools import product
from concurrent.futures import ProcessPoolExecutor
from flint import fq_default_ctx, nmod
from sympy import sqrt_mod
from verify_klein_center_lift import central_equations, echelon, field_paths
from verify_rank_fourteen_central_cut import central_sections, quotient_layers
from rank_fifteen_octic import finite_windows
from verify_rank_fourteen_geometric_reduction import Scalar, verify_central_profile

P = 101


def graded_profile(windows, unit):
    """Use centrality to impose the four quartics just at the initial position."""
    basis, image = field_paths(windows, unit)
    assert all(
        sum(len(basis[p, n, h]) for h in range(4)) == d
        for p in range(4)
        for n, d in enumerate((1, 4, 10, 20, 35, 56))
    )
    for h, d in enumerate((3, 2, 0, 0)):
        rows, width = central_equations(basis, image, h)
        assert width - len(echelon(rows)) == d
    sections = central_sections(basis, image, unit, ((1, 3, 5), (2, 7, 11)))

    def reduce(n, h, rows, size):
        piv = echelon(rows)
        free = [k for k in range(size) if k not in piv]
        columns = {v: i for i, v in enumerate(free)}
        return free, [
            {columns[k]: unit}
            if k in columns
            else {columns[i]: -v for i, v in piv[k].items() if i != k}
            for k in range(size)
        ]

    for start in range(4):
        for _ in quotient_layers(windows, sections, start, 13, unit, reduce):
            pass
    return True


def verify_residue(x0):
    t = nmod(x0, P)
    s = nmod(25, P) / 12
    k = (t + 1 / t) / 2
    roots = [sqrt_mod(int((k + v * s) ** 2 - 1), P, all_roots=True) for v in (1, -1)]
    if all(roots):
        unit = nmod(1, P)
        wrap = lambda x: x
        ordinates = [[nmod(y, P) for y in pair] for pair in roots]
    else:
        field = fq_default_ctx(P, 2, "a")
        unit = field.one()
        wrap = Scalar
        t = field(x0)
        s = 25 * unit / 12
        k = (t + 1 / t) / 2
        ordinates = []
        # Frobenius reverses all nonsplit square roots simultaneously.
        fixed = False
        for pair, v in zip(roots, (1, -1)):
            if pair:
                ordinates.append([field(y) for y in pair])
                continue
            y = ((k + s * v) ** 2 - 1).sqrt()
            assert y**101 == -y
            ordinates.append([y, -y] if fixed else [y])
            fixed = True
    c = 3 * unit / 4
    L = 5 * unit / 4
    H = (t - 1 / t) / 2
    n = H / L
    for A, B in product(*ordinates):
        r = k + (H * H + H * A) / s
        w = k + (-H * H + H * B) / s
        z = r * r
        assert not (
            z * (z * z - 1) * (w * w - 1) * (w * w - z) * (w * w * z - 1) == 0 * unit
        )
        for h in (0, 1):
            a, b, e, f = 1 + n, (k + c * n) / r, 1 - n, (k - c * n) / r
            if h:
                a, b, e, f = z * b, a, z * f, e
            windows = finite_windows((z, w / r, a, b, e, f), unit, wrap)
            try:
                assert graded_profile(windows, wrap(unit))
            except AssertionError:
                # A fixed pair of invariant quartics can fail; test alternative pairs on the same fiber.
                assert verify_central_profile(windows, wrap(unit))
    print("RATIONAL x", x0, len(ordinates[0]) * len(ordinates[1]), flush=True)


def representatives():
    s = 25 * pow(12, -1, P) % P
    bad = {0}
    reps = []
    for x in range(1, P):
        k = (x + pow(x, -1, P)) * pow(2, -1, P) % P
        if (
            k in {0, 1, P - 1, *[(j * s + q) % P for j in (-1, 1) for q in (-1, 1)]}
            or (k * k + s * k - 1) % P == 0
        ):
            bad.add(x)
        elif x < pow(x, -1, P):
            reps.append(x)
    assert bad == {0, 1, 6, 10, 17, 33, 49, 52, 68, 84, 91, 95, 100} and len(reps) == 44
    return reps


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=3) as pool:
        list(pool.map(verify_residue, representatives()))
