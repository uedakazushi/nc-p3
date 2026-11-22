"""All nonordinary residue disks at 101 on the genus-five level c=3/4."""

from flint import fq_default_ctx, fq_default_poly_ctx
from rank_fifteen_octic import finite_windows
from verify_rank_fourteen_geometric_reduction import Scalar, verify_central_profile
from verify_rank_fifteen_pencil_infinity import verify as verify_infinity
from functools import cache
from itertools import product
from flint import fmpq, fmpq_poly, ctx
from concurrent.futures import ProcessPoolExecutor
from verify_rank_fifteen_pencil_ramified_ends import minors, grassmann_chart
from verify_rank_fourteen_all_ordinary import (
    DIMENSIONS,
    gaussian_jets,
    saturated_projection,
)
from verify_rank_fourteen_central_cut import quotient_layers
from verify_klein_center_lift import field_paths, central_equations


def local_relations(Jet, k, sign_r, sign_w, color, *, level=None, parameter=None):
    def power(v, n):
        if not n:
            return Jet.constant(1)
        half = power(v, n // 2)
        square = half * half
        return square * v if n % 2 else square

    n = Jet(fmpq_poly([0, 1])) if parameter is None else parameter
    c = Jet.constant(fmpq(3, 4)) if level is None else level
    L2 = 1 + c * c

    def sqrt(a, x):
        v = x if isinstance(x, Jet) else Jet.constant(x)
        for _ in range(Jet.precision.bit_length()):
            v = (v + a / v) * Jet.constant(fmpq(1, 2))
        assert not (v * v - a)
        return v

    K = sqrt(1 + L2 * n * n, k)
    r = K + c * n * n + sign_r * n * sqrt(L2 + 2 * c * K + c * c * n * n, c + k)
    w = K - c * n * n + sign_w * n * sqrt(L2 - 2 * c * K + c * c * n * n, c - k)
    z = r * r
    la = w / r
    a, b, e, f = 1 + n, (K + c * n) / r, 1 - n, (K - c * n) / r
    if color:
        a, b, e, f = z * b, a, z * f, e
    vals = (z, la, a, b, e, f)
    unit = Jet.constant(1)

    @cache
    def entry(g):
        v = Jet.constant(0)
        for m, c in g.terms():
            term = Jet.constant(int(c))
            for a, n in zip(vals, m):
                if n:
                    term = term * power(a, n)
            v = v + term
        return v

    out = []
    loss = 0
    for side in minors():
        rules = []
        for h, group in side:
            ps = {key: entry(g) for key, g in group.items()}
            if not any(v.unit() for v in ps.values()):
                ps = {key: v.divide_by_parameter() for key, v in ps.items()}
                loss = 1
            pivot = next(key for key, v in ps.items() if v.unit())
            order = {pivot: ps[pivot], **ps}
            for row in grassmann_chart(order, unit):
                rules.append((h, {(a, a ^ h): v for a, v in row.items() if v}))
        out.append(rules)
    Jet.precision -= loss
    return out


def local_central_sections(
    relations,
    Jet,
    *,
    project=saturated_projection,
    invariant_pairs=((1, 3, 5), (2, 7, 11)),
):
    """Lift the five-dimensional center after saturating the actual path lattices."""
    unit = Jet.constant(1)
    widths = (1, 1, 4, 10, 20, 35)
    dims = [(), (1,) * 4, (2, 2, 3, 3), (5,) * 4, (9, 10, 8, 8), (14,) * 4]
    counter = 0
    losses = []

    def projection(rows):
        nonlocal counter
        # field_paths visits four positions and their four characters at each length.
        n, h = counter // 16 + 1, counter % 4
        counter += 1
        free, images, l = project(rows, widths[n], dims[n][h], Jet)
        losses.append(l)
        result = {
            j: {free[i]: -v for i, v in row.items() if v}
            for j, row in enumerate(images)
            if j not in free
        }
        if counter % 16 == 0:
            Jet.precision -= max(losses)
            losses.clear()
        return result

    basis, image = field_paths(relations, unit, projection)
    cuts = [[] for p in range(4)]
    for h, pairs in enumerate((invariant_pairs, ((1, 0), (0, 1)))):
        rows, size = central_equations(basis, image, h)
        free, images, l = project(rows, size, 3 - h, Jet)
        losses.append(l)
        for weights in pairs:
            coords = [
                sum((weights[i] * v for i, v in row.items()), Jet.constant(0))
                for row in images
            ]
            offset = 0
            for p in range(4):
                words = basis[p, 4, h]
                cuts[p].append(
                    (
                        h,
                        {
                            w: coords[offset + j]
                            for j, w in enumerate(words)
                            if coords[offset + j]
                        },
                    )
                )
                offset += len(words)
            assert offset == size
    Jet.precision -= max(losses)
    assert Jet.precision > 0
    return cuts


def run(job):
    k, sign_w, color = job
    ctx.cap = 32
    Jet = gaussian_jets(initial_precision=32)
    relations = local_relations(Jet, k, 1, sign_w, color)
    cuts = local_central_sections(relations, Jet)
    N = Jet.precision
    print("SECTIONS", job, N, flush=True)
    for p in range(4):
        Jet.precision = N
        losses = []

        def project(n, h, rows, size):
            free, images, l = saturated_projection(rows, size, DIMENSIONS[n][h], Jet)
            losses.append(l)
            return free, images

        for n in quotient_layers(relations, cuts, p, 13, Jet.constant(1), project):
            Jet.precision -= max(losses)
            losses.clear()
            assert Jet.precision > 0
        print("GENUS5CUT", job, p, Jet.precision, flush=True)


def complementary_ends(color):
    F = fq_default_ctx(101, 4, "a")
    c = F(3) / 4
    L = F(5) / 4
    s = c + 1 / c
    G = fq_default_poly_ctx(F)
    x = G.gen()
    roots = (x * x + s * x - 1).roots()
    assert len(roots) == 2

    def both(v):
        a = v.sqrt()
        assert not a.is_zero()
        return a, -a

    for (k, _), r in product(roots, both(F(-1))):
        for H, v in product(both(k * k - 1), both(4 * k * k - 1)):
            w = 2 * k + v
            n = H / L
            z = r * r
            a, b, e, f = 1 + n, (k + c * n) / r, 1 - n, (k - c * n) / r
            if color:
                a, b, e, f = z * b, a, z * f, e
            assert verify_central_profile(
                finite_windows((z, w / r, a, b, e, f), F.one(), Scalar), Scalar(F.one())
            )
    print("GOLDEN101", color, flush=True)
    for r, w in product(range(2), repeat=2):
        for v in (2, 50):
            verify_infinity(color, r, w, 101, v)
    print("INFINITY101", color, flush=True)


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=2) as pool:
        list(pool.map(run, product((1, -1), (1, -1), range(2))))
        list(pool.map(complementary_ends, range(2)))
