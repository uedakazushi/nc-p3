"""Certify all eight fixed branch windows at slopes (3,2) at one split prime."""

import sympy as s
from flint import nmod, nmod_mat
from verify_klein_elliptic_completion import coeffs, gram_twist
from verify_klein_center_lift import field_paths, central_equations
from verify_central_window_cutoff import quotient_dimensions
from verify_fixed_line_elliptic_curves import branch_polynomial

p = 307


def red(x, prime=p):
    return int(s.numer(x)) * pow(int(s.denom(x)), -1, prime) % prime


def matrix(A, prime=p):
    return nmod_mat(A.rows, A.cols, [red(x, prime) for x in A], prime)


def check(sig, aa):
    a = s.Integer(aa)
    z = s.Integer(3)
    l = s.Integer(2)
    A = s.Matrix([[a, 1], [sig * z, sig * a]])
    co, P, Q = coeffs(*A, z, l)
    x = -P[1] / (2 * P[0])
    y = -Q[1] / (2 * Q[0])
    f, F, g, G, e, E, h, H = co
    B = s.Matrix([[e + E * y, f + F * x], [g + G * x, h + H * y]])
    check_window(A, B, z, l, (sig, aa), p)


def check_window(A, B, z, l, label, prime):
    w, opp, rev = gram_twist(A, B, z, l)
    rel = []
    for j, M in enumerate((w, rev, w.T, rev.T)):
        rr, rank = matrix(M.T, prime).rref()
        assert rank == 6
        rules = []
        for ii in range(rank):
            row = {
                (k // 4, k % 4) if j % 2 == 0 else (k % 4, k // 4): nmod(
                    int(rr[ii, k]), prime
                )
                for k in range(16)
                if rr[ii, k]
            }
            h, jj = next(iter(row))
            rules.append((h ^ jj, row))
        rel.append(rules)
    basis, image = field_paths(rel, nmod(1, prime))
    counts = [
        sum(len(basis[i, n, h]) for h in range(4)) for i in range(4) for n in (4, 5)
    ]
    print(label, "counts", counts, "opp", matrix(opp, prime).rank(), flush=True)
    assert counts == [35, 56] * 4 and matrix(opp, prime).rank() == 16
    sections = []
    for h, choices in [
        (0, ((1, 3, 5), (2, 7, 11))),
        (1, ((1, 0), (0, 1))),
        (2, ()),
        (3, ()),
    ]:
        rows, n = central_equations(basis, image, h)
        m = nmod_mat(
            len(rows), n, [int(row.get(j, 0)) for row in rows for j in range(n)], prime
        )
        ker, dim = m.nullspace()
        print("ker", h, dim, flush=True)
        assert dim == (3, 2, 0, 0)[h]
        for choice in choices:
            vec = [sum(ker[i, j] * c for j, c in enumerate(choice)) for i in range(n)]
            start = 0
            section = []
            for i in range(4):
                words = basis[i, 4, h]
                section.append(
                    {
                        w: nmod(int(vec[start + j]), prime)
                        for j, w in enumerate(words)
                        if vec[start + j]
                    }
                )
                start += len(words)
            sections.append(section)
    windows = [[r for _, r in rs] for rs in rel]
    for i in range(4):
        ds = quotient_dimensions(windows, sections, i, nmod(1, prime))
        assert ds == [1, 4, 10, 20, 31, 40, 44, 40, 31, 20, 10, 4, 1, 0]
        print(label, "cut", i, ds, flush=True)


if __name__ == "__main__":
    for sig, roots in [(1, [102, 153, 298, 301]), (-1, [80, 119, 188, 227])]:
        t = s.Symbol("t")
        quartic = branch_polynomial(sig, t, s.Integer(3), s.Integer(2))
        assert [v for v in range(p) if quartic.eval(v) % p == 0] == roots
        for a in roots:
            check(sig, a)
