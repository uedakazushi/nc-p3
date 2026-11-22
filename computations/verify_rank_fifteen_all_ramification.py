"""Dual Fitting support and residual finite-field checks on the two ramification curves."""

from concurrent.futures import ProcessPoolExecutor
from itertools import product

from rank_fifteen_ramification import presentation
from rank_fourteen_dual import fitting_support, polynomial
from verify_rank_fifteen_rational_residues import boundary_at_101, verify_residue


def dual_rules(family):
    out = []
    for side in family:
        rules = []
        for h, row in side:
            assert all(not b and not d for a, b, c, d in row.values())
            rules.append(
                (
                    h,
                    {
                        word: (polynomial(a), polynomial(c))
                        for word, (a, b, c, d) in row.items()
                    },
                )
            )
        out.append(rules)
    return out


def support(color, start):
    rules, cuts = map(dual_rules, presentation(color))
    return fitting_support(rules, cuts, start, 10, radicand="-t^4+14*t^2-1")


def run():
    boundary = boundary_at_101()
    with ProcessPoolExecutor(max_workers=4) as pool:
        jobs = {
            (c, p): pool.submit(support, c, p) for c, p in product(range(2), repeat=2)
        }
        for color in range(2):
            supports = [jobs[color, p].result() for p in (0, 1)]
            assert all(len(group) == (29, 26)[color] for group in supports)
            bad = set.union(*supports)
            residual = sorted(
                (q for q in bad - boundary if len(q) > 2), key=lambda f: (len(f), f)
            )
            assert (
                sorted(len(g) - 1 for g in residual)
                == ([2, 2, 3, 3, 4, 4, 4, 6], [2, 2, 3, 3, 4, 4, 12])[color]
            )
            print(
                "Curve",
                color,
                "support",
                len(bad),
                "higher-degree ordinary factors",
                len(residual),
                flush=True,
            )
            checks = [pool.submit(verify_residue, color, q) for q in residual]
            for q, check in zip(residual, checks):
                check.result()
                print(color, q, "certified", flush=True)


if __name__ == "__main__":
    run()
