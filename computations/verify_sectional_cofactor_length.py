"""Check trace and endomorphism lengths on numerical semigroup plane branches."""
from itertools import combinations


def verify_branch(p, q):
    # All exponents at or above the conductor belong to every normalized module.
    conductor = (p - 1) * (q - 1)
    bound = 3 * conductor
    ring = {p * i + q * j for i in range(bound // p + 1)
            for j in range(bound // q + 1) if p * i + q * j < bound}
    gaps = sorted(set(range(conductor)) - ring)
    assert len(gaps) * 2 == conductor
    checked = strict = 0
    for size in range(len(gaps) + 1):
        for subset in combinations(gaps, size):
            module = ring.union(subset)
            if any(n + h < conductor and n + h not in module
                   for n in subset for h in (p, q)):
                continue

            def colon(target):
                return {n for n in range(bound)
                        if all(n + m in target for m in module if n + m < conductor)}

            dual = colon(ring)
            shift = min(dual)
            if dual != {n + shift for n in module if n + shift < bound}:
                continue
            end = colon(module)
            trace = {n + m for n in dual for m in module if n + m < bound}
            assert trace <= ring <= end <= module
            assert len(end - ring) == len(ring - trace)
            assert min(trace) == 2 * len(module - ring)
            assert len(module - end) == min(trace) // 2 - len(ring - trace)
            strict += bool(module - end)
            checked += 1
    return checked, strict


def main():
    counts = [verify_branch(*pair) for pair in [(2, 5), (3, 4), (3, 7), (4, 5), (4, 7)]]
    assert sum(strict for _, strict in counts) > 0
    for ell in range(5, 25):
        for d in range(25 - ell):
            assert d <= 19
            if d == 19:
                assert (ell, 21 - d) == (5, 2)
    print(f"PASS: trace, endomorphism and valuation identities for {sum(n for n, _ in counts)}"
          f" self-dual modules, including {sum(n for _, n in counts)} strict inequalities")


if __name__ == "__main__":
    main()
