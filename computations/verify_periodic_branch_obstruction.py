"""Check degree budgets and the periodic first relation count."""
from math import comb

assert comb(5, 2) == comb(5, 3) == 10
assert 4 * (comb(5, 2) + comb(5, 3)) // 2 == 40
for pattern in [(2, 3, 3, 3), (3, 2, 2, 3), (3, 3, 3, 3)]:
    higher = [3-pattern[(n-2) % 4] for n in range(2, 42)]
    assert sum(higher) == 10 * pattern.count(2)
print("PASS: complete-intersection degree budget and periodic first Tor counts")
