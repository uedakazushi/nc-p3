"""Check the graded generation bounds for the rank-fifteen central cut."""

from math import comb

from sympy import Poly, symbols


def verify():
    t = symbols("t")
    quotient = Poly((1 + t + t**2 + t**3) ** 4, t)
    known = [1, 4, 10, 20, 31, 40, 44, 40]
    assert known == [quotient.nth(n) for n in range(8)]
    for n in range(8):
        generated = sum(comb(j + 3, 3) * known[n - 4 * j] for j in range(n // 4 + 1))
        assert generated == comb(n + 3, 3)
    floors = [
        comb(n + 3, 3) - 4 * known[n - 4] - 10 * known[n - 8] for n in range(8, 12)
    ]
    assert floors == [31, 20, 10, 4] == [quotient.nth(n) for n in range(8, 12)]
    print("Uncut dimensions through seven:", [comb(n + 3, 3) for n in range(8)])
    print("Quotient lower bounds through eleven:", floors)


if __name__ == "__main__":
    verify()
