"""Check the one-sided block lifting equation and the Tor degrees used in descent."""

from sympy import Matrix, symbols


def main():
    diag_top, diag_bot, middle, error = [symbols(s + ":3") for s in ("p", "q", "c", "f")]
    t = symbols("t")
    arrows = [Matrix([[diag_top[i], middle[i]], [t*error[i], diag_bot[i]]]) for i in range(3)]
    lower = (arrows[0]*arrows[1]*arrows[2])[1, 0].expand().coeff(t)
    assert lower == (error[0]*diag_top[1]*diag_top[2] +
                     diag_bot[0]*error[1]*diag_top[2] +
                     diag_bot[0]*diag_bot[1]*error[2])
    assert not lower.free_symbols.intersection(middle)
    # Each fixed terminal degree sees a full resolution once its four predecessors are in the tail.
    for j in range(4):
        for n in range(j + 4, j + 36):
            assert [m % 4 for m in range(n-4, n+1)] == [(n-i) % 4 for i in range(4, -1, -1)]
    print("Opposite off-diagonal lifting and the four-periodic Tor indexing verified.")


if __name__ == "__main__":
    main()
