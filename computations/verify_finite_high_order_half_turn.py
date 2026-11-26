"""Check the off-mark graph and residual equations in the two unit coordinates."""
from sympy import Matrix, symbols


def main():
    x, y, s, t, n, m = symbols("x y s t n m")
    c, k = symbols("c k", nonzero=True)
    first, second = Matrix([1, x]), Matrix([1, y])
    lower = first*m - second*n
    skew = Matrix.hstack(first, second).det()
    assert lower[0] == m-n
    assert (lower[1].subs(m, n) + n*skew).expand() == 0
    d = s-c*t
    solved = {m: n, y: x+n*d/k}
    assert (s*m-c*n*t-k*skew).subs(solved).expand() == 0
    assert (lower[1].subs(solved) + n**2*d/k).expand() == 0
    assert (y-x).subs(solved).subs(s, c*t).expand() == 0
    print("PASS: off-mark linkage multiplier n^2 and the linear graph ideal")


if __name__ == "__main__":
    main()
