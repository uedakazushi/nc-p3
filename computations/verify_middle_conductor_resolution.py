"""Check the subresultant Hilbert--Burch degrees and the binary column syzygies."""
from itertools import combinations
from sympy import Matrix, Poly, Symbol, binomial, expand, simplify

t = Symbol("t")


def columns(f, g):
    rows = [[0]*6 for _ in range(7)]
    j = 0
    for coeff in (f,g):
        for i in range(8-len(coeff)):
            for k, val in enumerate(coeff):
                rows[i+k][j] = val
            j += 1
    assert j == 6
    return Matrix(rows)


def main():
    for f,g in (([1,2,3,5],[3,1,7,1,2,6]),([2,5,1,7,3],[1,0,8,2,3])):
        mat = columns(f,g)
        cof = Matrix([(-1)**i*mat.extract([j for j in range(7) if j != i],range(6)).det()
                      for i in range(7)])
        assert cof != Matrix.zeros(7,1) and mat.T*cof == Matrix.zeros(6,1)
        squarefree = Poly((t-1)*(t-2)*(t-3),t)
        fs = squarefree*Poly(t**(len(f)-4),t)
        gs = squarefree*Poly((t+1)**(len(g)-4),t)
        assert columns(list(reversed(fs.all_coeffs())),list(reversed(gs.all_coeffs()))).rank() == 4
    term = lambda n: expand((n+3)*(n+2)*(n+1)/6)
    assert simplify(term(t)-7*term(t-6)+6*term(t-7) - (21*t-49)) == 0
    print("PASS: both binary subresultants, six linear syzygies, double vanishing at common cubic factors, and conductor Hilbert polynomial")


if __name__ == "__main__":
    main()
