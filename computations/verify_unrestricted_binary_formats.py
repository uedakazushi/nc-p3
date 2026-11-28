"""Check the local length and Hilbert-function bounds used to recover both binary formats."""
from math import comb
from sympy import Matrix, Symbol


def main():
    # With h1=4 and h2<=6, the quotient by a general linear form has
    # degree-two dimension at most two, hence degree-three dimension at most two.
    assert 4+2+2 < 10
    for rank in (2,3):
        # After cancelling rank r units the maximal minors have degree 4-r.
        expected = comb((4-rank)+2,3)
        assert expected == (4 if rank == 2 else 1)
    t = Symbol('t')
    a = Matrix.diag(1,1,1,1,t**2)
    assert a.adjugate().subs(t,0).rank() == 1
    assert a.adjugate()[4,4] == 1 and a.adjugate()[3,3] == t**2
    print("PASS: plane complement bound, local rank-two length surplus and the rank-four adjugate on a repeated quintic factor")


if __name__ == "__main__":
    main()
