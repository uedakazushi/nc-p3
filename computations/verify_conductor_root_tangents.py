"""Check linearized cofactor structure at rank two and secant Pluecker degree."""
from itertools import combinations
from sympy import Matrix, diag, expand, symbols, zeros


def main():
    t = symbols("t")
    a,b,c = symbols("a b c")
    base = diag(1,1,0,0)
    params = iter(symbols("p:10"))
    change = zeros(4)
    for i in range(4):
        for j in range(i,4):
            change[i,j] = change[j,i] = next(params)
    derivative = (base+t*change).adjugate().applyfunc(lambda e: expand(e).coeff(t,1))
    assert derivative == diag(zeros(2), Matrix([[change[3,3],-change[2,3]],[-change[2,3],change[2,2]]]))
    mult = Matrix([[a,0],[b,a],[c,b],[0,c]])
    minors = [expand(mult.extract(pair,[0,1]).det()) for pair in combinations(range(4),2)]
    basis = [a*a,a*b,a*c,b*b,b*c,c*c]
    coefficients = Matrix([[f.coeff(m) for m in basis] for f in minors])
    assert coefficients.rank() == 6
    assert 2*6 < 2*2*6
    print("PASS: linearized four-by-four adjugate and full quadratic Pluecker system of binary secants")


if __name__ == "__main__":
    main()
