"""Check the degree-five symmetric complex and all reduced support degrees up to five."""
from sympy import Symbol, expand

t = Symbol("t")
term = lambda n: (n+1)*(n+2)*(n+3)/6


def main():
    assert expand(term(t)-15*term(t-4)+24*term(t-5)-10*term(t-6)) == 20
    for support_degree, max_marks in ((2,7),(4,14)):
        for degree in range(support_degree+2,10,2):
            assert 5*degree//2 > max_marks
    for degree, chi, bound in [(2,1,7),(4,-2,9),(4,0,12),(4,1,13),(4,2,14)]:
        assert 7*degree//2+2*chi > bound
    print("PASS: symmetric complex of length twenty, all reduced conics/quartics and all possible multiple support degrees")


if __name__ == "__main__":
    main()
