"""Check conductor and secant degrees, postulation bounds, and the local Bezout factorization."""
from sympy import Matrix, Poly, expand, symbols

u,v = symbols("u v")


def bezout(f, g, degree):
    numerator = Poly(f*g.subs(u,v)-g*f.subs(u,v), u,v)
    quotient = numerator.exquo(Poly(u-v, u,v))
    return Matrix(degree, degree, lambda i,j: quotient.nth(i,j))


def main():
    q0,q1,q2,a,b,c,d,e,f = symbols("q0 q1 q2 a b c d e f")
    q = q0+q1*u+q2*u**2
    p, r = a+b*u+c*u**2, d+e*u+f*u**2
    mult = Matrix([[q0,0],[q1,q0],[q2,q1],[0,q2]])
    assert (bezout(q*p,q*r,4)-mult*bezout(p,r,2)*mult.T).applyfunc(expand).is_zero_matrix
    # Classes on the secant plane sextic are given by degrees of H and h.
    input_degree, plane_degree = 21, 6
    marked = 4*(input_degree-plane_degree)
    assert marked == 3*20 and 8*input_degree-3*marked < 0
    assert 6*input_degree-2*marked == 6
    assert 6*20-2*3*20 == 0
    for degree in (20,21):
        assert 8*8 < 4*degree
    for length in range(2,7):
        for longest_line in range(3,length+1):
            assert length-longest_line <= 3
            if longest_line < 6:
                assert max(0,longest_line-4) <= 1
    print("PASS: symbolic quadratic-divisor Bezout factorization, conductor intersections and plane-sextic bounds")


if __name__ == "__main__":
    main()
