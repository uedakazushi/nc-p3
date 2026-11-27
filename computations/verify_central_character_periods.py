"""Check the index-four orbit/weight bound and the codimension-one Hilbert comparison."""
from itertools import product
from sympy import Rational, expand, symbols


def main():
    single_nonthin = []
    for h,m,d,weight in product(range(2,5),range(1,5),range(2,5),range(1,5)):
        if h % m == 0 and d >= 2*m and weight*d*(h//m) <= 4:
            single_nonthin.append((h,m,d,weight))
    assert single_nonthin == [(2,1,2,1), (2,2,4,1)]
    n,b = symbols('n b')
    def hilb(t):
        return Rational(2,3)*t**3 + t**2 + Rational(7,3)*t + 1
    def difference(poly):
        return expand(poly - 2*poly.subs(n,n-b) + poly.subs(n,n-2*b))
    assert difference(hilb(2*n+2)-hilb(2*n)) == 2*difference(hilb(2*n+1)-hilb(2*n))
    print('All nonthin orbit profiles at index four and the hyperplane-curve degree comparison verified.')


if __name__ == '__main__':
    main()
