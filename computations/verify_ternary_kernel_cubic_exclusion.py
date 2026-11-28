"""Verify the rank-one equations restricted to a horizontal cubic and after the row operation."""

from itertools import combinations_with_replacement

import sympy as sp


def main():
    a, b, l, m, f = sp.symbols("a b l m f")
    v = sp.Matrix(sp.symbols("v:3"))
    first = sp.Matrix(sp.symbols("p:3"))
    second = sp.Matrix(sp.symbols("q:3"))
    corners = 2*a*b*(v*v.T)
    column = (a*l+b*m)*v
    discriminant = 2*a*b*f-(a*l+b*m)**2
    assert (f*corners-column*column.T-discriminant*v*v.T).applyfunc(sp.expand).is_zero_matrix

    block = first*second.T + second*first.T
    reduced_column = l*first
    cubic = f*second-l**2*first/2
    assert (f*block-reduced_column*reduced_column.T
            -first*cubic.T-cubic*first.T).applyfunc(sp.expand).is_zero_matrix
    # Cubic support cycles whose components are smooth of degrees 1, 2, 3.
    cycles = [parts for count in range(1, 4) for parts in combinations_with_replacement((1, 2, 3), count)
              if sum(parts) == 3]
    assert cycles == [(3,), (1, 2), (1, 1, 1)]
    assert 10 > 2*3+2 and 10-3 > 2*2+2 and 3*3 < 10
    print("Horizontal determinant and transformed cubic equations verified for the three cubic cycles.")


if __name__ == "__main__":
    main()
