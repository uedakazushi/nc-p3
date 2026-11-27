"""Euler determinant recurrence for four indexed polynomial-growth lattices."""
from math import comb
import sympy as sp


def verify():
    x, m, v = sp.symbols("x m v")
    numerator = v*x-12*m*x**2+(16*m-v)*x**3-4*m*x**4
    assert sp.expand(numerator.diff(x).subs(x,1)) == 8*m-2*v
    assert sp.cancel(numerator.subs(v,4*m)/(1-x)**8 - 4*m*x/(1-x)**5) == 0
    # Exercise the four distinct rows with nonsymmetric formal valuations.
    letters = [9,-5,3,1]
    top = sum(letters)//4
    relations = [7,-3,6*top-7,6*top+3]
    valuations = [[0] for _ in letters]
    ranks = [1,4,6,4,1]
    for length in range(1,31):
        current = []
        for i in range(4):
            syzygies = [0,letters[i],relations[i],4*top-letters[(i+3)%4],top]
            value = sum((-1)**(j+1)*(ranks[j]*valuations[(i+j)%4][length-j]
                                      + comb(length-j+3,3)*syzygies[j])
                        for j in range(1,min(length,4)+1))
            current.append(value)
        for row,value in zip(valuations,current):
            row.append(value)
        assert sum(current) == 4*top*comb(length+3,4)


if __name__ == '__main__':
    verify()
    print("PASS: determinant growth obstruction, generating series and four-start Euler recurrence")
