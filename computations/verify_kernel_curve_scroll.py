"""Compute the tautological triple-divisor class on an ordered cube of a curve."""
from sympy import Poly, expand, factor, symbols

points = symbols("h0:3")
diagonals = symbols("d01 d02 d12")
degree, genus, parameter = symbols("d g t")
variables = (*points, *diagonals)


def diagonal_integral(exponents):
    components = [{i} for i in range(3)]
    marks = [i for i in range(3) for _ in range(exponents[i])]
    excess = 0
    for (a, b), power in zip(((0, 1), (0, 2), (1, 2)), exponents[3:]):
        for _ in range(power):
            left = next(group for group in components if a in group)
            right = next(group for group in components if b in group)
            if left is right:
                marks.append(a)
                excess += 1
            else:
                left.update(right)
                components.remove(right)
    if all(sum(i in group for i in marks) == 1 for group in components):
        return (2 - 2 * genus)**excess
    return 0


total = 1
for point, divisor in zip(points, (0, diagonals[0], diagonals[1] + diagonals[2])):
    total *= (1 + (degree * point - 2 * divisor) * parameter
              + (divisor**2 - degree * point * divisor) * parameter**2)
third = Poly(expand(total).coeff(parameter, 3), *variables)
integral = sum(coefficient * diagonal_integral(powers) for powers, coefficient in third.terms()) / 6
formula = (degree - 4) * ((degree - 2) * (degree - 3) - 6 * genus) / 6
assert factor(integral - formula) == 0
assert factor(integral.subs(degree, 8) - (20 - 4 * genus)) == 0
print("PASS: ordered-cube diagonal intersections give 20 - 4g for a degree-eight ruling")
