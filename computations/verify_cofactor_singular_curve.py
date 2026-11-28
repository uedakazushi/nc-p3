"""Check the Hilbert functions, punctual lengths and incidence budgets in the proof."""
from math import comb
from sympy import Matrix, expand, symbols

directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1),
              (1, 2, 3), (1, 3, 2)]
for count in range(2, 7):
    ranks = []
    for degree in range(5):
        rows = [[x**i * y**j * z**(degree-i-j)
                 for i in range(degree+1) for j in range(degree+1-i)]
                for x, y, z in directions[:count]]
        ranks.append(Matrix(rows).rank())
    assert ranks == [1, min(3, count), count, count, count]
    delta = sum(count-rank for rank in ranks)
    delta_w = sum(count-rank for rank in ranks[1:])
    assert delta_w == delta-(count-1)
    local_length = 3-count if count < 4 else comb(6, 3)-sum(ranks[:4])
    assert 4-local_length+delta+delta_w <= 2*count
s, t = symbols("s t")
quad = Matrix([s*s, s*t, t*t])
relations = Matrix([[-t,s,0],[0,-t,s]])
assert relations*quad == Matrix([0,0])
assert relations.row(0).cross(relations.row(1)).T == quad
m, marked_on_curve, total_rank_one, beta, d, chi = symbols("m m_C s beta d chi")
assert expand(4*(20-marked_on_curve)+2*total_rank_one+beta+6*d+2*chi-80) == expand(
    6*d+2*chi+beta+2*(total_rank_one-marked_on_curve)-2*marked_on_curve)
print("PASS: direction evaluations, exterior-square defects and marked incidence budgets")
