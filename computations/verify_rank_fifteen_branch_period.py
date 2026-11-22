"""Exclude period one on the entire ordinary elliptic ramification curve."""

from flint import fmpz_poly
from rank_fifteen_ramification import BranchElement, T, Y
from rank_fourteen_function_field import irreducible_factors
from verify_rank_fifteen_coordinate_collisions import fifth_contraction
from verify_rank_fifteen_covers import frame_on_branch
from verify_rank_fifteen_endpoint_orbits import multiply, trace, trace_blocks


def moments(blocks):
    """H preserves the xor c of the two letters on which it acts."""
    second = third = 0
    for c in range(4):
        matrix = [[blocks[x ^ y][x][y ^ c] for y in range(4)] for x in range(4)]
        square = multiply(matrix, matrix)
        second += trace(square)
        third += trace(multiply(square, matrix))
    return second, third


def verify():
    blocks, _ = trace_blocks(frame_on_branch(T, Y, 0))
    rotated = [
        [[block[u][i ^ h] for u in range(4)] for i in range(4)]
        for h, block in enumerate(blocks)
    ]
    second, third = moments(blocks)
    assert moments(rotated) == (second, third)
    d = fifth_contraction(blocks, BranchElement)
    e = fifth_contraction(rotated, BranchElement)
    eq2 = (second * (d**2 - e**2)).norm_numerator()
    eq3 = (third * (d**3 - e**3)).norm_numerator()
    assert eq2 and eq3
    common = irreducible_factors(eq2.gcd(eq3))
    t = fmpz_poly([0, 1])
    expected = (
        t,
        t - 1,
        t + 1,
        t * t + 1,
        t * t - 4 * t + 1,
        t * t + 4 * t + 1,
        t**4 + 10 * t * t + 1,
    )
    assert common == {tuple(map(int, g)) for g in expected}
    print(
        "The two cyclic period-one obstructions have only seven boundary factors in common."
    )


if __name__ == "__main__":
    verify()
