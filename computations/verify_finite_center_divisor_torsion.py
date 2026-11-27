"""Integral augmentation lattice for the four periodic vertex determinants."""
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form


def main():
    relations = Matrix([[4,6,4], [0,4,6], [0,0,4]])
    invariants = smith_normal_form(relations, domain=ZZ)
    assert sorted(abs(invariants[i,i]) for i in range(3)) == [2,2,16]
    assert abs(relations.det()) == 64
    assert (16*relations.inv()).applyfunc(lambda x: int(x)) == 16*relations.inv()
    assert any(x.q != 1 for x in 8*relations.inv())
    print("The determinant differences have universal lattice invariants 2, 2, 16.")


if __name__ == "__main__":
    main()
