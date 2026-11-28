"""Check numerical splitting bounds and the connecting-map description of the second format."""
from itertools import combinations

import numpy as np

from verify_degree_twenty_hankel_pairs import field, forms, hankel, nullrows, quartic_basis
from verify_quadric_biform_window import PRIME


def main():
    allowed = []
    for small in range(1,5):
        for a in (2,3,4):
            for b in range(-20,21):
                quotient_fiber = 4-a
                quotient_min_section = quotient_fiber*small-b
                if quotient_min_section < 0:
                    continue
                defect = 20-(8*a*(4-a)+(4-2*a)*b)
                if defect < 0:
                    continue
                sections = sum(max(0,1+j*(8-small)+(quotient_fiber-j)*small-b)
                               for j in range(quotient_fiber+1))
                if sections >= 6:
                    allowed.append((small,a,b,defect))
    assert {(a,b,length) for _,a,b,length in allowed} == {(3,2,0),(4,-5,0)}
    v = forms()
    transfer = quartic_basis(v,1,0).inv()*quartic_basis(v,0,1)
    t = np.array(transfer.tolist(),dtype=int)
    equations = []
    brackets = []
    for basis_vector in np.eye(9,dtype=int):
        form = hankel(basis_vector)
        transformed = t.T @ form @ t % PRIME
        equations.append((transformed[1:,:-1]-transformed[:-1,1:]).flatten())
        brackets.append([transformed[i,j+1]-transformed[i+1,j]
                         for i,j in combinations(range(4),2)])
    assert field(np.array(equations).T).rank() == field(np.array(brackets).T).rank() == 5
    for eta in nullrows(field(np.array(equations).T)):
        g = hankel(eta)
        psi = (t.T @ g @ t) % PRIME
        moments = [int(psi[min(i,4),max(i-4,0)]) for i in range(9)]
        for i in range(4):
            row = np.array(moments[i:i+6])
            assert np.array_equal((t[:,i+1].T @ g @ t) % PRIME, row[1:])
            assert np.array_equal((t[:,i].T @ g @ t) % PRIME, row[:-1])
    print("PASS: the two splitting types, zero residual lengths, and the common Hankel column formula")


if __name__ == "__main__":
    main()
