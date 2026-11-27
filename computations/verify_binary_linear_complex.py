"""Single-mark reconstruction, generic unmarked uniqueness, and the full opposite tangent."""
from itertools import combinations, product

import numpy as np
import sympy as sp

from verify_binary_component_crossing import PRIME, mod_kernel, reduced
from verify_binary_incidence_sextic import convolution
from verify_binary_opposite_determinant import matrix, multiply, opposite, parity
from verify_binary_secant_obstructions import rank


def main():
    lift = sp.Matrix(6,4,[
        sp.Rational(808,243),sp.Rational(22111,729),sp.Rational(5344,729),sp.Rational(10933,729),
        sp.Rational(-4613,243),-1,-9,5,sp.Rational(-1543,81),-10,6,0,
        sp.Rational(-50,27),-4,4,5,sp.Rational(23,9),7,1,9,
        sp.Rational(-5,3),-6,7,-2])
    pairs = list(combinations(range(4),2))
    plucker = sp.zeros(9,6)
    for k,(i,j) in enumerate(pairs):
        for h in range(6):
            plucker[h+i,k] += lift[h,j]
            plucker[h+j,k] -= lift[h,i]
    omega = sp.Matrix([1,2,-3,1,2,1])
    assert plucker*omega == sp.zeros(9,1) and rank(plucker)==5
    assert omega[0]*omega[5]-omega[1]*omega[4]+omega[2]*omega[3] != 0
    octavics = sp.Matrix.hstack(*plucker.T.nullspace())
    hankels = [sp.Matrix(6,4,lambda i,j:octavics[i+j,h]) for h in range(4)]
    maps = [sp.Matrix(6,4,lambda i,j:int(i==j+k)) for k in range(3)] + [lift]
    w = sp.Matrix(16,16,lambda u,v:(maps[u%4][:,u//4].T*hankels[v%4][:,v//4])[0])
    assert matrix(opposite(reduced(w))).rank() == 10

    # The six-plane projective rank-two locus has only one rational quintic.
    coords = sp.symbols("a b c")
    lift_mod = sp.Matrix(reduced(lift).tolist())
    form = convolution(3,coords).col_join(convolution(1,coords)*lift_mod.T).det(method='domain-ge')
    gradient = [sp.diff(form,z) for z in coords]
    for i,z in enumerate(coords):
        patch = {q:0 for q in coords[:i]}
        patch[z] = 1
        equations = [q.subs(patch) for q in gradient]
        if i<2:
            assert list(sp.groebner(equations,*coords[i+1:],modulus=PRIME,order="grevlex")) == [1]
        else:
            assert any(int(c)%PRIME for c in equations)

    # The normal to the one-symmetry divisor lifts all six opposite kernel rows.
    base = reduced(w)
    left,right = mod_kernel(base.T),mod_kernel(base)
    other_kernel = mod_kernel(opposite(base))
    basis = []
    for i,j,k,l in product(range(4),repeat=4):
        if i>=k:
            continue
        entry = np.zeros((16,16),dtype=int)
        entry[4*i+j,4*k+l],entry[4*k+j,4*i+l] = 1,-1
        basis.append(entry)
    basis = np.array(basis)
    jac = multiply(left.T,basis,right).reshape(len(basis),100).T
    directions = mod_kernel(jac).T@basis.reshape(-1,256)
    assert len(directions) == 16
    ranks = [matrix(multiply(parity(-1).T,opposite(direction),other_kernel)).rank()
             for direction in directions.reshape(-1,16,16)%PRIME]
    assert max(ranks) == 6
    print("Nondegenerate linear complex, unique generic binary mark, opposite rank ten and six independent normal rows verified.")


if __name__ == '__main__':
    main()
