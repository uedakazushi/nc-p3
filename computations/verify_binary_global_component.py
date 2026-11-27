"""A simple rational marked point with smooth ternary discriminant on the binary contact family."""
from itertools import combinations

import numpy as np
import sympy as sp

from verify_binary_component_crossing import PRIME, reduced
from verify_binary_opposite_determinant import matrix, pivots
from verify_twenty_mixed_bracket_invariants import marked_tensor


def main():
    L = sp.Matrix(6,4,[-150,-434,-206,-282,151,-86,-1,-236,198,205,-9,5,
                           58,63,-10,6,21,14,0,-4,12,16,4,5])
    q = sp.Matrix([1,-1,0,1])
    p = sp.Matrix([2,3,1])
    mult = sp.Matrix(6,3,lambda i,j:q[i-j] if 0<=i-j<4 else 0)
    assert L*q == mult*p
    maps = [sp.Matrix(6,4,lambda i,j:int(i==j+k)) for k in range(3)] + [L]
    # Fix the extra coordinate of U and the leading coordinate of the monic cubic.
    jac = mult.row_join((L-sum((p[k]*maps[k] for k in range(3)),sp.zeros(6,4)))[:,:3])
    assert jac.det() != 0
    plucker = sp.zeros(9,6)
    for k,(i,j) in enumerate(combinations(range(4),2)):
        for h in range(6):
            plucker[h+i,k] += L[h,j]
            plucker[h+j,k] -= L[h,i]
    assert plucker*sp.Matrix([1,2,-3,1,2,1]) == sp.zeros(9,1)
    S = sp.Matrix.hstack(*plucker.T.nullspace())
    hankels = [sp.Matrix(6,4,lambda i,j:S[i+j,k]) for k in range(4)]
    kernel = sp.Matrix.hstack(*(h*q for h in hankels)).nullspace()
    assert len(kernel) == 1
    points = [-p.col_join(sp.Matrix([-1])),q,kernel[0],q]
    frames = np.array([reduced(v).flatten() for v in points],dtype=object)
    w = sp.Matrix(16,16,lambda a,b:(maps[a%4][:,a//4].T*hankels[b%4][:,b//4])[0])
    rotated = reduced(w).reshape(4,4,4,4).transpose(1,2,3,0)
    t = np.asarray(marked_tensor(rotated.astype(object),frames)%PRIME,dtype=int)
    G, F = t[1:,1:,0,1:],t[0,1:,1:,1:].transpose(1,0,2)
    mixed = (G[:,None]@F[None,:]-F[None,:]@G[:,None])%PRIME
    entries = pivots(mixed.reshape(9,9))
    assert len(entries) == 3
    coords = sp.symbols("x y z")
    M = sp.Matrix(3,3,lambda i,j:sum(coords[k]*int(mixed.reshape(3,3,9)[i,j,v])
                                      for k,v in enumerate(entries)))
    cubic = sp.Poly(M.det(),*coords,modulus=PRIME).as_expr()
    for i,z in enumerate(coords):
        patch = {q:0 for q in coords[:i]}
        patch[z] = 1
        eqs = [sp.diff(cubic,v).subs(patch) for v in coords]
        if i<2:
            assert list(sp.groebner(eqs,*coords[i+1:],modulus=PRIME)) == [1]
        else:
            assert any(int(v)%PRIME for v in eqs)
    # Compute the full two-cut tangent as the sum of the two opposite-swap signs.
    from verify_binary_component_crossing import mod_kernel
    from verify_binary_opposite_determinant import multiply
    original = reduced(w)
    left,right = mod_kernel(original.T),mod_kernel(original)
    for sign,expected in ((1,60),(-1,16)):
        basis = []
        for i,j,k,l in np.ndindex(4,4,4,4):
            if i>k or (i==k and sign==-1):
                continue
            a = np.zeros((16,16),dtype=int)
            a[4*i+j,4*k+l],a[4*k+j,4*i+l] = 1,sign
            basis.append(a)
        normal = multiply(left.T,np.array(basis),right).reshape(len(basis),100)
        assert len(basis)-matrix(normal).rank() == expected
    print("Simple rational mark, smooth mixed plane cubic, and unrestricted tangent dimension 76 verified.")


if __name__ == '__main__':
    main()
