"""Recover the four canonical cross contractions and the repeated structure's extra map."""
from itertools import combinations, product

import numpy as np
from flint import nmod_mat

from verify_quadric_biform_window import PRIME, forms
from verify_quadric_cone_middle import binary_coordinates, cone_edges, field
from verify_quadric_ruling_identifiability import quartic_basis


def nullrows(mat):
    cols, count = mat.nullspace()
    return np.array([[int(cols[i,j]) for i in range(cols.nrows())] for j in range(count)],dtype=int)


def hankel(functional):
    return np.array([[functional[i+j] for j in range(5)] for i in range(5)],dtype=int)


def distinct():
    v = forms()
    transfer = quartic_basis(v,1,0).inv()*quartic_basis(v,0,1)
    t = np.array(transfer.tolist(),dtype=int)
    constraints = []
    for direction in np.eye(9,dtype=int):
        a = t.T @ hankel(direction) @ t % PRIME
        constraints.append((a[1:,:-1]-a[:-1,1:]).flatten())
    common = nullrows(field(np.array(constraints).T))
    assert len(common) == 4
    rows = []
    for i,j in product(range(2),repeat=2):
        rows.append(np.concatenate([(t[:,i:i+4].T @ hankel(g)[:,j:j+4]).flatten()
                                    for g in common]))
    assert field(np.array(rows)).rank() == 4


def repeated():
    _,edge = cone_edges()
    _,d = binary_coordinates(edge)
    brackets = np.zeros((9,6),dtype=int)
    for k,(i,j) in enumerate(combinations(range(4),2)):
        for h in range(6):
            brackets[h+i,k] += int(d[h,j])
            brackets[h+j,k] -= int(d[h,i])
    common = nullrows(field(brackets).transpose())
    eq = [[int(k == i)*g[j+h]-int(k == j)*g[i+h]
           for k in range(4) for h in range(6)]
          for g in common for i,j in combinations(range(4),2)]
    assert field(np.array(eq)).rank() == 20
    generators = np.array([[int(h == k+r) for k in range(4) for h in range(6)]
                           for r in range(3)]+[[int(d[h,k]) for k in range(4) for h in range(6)]])
    assert field(generators).rank() == 4
    assert (field(np.array(eq)) * field(generators).transpose()).rank() == 0
    print("PASS: four independent cross contractions; equal structures recover exactly the three multipliers and the weighted-cone operator")


if __name__ == "__main__":
    distinct()
    repeated()
