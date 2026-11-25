"""Stabilizers, normalization branches and rank-two pencils in the order-four half-turn chart."""

from itertools import permutations, product

from flint import nmod_mat

from verify_walsh_full_component import K, P, matrices


def sign(x, character):
    return (-1) ** ((x & character).bit_count())


def tensor_product(a, b):
    return nmod_mat([[int(a[i,k]*b[j,l]) for k,l in product(K,repeat=2)]
                     for i,j in product(K,repeat=2)],P)


def submatrix(matrix, rows, cols):
    return nmod_mat([[int(matrix[i,j]) for j in cols] for i in rows],P)


def half_turn(x, even):
    return nmod_mat([[int(j ^ x == i)*(-1)**(int(j == 3) + int(even)*(x >> 1)*(j & 1))
                      for j in K] for i in K],P)


def determinant_frame(matrix):
    rows,n = matrix.rref()
    assert n == 2
    roots = []
    for k in (*range(P), None):
        entries = [rows[1,j] if k is None else rows[0,j]+k*rows[1,j] for j in K]
        point = nmod_mat(2,2,entries,P)
        if not point.det():
            roots.append(point)
    assert len(roots) == 2
    us,vs = [],[]
    for point in roots:
        i,j = next((i,j) for i,j in product(range(2),repeat=2) if point[i,j])
        us.append([int(point[k,j]) for k in range(2)])
        vs.append([int(point[i,k]) for k in range(2)])
    return nmod_mat(us,P).transpose().inv(), nmod_mat(vs,P).transpose().inv()


def pencil_chart(x):
    first,_ = matrices()
    a,b = half_turn(x,True),half_turn(x,False)
    g = tensor_product(a,b)
    assert first == g*first.transpose()*g.transpose()
    m = first*g.transpose()
    involutions = [a*a,b*b]
    spaces = [[[i for i in K if int(d[i,i]) == (-1)**parity % P] for parity in range(2)]
              for d in involutions]

    def indices(u,v):
        return [4*i+j for i,j in product(spaces[0][u],spaces[1][v])]

    changes = [nmod_mat(4,4,[0]*16,P) for _ in range(2)]
    for u in range(2):
        idx = indices(u,1-u)
        de,df = determinant_frame(submatrix(m,idx,idx))
        for change, labels, frame in zip(changes,(spaces[0][u],spaces[1][1-u]),(de,df),strict=True):
            for i,j in product(range(2),repeat=2):
                change[labels[i],labels[j]] = frame[i,j]
    move = tensor_product(*changes)
    m = move*m*move.transpose()
    n = nmod_mat([[m[4*k+j,4*i+l]*involutions[0][k,k] for k,l in product(K,repeat=2)]
                  for i,j in product(K,repeat=2)],P)
    for cut in (m,n):
        assert (cut+cut.transpose()).rank() == 4
        assert (cut-cut.transpose()).rank() == 2
    pos,neg = indices(0,0),indices(1,1)
    return submatrix(m,pos,pos),submatrix(m,neg,neg),submatrix(m,pos,neg),submatrix(n,pos,neg)


def characters():
    elements = list(product(K,repeat=2))
    sections = {(h,c) for h in (1,2,3) for c in K}
    def stabilizer(active):
        return {(x,z) for x,z in elements if all(sign(x,c)*sign(z,h) == 1 for h,c in active)}
    assert stabilizer(sections) == {(0,0)}
    for h in (1,2,3):
        flat = ((h & 1) << 1) + (h >> 1)
        assert stabilizer({(h,c) for c in K}) == {(0,0),(0,flat)}
        reps = [c for c in K if c < c ^ flat]
        signs = {(sign(x,reps[0])*sign(z,h),sign(x,reps[1])*sign(z,h))
                 for x,z in elements if x in (0,h)}
        assert len(signs) == 4
    assert len({tuple(p) for p in permutations((1,2,3))}) == 6


if __name__ == "__main__":
    characters()
    for x in (1,2,3):
        s,t,h,j = pencil_chart(x)
        assert all(a.det() for a in (s,t,h,j))
    print("PASS: three order-four lifts, transverse determinant pencils and generic stabilizers")
