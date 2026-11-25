"""An unrestricted 76-coordinate etale chart at the regular residue window in characteristic 53."""

from itertools import product

from flint import nmod_mat

from verify_half_period_components import submatrix
from verify_walsh_full_component import P, matrices, point_certificate


def pivots(matrix):
    rows, rank = matrix.rref()
    return [next(j for j in range(rows.ncols()) if rows[i,j]) for i in range(rank)]


def cut(values, side, modulus):
    def index(row, col):
        i,j = divmod(row,4)
        k,l = divmod(col,4)
        if side:
            i,k = k,i
        return 64*i+16*j+4*k+l
    entries = [[values[index(r,c)] for c in range(16)] for r in range(16)]
    return nmod_mat(entries,modulus)


def extract(matrix, rows, cols, modulus):
    return nmod_mat([[int(matrix[i,j]) for j in cols] for i in rows],modulus)


def schur(matrix, frame, modulus):
    rows,cols,rest_rows,rest_cols,base_inverse = frame
    core = extract(matrix,rows,cols,modulus)
    # Newton inversion of the unit minor uses field inversion only at the residue.
    inverse = nmod_mat(6,6,[int(x) for x in base_inverse.entries()],modulus)
    eye = nmod_mat([[int(i == j) for j in range(6)] for i in range(6)],modulus)
    accuracy = P
    while accuracy < modulus:
        inverse = inverse*(2*eye-core*inverse)
        accuracy *= accuracy
    assert core*inverse == eye
    lower = extract(matrix,rest_rows,cols,modulus)
    upper = extract(matrix,rows,rest_cols,modulus)
    return extract(matrix,rest_rows,rest_cols,modulus)-lower*inverse*upper


def verify_chart():
    first,second = matrices()
    point_certificate(first,second)
    values = [int(x) for x in first.entries()]
    frames, gradients = [],[]
    for side, matrix in enumerate((first,second)):
        rows,cols = pivots(matrix.transpose()),pivots(matrix)
        rr = [r for r in range(16) if r not in rows]
        cc = [c for c in range(16) if c not in cols]
        inv = submatrix(matrix,rows,cols).inv()
        frames.append((rows,cols,rr,cc,inv))
        # dSchur = L dM R with each outer null frame in the original word order.
        lower = -submatrix(matrix,rr,cols)*inv
        upper = -inv*submatrix(matrix,rows,cc)
        left, right = nmod_mat(10,16,P), nmod_mat(16,10,P)
        for i in range(10):
            left[i,rr[i]],right[cc[i],i] = 1,1
            for j in range(6):
                left[i,rows[j]],right[cols[j],i] = lower[i,j],upper[j,i]
        for a,b in product(range(10),repeat=2):
            gradients.append([left[a,4*(k if side else i)+j]*right[4*(i if side else k)+l,b]
                              for i,j,k,l in product(range(4),repeat=4)])
    jac = nmod_mat(gradients,P)
    equations,unknowns = pivots(jac.transpose()),pivots(jac)
    assert len(equations) == len(unknowns) == 180
    inverse = submatrix(jac,equations,unknowns).inv()
    free = [i for i in range(256) if i not in unknowns]
    assert len(free) == 76
    signed = nmod_mat([[(-1)**(i.bit_count()+j.bit_count())*first[i,4*(i^j)+j]**2
                        for j in range(3)] for i in range(3)],P)
    # first[0*4+i, (i xor j)*4+j] is the coefficient m_ij.
    assert signed.rank() == 2
    assert all(sum(signed[i,j] for j in range(3)) == 0 for i in range(3))
    assert all(sum(signed[i,j] for i in range(3)) == 0 for j in range(3))

    def residual(modulus):
        return [int(x) for side,frame in enumerate(frames)
                for x in schur(cut(values,side,modulus),frame,modulus).entries()]
    assert not any(residual(P))
    for power in (1,2):
        accuracy = P**power
        modulus = P*accuracy
        for i,word in enumerate(free):
            values[word] += accuracy*((i+1)**power % P)
        error = residual(modulus)
        assert all(x % accuracy == 0 for x in error)
        correction = -inverse*nmod_mat([[error[i]//accuracy] for i in equations],P)
        for i,word in enumerate(unknowns):
            values[word] += accuracy*int(correction[i,0])
        assert not any(residual(modulus))
    return len(free)


if __name__ == "__main__":
    print(f"PASS: {verify_chart()} free coefficients; both full Schur complements vanish modulo 53^3")
