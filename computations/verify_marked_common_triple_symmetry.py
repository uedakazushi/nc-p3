"""Check the trace-form radical and the dimension of symmetric star spaces."""
from sympy import Matrix, eye, symbols


def cross(v):
    a, b, c = v
    return Matrix([[0, -c, b], [c, 0, -a], [-b, a, 0]])


def symmetric_basis():
    basis = [Matrix.diag(1, 0, -1), Matrix.diag(0, 1, -1)]
    for i, j in ((0, 1), (0, 2), (1, 2)):
        item = Matrix.zeros(3)
        item[i, j] = item[j, i] = 1
        basis.append(item)
    return basis


s = symmetric_basis()
for n in (Matrix([0, 0, 1]), Matrix([1, 2, 3])):
    k = cross(n)
    form = Matrix(5, 5, lambda i, j: (k*(s[i]*s[j]-s[j]*s[i])).trace())
    radical = k*k - (k*k).trace()*eye(3)/3
    assert form.rank() == 4 and k*radical == radical*k

variables = symbols("a:6")
a,b,c,d,e,f = variables
m = Matrix([[a,b,c], [b,d,e], [c,e,f]])
# alpha is the first coordinate, beta first or second.
for beta, expected in ((Matrix([1,0,0]), 3), (Matrix([0,1,0]), 2)):
    conditions = list(beta.cross(m[:,1])) + list(beta.cross(m[:,2]))
    assert 6 - Matrix(conditions).jacobian(variables).rank() == expected
print("PASS: symplectic quotient dimension four and both symmetric-image dimensions")
