"""Polarized pencil conic and all prescribed brackets on a defect chart."""
import sympy as sp


def comm(a, b):
    return a*b-b*a


def anti(a, b):
    return a*b+b*a


def zero_trace(a):
    return a-a.trace()*sp.eye(3)/3


def inverse_anti(k, d):
    return ((k*k*d).trace()*sp.eye(3)+(k*d).trace()*k-2*k*d*k)/(2*k.det())


def check():
    # A nondiagonalizable invertible K exercises the polynomial extension.
    k = sp.Matrix([[1, 1, 0], [0, 1, 0], [0, 0, -2]])
    values = sp.symbols("t:8")
    t = sp.Matrix(3, 3, [*values, -values[0]-values[4]])
    assert (anti(k, inverse_anti(k, t))-t).applyfunc(sp.expand) == sp.zeros(3)
    f = (k*t).trace()**2-2*((k*t)**2).trace()
    polar = (k*t).trace()**2-(k*k).trace()*(t*t).trace()-2*(comm(k, t)**2).trace()
    assert sp.expand(3*f-polar) == 0

    x = sp.diag(1, 2, -3)
    y = sp.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    z = zero_trace(anti(x, y))
    k, d, ell = comm(x, y), comm(x, z), comm(y, z)
    assert k.det() != 0
    assert sp.Matrix.hstack(*(v.reshape(9, 1) for v in (k, d, ell))).rank() == 3
    u, v = inverse_anti(k, d), inverse_anti(k, ell)
    assert -comm(u, v) == k
    h = sp.Symbol("h")
    for beta in (1, -1):
        a = beta*zero_trace(u)
        b = -beta*zero_trace(v)+h*a
        c = zero_trace(anti(a, b)/beta+2*u.trace()*b/3-2*v.trace()*a/3-h*a*a/beta)
        for found, target in zip([comm(a, b), comm(a, c), comm(b, c)], [k, d, ell]):
            assert (found-target).applyfunc(sp.expand) == sp.zeros(3)


if __name__ == "__main__":
    check()
    print("PASS: Sylvester inverse at a Jordan block, pencil conic polarization, and both defect families")
