"""Intrinsic diagonal endomorphisms on the mixed transverse curve, over Q(t)."""
import numpy as np
import sympy as sp

from verify_mixed_component_smoothing import along

t = sp.symbols('t')


def tensor_pair():
    base, plus, minus = [along(k) for k in (0, 1, -1)]
    w = base + (plus-minus)*t/2 + (plus+minus-2*base)*t*t/2
    m = sp.Matrix(w.transpose(0, 2, 1, 3).reshape(16, 16))
    # Total color makes the opposite flattening a sum of four square blocks.
    inv = sp.zeros(16)
    for h in range(4):
        ids = [4*i + (i ^ h) for i in range(4)]
        block = m.extract(ids, ids).inv().applyfunc(sp.factor)
        for i, row in enumerate(ids):
            for j, col in enumerate(ids):
                inv[row, col] = block[i, j]
    assert (m*inv-sp.eye(16)).applyfunc(sp.cancel).is_zero_matrix
    dual = np.array(inv.T.tolist(), dtype=object).reshape(4, 4, 4, 4).transpose(0, 2, 1, 3)
    return w, dual


def transfer(w, dual, i, j, h):
    order = (i, j, *[k for k in range(4) if k not in (i, j)])
    x, y = w.transpose(order), dual.transpose(order)
    block = sp.Matrix(4, 4, lambda a, b: sum((x[a, b]*y[a ^ h, b ^ h]).flat))
    flip = [a ^ h for a in range(4)]
    composition = block*block.extract(flip, flip).T
    return ((composition-16*sp.eye(4)/9)*144*t*t*(1+t*t)**2).applyfunc(sp.cancel)


def main():
    w, dual = tensor_pair()
    minors = []
    for i in (0, 1):
        neighbors = [(i+1) % 4, (i-1) % 4]
        for j in neighbors:
            diag = transfer(w, dual, i, j, 0)
            assert diag == sp.ones(4)*sp.expand(512*t*t*(1+t*t)**2)
        block = sp.Matrix.vstack(*(transfer(w, dual, i, j, 1) for j in neighbors))
        for rows in [(0, 1, 2, 3), (0, 1, 2, 4)]:
            minors.append(block.extract(rows, range(4)).det(method='domain-ge'))
    assert sp.Poly(sp.gcd_list(minors), t).monic().as_expr() == sp.expand((t*t-1)**4*(t*t+1)**4)

    labels = [(0, 1), (1, 0), (2, 3), (3, 2)]
    u, v = t+1, t-1
    expected = [sp.Matrix([[1, 0, u, -v], [0, 1, v, u]]),
                sp.Matrix([[1, 0, v, u], [0, 1, -u, v]])]
    for n in range(4):
        arr = w.transpose(tuple((j+n) % 4 for j in range(4)))
        rows = sp.Matrix([[arr[a, b, c, d] for a, b in labels] for c, d in labels[:2]])
        target = expected[n % 2].copy()
        if n >= 2:
            target[:, 2:] = -target[:, 2:]
        assert (rows-target).applyfunc(sp.expand).is_zero_matrix
    print("PASS: diagonal block, four polynomial minors and gcd, four edge relation spaces")


if __name__ == '__main__':
    main()
