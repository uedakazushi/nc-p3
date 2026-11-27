"""The five fixed-locus components at the period-one no-flip surface collision."""
from itertools import product

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from verify_split_surface_lifts import U, V, WORDS, split_tensor

I = sp.I
T = 1+sp.sqrt(2)
FIELD = sp.QQ.algebraic_field(I, sp.sqrt(2))


def scalar(value):
    return FIELD.from_sympy(sp.sympify(value))


def matrix(rows):
    return DomainMatrix(rows, (len(rows), len(rows[0])), FIELD)


def phase_product(values):
    answer = FIELD.one
    for value in values:
        answer *= value
    return answer


PHASE = tuple(map(scalar, [I,1,I*T**2,-T**2, -I,-1,I,1, -I/T**2,-1/T**2,-I,1, I,1,I,1]))


def shift(q):
    rotated = q[1:]+q[:1]
    factor = phase_product(PHASE[4*j+a] for j,a in enumerate(q))*scalar(T)**(sum(a<2 for a in q)-2)
    return tuple(a ^ int(a>=2) for a in rotated), factor


def pauli(q, flip, sign):
    n = sum(a<2 for a in q)
    factor = FIELD.one*(-1)**(sum(q)*sign)
    if flip:
        unit = scalar(I)
        factor *= unit**(n-2)*phase_product((unit,unit,-unit,-unit)[j] for j,a in enumerate(q) if a<2)
    return tuple(a^flip for a in q), factor


def representations():
    w = {q:scalar(sp.sympify(a).subs({U:T,V:T})) for q,a in split_tensor(sp.eye(2)).items() if a}
    for q,a in w.items():
        target, phase = shift(q)
        assert phase*a == w[target]
    pairs = tuple(product(range(4), repeat=2))
    cuts = []
    for i in (0, 1):
        def value(q):
            return w.get(q[-i:]+q[:-i] if i else q,FIELD.zero)
        m = matrix([[value(a+b) for b in pairs] for a in pairs])
        cuts.append((m.transpose().nullspace().to_list(),m.nullspace().to_list()))
    result = []
    for n in (3, 4):
        words = [q for q in WORDS if sum(j<2 for j in q)==n]
        pos = {q:j for j,q in enumerate(words)}
        rows = []
        for i,(left,right) in enumerate(cuts):
            for a,b in product(left,right):
                row = [a[4*q[i]+q[(i+1)%4]]*b[4*q[(i+2)%4]+q[(i+3)%4]] for q in words]
                if any(row):
                    rows.append(row)
        frame = matrix(rows).nullspace().transpose()
        d = 0
        if n == 3:
            shears = matrix([[w.get(q[:j]+(v,)+q[j+1:],FIELD.zero) if q[j]==a else FIELD.zero
                              for j,a,v in product(range(4),range(2),range(2,4))] for q in words])
            frame = shears.hstack(frame)
            frame = frame.extract(range(len(words)),frame.rref()[1])
            d = 16
        rs = frame.transpose().rref()[1]
        lift = frame.extract(rs,range(frame.shape[1])).inv()

        def action(fn):
            basis = frame.extract(range(len(words)),range(d,frame.shape[1])).to_list()
            out = [None]*len(words)
            for q,vector in zip(words,basis):
                target,c = fn(q)
                out[pos[target]] = [c*x for x in vector]
            out = matrix(out)
            coords = lift*out.extract(rs,range(out.shape[1]))
            assert frame*coords == out
            return coords.extract(range(d,frame.shape[1]),range(out.shape[1])).to_Matrix().applyfunc(sp.simplify)
        s,x,z = action(shift), action(lambda q:pauli(q,1,0)), action(lambda q:pauli(q,0,1))
        assert sp.simplify(s**4) == sp.eye(s.rows)
        result.append((s,x,z))
    return result


def eig(matrix, value):
    vectors = (matrix-value*sp.eye(matrix.rows)).nullspace()
    return sp.Matrix.hstack(*vectors) if vectors else sp.zeros(matrix.rows,0)


def verify():
    (sb,xb,zb),(sc,xc,zc) = representations()
    group = [(sp.eye(4),sp.eye(3)), (zb,zc), (xb,xc), (xb*zb,xc*zc)]
    charts = []
    for j,(b,c) in enumerate(group):
        for phase in (1,-1,I,-I):
            mb,mc = eig(b*sb,phase),eig(c*sc,phase**2)
            assert mb.cols == 1
            charts.append((j,str(phase),mb,mc))
            pure = eig(c*sc,phase)
            if pure.cols:
                charts.append((j,'c'+str(phase),sp.zeros(4,0),pure))
    charts.sort(key=lambda chart: -(chart[2].cols+chart[3].cols))
    reps = []
    def contained(a,b):
        return b.row_join(a).rank() == b.cols
    for j,p,b,c in charts:
        if not any(contained(b,g*b0) and contained(c,h*c0)
                   for _,_,b0,c0 in reps for g,h in group):
            reps.append((j,p,b,c))
    assert [(j,p,b.cols,c.cols) for j,p,b,c in reps] == [
        (0,'1',1,2),(0,'I',1,1),(2,'1',1,1),(2,'I',1,0),(2,'cI',0,1)]


if __name__ == '__main__':
    verify()
    print('The collision shift, both Pauli actions and all five period-one components verified.')
