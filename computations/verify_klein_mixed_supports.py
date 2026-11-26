"""Exact support quotients, low-degree obstructions and regular Walsh limits."""
from collections import defaultdict
from itertools import permutations, product
from operator import xor
from functools import reduce

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

AXES = range(3)
LABELS = range(4)
QUADRUPLES = list(product(LABELS, repeat=4))
MASKS = [v for v in product((0, 1), repeat=4) if v[0]*v[3] or v[1]*v[2]]
# A six marks q- and p+, and a nine marks q+ and p-.
EXPECTED = [
    ((6,6,6),(6,6,6)), ((6,6,6),(9,9,9)), ((6,6,9),(6,6,9)),
    ((6,6,9),(6,6,13)), ((6,6,9),(6,6,15)), ((6,6,13),(6,6,13)),
    ((6,6,13),(6,6,15)), ((6,6,15),(6,6,15)), ((6,14,9),(6,14,9)),
    ((6,14,9),(6,14,13)), ((6,14,13),(6,14,9)), ((6,14,13),(6,14,13)),
    ((6,9,15),(6,9,15)), ((6,13,13),(6,13,11)), ((6,13,11),(14,13,13)),
    ((6,13,15),(6,13,15)), ((6,15,15),(6,15,15)), ((14,13,15),(14,13,15)),
    ((15,15,15),(15,15,15)),
]


def off_diagonals():
    for s, i, step in product((0, 1), AXES, (1, -1)):
        j = (i+step) % 3
        lhs = [(1,i,s), (1,j,s+2)]
        if s == 0 and step == 1:
            rhs = [(0,j,1), (0,i,3)]
        elif s == 0:
            rhs = [(0,i,0), (0,j,2)]
        elif step == 1:
            rhs = [(0,i,1), (0,j,3)]
        else:
            rhs = [(0,j,0), (0,i,2)]
        yield lhs, rhs


def support_orbits():
    by_image = defaultdict(list)
    def flags(side, terms):
        return tuple(side[i][s]*side[j][t] for (_,i,s), (_,j,t) in terms)
    left, right = zip(*off_diagonals())
    for side in product(MASKS, repeat=3):
        by_image[flags(side, left)].append(side)
    pairs = {(q,g) for q in product(MASKS, repeat=3) for g in by_image[flags(q, right)]}
    assert len(pairs) == 293
    representatives = []
    while pairs:
        pair = min(pairs)
        representatives.append(pair)
        orbit = set()
        for _ in LABELS:
            for pi in permutations(AXES):
                parity = sum(a > b for a,b in ((pi[0],pi[1]),(pi[0],pi[2]),(pi[1],pi[2]))) % 2
                orbit.add(tuple(tuple(tuple(side[i][s^parity] for s in LABELS) for i in pi) for side in pair))
            q,g = pair
            pair = (tuple((c,d,a,b) for a,b,c,d in g), tuple((b,a,d,c) for a,b,c,d in q))
        assert orbit <= pairs
        pairs -= orbit
    codes = [tuple(tuple(sum(bit << j for j, bit in enumerate(v)) for v in side) for side in pair)
             for pair in representatives]
    assert codes == EXPECTED
    return representatives


def character(row, values):
    return sp.factor(sp.prod(v**int(e) for v,e in zip(values, row)))


def smith(matrix):
    diagonal, left, right = [m.to_Matrix() for m in smith_normal_decomp(DomainMatrix.from_Matrix(matrix))]
    assert left*matrix*right == diagonal
    return diagonal, left, right


def solve_support(pair, number):
    atoms = [(h,i,s) for h,i,s in product(range(2), AXES, LABELS) if pair[h][i][s]]
    rows, targets, products = [], [], []
    def insert(positive, negative, target):
        row = [0]*len(atoms)
        for atom in positive:
            row[atoms.index(atom)] += 1
        for atom in negative:
            row[atoms.index(atom)] -= 1
        rows.append(row)
        targets.append(sp.sympify(target))
    for h,side in enumerate(pair):
        beta = []
        for i, (a,b,c,d) in enumerate(side):
            v = sp.Symbol(f'b{h}{i}') if a*b*c*d else sp.Integer(0 if a*d else -1)
            beta.append(v)
            if a*d:
                insert([(h,i,0),(h,i,3)], [], 1+v)
            if b*c:
                insert([(h,i,1),(h,i,2)], [], v)
        products.append(beta)
    for lhs,rhs in off_diagonals():
        if all(x in atoms for x in lhs):
            insert(lhs, rhs, -1)
    matrix = sp.Matrix(rows)
    diagonal,left,right = smith(matrix)
    rank = matrix.rank()
    assert list(diagonal.diagonal()[:rank]) == [1]*rank
    equations = [sp.fraction(sp.factor(character(row, targets)-1))[0] for row in left.tolist()[rank:]]
    beta,gamma = products
    variables = sorted(set().union(*(v.free_symbols for v in beta+gamma)), key=str)
    # Compare ideals on the open set where every listed coefficient is nonzero.
    localization = sp.Symbol('inv')
    denominator = sp.prod(v*(v+1) for v in variables)
    def ideal(gens):
        return list(sp.groebner([*gens, localization*denominator-1], localization, *variables))
    if number in (5,7):
        assert ideal(equations) == [1]
        return None
    equalities = [] if number == 2 else [v-w for v,w in zip(beta,gamma)]
    assert ideal(equations) == ideal(equalities)
    substitutions = {w:v for v,w in zip(beta,gamma) if w.is_Symbol}
    targets = [v.subs(substitutions) for v in targets]
    gauge = sp.zeros(len(atoms),6)
    for k,(h,i,s) in enumerate(atoms):
        sign = (-1)**(s%2)
        gauge[k,h*3+(i+1)%3],gauge[k,h*3+(i+2)%3] = sign,-sign
        gauge[k,(1-h)*3+i] = 1 if s<2 else -1
    assert matrix*gauge == sp.zeros(matrix.rows,6)
    # In kernel Smith coordinates, quotient the diagonal reparametrizations.
    gauge_coords = (right.inv()*gauge)[rank:,:]
    gd, gl, _ = smith(gauge_coords)
    gauge_rank = sum(entry != 0 for entry in gd.diagonal())
    dim = len(atoms)-rank-gauge_rank
    assert dim == int(number in (12,15,18))
    free = [sp.Integer(1)]*gauge_rank + [sp.Symbol('T')]*dim
    kernel_point = [character(row, free) for row in gl.inv().tolist()]
    base = [character(row,targets) for row in left.tolist()[:rank]] + kernel_point
    return {x:character(row,base) for x,row in zip(atoms,right.tolist())}


def tensor(solution):
    nets = []
    for h,s in ((1,0),(1,2),(0,0),(0,2)):
        side = []
        for i in AXES:
            j,k = (i+1)%3,(i+2)%3
            m = sp.zeros(3)
            m[j,k],m[k,j] = (solution.get((h,i,t),0) for t in (s,s+1))
            side.append(m)
        nets.append(side)
    g,f,q,p = nets
    def entry(i,j,k,l):
        if i==k==0 and j==l and j or j==l==0 and i==k and i:
            return sp.Integer(1)
        if (i,j,k,l).count(0)>=2:
            return sp.Integer(0)
        if i*j*k*l:
            return sp.factor((g[i-1]*f[k-1])[j-1,l-1]+(q[j-1]*p[l-1])[i-1,k-1])
        if not i:
            return f[k-1][j-1,l-1]
        if not j:
            return p[l-1][i-1,k-1]
        if not k:
            return g[i-1][j-1,l-1]
        return q[j-1][i-1,k-1]
    return {key:value for key in QUADRUPLES if (value:=entry(*key)) != 0}


def dimension(w, start, length):
    rows_by_color = defaultdict(list)
    for pos in range(length-1):
        cut = (start+pos)%4
        for c in (1,2,3):
            for end in ((0,c),(c,0)):
                relation = {}
                for letters in product(LABELS, repeat=2):
                    key = letters+end
                    key = key[-cut:]+key[:-cut] if cut else key
                    if key in w:
                        relation[letters] = w[key]
                for a in product(LABELS, repeat=pos):
                    for b in product(LABELS, repeat=length-pos-2):
                        color = reduce(xor,a+b,c)
                        rows_by_color[color].append({a+k+b:v for k,v in relation.items()})
    dim = 4**length
    for color,rows in rows_by_color.items():
        columns = [v for v in product(LABELS, repeat=length) if reduce(xor,v)==color]
        m = sp.Matrix([[row.get(v,0) for v in columns] for row in rows])
        dim -= DomainMatrix.from_Matrix(m).to_field().rank()
    return dim


def walsh(w, number):
    rows, values = [], []
    factors = set()
    for (i,j,k,l),v in w.items():
        b,e = i^j,i^l
        assert k == i^b^e and b and e
        assert all((z,z^b,z^b^e,z^e) in w for z in LABELS)
        row = [0]*25
        for pos,label in enumerate((i,j,k,l)):
            row[4*pos+label] = 1
        row[16+3*(b-1)+e-1] = -1
        rows.append(row)
        values.append((-1)**((i>>1)*((b^e)&1))/v)
        factors.update(f for f,_ in sp.factor_list(sp.fraction(v)[0])[1])
    matrix = sp.Matrix(rows)
    d,left,_ = smith(matrix)
    rank = sum(entry != 0 for entry in d.diagonal())
    assert all(character(row,values) == 1 for row in left.tolist()[rank:])
    assert (0,3,1,2) not in w  # missing (b,e)=(3,2)
    z,y = sp.symbols('b01 b02')
    allowed = {z,y,z+1,y+1}
    if number == 17:
        assert factors-allowed == {z+y+1}
        assert {k for k,v in w.items() if sp.rem(sp.fraction(v)[0],z+y+1,z)==0} == {(2,3,2,3),(3,2,3,2)}
    else:
        assert factors <= allowed


def main():
    failing = {3:(0,3,21),4:(0,3,21),6:(0,3,21),9:(0,3,21),
               10:(0,3,21),11:(0,3,21),12:(0,3,21),13:(0,4,38),
               14:(0,3,21),15:(3,3,21),16:(0,4,38),18:(0,4,38)}
    for n,pair in enumerate(support_orbits()[:-1],1):
        solution = solve_support(pair,n)
        if solution is None:
            continue
        w = tensor(solution)
        if n in failing:
            start,length,expected = failing[n]
            assert dimension(w,start,length) == expected
        elif n in (1,8,17):
            walsh(w,n)
        else:
            assert n == 2 and dimension(w,0,5) == 60
    print('PASS: 293 supports, 19 orbits, saturated character ideals, exact Hilbert obstructions and Walsh certificates')


if __name__ == '__main__':
    main()
