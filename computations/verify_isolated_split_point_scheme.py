"""Local line pencils and embedded points over the exceptional number field."""
from itertools import product

import sympy as sp
from sympy.polys.agca.extensions import MonogenicFiniteExtension
from sympy.polys.matrices import DomainMatrix

from verify_isolated_split_component import exceptional_maps, polynomials, spectrum_at
from verify_split_inverse_windows import BITS, PAIRS, bell, borders, positions


def tensor_over_field(field, d):
    p = {q: field.convert(value) for q, value in bell(spectrum_at(d)).items()}
    maps = exceptional_maps()
    targets = []
    for i in (0, 1):
        cut = DomainMatrix.from_list([[p[positions(j,k,i)] for k in range(4)] for j in range(4)], field)
        inv = cut.inv().transpose().to_list()
        dual = {positions(j,k,i):inv[j][k] for j,k in product(range(4), repeat=2)}
        targets.append({q:sum((value*field.convert(sp.prod(maps[i][j][q[PAIRS[i][j]],a[j]]
                                                           for j in range(4)))
                               for a,value in dual.items()), field.zero) for q in BITS})
    pivot = next(q for q in BITS if targets[0][q])
    scale = targets[0][pivot]/targets[1][pivot]
    assert all(targets[0][q] == scale*targets[1][q] for q in BITS)
    colors = (0,1,0,1)
    def word(q, pattern):
        return tuple(j+2*c for j,c in zip(q,pattern))
    w = {word(q, colors):value for q,value in p.items()}
    w.update({word(q, [1-c for c in colors]):value for q,value in targets[0].items()})
    for i,q in product(range(2),BITS):
        values = borders(maps[i],PAIRS[i],i,q)
        for k in (0,1):
            pattern = [c ^ (int(j not in {i,i+1}) ^ k) for j,c in enumerate(colors)]
            w[word(q,pattern)] = field.convert(values[k])*(scale if i == 1 and k == 0 else field.one)
    return w


def local_pencil(field, normal_relation, tangent_relation, nt, tn):
    ring = field.poly_ring('y')
    y = ring.gens[0]
    constant = ring.ring.ground_new
    skew = DomainMatrix.from_list([[0,-1],[1,0]],field)
    rows = []
    for first, second in zip(nt, tn):
        # n=(1,y). The unknowns are the two tangential target then source coordinates.
        matrices = [first,(second*skew*normal_relation.transpose()).transpose()]
        rows.append([constant(m[0][j])+y*constant(m[1][j])
                     for a in matrices for m in [a.to_list()] for j in (0,1)])
    matrix = DomainMatrix.from_list(rows,ring)
    det = matrix.det()
    assert det.degree() == 4 and det.gcd(det.diff(y)).degree() == 0
    adj = matrix.adjugate().to_list()
    coeff = tangent_relation.to_list()
    common = det
    for j in range(4):
        value = sum((adj[2+k][j]*constant(coeff[k][l])*adj[l][j]
                     for k,l in product(range(2),repeat=2)),ring.zero)
        common = common.gcd(value)
    assert common.degree() == 3


def verify():
    d, modulus = polynomials()
    polynomial = sp.Poly(modulus,d,domain=sp.QQ_I)
    assert polynomial.is_irreducible
    field = MonogenicFiniteExtension(polynomial)
    # Polynomial Euclidean division requires quotient in the residue field.
    field.quo = field.exquo = lambda a,b: a/b
    tensor = tensor_over_field(field,d)
    def column(i,j,k,l):
        return DomainMatrix.from_list([[tensor.get((a,b,i,j),field.zero)
                                        for b in range(l,l+2)] for a in range(k,k+2)],field)
    def pure_relation(k,l):
        candidates = [column(i,j,k,k) for i,j in product(range(l,l+2), repeat=2)]
        return next(m for m in candidates if not m.is_zero_matrix)
    nn, tt = pure_relation(0,2), pure_relation(2,0)
    pairs = list(product(range(2),range(2,4)))
    nt, tn = ([column(i,j,k,l) for i,j in pairs] for k,l in ((0,2),(2,0)))
    local_pencil(field,nn,tt,nt,tn)
    local_pencil(field,tt,nn,tn,nt)
    print("Both quartics are squarefree with cubic transverse factors over the defining field.")


if __name__ == '__main__':
    verify()
