"""Character strata, attracting obstructions and normalized half-shift at reflected collisions."""
from itertools import product
from subprocess import run

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from verify_isolated_split_component import WORDS
from verify_split_repeated_surface import C, S, surface, tangent_blocks
from verify_split_surface_all_fibers import ELIMINATE

AXES = {(2,0),(-2,0),(0,2),(0,-2)}
DIAGONALS = {(a,b) for a,b in product((-2,2),repeat=2)}
MIXED = {(a,b) for a,b in product((-3,-1,1,3),repeat=2) if abs(a)+abs(b)==4}
LOCALS = {'s':S,'c':C}


def character(q):
    return q.count(0)-q.count(2), q.count(1)-q.count(3)


def orbit_rows(w, words):
    return [[w.get(q[:i]+(b,)+q[i+1:],0) if q[i]==a else 0
             for i,a,b in product(range(4),repeat=3)] for q in words]


def uniform_ranks():
    source = ELIMINATE.split('poly base=')[0].replace('ring r=0,(u,v,z),dp;','ring r=0,(s,c,z),dp;')
    source += 'poly base=s*c*(s^4-1);'
    script = [source]
    checks = 0
    def certify(rows, rank, locus):
        nonlocal checks
        mat = DomainMatrix.from_Matrix(sp.Matrix(rows))
        assert mat.to_field().rank() == rank
        if rank:
            vals = ','.join(str(x).replace('**','^') for row in rows for x in row)
            script.append(f'matrix a[{len(rows)}][{len(rows[0])}]={vals};'
                          f'rankBound(a,{rank},{locus});kill a;')
            checks += 1
    for exceptional in (False,True):
        w = surface(S,sp.Integer(1) if exceptional else C)
        for h, words, raw in tangent_blocks(w):
            rows = [[sp.sympify(v,locals=LOCALS) for v in row] for row in raw]
            if h==(0,0):
                orbit = 11
                transverse = 2
            else:
                orbit = sum(tuple(x-y for x,y in zip(character((a,)),character((b,))))==h
                            for a,b in product(range(4),repeat=2))*4
                transverse = int(h in AXES or exceptional and h in DIAGONALS | MIXED)
            locus = 'base' if exceptional else 'base*(c^2-1)'
            certify(rows,len(words)-orbit-transverse,locus)
            if not exceptional and orbit:
                orows = orbit_rows(w,words)
                cols = [j for j in range(64) if any(row[j] for row in orows)]
                clean = []
                for row in orows:
                    vals = [sp.factor(row[j]) for j in cols]
                    den = sp.lcm([sp.denom(v) for v in vals])
                    assert all(sp.denom(sp.cancel(S*C*(S**4-1)/f)).is_number
                               for f,_ in sp.factor_list(den)[1])
                    clean.append([sp.factor(v*den) for v in vals])
                certify(clean,orbit,'base')
    out = run(['Singular','-q'],input='\n'.join(script+['quit;']),text=True,
              capture_output=True,check=True).stdout
    assert '?' not in out and 'FAILED' not in out and out.count('VERIFIED')==checks, out


class CharacterEquations:
    def __init__(self, w, field):
        self.field = field
        self.w = {q:field.from_sympy(sp.sympify(v)) for q,v in w.items()}
        self.cuts = []
        for i in (0,1):
            m = self.cut(self.w,i)
            rows,cols = m.transpose().rref()[1],m.rref()[1]
            inverse = m.extract(rows,cols).inv().to_list()
            g = DomainMatrix.zeros((16,16),field).to_dense()
            for a,b in product(range(6),repeat=2):
                g[cols[a],rows[b]] = inverse[a][b]
            self.cuts.append((m.transpose().nullspace(divide_last=True),
                              m.nullspace(divide_last=True).transpose(),g))

    def matrix(self, rows, n):
        return DomainMatrix(rows,(len(rows),n),self.field)

    def cut(self, tensor, i):
        return self.matrix([[tensor.get(q[-i:]+q[:-i] if i else q,self.field.zero)
                             for q in WORDS[j:j+16]] for j in range(0,256,16)],16)

    def linear(self,h):
        words = [q for q in WORDS if character(q)==h]
        rows = []
        for i,(l,r,_) in enumerate(self.cuts):
            for a,b in product(l.to_list(),r.transpose().to_list()):
                rows.append([a[4*q[i]+q[(i+1)%4]]*b[4*q[(i+2)%4]+q[(i+3)%4]] for q in words])
        return words,self.matrix(rows,len(words))

    def tangent(self,h):
        words,m = self.linear(h)
        orbit = self.matrix([[self.w.get(q[:i]+(b,)+q[i+1:],self.field.zero) if q[i]==a else self.field.zero
                              for i,a,b in product(range(4),repeat=3)] for q in words],64)
        ker = m.nullspace(divide_last=True).transpose()
        all_columns = orbit.hstack(ker)
        pivots = all_columns.rref()[1]
        assert sum(j>=64 for j in pivots)==1
        frame = all_columns.extract(range(len(words)),pivots)
        return words,frame


def obstructions():
    eq = CharacterEquations(surface(S,C),sp.QQ.frac_field(S,C))
    directions = {}
    for h in ((2,0),(0,2),(0,-2)):
        words,frame = eq.tangent(h)
        b = dict(zip(words,[row[-1] for row in frame.to_list()]))
        directions[h] = [eq.cut(b,i) for i in (0,1)]
    for a,b in [((2,0),(2,0)),((0,2),(0,2)),((2,0),(0,2)),((2,0),(0,-2))]:
        h = tuple(x+y for x,y in zip(a,b))
        _,linear = eq.linear(h)
        rhs = []
        for i,(l,r,g) in enumerate(eq.cuts):
            x,y = directions[a][i],directions[b][i]
            block = l*(x*g*y+y*g*x)*r
            rhs.extend([v] for row in block.to_list() for v in row)
        assert linear.hstack(eq.matrix(rhs,1)).rank() == linear.rank()


def half_shift():
    f = sp.Rational(3,4)
    w = surface(sp.Integer(2),sp.Integer(1))
    def phase(q):
        return sp.prod(f if i<2 else 1/f for i,a in enumerate(q) if a<2)
    assert all(w[q[2:]+q[:2]]==sp.cancel(phase(q)*v) for q,v in w.items())
    incidence = sp.Matrix([[int(q[i]==a) for i,a in product(range(4),repeat=2)] for q in w])
    d,_,_ = smith_normal_decomp(DomainMatrix.from_Matrix(incidence).convert_to(sp.ZZ))
    assert list(d.to_Matrix().diagonal()) == [1]*11+[0]*5
    eq = CharacterEquations(w,sp.QQ)
    negatives = {(-3,1),(-2,0),(-1,3),(1,3),(2,0),(3,1)}
    for h in AXES | MIXED | DIAGONALS:
        words,frame = eq.tangent(h)
        vector = [row[-1] for row in frame.to_list()]
        shifted = eq.matrix([[sp.QQ.from_sympy(sp.sympify(phase(q[2:]+q[:2])))*
                              vector[words.index(q[2:]+q[:2])]] for q in words],1)
        rows = frame.transpose().rref()[1]
        coords = frame.extract(rows,range(frame.shape[1])).inv()*shifted.extract(rows,[0])
        assert frame*coords == shifted
        assert coords.to_list()[-1][0] == (-1 if h in negatives else 1)


if __name__ == '__main__':
    uniform_ranks()
    obstructions()
    half_shift()
    print('Reflected collision weights, quadratic solvability and every corrected half-shift sign verified.')
