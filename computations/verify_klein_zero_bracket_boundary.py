"""Exponential monomial quotients and a normal sequence for the single sparse survivor."""
from collections import Counter, defaultdict
from itertools import permutations, product

import sympy as sp
from verify_klein_mixed_supports import off_diagonals, tensor


def rotate(pair):
    q,g = pair
    return tuple((c,d,a,b) for a,b,c,d in g), tuple((b,a,d,c) for a,b,c,d in q)


def has_free_cycle(pair):
    reachable = set()
    for _ in range(4):
        q,g = pair
        for j,k in permutations(range(3),2):
            # Original edge column *i reads (G_j)[k,i], and i* reads (Q_k)[j,i].
            if not g[j][int(k!=(j+1)%3)] and not q[k][int(j!=(k+1)%3)]:
                reachable.add((j,k))
        pair = rotate(pair)
    for _ in range(3):
        reachable |= {(a,d) for a,b in reachable for c,d in reachable if b==c}
    return any(a==b for a,b in reachable)


def supports():
    masks = [v for v in product((0,1),repeat=4) if v[0]*v[3]==v[1]*v[2]]
    left,right = zip(*off_diagonals())
    def signature(side,terms):
        return tuple(side[i][s]*side[j][t] for (_,i,s),(_,j,t) in terms)
    sides = list(product(masks, repeat=3))
    images = defaultdict(list)
    for side in sides:
        images[signature(side,left)].append(side)
    unseen = {(q,g) for q in sides for g in images[signature(q,right)]}
    assert len(unseen)==47492
    survivors = []
    counts = Counter()
    while unseen:
        pair = min(unseen)
        fails = has_free_cycle(pair)
        counts[fails] += 1
        if not fails:
            survivors.append(tuple(tuple(sum(bit<<s for s,bit in enumerate(v)) for v in side) for side in pair))
        orbit = set()
        for _ in range(4):
            for pi in permutations(range(3)):
                odd = sum(pi[i]>pi[j] for i,j in ((0,1),(0,2),(1,2))) % 2
                orbit.add(tuple(tuple(tuple(side[i][s^odd] for s in range(4)) for i in pi) for side in pair))
            pair = rotate(pair)
        assert orbit <= unseen
        unseen -= orbit
    assert counts == {True:2036,False:2}
    assert survivors == [((0,15,15),)*2, ((15,15,15),)*2]


def boundary(x,b):
    entries = {}
    for i,(m,n) in enumerate(((0,0),(1,b),(1,x/b))):
        for j,v in enumerate((m,n,-m,-n)):
            entries[1,i,j] = sp.sympify(v)
    for i,r in enumerate((0,1,x)):
        for j,v in enumerate((int(i!=0),int(i!=0),r,r)):
            entries[0,i,j] = sp.sympify(v)
    return tensor(entries)


def square_data(w,p):
    squares = sp.zeros(4,2)
    squares[2:,:] = sp.eye(2)
    cross = {}
    for color,direction in product(range(1,4), range(2)):
        end = (color,0) if direction else (0,color)
        relation = {}
        for a,b in product(range(4),repeat=2):
            key = (a,b)+end
            key = key[-p:]+key[:-p] if p else key
            if key in w:
                relation[a^(p%2),b^((p+1)%2)] = w[key]
        i,j = sorted(next(iter(relation)))
        if i==j:
            i = 0 if (0,0) in relation else 1
            squares[i,:] = sp.Matrix([[-relation[k,k]/relation[i,i] for k in (2,3)]])
        else:
            cross[i,j] = -relation[i,j]/relation[j,i]
    return squares,cross


def diamonds():
    x,b = sp.symbols('x b',nonzero=True)
    w = boundary(x,b)
    data = [square_data(w,p) for p in range(4)]
    moves = {}
    for p,(sq,c) in enumerate(data):
        sq1,c1 = data[(p+1)%4]
        for i in (0,1):
            moves[p,i] = sp.diag(*[1/(c[i,j]*c1[i,j]) for j in (2,3)])
        for j in (2,3):
            # Its own cube and the square of letter 0 give two independent rows, even at x=b.
            source = sq1[j,:].col_join(sq1[0,:])
            target = sq[j,:].col_join(c[0,j]*c1[0,j]*sq[0,:])
            m = (source.inv()*target).applyfunc(sp.factor)
            assert m.is_diagonal()
            for a in (0,1):
                assert (sq1[a,:]*m-c[a,j]*c1[a,j]*sq[a,:]).applyfunc(sp.factor)==sp.zeros(1,2)
            moves[p,j] = m
            for item in (source.det(),m.det()):
                assert all(f in (x,b) for f,_ in sp.factor_list(sp.fraction(sp.factor(item))[0])[1])
    def rule(p,i,j):
        sq,c = data[p]
        if i==j and i<4:
            return {(4+k,):v for k,v in enumerate(sq.row(i)) if v}
        if 2<=i<4 and j<2:
            return {(j,i):c[j,i]}
        if i<4<=j:
            return {(j,i):moves[p,i][j-4,j-4]}
        if i==5 and j==4:
            return {(4,5):moves[p,3][0,0]*moves[(p+1)%4,3][0,0]}
        return {}
    def weight(word):
        return sum(1 if i<4 else 2 for i in word)
    def normal(p,word,first=None):
        for h in range(len(word)-1):
            rhs = rule((p+weight(word[:h]))%4,*word[h:h+2])
            if rhs and (first is None or first==h):
                terms = defaultdict(int)
                for r,f in rhs.items():
                    for k,v in normal(p,word[:h]+r+word[h+2:]).items():
                        terms[k] += f*v
                return {k:z for k,v in terms.items() if (z:=sp.factor(v))!=0}
        return {word:1}
    for p in range(4):
        overlaps = [(i,j,k) for i,j,k in product(range(6),repeat=3)
                    if rule(p,i,j) and rule((p+weight((i,)))%4,j,k)]
        assert len(overlaps)==32
        for word in overlaps:
            assert normal(p,word,0)==normal(p,word,1)


if __name__ == '__main__':
    supports()
    diamonds()
    print('PASS: 47492 supports, 2036 exponential orbits, all 128 Laurent diamond identities')
