"""Rational coordinates and the simultaneous sign on the ordered quartet chart."""
import sympy as s

import verify_involution_quartet_component as quartet


def check():
    M, choices, obstructions = quartet.matrices()
    quartet.setup(obstructions)
    n, d = s.fraction(s.factor(quartet.delta.as_expr()))
    X = quartet.xx
    A, B = (s.Poly(f, X) for f in (n, d))
    assert A.degree() == B.degree() == 1
    resultant = A.nth(1)*B.nth(0)-A.nth(0)*B.nth(1)
    assert resultant.subs(dict(zip((quartet.x, quartet.u, quartet.B, quartet.F), (2, 3, 5, 7)))) != 0
    # Evaluate the split-conic constants and the denominator of the inverse at t.
    T, p = quartet.sample()
    raw = {name: int(T[tuple(map(int, word))]) for word, name in quartet.renames.items()}
    def ratio(n, d):
        return n*pow(d, -1, p) % p
    coordinates = [ratio(raw[l]*raw[r], raw[k]*raw['P'])
                   for l, r, k in [('l','m','E'), ('i','o','G'), ('i','l','B'), ('m','o','D')]]
    coordinates += [ratio(raw['B'], raw['A']), ratio(raw['F'], raw['E']), ratio(raw['k'], raw['l'])]
    substitution = dict(zip((quartet.x, X, quartet.u, quartet.uu, quartet.B, quartet.F, quartet.k),
                            coordinates))
    assert quartet.residue((A.nth(1)-quartet.k**2*B.nth(1)).subs(substitution), p)
    free = {v: 0 for v, expression in choices.items() if v == expression}
    pos = M[0]
    positive = quartet.norm(pos[2:, 2:]-pos[2:, :2]*pos[:2, :2].inv()*pos[:2, 2:])
    for i, (v0, w0) in zip((0, 2), (('w0312', 'w0313'), ('w2130', 'w2131'))):
        block = positive[i:i+2, i:i+2]
        for entry, var, parity in zip(block, (v0,w0,w0,v0), (1,-1,-1,1)):
            c = entry.subs(free)
            assert s.factor(entry-c-2*s.Symbol(var)) == 0
            assert s.factor(c.xreplace({quartet.k: -quartet.k})-parity*c) == 0
    base = M[0].subs(substitution).subs(free)
    schur = base[2:, 2:]-base[2:, :2]*base[:2, :2].inv()*base[:2, 2:]
    constants = []
    for i in (0, 2):
        a, b, c, d = schur[i:i+2, i:i+2]
        constants.append(quartet.residue((a-d)**2-(b-c)**2, p))
    assert constants == [2218, 1246]
    def sign(word):
        return (-1)**sum(int(i) % 2 + int(i)//2 for i in word)
    named = {k: sign(word) for word, k in quartet.renames.items()}
    assert all(named[k] == 1 for k in ('P', 'Q', 'A', 'E', 'l', 'B', 'F', 'm'))
    assert named['k'] == -1
    assert [sign(str(f)[1:]) for f in sorted(free, key=str)] == [1, -1, 1, -1]


if __name__ == '__main__':
    check()
    print('PASS: inverse fractional-linear coordinate, both split conics and the residual sign')
