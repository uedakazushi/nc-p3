"""Filtered central fibers, a transverse plane quartic, and the regular trace pairing."""

from functools import cache
from itertools import product
from math import gcd
from functools import reduce

from flint import nmod_mat, nmod_poly
import sympy as sp

from verify_klein_elliptic_completion import rational_example, four_relations
from verify_klein_period_four_center import K, paths, commuting_space


def integral_quartics(relations, prime):
    basis, image = paths(relations, 5)

    def residue(x):
        return int(x.p) * pow(int(x.q), -1, prime) % prime

    rules = [
        [
            ({divmod(i, 4): residue(x) for i, x in enumerate(row) if x}, None)
            for _, row in rs
        ]
        for rs in relations
    ]
    quartics = []
    for h, choices in [(0, [(1, 3, 5), (2, 7, 11), (1, 0, 0)]), (1, [(1, 0), (0, 1)])]:
        space = commuting_space(basis, image, h, 0)
        for coeff in choices:
            section = space * sp.Matrix(coeff)
            section *= sp.ilcm(*(x.q for x in section))
            section /= reduce(gcd, map(int, section))
            pieces = []
            start = 0
            for i in K:
                words = basis[i, 4, h]
                pieces.append(
                    {
                        w: residue(v)
                        for w, v in zip(words, section[start : start + len(words)])
                        if v
                    }
                )
                start += len(words)
            quartics.append(pieces)
    return rules, [quartics[i] for i in [0, 1, 3, 4, 2]]


class CentralFiber:
    """Filtered spaces of words of a given length modulo four, at fixed parameters."""

    def __init__(self, rules, quartics, parameters, prime):
        self.prime = prime
        self.prefix = cache(self.prefix)
        self.basis = {(p, 0): [()] for p in K}
        self.unit = {(p, 0): nmod_mat([[1]], prime) for p in K}
        self.left = {}
        relations = [list(rs) for rs in rules]
        for section, value in zip(quartics, parameters):
            for p in K:
                relations[p].append((section[p], value))
        for n in range(1, 17):
            for p in K:
                self.reduce_degree(p, n, relations[p])
        hilbert = [1, 4, 10, 20, 31, 40, 44, 40, 31, 20, 10, 4, 1]
        for n in range(17):
            expected = sum(hilbert[d] for d in range(n % 4, min(n, 12) + 1, 4))
            assert all(len(self.basis[p, n]) == expected for p in K)

    def matrix(self, rows):
        return nmod_mat(rows, self.prime)

    def pivot_columns(self, constraint):
        reduced, rank = constraint.rref()
        pivots = [
            next(j for j in range(constraint.ncols()) if reduced[i, j])
            for i in range(rank)
        ]
        return reduced, pivots

    def invert(self, matrix):
        return matrix.inv()

    def prefix(self, p, n, word):
        if not word:
            dim = len(self.basis[p, n])
            return self.matrix([[int(i == j) for j in range(dim)] for i in range(dim)])
        return self.left[p, n, word[0]] * self.prefix((p + 1) % 4, n - 1, word[1:])

    def evaluate(self, p, n, word):
        return self.prefix(p, n, word) * self.unit[(p + len(word)) % 4, n - len(word)]

    def reduce_degree(self, p, n, rules):
        prime = self.prime
        after = (p + 1) % 4
        old = self.basis[after, n - 1]
        size = len(old)
        words = [(j,) + w for j in K for w in old]
        if n % 4 == 0:
            words.append(())
        rows = []
        if n >= 4:
            lower = self.basis[p, n - 4]
            embeddings = []
            for w in lower:
                vec = [0] * len(words)
                if w:
                    start = w[0] * size
                    vec[start : start + size] = self.evaluate(
                        after, n - 1, w[1:]
                    ).entries()
                else:
                    vec[-1] = 1
                embeddings.append(vec)
            include = self.matrix(embeddings)
        for rule, constant in rules:
            length = len(next(iter(rule)))
            if length > n:
                continue
            tails = len(self.basis[(p + length) % 4, n - length])
            block = [nmod_mat(size, tails, prime) for _ in K]
            for term, coeff in rule.items():
                block[term[0]] += coeff * self.prefix(after, n - 1, term[1:])
            terms = [row for b in block for row in b.table()]
            if n % 4 == 0:
                terms.append([0] * tails)
            equations = self.matrix(terms).transpose()
            if constant is not None:
                equations -= constant * include
            rows.extend(equations.table())
        constraint = self.matrix(rows) if rows else nmod_mat(0, len(words), prime)
        reduced, pivots = self.pivot_columns(constraint)
        free = [j for j in range(len(words)) if j not in pivots]
        self.basis[p, n] = [words[j] for j in free]
        projection = []
        for j in free:
            row = [0] * len(words)
            row[j] = 1
            for i, k in enumerate(pivots):
                row[k] = -int(reduced[i, j]) % prime
            projection.append(row)
        for i in K:
            self.left[p, n, i] = self.matrix(
                [row[i * size : (i + 1) * size] for row in projection]
            )
        if n % 4 == 0:
            self.unit[p, n] = self.matrix([[row[-1]] for row in projection])

    def fifth_action(self, section):
        words = self.basis[0, 12]
        include = self.matrix(
            [self.evaluate(0, 16, w).entries() for w in words]
        ).transpose()
        action = nmod_mat(64, 64, self.prime)
        for t, c in section[0].items():
            action += c * self.prefix(0, 16, t)
        return self.invert(include) * action

    def center_dimension(self):
        rows = []
        for i, j in product(K, repeat=2):
            k = (i + 1) % 4
            left = self.matrix(
                [self.evaluate(i, 13, w + (j,)).entries() for w in self.basis[i, 12]]
            ).transpose()
            right = self.left[i, 13, j]
            for v, w in zip(left.table(), right.table()):
                row = [0] * 256
                row[64 * i : 64 * (i + 1)] = v
                row[64 * k : 64 * (k + 1)] = [-int(x) for x in w]
                rows.append(row)
        return 256 - self.matrix(rows).rank()

    def loop_operators(self, words):
        # All four spaces end at vertex zero; stabilized dimensions are 64.
        action = {}
        levels = [12, 11, 10, 9]
        for i, j in product(K, repeat=2):
            n = levels[(i + 1) % 4] + 1
            mat = self.left[i, n, j]
            if n != levels[i]:
                embed = self.matrix(
                    [self.evaluate(i, n, w).entries() for w in self.basis[i, levels[i]]]
                ).transpose()
                mat = embed.inv() * mat
            action[i, j] = mat
        operators = []
        for w in words:
            op = self.prefix(0, 12, ())
            for k, j in enumerate(w):
                op = op * action[k % 4, j]
            operators.append(op)
        return operators

    def trace_rank(self):
        operators = self.loop_operators(self.basis[0, 12])
        a = self.matrix([op.entries() for op in operators])
        b = self.matrix([op.transpose().entries() for op in operators])
        return (a * b.transpose()).rank()


def plane_section(rules, quartics, prime):
    samples = []
    for t in range(5):
        params = (1, t, 1 + 2 * t, 3 - t)
        f = CentralFiber(rules, quartics, params, prime)
        action = f.fifth_action(quartics[4])
        polynomial = action.minpoly()
        # Faithfulness of the quartic powers holds over every central fiber once
        # their independence modulo the four parameters is certified.
        assert polynomial.degree() == 4 and polynomial**16 == action.charpoly()
        samples.append([int(polynomial[4 - i]) for i in range(5)])
        f.prefix.cache_clear()
    interpolate = nmod_mat(
        [[j**i % prime for i in range(5)] for j in range(5)], prime
    ).inv()
    coeff = interpolate * nmod_mat(samples, prime)
    s, t, z = sp.symbols("s t z")
    assert all(not coeff[j, i] for i in range(5) for j in range(i + 1, 5))
    form = sum(
        int(coeff[j, i]) * s ** (i - j) * t**j * z ** (4 - i)
        for i in range(5)
        for j in range(i + 1)
    )
    variables = (s, t, z)
    for i, v in enumerate(variables):
        sub = {**dict.fromkeys(variables[:i], 0), v: 1}
        equations = [sp.diff(form, u).subs(sub) for u in variables]
        if i == 2:
            assert any(int(e) % prime for e in equations)
        else:
            assert list(sp.groebner(equations, *variables[i + 1 :], modulus=prime)) == [
                1
            ]
    print("Smooth plane section of the monic quartic, including the chart at infinity.")


def main():
    prime = 1009
    w, _, n = rational_example()
    rules, z = integral_quartics(four_relations(w, n), prime)
    zero = CentralFiber(rules, z, (0, 0, 0, 0), prime)
    assert zero.fifth_action(z[4]).minpoly() == nmod_poly([0, 0, 0, 0, 1], prime)
    zero.prefix.cache_clear()
    fiber = CentralFiber(rules, z, (1, 2, 3, 4), prime)
    action = fiber.fifth_action(z[4])
    f = nmod_poly([727, 542, 601, 421, 1], prime)
    assert (
        action.minpoly() == f
        and action.charpoly() == f**16
        and f.gcd(f.derivative()).degree() == 0
    )
    assert fiber.center_dimension() == 4 and fiber.trace_rank() == 64
    fiber.prefix.cache_clear()
    print(
        "Four central powers independent at zero; center dimension 4, corner trace rank 64 at (1,2,3,4)."
    )
    plane_section(rules, z, prime)


if __name__ == "__main__":
    main()
