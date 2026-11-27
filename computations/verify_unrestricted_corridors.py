"""Enumerate necessary corridor data, allowing every finite Kronecker pencil.

Indices of a degree tuple are zero-based; the defect divisors use paper
indices d_1,...,d_m. A center key records its direction, phase and number.
The final geometric exclusions are also checked on every surviving row.
"""
from collections import defaultdict
from dataclasses import dataclass
from itertools import product


def centers(phase_counts, direction):
    return [(direction, phase, number) for phase, count in enumerate(phase_counts)
            for number in range(count)]


@dataclass
class Path:
    start: tuple
    end: tuple
    degrees: tuple
    budget: tuple
    restrictions: dict
    targets: dict
    neighbors: list


def differences(degrees):
    values = (0, 0, *degrees, 0, 0)
    return [values[j+3] - 3*values[j+2] + 3*values[j+1] - values[j]
            for j in range(len(degrees)+1)]


def possible_paths(start, end, center_groups, counts):
    phase = start[1]
    end_phase = end[1]
    n = sum(counts[1])
    total = 28 - 2*(n-1)
    for length in range(2, min(20, total)+1):
        if (phase+length+1) % 4 != end_phase:
            continue
        phases = [(phase+k+1) % 4 for k in range(length)]
        allowed = [set(range(1, 8)) for _ in phases]
        allowed[0] = {1}
        allowed[-1] = {1}
        next_groups = [(start, 1, center_groups[1, (phase+1) % 4]),
                       (end, length-2, center_groups[-1, (end_phase-1) % 4])]
        for _, index, neighbors in next_groups:
            allowed[index] &= {1, 2, 3} if neighbors else {3}
            if length > 2:
                allowed[index].discard(1)
        moves = []
        for k, source_phase in enumerate(phases[:-1]):
            visits = [(1, source_phase-1, k+1, k == 0, start),
                      (-1, source_phase+2, k, k == length-2, end)]
            for direction, origin, target, first, own in visits:
                for center in center_groups[direction, origin % 4]:
                    if first and center == own:
                        continue
                    allowed[target] &= {1, 2, 3}
                    moves.append((center, target, direction, first))
        if any(not choices for choices in allowed) or sum(map(min, allowed)) > total:
            continue

        def sign_ok(j, delta):
            if delta > 0:
                return j >= 2 and counts[-1][(phase+j) % 4]
            if delta < 0:
                return j+1 < length and counts[1][(phase+j+1) % 4]
            return True

        def degree_rows(prefix, used, remaining):
            k = len(prefix)
            if k == length:
                yield prefix, tuple(used)
                return
            slot = phases[k]
            for degree in allowed[k]:
                if used[slot]+degree > 7 or degree+length-k-1 > remaining:
                    continue
                if k >= 1:
                    previous = (0, 0, *prefix)
                    delta = degree-3*previous[-1]+3*previous[-2]-previous[-3]
                    if not sign_ok(k-1, delta):
                        continue
                used[slot] += degree
                yield from degree_rows((*prefix, degree), used, remaining-degree)
                used[slot] -= degree

        for degrees, budget in degree_rows((), [0]*4, total):
            if not all(sign_ok(j, delta) for j, delta in enumerate(differences(degrees))):
                continue
            embedded = {k for k, d in enumerate(degrees) if d <= 2} | {1, length-2}
            embedded.update(k for _, k, _, _ in moves if degrees[k] == 3)
            restrictions = defaultdict(lambda: set(range(4)))
            targets = defaultdict(list)
            for center, k, direction, first in moves:
                degree = degrees[k]
                predecessor_embedded = k-2*direction in embedded
                source_degree = degrees[k-direction]
                if predecessor_embedded and source_degree > 2:
                    break
                if first or degree != 1 or (predecessor_embedded and source_degree > 1):
                    restrictions[center] &= {degree}
                targets[center].append(degree)
            else:
                neighbor_tests = []
                for center, index, neighbors in next_groups:
                    if degrees[index] == 2:
                        required = degrees[2] if center == start else degrees[-3]
                        neighbor_tests.append((neighbors, required))
                    # -1 consumes an exceptional line for a next-phase center
                    # outside the initial line, distinct from any moving target.
                    off_line = max(0, len(neighbors) - (3-degrees[index]))
                    targets[center].extend([-1]*off_line)
                if all(restrictions.values()):
                    yield Path(start, end, degrees, budget, restrictions, targets, neighbor_tests)


def target_capacity(labels, degree):
    return (all(labels.count(d) <= 1 for d in (2, 3))
            and labels.count(-1) <= 3-degree
            and labels.count(-1)+labels.count(1) <= 3-degree+int(degree == 1))


def contact_bounds(paths, types, counts):
    short_only = {(direction, phase): all(len(p.degrees) == 2 and
                  max(types[p.start], types[p.end]) > 1 for p in paths
                  if (p.start if direction == 1 else p.end)[1] == phase)
                  for direction, phase in product((1, -1), range(4))}
    for p in paths:
        phase, end_phase, degrees = p.start[1], p.end[1], p.degrees
        m = len(degrees)
        embedded = {k for k, d in enumerate(degrees, 1) if d <= 2} | {2, m-1}
        lower_b2 = max(0, counts[-1][(phase+2) % 4] - int(types[p.start] == 0))
        lower_fpen = max(0, counts[1][(end_phase-2) % 4] - int(types[p.end] == 0))
        for j, delta in enumerate(differences(degrees)):
            slot, next_slot = (phase+j) % 4, (phase+j+1) % 4
            ub = 3*degrees[j-1] if j >= 2 and counts[-1][slot] else 0
            uf = 3*degrees[j] if j+1 < m and counts[1][next_slot] else 0
            if j in embedded and short_only[-1, slot]:
                ub = min(ub, counts[-1][slot])
            if j+1 in embedded and short_only[1, next_slot]:
                uf = min(uf, counts[1][next_slot])
            lb = lower_b2 if j == 2 else 0
            lf = lower_fpen if j+1 == m-1 else 0
            if lb > ub or lf > uf or not lb-uf <= delta <= ub-lf:
                return False
    return True


def enumerate_graph(forward, backward):
    counts = {1: forward, -1: backward}
    starts, ends = centers(forward, 1), centers(backward, -1)
    groups = {(direction, phase): [c for c in starts+ends if c[:2] == (direction, phase)]
              for direction, phase in product((1, -1), range(4))}
    types = {c: set(range(int(len(groups[c[:2]]) == 2), 4)) for c in starts+ends}
    alternatives = {}
    for c in starts:
        alternatives[c] = [p for z in ends for p in possible_paths(c, z, groups, counts)]
        if not alternatives[c]:
            return

    def matching(paths, ends_used, budget, possible, labels):
        if len(paths) == len(starts):
            for values in product(*(possible[c] for c in starts+ends)):
                assigned = dict(zip(starts+ends, values))
                if any(not any(assigned[n] == e for n in ns)
                       for p in paths for ns, e in p.neighbors):
                    continue
                if contact_bounds(paths, assigned, counts):
                    yield tuple((p.start[1], p.degrees, p.end[1], assigned[p.start], assigned[p.end])
                                for p in paths)
            return
        for p in alternatives[starts[len(paths)]]:
            if p.end in ends_used:
                continue
            new_budget = tuple(a+b for a, b in zip(budget, p.budget))
            if max(new_budget) > 7:
                continue
            feasible, combined = possible.copy(), labels.copy()
            for center in p.targets.keys() | p.restrictions.keys():
                used = labels.get(center, ()) + tuple(p.targets.get(center, []))
                combined[center] = used
                feasible[center] = {e for e in possible[center]
                                    if e in p.restrictions.get(center, range(4))
                                    and target_capacity(used, e)}
            if all(feasible.values()):
                yield from matching(paths+[p], ends_used | {p.end}, new_budget, feasible, combined)

    yield from matching([], set(), (0,)*4, types, {})


def incidence_possible(row, forward, backward):
    for direction, counts in ((1, forward), (-1, backward)):
        for phase in range(4):
            paths = [d if direction == 1 else d[::-1] for p, d, z, _, _ in row
                     if (p if direction == 1 else z) == phase]
            if sum(len(d) == 2 for d in paths) > 1:
                return False
            first_hits = sum(d[1] < 3 for d in paths)
            second_hits = sum(3*d[1]-3-(d[2] if len(d) > 2 else 0) > 0 for d in paths)
            if first_hits > counts[(phase+direction) % 4]:
                return False
            if second_hits > counts[(phase+2*direction) % 4]:
                return False
    return True


def paired_cubic_obstruction(row):
    # An initial second cubic and a terminal penultimate cubic lie at one
    # phase. If the opposite moving maps demand the primary cubics of the
    # two endpoint pencils, both contact divisors lie in their intersection.
    for p, d, _, e, _ in row:
        for _, r, z, _, h in row:
            if len(d) < 4 or len(r) < 4 or (p+2) % 4 != (z-2) % 4:
                continue
            if d[1] != 3 or r[-2] != 3 or e != 3 or h != 3:
                continue
            b2 = d[3]-3*d[2]+8
            fpen = r[-4]-3*r[-3]+8
            if b2+fpen > 5:
                return True
    return False


def shape(row):
    return min(tuple(sorted(((p+k) % 4, d, (z+k) % 4) for p, d, z, _, _ in row))
               for k in range(4))


def main():
    distributions = [s for s in product(range(3), repeat=4) if s.count(0) <= 2]
    survivors = {}
    for f, b in product(distributions, repeat=2):
        if sum(f) != sum(b):
            continue
        if any((f[k:]+f[:k], b[k:]+b[:k]) < (f, b) for k in range(1, 4)):
            continue
        for row in enumerate_graph(f, b):
            if not incidence_possible(row, f, b) or paired_cubic_obstruction(row):
                continue
            degrees = [d for _, d, _, _, _ in row]
            if all(d == (1, 3, 1) for d in degrees):
                if len(row) == 2:
                    assert sorted(f) == [0, 0, 1, 1] and f == b and f[0] == f[2]
                    continue
                assert len(row) == 4 and all(e == h == 3 for _, _, _, e, h in row)
                continue
            if all(d == (1, 1) for d in degrees):
                assert f == b == (1, 1, 1, 1)
                continue
            survivors.setdefault(shape(row), row)
    representatives = [
        [(1, (1,1)), (2, (1,1)), (3, (1,3,1))],
        [(1, (1,2,1)), (2, (1,2,3,2,1)), (3, (1,3,1))],
        [(k, (1,2,1)) for k in range(4)],
        [(k, (1,2,2,1) if k % 2 == 0 else (1,1)) for k in range(4)],
    ]
    expected = {shape([(p, d, (p+len(d)+1) % 4, None, None) for p, d in rep])
                for rep in representatives}
    assert set(survivors) == expected, (set(survivors)-expected, expected-set(survivors))
    for configuration in sorted(survivors):
        print(configuration)
    print('PASS: four unrestricted corridor configurations up to notation rotation')


if __name__ == '__main__':
    main()
