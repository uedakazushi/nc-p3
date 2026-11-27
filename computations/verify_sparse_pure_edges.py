"""Enumerate finite center-pencil defects when every edge is simple or all pure points lie on one edge."""
from collections import Counter
from itertools import product

from verify_single_pure_corridors import centers


def capacity(loss, degree):
    return loss + int(loss == 2) if degree == 1 else 1


def path_choices(source, target, f, b, fl, bl, heavy):
    start, end = source[0], target[0]
    # A phase with an ordinary center permits at most two moves. Otherwise
    # a simple exceptional center permits at most three.
    limit = 12 if sum(f) > len(fl) or sum(b) > len(bl) else 16
    if source in fl and (start+1) % 4 not in heavy:
        paired = ((start-1) % 4, 0)
        if paired in bl and target != paired:
            return
    for m in range(2, limit+1):
        if (start+m+1) % 4 != end:
            continue
        if source in fl and target in bl and end == (start-1) % 4 and (start+1) % 4 not in heavy and m != 2:
            continue
        choices = [set(range(1, 8)) for _ in range(m)]
        choices[0] = choices[-1] = {1}
        for p, counts, loss, k, step in [(source, f, fl, 1, 1), (target, b, bl, m-2, -1)]:
            neighbors = counts[(p[0]+step) % 4]
            contact = max(0, neighbors-loss.get(p, 0))
            choices[k] &= set(range(1, 4-contact)) if neighbors else {3}
            if counts[p[0]] == 2:
                choices[k] &= {3-loss.get((p[0], 1-p[1]), 0)}
        if m > 2:
            choices[1] &= {2, 3}
            choices[-2] &= {2, 3}
        labels = []
        for j in range(m-1):
            phase = (start+j+1) % 4
            for direction, counts, losses, mark, h, k, first in [
                    (1, f, fl, source, (phase-1) % 4, j+1, j == 0),
                    (-1, b, bl, target, (phase+2) % 4, j, j == m-2)]:
                for n in range(counts[h]):
                    p = (h, n)
                    if first and p == mark:
                        continue
                    loss = losses.get(p, 0)
                    choices[k] &= ({1, 3-loss} - {0}) if loss else {3}
                    labels.append((direction, p, k, loss))

        for path in product(*choices):
            d = (0, 0, *path, 0, 0)
            for j in range(m+1):
                delta = d[j+3]-3*d[j+2]+3*d[j+1]-d[j]
                if delta > 0 and (j < 2 or not b[(start+j) % 4]):
                    break
                if delta < 0 and (j+1 >= m or not f[(start+j+1) % 4]):
                    break
            else:
                embedded = {j for j, degree in enumerate(path) if degree <= 2}
                embedded.update((1, m-2))
                embedded.update(k for _, _, k, _ in labels if path[k] == 3)
                for direction, _, k, loss in labels:
                    if k-2*direction in embedded:
                        bound = 2 if path[k] == 3-loss else 1
                        if path[k-direction] > bound:
                            break
                else:
                    used = Counter((direction, p, loss, path[k]) for direction, p, k, loss in labels)
                    budget = tuple(sum(path[k] for k in range(m) if (start+k+1) % 4 == h) for h in range(4))
                    if max(budget) <= 7 and all(n <= capacity(loss, degree) for (_, _, loss, degree), n in used.items()):
                        yield target, path, budget, used


def match(choices, used_targets, budget, used_components):
    if not choices:
        yield ()
        return
    for target, path, add, labels in choices[0]:
        new_budget = tuple(x+y for x, y in zip(budget, add))
        if target in used_targets or max(new_budget) > 7:
            continue
        used = used_components+labels
        if any(n > capacity(loss, degree) for (_, _, loss, degree), n in used.items()):
            continue
        for tail in match(choices[1:], used_targets | {target}, new_budget, used):
            yield (path,)+tail


def phase_profiles(count, phase, maximum):
    # A pair of centers has e >= 1, hence loss at most two.
    if count == 2:
        for a in range(min(maximum, 2)+1):
            for b in range(a+1):
                if a+b <= maximum:
                    yield {(phase, j): n for j, n in enumerate((a, b)) if n}
    else:
        for n in range(min(3, maximum)*count+1):
            yield {(phase, 0): n} if n else {}


def profiles(counts, caps, step):
    each = [phase_profiles(n, i, caps[(i+step) % 4]) for i, n in enumerate(counts)]
    for pieces in product(*each):
        yield {p: v for piece in pieces for p, v in piece.items()}


def configurations(caps):
    heavy = {i for i, v in enumerate(caps) if v > 1}
    counts = [cs for cs in product(range(3), repeat=4) if cs.count(0) <= 2]
    for f, b in product(counts, repeat=2):
        if sum(f) != sum(b):
            continue
        sources, targets = centers(f), centers(b)
        for fl, bl in product(profiles(f, caps, 1), profiles(b, caps, -2)):
            if not (fl or bl):
                continue
            if not heavy and (0, 0) not in fl:
                continue  # cyclic notation and reversing all paths
            opts = []
            for p in sources:
                possibilities = [option for q in targets for option in path_choices(p, q, f, b, fl, bl, heavy)]
                if not possibilities:
                    break
                opts.append(possibilities)
            else:
                for paths in match(opts, set(), (0,)*4, Counter()):
                    short = {p[0] for p, d in zip(sources, paths) if d == (1, 1)}
                    rev = {(i+3) % 4 for i in short}
                    if any((i+1) % 4 in short and f[(i+1) % 4] == 1 and (i+2) % 4 not in heavy for i in short):
                        continue
                    if any((i-1) % 4 in rev and b[(i-1) % 4] == 1 and (i-3) % 4 not in heavy for i in rev):
                        continue
                    cubic = Counter(p[0] for p, d in zip(sources, paths) if d == (1, 3, 1))
                    if any(n > 1 and (f[(i+2) % 4] == 1 or b[(i+2) % 4] == 1) for i, n in cubic.items()):
                        continue
                    yield f, b, fl, bl, paths


def verify():
    found = set()
    for caps in ((1, 1, 1, 1), (0, 4, 0, 0)):
        for f, b, fl, bl, paths in configurations(caps):
            assert f == b and f in ((1, 0, 1, 0), (0, 1, 0, 1))
            if paths == ((1, 3, 1),)*2:
                found.add((caps, 3))
                continue
            assert caps == (1, 1, 1, 1)
            for p, path in zip(centers(f), paths):
                if len(path) == 7:
                    assert p not in fl and p not in bl
                    a, m, c = path[2:5]
                    assert path == (1, 3, a, m, c, 3, 1)
                    assert all(n in (1, 2) for n in (a, m, c))
                    assert (m-3*a+8)+(m-3*c+8) >= 6
                else:
                    assert path == (1, 3, 1)
            found.add((caps, 7))
    assert found == {((1, 1, 1, 1), 3), ((1, 1, 1, 1), 7), ((0, 4, 0, 0), 3)}
    print("PASS: opposite single-center phases; long cubic contact is at least six")


if __name__ == '__main__':
    verify()
