"""Exhaustive phase and degree test when the whole pure locus is a single simple point."""
from collections import Counter
from itertools import product


def centers(counts):
    return [(i, j) for i, count in enumerate(counts) for j in range(count)]


def path_options(source, target, forward, backward, bad_source, bad_target):
    i, j = source[0], target[0]
    for length in range(2, 13):
        if (i+length+1) % 4 != j:
            continue
        # If both exceptional pencils occur they share the unique pure point.
        if source == bad_source and bad_target is not None and (target != bad_target or length != 2):
            continue
        choices = [set(range(1, 8)) for _ in range(length)]
        choices[0] = choices[-1] = {1}
        for mark, counts, bad, place, step in [(source, forward, bad_source, 1, 1),
                                               (target, backward, bad_target, length-2, -1)]:
            neighbors = counts[(mark[0]+step) % 4]
            min_contact = max(0, neighbors - int(mark == bad))
            choices[place] &= set(range(1, 4-min_contact)) if neighbors else {3}
            if counts[mark[0]] == 2:
                other = (mark[0], 1-mark[1])
                choices[place] &= {2 if other == bad else 3}
        if length > 2:
            choices[1] &= {2, 3}
            choices[-2] &= {2, 3}

        # Noninitial transitions are labelled by their forced target component;
        # the backward test labels their forced source component.
        labels = []
        for k in range(length-1):
            h = (i+k+1) % 4
            for direction, counts, mark, bad, phase, index, own in [
                    (1, forward, source, bad_source, (h-1) % 4, k+1, k == 0),
                    (-1, backward, target, bad_target, (h+2) % 4, k, k == length-2)]:
                for n in range(counts[phase]):
                    center = (phase, n)
                    if own and center == mark:
                        continue
                    choices[index] &= {1, 2} if center == bad else {3}
                    labels.append((direction, center, index))

        for path in product(*choices):
            extended = (0, 0, *path, 0, 0)
            for k in range(length+1):
                a, b, c, d = extended[k:k+4]
                delta = d-3*c+3*b-a
                if delta > 0 and (k < 2 or not backward[(i+k) % 4]):
                    break
                if delta < 0 and (k+1 >= length or not forward[(i+k+1) % 4]):
                    break
            else:
                used = frozenset((direction, center, path[index]) for direction, center, index in labels)
                budget = [0]*4
                for k, degree in enumerate(path, 1):
                    budget[(i+k) % 4] += degree
                if len(used) == len(labels) and max(budget) <= 7:
                    yield target, path, budget, used


def match(options, occupied=frozenset(), budget=(0, 0, 0, 0), used=frozenset()):
    if not options:
        yield ()
        return
    for target, path, increment, labels in options[0]:
        total = tuple(a+b for a, b in zip(budget, increment))
        if target in occupied or labels & used or max(total) > 7:
            continue
        for rest in match(options[1:], occupied | {target}, total, used | labels):
            yield (path,)+rest


def candidates():
    counts = [cs for cs in product(range(3), repeat=4) if cs.count(0) <= 2]
    for forward, backward in product(counts, repeat=2):
        if sum(forward) != sum(backward):
            continue
        sources, targets = centers(forward), centers(backward)
        for bad_source in [None]+([(0, 0)] if forward[0] else []):
            for bad_target in [None]+([(3, 0)] if backward[3] else []):
                if bad_source is None and bad_target is None:
                    continue
                options = []
                for source in sources:
                    available = [item for target in targets for item in
                                 path_options(source, target, forward, backward, bad_source, bad_target)]
                    if not available:
                        break
                    options.append(available)
                else:
                    for paths in match(options):
                        short = {p[0] for p, d in zip(sources, paths) if d == (1, 1)}
                        if any(forward[(i+1) % 4] == 1 and (i+1) % 4 in short for i in short):
                            continue
                        cubic = Counter(p[0] for p, d in zip(sources, paths) if d == (1, 3, 1))
                        if any(n > 1 and forward[(i+2) % 4] == 1 for i, n in cubic.items()):
                            continue
                        yield forward, backward, bad_source, bad_target, paths


if __name__ == '__main__':
    even, odd = (1, 0, 1, 0), (0, 1, 0, 1)
    pair = ((1, 3, 1),)*2
    assert set(candidates()) == {(even, even, (0, 0), None, pair),
                                 (odd, odd, None, (3, 0), pair)}
    print("PASS: only the two labelled alternating cubic pairs survive the integer tests")
