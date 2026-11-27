"""Finite enumeration of necessary integer path degrees for pure-free branches."""
from collections import Counter
from itertools import product


def paths(start, length, forward, backward):
    if length == 2:
        if (start+1) % 4 in forward and (start+2) % 4 in backward:
            yield (1, 1)
        return
    choices = [set(range(1, 8)) for _ in range(length)]
    choices[0] = choices[-1] = {1}
    choices[1] &= {2 if (start+1) % 4 in forward else 3}
    choices[-2] &= {2 if (start+length) % 4 in backward else 3}
    for j in range(length-1):
        phase = (start+j+1) % 4
        if (phase-1) % 4 in forward and j != 0:
            choices[j+1] &= {3}
        if (phase+2) % 4 in backward and j != length-2:
            choices[j] &= {3}
    for sequence in product(*choices):
        extended = (0, 0, *sequence, 0, 0)
        for j in range(length+1):
            a, b, c, d = extended[j:j+4]
            delta = d-3*c+3*b-a
            if delta > 0 and (j < 2 or (start+j) % 4 not in backward):
                break
            if delta < 0 and (j+1 >= length or (start+j+1) % 4 not in forward):
                break
        else:
            yield sequence


def candidates():
    for fm, bm in product(range(16), repeat=2):
        forward = tuple(i for i in range(4) if fm & (1 << i))
        backward = tuple(i for i in range(4) if bm & (1 << i))
        if not {0, 1} <= set(forward) or not {2, 3} <= set(backward):
            continue
        if len(forward) != len(backward):
            continue
        for lengths in product((2,), *(range(2, 13) for _ in forward[1:])):
            if sorted((i+m+1) % 4 for i, m in zip(forward, lengths)) != list(backward):
                continue
            moves = Counter((i+j) % 4 for i, m in zip(forward, lengths) for j in range(1, m))
            if any(moves[h] > 2 for h in range(4)
                   if (h-1) % 4 in forward or (h+2) % 4 in backward):
                continue
            for ds in product(*(paths(i, m, forward, backward) for i, m in zip(forward, lengths))):
                degrees = Counter()
                for i, sequence in zip(forward, ds):
                    for j, degree in enumerate(sequence, 1):
                        degrees[(i+j) % 4] += degree
                if max(degrees.values()) <= 7:
                    yield forward, backward, ds


if __name__ == '__main__':
    two, short, long = (1, 1), (1, 3, 1), (1, 3, 3, 3, 3, 3, 1)
    all_phases = tuple(range(4))
    rows = {((0, 1, 2), (0, 2, 3), (two, two, d)) for d in (short, long)}
    rows |= {((0, 1, 3), (1, 2, 3), (two, d, two)) for d in (short, long)}
    rows.add((all_phases, all_phases, (two,)*4))
    assert set(candidates()) == rows
    print("PASS: five necessary rows before the intersection bound, two up to cyclic notation after it")
