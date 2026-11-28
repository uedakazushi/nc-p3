"""Dimension-four profiles under finite-edge ranks and generic extension inequalities."""

from itertools import product


def row_space_bound(size):
    if size == 4:
        return 6
    if size == 5:
        return 10
    full, extra = divmod(size, 4)
    relation_rank = 6*full + (0, 3, 5, 6)[extra]
    return min(3*(size-1), 4*size-relation_rank)


def extension_inequalities(d, a, b, h, euler):
    occupied = [x*y for x, y in zip(a, b)]
    if sum(occupied) > euler:
        return False
    defects = [x+y+z for x, y, z in zip(a, b, h)]
    # A universal extension at one vertex on either side; -1 means none.
    for i, k in product(range(-1, 4), repeat=2):
        if i == k and i >= 0:
            continue
        source = a[i] if i >= 0 else 0
        target = b[k] if k >= 0 else 0
        lower = occupied.copy()
        if source:
            j = (i+1) % 4
            lower[j] = (a[j]+max(0, 4*source-defects[j]))*b[j]
        if target:
            j = (k-1) % 4
            rows = a[j]
            if source and j == (i+1) % 4:
                rows += max(0, 4*source-defects[j])
            lower[j] = rows*(b[j]+max(0, 4*target-defects[j]))
        if source:
            lower[i] = 0
        if target:
            lower[k] = 0
        omit = min(lower[j] for j in range(4)
                   if (not source or j != i) and (not target or j != k))
        upper = euler
        if source:
            upper += source*(h[i]+a[i])
        if target:
            upper += target*(h[k]+b[k])
        if source and target and (k-i) % 4 == 2:
            upper += 6*source*target
        if sum(lower)-omit > upper:
            return False
    # The effective adjacent-arrow product in the six-dimensional middle space.
    for i in range(4):
        j, k, middle = (i+1) % 4, (i+2) % 4, (i+3) % 4
        if d[middle] != 6:
            continue
        source, target = a[i], b[k]
        p, length = d[i]+source, d[k]+target
        rank_a = max(0, p-8) if length > 10 else 0
        rank_b = max(0, length-8) if p > 10 else 0
        u = a[j]+max(0, 4*source-defects[j])
        v = b[j]+max(0, 4*target-defects[j])
        upper = (euler + source*(h[i]+a[i]) + target*(h[k]+b[k])
                 + 6*source*target - max(source*rank_a, target*rank_b))
        if u*v > upper:
            return False
    return True


def survives(d):
    h = [2*d[i]-4*d[i-1]-4*d[(i+1) % 4]+6*d[(i+2) % 4] for i in range(4)]
    euler = sum(x*y for x, y in zip(d, h))
    if euler < 0:
        return False
    options = []
    for i in range(4):
        left = row_space_bound(d[i-1])-d[i]
        right = row_space_bound(d[(i+1) % 4])-d[i]
        options.append([(a, b) for a in range(left+1)
                        for b in range(max(0, -h[i]-a), right+1) if a*b <= euler])
    for profile in product(*options):
        a, b = zip(*profile)
        if extension_inequalities(d, a, b, h, euler):
            return True
    return False


def main():
    profiles = set()
    for x in range(4, row_space_bound(4)+1):
        for z in range(x, row_space_bound(4)+1):
            for y in range(4, min(row_space_bound(x), row_space_bound(z))+1):
                if survives((4, x, y, z)):
                    profiles.add((4, x, y, z))
    assert profiles == {(4, 4, 4, 4)}, profiles
    print("The unique allowed dimension vector with minimum four is (4, 4, 4, 4).")


if __name__ == "__main__":
    main()
