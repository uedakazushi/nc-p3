"""Check the stability weights and the profiles in the four-vertex compactification."""
from itertools import product


def main():
    theta = (1, 3, 9, -13)

    def weight(d):
        return sum(a*b for a, b in zip(theta, d))

    assert [d for d in product(range(3), repeat=4) if not weight(d)] == [
        (0,)*4, (1,)*4, (2,)*4]
    assert all(sum(theta[:j]) > 0 for j in range(1, 4))
    assert all(-sum(theta[-j:]) > 0 for j in range(1, 4))
    for d in product(range(3), repeat=4):
        if d in [(0,)*4, (2,)*4]:
            continue
        if all(not d[i] or d[i+1] == 2 for i in range(3)):
            assert weight(d) < 0
    print('PASS: stable cut profiles, proper-prefix signs and equal-slope factors')


if __name__ == '__main__':
    main()
