"""Check the four-corner Hilbert functions of the canonical genus-three section."""
from math import comb

from verify_successively_central_cutoff import main as check_central_sequence


def paths(n):
    return comb(n+3, 3) if n >= 0 else 0


def main():
    check_central_sequence()
    hilbert = [paths(n)-2*paths(n-4)+paths(n-8) for n in range(41)]
    assert hilbert[:5] == [1, 4, 10, 20, 33]
    assert all(hilbert[n] == 16*n-32 for n in range(5, len(hilbert)))
    assert all(hilbert[n]-hilbert[4-n] == 16*n-32 for n in range(5))
    assert sum((-1)**j*comb(4, j)*paths(13-4*j) for j in range(5)) == 0
    print('PASS: genus-three Riemann--Roch, duality and finite section quotient')


if __name__ == '__main__':
    main()
