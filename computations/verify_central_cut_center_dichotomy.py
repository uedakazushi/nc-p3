"""Socle-three Hilbert functions and the two generic vertex ranks of finite central windows."""
from itertools import product


def main():
    choices = []
    for degree, rank, middle in product(range(1,65), range(2,9), range(32)):
        if degree*rank**2 == 64 and degree == 2 + 2*middle:
            choices.append((degree, rank, 4+middle))
    assert choices == [(4,4,5), (16,2,11)]
    assert all((-1)**((i+1)%4) == -(-1)**i for i in range(4))
    assert 4*4 - 4 == 12
    print("Two center Hilbert functions, parity conjugation, and relative dual degree twelve.")


if __name__ == "__main__":
    main()
