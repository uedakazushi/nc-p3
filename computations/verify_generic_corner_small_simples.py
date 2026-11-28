"""Enumerate dimension profiles with generic separation after deleting one vertex."""

from verify_opposite_four_matrix_corners import surviving_extensions


def main():
    for size in (2, 3):
        survivors = set()
        for x in range(2, 3*(size-1)+1):
            for z in range(x, 3*(size-1)+1):
                for y in range(2, 3*(x-1)+1):
                    d = (size, x, y, z)
                    if surviving_extensions(d, set(), separated=True):
                        survivors.add(d)
        assert survivors == {(size,)*4}, survivors
    print("All profiles with a dimension two or three are balanced under generic corner separation.")


if __name__ == "__main__":
    main()
