"""Binary cubic polar pencil, and a free pair seen on cubic words."""
from itertools import product
from sympy import Matrix, symbols


def main():
    a, b = symbols('a b')
    polar = Matrix([[a,b,0,0], [0,a,b,0], [0,0,a,b]])
    minors = [polar[:, [i for i in range(4) if i != omitted]].det()
              for omitted in range(4)]
    assert minors == [b**3, a*b**2, a**2*b, a**3]
    # Common row (1,1) evaluates to 1 on both chosen cubic end images.
    row = Matrix([[1,1]])
    coefficients = [Matrix([1,0])*row, Matrix([0,1])*row]
    for degree in range(1, 7):
        words = {}
        for word in product(range(2), repeat=degree):
            value = Matrix.eye(2)
            for letter in word:
                value *= coefficients[letter]
            words[word] = row * value
        assert len(words) == 2**degree and all(value == row for value in words.values())
    print('Polar pencil maximal minors and nonzero coefficients of distinct free words through length six verified.')


if __name__ == '__main__':
    main()
