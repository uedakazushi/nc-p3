"""Right multiplication after specializing central quartics to multiples of one variable."""

from verify_klein_center_lift import add_multiple


def pivot_words(rows, length):
    """Reduce only the columns of uncut words of the current length."""
    pivots = {}
    residuals = []
    for source in rows:
        row = source.copy()
        for column, pivot in pivots.items():
            if column in row:
                add_multiple(row, pivot, -row[column])
        words = [w for w in row if len(w) == length]
        if words:
            column = min(words)
            pivots[column] = {w: z / row[column] for w, z in row.items()}
        elif row:
            residuals.append(row)
    reduced = {}
    for column in reversed(pivots):
        row = pivots[column].copy()
        for w in list(row):
            if w in reduced:
                add_multiple(row, reduced[w], -row[w])
        reduced[column] = row
    return reduced, residuals


def independent_sources(rows):
    """Propagate original rows, retaining an independent set over the input field."""
    selected = []
    pivots = {}
    for source in rows:
        row = source.copy()
        for column, pivot in pivots.items():
            if column in row:
                add_multiple(row, pivot, -row[column])
        if row:
            column = min(row)
            pivots[column] = {w: z / row[column] for w, z in row.items()}
            selected.append(source)
    return selected


def central_border(
    relations,
    centers,
    direction,
    one,
    *,
    limit=13,
    reduce_words=pivot_words,
    select=independent_sources,
):
    """Compute one start of the four-cut using the five-cut's normal words.

    A dictionary key w at length n denotes z**((n-len(w))/4) times w.
    relations(n,h) gives the quadratic rows for each incoming character.
    Both ordinary arrows and the fixed word-parity normalization are allowed.
    """
    assert len(direction) == len(centers) == 5 and any(direction)
    normal = {(0, h): ([()] if h == 0 else []) for h in range(4)}
    arrows = {}
    relations_over_center = {}

    def multiply(row, letter):
        out = {}
        for word, coefficient in row.items():
            key = (word, letter)
            if key in arrows:
                add_multiple(out, arrows[key], coefficient)
            else:
                add_multiple(out, {word + (letter,): one}, coefficient)
        return out

    def centered_word(word):
        row = {(): one}
        for letter in word:
            row = multiply(row, letter)
        return row

    for n in range(1, limit + 1):
        dimensions = []
        for h in range(4):
            rows = []
            if n >= 2:
                for color, relation in relations(n, h):
                    for w in normal[n - 2, h ^ color]:
                        row = {}
                        for (a, b), coefficient in relation.items():
                            term = multiply(multiply({w: one}, a), b)
                            add_multiple(row, term, coefficient)
                        rows.append(row)
            if (n, h) == (4, 0):
                for central, value in zip(centers, direction, strict=True):
                    row = {}
                    for w, coefficient in central.items():
                        add_multiple(row, centered_word(w), coefficient)
                    add_multiple(row, {(): one}, -value)
                    rows.append(row)
            words = [w + (a,) for a in range(4) for w in normal[n - 1, h ^ a]]
            pivots, residuals = reduce_words(rows, n)
            normal[n, h] = [w for w in words if w not in pivots]
            for w in words:
                if w in pivots:
                    projection = {v: -z for v, z in pivots[w].items() if v != w}
                else:
                    projection = {w: one}
                arrows[w[:-1], w[-1]] = projection
            if n >= 8:
                assert not normal[n, h], "five-center quotient survives in length eight"
            for a in range(4):
                residuals.extend(
                    multiply(row, a)
                    for row in relations_over_center.get((n - 1, h ^ a), [])
                )
            residuals.extend(relations_over_center.get((n - 4, h), []))
            basis = select(residuals)
            relations_over_center[n, h] = basis
            free_dimension = sum(len(normal[d, h]) for d in range(n % 4, n + 1, 4))
            dimensions.append(free_dimension - len(basis))
        yield n, dimensions
