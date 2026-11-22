"""Compute the generic dual paths and the exact Fitting support over F_101[t] with Singular."""

import os
import re
import subprocess

from flint import nmod_poly
from rank_fourteen_function_field import ONE

PRELUDE = """
option(redSB);
LIB "matrix.lib";
proc kernelstep(matrix A, matrix Y, string label) {
    module columnBasis=std(module(A));
    int r=size(columnBasis);
    matrix C[nrows(A)][r]=columnBasis;
    module rowBasis=std(module(transpose(C)));
    if(size(rowBasis)!=r){ERROR("unexpected number of row generators");}
    matrix D[r][r]=rowBasis;
    poly defect=det(D);
    print("FITTING "+label+" "+string(defect));
    module G=std(syz(module(A)));
    print("DIM "+label+" "+string(size(G)));
    if(size(G)==0){matrix dummy[1][1];return(list(dummy,dummy));}
    matrix J[ncols(A)][size(G)]=G;
    matrix U=Y*J;
    matrix M=lift(G,module(U));
    if(size(ideal(J*M-U))!=0){ERROR("coefficient action does not lift");}
    return(list(J,M));
}
proc actmat(poly a,poly b,matrix Y){return(a*unitmat(nrows(Y))+b*Y);}
proc part(matrix J,int first,int len){return(submat(J,first..(first+len-1),1..ncols(J)));}
int k,m;
poly cubic=2*(t-1)*(t+1)*(t+3);
matrix M0h0[2][2]=0,cubic,1,0;
"""


def polynomial(value):
    return (
        "(" + ("+".join(f"({c})*t^{j}" for j, c in enumerate(value) if c) or "0") + ")"
    )


def integral_rules(rules):
    result = []
    for h, row in rules:
        common = ONE
        for value in row.values():
            common *= value.denominator // common.gcd(value.denominator)
        entries = {}
        for word, value in row.items():
            scale = common // value.denominator
            a, b, c, d = (polynomial(x * scale) for x in value.coordinates)
            entries[word] = (f"({a}+i*{b})", f"({c}+i*{d})")
        result.append((h, entries))
    return result


def program(windows, sections, start, imaginary, prime=101):
    lines = [f"ring S={prime},(t),(dp,C); number i={imaginary};", PRELUDE]
    emit = lines.append

    def action(n, h):
        return f"M{n}h{h}"

    def arrow(n, h, a):
        return f"L{n}h{h}a{a}"

    def dimension(n, h):
        if n == 12 and h:
            return "0"
        return f"ncols({action(n, h)})"

    def total(dimensions):
        return "+".join(dimensions) or "0"

    for h in range(4):
        emit(f"matrix {action(1, h)}=M0h0; matrix {arrow(1, h, h)}=unitmat(2);")
    for n in range(2, 14):
        for h in range(4):
            widths = [dimension(n - 1, h ^ b) for b in range(4)]
            rows = []
            for length, rules in (
                (2, windows[(start + n - 2) % 4]),
                (4, sections[start] if n == 4 else []),
            ):
                prefix = n - length
                for g, row in rules:
                    char = h ^ g
                    if prefix == 0 and char:
                        continue
                    blocks = [[] for _ in range(4)]
                    for word, (u, v) in row.items():
                        if widths[word[-1]] == "0":
                            continue
                        parts = [f"actmat({u},{v},{action(prefix, char)})"]
                        current = char
                        for pos, a in enumerate(word[:-1], prefix + 1):
                            current ^= a
                            parts.append(arrow(pos, current, a))
                        blocks[word[-1]].append("*".join(parts))
                    rows.append((dimension(prefix, char), blocks))
            heights = [height for height, _ in rows]
            emit(f"matrix A[{total(heights)}][{total(widths)}];")
            for j, (_, blocks) in enumerate(rows):
                for b, terms in enumerate(blocks):
                    if terms:
                        emit("matrix V=" + "+".join(terms) + ";")
                        emit(
                            "for(k=1;k<=nrows(V);k++){for(m=1;m<=ncols(V);m++){"
                            f"A[({total(heights[:j])})+k,({total(widths[:b])})+m]=V[k,m];"
                            + "}} kill V;"
                        )
            previous = [action(n - 1, h ^ b) for b in range(4) if widths[b] != "0"]
            direct_sum = (
                previous[0]
                if len(previous) == 1
                else "dsum(" + ",".join(previous) + ")"
            )
            emit(f'list rec=kernelstep(A,{direct_sum},"{n} {h}");')
            if n < 12 or (n == 12 and h == 0):
                emit(f"matrix {action(n, h)}=rec[2];")
                for b in range(4):
                    emit(
                        f"matrix {arrow(n, h, b)}=part(rec[1],1+({total(widths[:b])}),{widths[b]});"
                    )
            emit("kill rec,A;")
    return "\n".join(lines + ["quit;"])


def parse_polynomial(expression, prime=101):
    coefficients = {}
    terms = re.findall(r"[+-]?(?:[0-9]*t[0-9]*|[0-9]+)", expression)
    assert "".join(terms) == expression
    for term in terms:
        sign = -1 if term.startswith("-") else 1
        coeff, x, exp = term.lstrip("+-").partition("t")
        degree = int(exp or 1) if x else 0
        assert degree not in coefficients
        coefficients[degree] = sign * int(coeff or 1)
    return nmod_poly(
        [coefficients.get(j, 0) for j in range(max(coefficients) + 1)], prime
    )


def fitting_support(windows, sections, start, imaginary, prime=101):
    code = program(windows, sections, start, imaginary, prime)
    run = subprocess.run(
        [os.environ.get("SINGULAR_BIN", "Singular"), "-q"],
        input=code,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=True,
    )
    if re.search(r"^\s*\?", run.stdout, re.MULTILINE):
        raise ValueError(run.stdout[-6000:])
    bad = set()
    sizes = {}
    seen = set()
    for line in run.stdout.splitlines():
        if line.startswith(("DIM ", "FITTING ")):
            tag, n, h, datum = line.split()
            key = int(n), int(h)
            if tag == "DIM":
                sizes[key] = int(datum)
            else:
                assert key not in seen
                seen.add(key)
                determinant = parse_polynomial(datum, prime)
                assert determinant
                bad.update(tuple(map(int, f)) for f, e in determinant.factor()[1])
    assert seen == set(sizes) == {(n, h) for n in range(2, 14) for h in range(4)}
    profile = (10, 20, 31, 40, 44, 40, 31, 20, 10, 4, 1, 0)
    for n, expected in enumerate(profile, 2):
        assert sum(sizes[n, h] for h in range(4)) == 2 * expected
    assert [sizes[12, h] for h in range(4)] == [2, 0, 0, 0]
    return bad
