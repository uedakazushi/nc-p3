"""Exact rational tools for four-window top tensors; no Koszul certification."""
from fractions import Fraction as F
from verify_quartic_obstruction import (exact_rank, insertion_columns,
                                        relation_matrix)


def orthogonal_windows(windows):
    result=[]
    for window in windows:
        r=relation_matrix(window)
        result.append([{divmod(j,4):F(x) for j,x in enumerate(v) if x}
                       for v in r.T.nullspace()])
    return result


def nullspace_sparse(rows, size):
    """Nullspace of sparse row equations over Q by triangular back substitution."""
    pivots={}
    for row in rows:
        row=dict(row)
        while row:
            i=min(row)
            if i not in pivots:
                a=row[i];pivots[i]={j:x/a for j,x in row.items()};break
            a=row[i]
            for j,x in pivots[i].items():
                val=row.get(j,F(0))-a*x
                if val:row[j]=val
                else:row.pop(j,None)
    free=sorted(set(range(size))-pivots.keys());out=[]
    for j in free:
        v={j:F(1)}
        for i,row in sorted(pivots.items(),reverse=True):
            val=-sum(x*v.get(k,F(0)) for k,x in row.items() if k!=i)
            if val:v[i]=val
        out.append(v)
    return out


def digits(i,n=4):
    return tuple((i//4**(n-1-j))%4 for j in range(n))


def flatten_columns(tensor, left):
    right=tuple(i for i in range(4) if i not in left)
    cols={}
    for i,x in tensor.items():
        word=digits(i);a=b=0
        for j in left:a=4*a+word[j]
        for j in right:b=4*b+word[j]
        cols.setdefault(b,{})[a]=x
    return list(cols.values())



def flatten(tensor, left):
    return exact_rank(flatten_columns(tensor, left))


def normalized(tensor):
    scale = tensor[min(tensor)]
    return {i: x / scale for i, x in tensor.items() if x}


def top_tensor(windows, start=0):
    dual = orthogonal_windows(windows)
    basis = nullspace_sparse(insertion_columns(dual, start, 4), 256)
    if len(basis) != 1:
        raise ValueError(f"top intersection dimension {len(basis)} != 1")
    return normalized(basis[0])


def independent_columns(columns):
    pivots = {}
    for column in columns:
        col = dict(column)
        while col:
            i = min(col)
            if i not in pivots:
                scale = col[i]
                pivots[i] = {j: x / scale for j, x in col.items()}
                break
            scale = col[i]
            for j, x in pivots[i].items():
                value = col.get(j, F(0)) - scale * x
                if value:
                    col[j] = value
                else:
                    col.pop(j, None)
    return list(pivots.values())


def reconstructed_windows(tensor):
    return [[{divmod(j, 4): x for j, x in col.items()}
             for col in independent_columns(
                 flatten_columns(tensor, (i, (i + 1) % 4)))]
            for i in range(4)]


def rotate(tensor, start):
    result = {}
    for index, value in tensor.items():
        word = digits(index)
        shifted = word[start:] + word[:start]
        new = 0
        for letter in shifted:
            new = 4 * new + letter
        result[new] = value
    return normalized(result)
