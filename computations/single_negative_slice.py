"""Rank factorization for a rational tensor alternating in its first two slots."""
import sympy as s


def linear_complex(w):
    assert all(w.get(64*a+16*b+4*c+d,0)==-w.get(64*b+16*a+4*c+d,0)
               for a in range(4) for b in range(4) for c in range(4) for d in range(4))
    matrix=s.Matrix(16,16,lambda i,j:w.get(64*(i//4)+16*(j//4)+4*(i%4)+j%4,0))
    columns=matrix.rref()[1];assert len(columns)==6
    first=matrix[:,list(columns)];rows=first.T.rref()[1];assert len(rows)==6
    second=first[list(rows),:].inv()*matrix[list(rows),:]
    assert first*second==matrix
    x=s.symbols('x0:4');y=s.symbols('y0:4')
    P=s.Matrix(4,6,lambda c,r:sum(x[a]*first[4*a+c,r] for a in range(4)))
    Q=s.Matrix(6,4,lambda r,d:sum(x[b]*second[r,4*b+d] for b in range(4)))
    assert (P*Q).applyfunc(s.expand)==s.zeros(4)
    assert (P*Q.subs(dict(zip(x,y)))+P.subs(dict(zip(x,y)))*Q).applyfunc(s.expand)==s.zeros(4)
    return P,Q,{'pivot_columns':list(columns),'pivot_rows':list(rows)}
