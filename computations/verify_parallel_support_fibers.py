"""All regularity decisions on seven sparse supports of the two invertible rectangles."""
from itertools import product

from verify_parallel_excluded_shapes import path_dimension
from verify_parallel_invertible_rectangles import (
    A,B,C,D,K,T,X,BETA,GAMMA,bit,candidates,character,completed,matrix,quadric,s,smith,template,top)
from verify_klein_three_lines import twist


def vertical(shape,t):
    a,y,d = s.symbols('a y d')
    if shape==76:
        w = completed(shape,t,1,1,A,B,C,D)
        det = A*C-B*D
        forms = [A,B,A,D,B,C,D,C,C-D,A-B,C-D,C,A-B,B,C,B]
    else:
        w = completed(shape,t,1,1,t*a,1,y,t*d)
        det = a*y-d
        cross = a+y-d-1
        forms = [a,1,a,d,1,y,d,y,cross,a-1,cross,y-1,a-1,1,y-1,1]
    for key,x in template(shape).items():
        if x in X:
            i = X.index(x)
            assert s.factor(w[key]*(det if i>=8 else 1)/forms[i]) in (-1,1)


def horizontal(shape):
    a,y,b,c,t = s.symbols('a y b c t')
    m = b+c-b*c
    first = a if shape==76 else t*a
    last = s.S.One if shape==76 else y
    d = s.solve(quadric(shape,t,b,c,first,1,last,D),D)[0]
    vals = top(shape,t,b,c,first,1,last,d)
    w = template(shape)
    subs = dict(zip(X[:8],vals))|{T:t}
    for color in (0,1):
        src = matrix(w,0,color)
        mat = src.subs(subs)
        top_block = mat[:2,:2]
        adj = mat[2:,:2]*top_block.adjugate()*mat[:2,2:]
        for x,v in zip(src[2:,2:],adj):
            subs[x] = s.cancel(v/top_block.det())
    def reduce(f):
        return s.factor(s.rem(s.expand(f),t*t-m,t))
    E = (a*(b+c)-b*c) if shape==76 else (a*y*(b+c)-b*c*(a+y)+b*b*(c-1))
    allowed = s.factor_list(t*b*c*m*E)[1]
    allowed = {p for p,_ in allowed}
    def is_unit(f):
        return f!=0 and all(p in allowed for p,_ in s.factor_list(f)[1])

    if shape==76:
        f = c*a+b*(1-c)
        forms = (f,1,a,b-a,1,1,a-1,1,a,a-b,f,1,a-1,1,1,1)
    else:
        n = -c*(a-b)*(y-b)
        forms = (c*a+b*(1-c),1,a,n+b*b,1,y,c*(a-1)*(y-1)-m,c*y+b*(1-c),
                 a*y-b,a-b,n+b*(1+b-a-y),y-1,a-1,1,y-b,1)
    for x,form in zip(X,forms):
        n,den = map(reduce,s.fraction(subs[x]))
        q = s.factor(n/form)
        assert is_unit(q) and is_unit(den),(x,q,den)
    delta = s.cancel(matrix(w,0,0).subs(subs)[:2,:2].det()/E)
    n,den = map(reduce,s.fraction(delta))
    assert is_unit(n) and is_unit(den)
    if shape==112:
        n3,n6,n10 = (forms[i] for i in (3,6,10))
        assert s.factor((n10-n3)+b*(a+y-1))==0
        assert s.factor(n10.subs(a,b/y)-b*(y-b)*(c-1)*(y-1)/y)==0
        assert s.factor(n6.subs({c:1,a:b/y})-(b-b/y-y))==0


def sign_change(shape):
    w = completed(shape,T,BETA,GAMMA,A,B,C,D)
    subs = {T:-T} if shape==76 else {T:-T,A:-A,D:-D}
    ratios = [s.factor(v.xreplace(subs)/v)*twist(k) for k,v in w.items()]
    assert set(ratios)<={1,-1}
    gauge = s.Matrix([[int(k[p]==i) for p,i in product(K,repeat=2)] for k in w])
    diag,left,_ = smith(gauge)
    rank = sum(bool(x) for x in diag.diagonal())
    assert all(character(row,ratios)==1 for row in left.tolist()[rank:])


def check():
    empty = {82:(9,12),118:(3,6,9,11,12,14),124:(9,11,12,14),
             129:(5,7,9,12),134:(3,6,9,12),139:(9,12),142:(3,6,8,10),
             143:(3,8,10),144:(3,10),145:(8,10),146:(10,),147:(),85:()}
    for index,required in empty.items():
        w = template(76 if index<86 else 112)
        actual = sorted(X.index(x) for k,x in w.items() if x in X and not(candidates(1)[index]&bit(k)))
        assert tuple(actual)==required
    for shape in (76,112):
        horizontal(shape)
        sign_change(shape)
        for t in (1,-1):
            vertical(shape,t)
    c,y = s.symbols('c y')
    for t in (1,-1):
        w = completed(112,t,1,c,t,1,y,t)
        assert [path_dimension(w,1,7,h) for h in K]==[34,33,34,33]
    print('PASS: horizontal and vertical supports, cocycle flip, two generic septic obstructions')


if __name__=='__main__':
    check()
