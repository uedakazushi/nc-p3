"""A common fourth-degree obstruction on the horizontal anti-diagonal."""
from verify_parallel_excluded_shapes import path_dimension
from verify_parallel_invertible_rectangles import K,completed,s


def check():
    b,x,y = s.symbols('b x y')
    for t in (b,-b):
        c = -b
        for shape in (76,112):
            A,Y,D = ((x,1,c*(1-x)/t**2) if shape==76 else
                     (t*x,y,(c*x*(1-y)+b*(1-c)+c*y)/t))
            raw = completed(shape,t,b,c,A,1,Y,D)
            assert [path_dimension(raw,0,4,h) for h in K]==[10,10,8,8]


if __name__=='__main__':
    check()
    print('PASS: H4 equals 36 generically on all four anti-diagonal loci')
