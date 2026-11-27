"""Two central-cut charts certify every point of the three-dimensional rectangle."""
from flint_rational_function import rational_function_field
from verify_parallel_third_rectangle import s,third_rectangle
from verify_rectangular_surface_cut_charts import chart,central_cut_data


def check():
    a,b,h = s.symbols('a y h')
    field = rational_function_field(('a','y','h'))
    factor = field.from_sympy(a*b*h*(b-1)*(h-1)).n
    base = third_rectangle(a,h/a,b)
    assert chart(2,central_cut_data(base,field,factor))==[]
    shifts = (2,2,0,0)
    moved = {tuple(x^d for x,d in zip(k,shifts)):f for k,f in base.items()}
    singular = chart(1,central_cut_data(moved,field,factor))
    assert singular==['a*y*h + a*y - 6*y*h + 4*h']
    f = a*b*h+a*b-6*b*h+4*h
    other = s.cancel((b-1)*f.xreplace({a:-a,b:b/(b-1)}))
    assert s.factor(f+other+8*h*b)==0
    print('PASS: uniform even cut and two complementary odd cuts on the entire third rectangle')


if __name__=='__main__':
    check()
