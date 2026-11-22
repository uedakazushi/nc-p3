"""Compare all four marked point configurations at the rational full-rank window."""
import sympy as sp
from verify_klein_elliptic_completion import rational_example
from verify_klein_generic_isomorphisms import kernel,framed,number_of_maps
from verify_klein_two_doublets import modular_points


def partner(points,matrix,transpose_relations,prime):
    def residue(x):
        return int(x.p)*pow(int(x.q),-1,prime)%prime
    relations=[sp.Matrix(4,4,col).applyfunc(residue) for col in matrix.columnspace()]
    if transpose_relations:
        relations=[r.T for r in relations]
    result=[kernel(sp.Matrix.vstack(*(sp.Matrix(1,4,q)*r for r in relations)),prime) for q in points]
    assert len(set(result))==len(result)==20
    return result


if __name__=='__main__':
    prime=2657
    w,_,n=rational_example()
    first=modular_points(w,prime)
    second=partner(first,w,False,prime)
    third=partner(second,n,True,prime)
    fourth=partner(third,w.T,False,prime)
    assert partner(fourth,n.T,True,prime)==first
    flip=[4*j+i for i in range(4) for j in range(4)]
    for points,mat in [(second,n[flip,:]),(third,w.T),(fourth,n.T[flip,:])]:
        assert set(modular_points(mat,prime))==set(points)
    configs=[framed(q,prime) for q in [first,second,third,fourth]]
    assert [[number_of_maps(x,y,prime) for y in configs] for x in configs]==(4*sp.eye(4)).tolist()
    print('Four reduced configurations, identity return, and exactly four automorphisms of each.')
