#!/usr/bin/env python3
"""Finite-Frobenius boundary of the symmetric chart with first defect in degree five."""
import argparse
import json
from pathlib import Path
from itertools import combinations_with_replacement,product
from fractions import Fraction as F
import sympy as s
from frobenius_tensor import (digits,flatten,reconstructed_windows,orthogonal_windows,
                             insertion_columns,nullspace_sparse,exact_rank)
from tensor_symmetrization import symmetrizer_space,symmetrize
from verify_quartic_obstruction import quintic_d3_rank,relation_matrix


def cubic_frame():
    t=s.symbols('t');monomials=list(combinations_with_replacement(range(4),3))
    lookup={m:i for i,m in enumerate(monomials)};rows=[]
    for a in range(4):
        for i in range(3):
            row=[0]*20
            row[lookup[tuple(sorted((i,i,a)))]]+=1
            row[lookup[tuple(sorted((i,3,a)))]]-=1
            rows.append(row)
        row=[0]*20
        for i in range(3):row[lookup[tuple(sorted((i,3,a)))]]+=1
        row[lookup[tuple(sorted((3,3,a)))]]+=t
        rows.append(row)
    equations=s.Matrix(rows)
    free=[(0,1,2),(1,3,3),(2,3,3),(3,3,3)]
    normalization=s.Matrix(4,20,lambda i,j:int(monomials[j]==free[i]))
    square=equations.col_join(normalization)
    assert s.expand(square.det())==-2
    basis=square.inv()*s.zeros(16,4).col_join(s.eye(4))
    assert equations*basis==s.zeros(16,4)
    assert normalization*basis==s.eye(4)
    return t,monomials,basis


def verify_involution():
    t,monomials,basis=cubic_frame();lookup={m:i for i,m in enumerate(monomials)}
    g=s.Matrix([[-1,0,0,1],[0,-1,0,1],[0,0,-1,1],[0,0,0,1]])
    changed=s.zeros(20,4)
    for row,(a,b,c) in enumerate(monomials):
        for u,v,w in product(range(4),repeat=3):
            coefficient=g[a,u]*g[b,v]*g[c,w]
            if coefficient:
                for j in range(4):changed[row,j]+=coefficient*basis[lookup[tuple(sorted((u,v,w)))],j]
    free=[(0,1,2),(1,3,3),(2,3,3),(3,3,3)]
    frame=changed[[lookup[m] for m in free],:]
    assert frame==s.Matrix([[-1,0,0,(t*t+3*t+2)/2],[0,-1,0,1],[0,0,-1,1],[0,0,0,1]])
    assert frame.det()==-1
    assert (changed-basis.subs(t,-3-t)*frame).applyfunc(s.expand)==s.zeros(20,4)


def tensor_at(parameter):
    t,monomials,basis=cubic_frame();lookup={m:i for i,m in enumerate(monomials)}
    w={}
    for encoded in range(256):
        a,b,j,c=digits(encoded)
        coefficient=F(basis[lookup[tuple(sorted((a,b,c)))],j].subs(t,parameter))
        if coefficient:w[encoded]=coefficient
    return w


def square_relation_cycles(windows):
    def quotient(start,degree):
        basis=nullspace_sparse(insertion_columns(windows,start,degree),4**degree)
        return s.Matrix(len(basis),4**degree,lambda i,j:basis[i].get(j,0))
    def section(quotient):
        pivots=quotient.rref()[1]
        result=s.zeros(quotient.cols,quotient.rows)
        inverse=quotient[:,list(pivots)].inv()
        for row,pivot in enumerate(pivots):result[pivot,:]=inverse[row,:]
        assert quotient*result==s.eye(quotient.rows)
        return result
    q3=quotient(2,3);s3=section(q3);q4=quotient(1,4)
    multiplication=q4[:,192:256]*s3
    annihilator=multiplication.nullspace()
    assert multiplication.rank()==14 and len(annihilator)==6
    square={ (3,3):F(1) }
    assert windows[0][5]==square
    d2=s.kronecker_product(s.eye(4),q4)*s.kronecker_product(relation_matrix(windows[0]),s3)
    first=s.kronecker_product(relation_matrix(windows[0]),s.eye(4))
    second=s.kronecker_product(s.eye(4),relation_matrix(windows[1]))
    syzygies=first.row_join(-second).nullspace()
    coordinates=s.Matrix.hstack(*(v[:24,:] for v in syzygies))
    q2=quotient(3,2);s2=section(q2)
    d3=s.kronecker_product(s.eye(6),q3)*s.kronecker_product(coordinates,s2)
    cycles=s.zeros(120,6)
    for j,col in enumerate(annihilator):cycles[100:120,j]=col
    assert d2*cycles==s.zeros(140,6) and d2*d3==s.zeros(140,40)
    assert d2.rank()==77 and d3.rank()==36
    assert d3.row_join(cycles).rank()==42
    # A mixed cycle completes the seven-dimensional quintic H2.  Relation
    # tensors here are specified in monomial coordinates, independently of
    # which six relation vectors reconstructed_windows chooses.
    relation_basis = relation_matrix(windows[0])
    def relation_coordinates(pairs):
        v = s.Matrix(16, 1, lambda i, _: pairs.get(divmod(i, 4), 0))
        solution, params = relation_basis.gauss_jordan_solve(v)
        assert params.rows == 0
        return solution
    def tensor_cycle(pairs, tail):
        j = 16 * tail[0] + 4 * tail[1] + tail[2]
        return s.kronecker_product(relation_coordinates(pairs), q3[:, j])
    r01 = {(0, 1): 1, (1, 0): 1}
    r13 = {(1, 1): 1, (1, 3): 1, (2, 2): -1, (2, 3): -1,
           (3, 1): 1, (3, 2): -1}
    mixed = (2 * tensor_cycle(square, (3, 3, 2))
             + 4 * tensor_cycle(square, (1, 3, 2))
             + 12 * tensor_cycle(square, (2, 2, 1))
             + 8 * tensor_cycle(square, (2, 3, 1))
             + 4 * tensor_cycle(square, (3, 3, 1))
             - 3 * tensor_cycle(r01, (2, 1, 3))
             - 6 * tensor_cycle(r13, (1, 2, 3)))
    assert d2 * mixed == s.zeros(140, 1)
    full = d3.row_join(cycles).row_join(mixed)
    assert exact_rank(({i: F(full[i, j]) for i in range(120) if full[i, j]}
                       for j in range(full.cols))) == 43
    # Since dim H2=7, the six square classes and this mixed class form a basis.
    # One especially small nonzero class is square relation tensor word 303.
    word=s.zeros(64,1);word[51]=1
    assert multiplication*q3*word==s.zeros(35,1)
    simple=s.zeros(120,1);simple[100:120,0]=q3*word
    assert d3.row_join(simple).rank()==37
    return {'left_x3_multiplication_rank':14,'annihilator_dimension':6,
            'd2_rank':77,'d3_rank':36,'d3_plus_square_cycles_rank':42,
            'simple_nonboundary_cycle':'(x3 tensor x3) tensor [303]',
            'remaining_H2_dimension_after_square_cycles':1,
            'mixed_cycle': '2*[33]x[332]+4*[33]x[132]+12*[33]x[221]+8*[33]x[231]+4*[33]x[331]-3*r01x[213]-6*r13x[123]',
            'd2_mixed_zero':True,'d3_plus_square_plus_mixed_rank':43,
            'remaining_H2_after_all_seven_cycles':0}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    verify_involution()
    w=tensor_at(0)
    assert len(w)==83
    assert [flatten(w,(i,)) for i in range(4)]==[4]*4
    assert [flatten(w,(0,j)) for j in (1,2,3)]==[6]*3
    windows=reconstructed_windows(w);dual=orthogonal_windows(windows)
    intersections=[[len(nullspace_sparse(insertion_columns(dual,i,n),4**n)) for n in (3,4)] for i in range(4)]
    assert intersections==[[4,1]]*4
    assert all(len(nullspace_sparse(insertion_columns(dual,i,5),1024))==0 for i in range(4))
    spaces=[symmetrizer_space(w,j) for j in (1,2,3)]
    assert spaces[0]==spaces[2]==[s.eye(4)]
    assert spaces[1]==[s.diag(0,0,0,1)]
    assert symmetrize(w)[0] is None
    assert symmetrize(w,-1)[0] is None
    quartic=[256-exact_rank(insertion_columns(windows,i,4)) for i in range(4)]
    quintic=[1024-exact_rank(insertion_columns(windows,i,5)) for i in range(4)]
    assert quartic==[35]*4 and quintic==[63,56,56,56]
    d3_ranks=[quintic_d3_rank(windows,i) for i in range(4)]
    assert d3_ranks==[36]*4
    square_cycles=square_relation_cycles(windows)
    generic=tensor_at(1)
    assert symmetrize(generic)[0] is not None
    bad=tensor_at(-1);bad_windows=reconstructed_windows(bad)
    bad_dual=orthogonal_windows(bad_windows)
    bad_intersections=[[len(nullspace_sparse(insertion_columns(bad_dual,i,n),4**n)) for n in (3,4)] for i in range(4)]
    assert bad_intersections==[[6,12],[12,12],[6,4],[4,4]]
    certificate={'result':'PASS','arithmetic':'exact rational polynomial and linear algebra',
        'quadrics':['e0^2-e0*e3','e1^2-e1*e3','e2^2-e2*e3','e3*(e0+e1+e2+t*e3)'],
        'cubic_normalization_coordinates':['012','133','233','333'],
        'cubic_normalization_determinant':-2,
        'symmetric_CI_parameters':'t not in {0,-1,-2,-3}',
        'whole_family_AS_parameters':'exactly t not in {0,-1,-2,-3}',
        'whole_family_finite_Frobenius_parameters':'exactly t not in {-1,-2}',
        'indexed_involution':'t -> -3-t; polynomial frame determinant -1',
        'intersection_dimensions_at_t_minus_one':bad_intersections,
        'special_tensor_terms':83,'singleton_ranks':[4]*4,'pair_cut_ranks':[6]*3,
        'intersection_dimensions_W3_W4':intersections,'dual_degree_five_dimensions':[0]*4,
        'symmetrizers':['k I','k E33','k I'],'alternating_membership':False,
        'primal_degree_four_dimensions':quartic,'primal_degree_five_dimensions':quintic,
        'quintic_d3_ranks':d3_ranks,'quintic_H2_dimensions':[7,0,0,0],
        'quintic_H3_dimensions':[0]*4,
        'square_relation_cycles':square_cycles,
        'primal_relation_insertion_rank_degree_five':[1024-d for d in quintic],
        'conclusions':['The symmetric chart is not relatively closed inside the finite-Frobenius low-rank candidate locus.',
                       'Finite Frobenius plus all four quartic Koszul filters does not imply Koszulity.',
                       'The special point is not a standard AS object and is not a LOW-UNCLASSIFIED AS witness.'],
        'limits':['No full AS classification or exhaustive analysis of the symmetric boundary.','Generic CI fibers are regular by the existing theorem; special-fiber non-Koszulity is certified by a degree-five Hilbert defect.']}
    if args.output:args.output.write_text(json.dumps(certificate,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(certificate,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
