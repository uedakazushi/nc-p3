#!/usr/bin/env python3
"""Exact linear-complex control at the alternating determinant tensor."""
import argparse,json
from itertools import permutations,combinations
from pathlib import Path
import sympy as s
from single_negative_slice import linear_complex


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args()
    w={}
    for word in permutations(range(4)):
        a,b,c,d=word;inversions=sum(word[i]>word[j] for i in range(4) for j in range(i+1,4))
        w[64*a+16*b+4*c+d]=(-1)**inversions
    P,Q,pivots=linear_complex(w);x=s.symbols('x0:4')
    assert all(s.expand(P[:,list(cols)].det())==0 for cols in combinations(range(6),4))
    assert all(s.expand(Q[list(rows),:].det())==0 for rows in combinations(range(6),4))
    sample=dict(zip(x,(1,2,3,5)));assert P.subs(sample).rank()==Q.subs(sample).rank()==3
    report={'result':'PASS','arithmetic':'exact QQ polynomial identities','tensor':'alternating volume in dimension4',
            'opposite_rank':6,'factorization':pivots,'P_size':[4,6],'Q_size':[6,4],
            'P_x_Q_y_plus_P_y_Q_x_zero':True,'all_four_minors_zero':True,
            'generic_ranks':[3,3],'rank_three_witness':[1,2,3,5],
            'P':[[str(P[i,j]) for j in range(6)] for i in range(4)],
            'Q':[[str(Q[i,j]) for j in range(4)] for i in range(6)],
            'limits':['This checks the polynomial-class control; it does not prove that every rank-(3,3) slice is polynomial.']}
    if args.output:args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
