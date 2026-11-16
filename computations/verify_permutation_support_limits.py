#!/usr/bin/env python3
"""Exhaust all coordinate faces of Birkhoff B4: concise limits are trivial."""
import argparse
import json
from itertools import permutations
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    perms=list(permutations(range(4)))
    edges=[sum(1<<(4*i+p[i]) for i in range(4)) for p in perms]
    families=set()
    for graph in range(1<<16):
        family=sum(1<<j for j,matching in enumerate(edges) if matching&graph==matching)
        if family:families.add(family)
    histogram={};concise=[]
    for family in families:
        allowed=0
        for j,matching in enumerate(edges):
            if family>>j&1:allowed|=matching
        ranks=[sum((allowed>>(4*i+a))&1 for a in range(4)) for i in range(4)]
        # Distinct singleton slices have disjoint supports: an assignment of
        # the other three labels uniquely determines the missing fourth label.
        minimum=min(ranks);histogram[minimum]=histogram.get(minimum,0)+1
        if minimum==4:concise.append(family)
    assert concise==[(1<<24)-1]
    certificate={'result':'PASS','arithmetic':'exhaustive finite combinatorics, exact integers',
                 'bipartite_graphs_examined':1<<16,'nonempty_Birkhoff_B4_faces':len(families),
                 'faces_by_minimum_singleton_rank':histogram,
                 'all_singletons_rank_four_faces':1,
                 'application':'For any tensor supported exactly on all 24 permutations with nonzero coefficients, an independent diagonal one-parameter limit preserving all singleton ranks is the original projective tensor.',
                 'limits':['Fixed-coordinate diagonal one-parameter limits only; arbitrary basis changes or moving source tensors are not covered.','No relative closedness theorem for the symmetric CI union.']}
    if args.output:args.output.write_text(json.dumps(certificate,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(certificate,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
