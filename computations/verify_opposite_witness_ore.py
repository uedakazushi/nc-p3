#!/usr/bin/env python3
"""Opposite-positive witnesses on the already classified signed Ore chart."""
import argparse,json
from itertools import product
from pathlib import Path
from fractions import Fraction as F
from verify_indexed_anticommutative_ore import H,G,normalizer,dual_rules,index,window
from tensor_symmetrization import symmetrizer_space,invertible_witness
from frobenius_tensor import flatten,top_tensor
import sympy as s


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args()
    nf=normalizer(dual_rules)
    tensor={word:s.expand(nf(0,word).get((0,1,2,3),0)) for word in product(range(4),repeat=4)}
    assert all(tensor[a,b,c,d]==tensor[c,b,a,d] for a,b,c,d in tensor)
    assert all(tensor[a,b,c,d]==tensor[a,d,c,b] for a,b,c,d in tensor)
    fixtures=[]
    for h,g in [([0]*9,[1,0,0]),(list(range(1,10)),[1,2,3])]:
        sub=dict(zip(H+G,h+g));w={index(word):F(value.subs(sub)) for word,value in tensor.items() if value.subs(sub)}
        independent=top_tensor([window(h,g,i) for i in range(4)]);scale=independent[27]
        assert {i:x/scale for i,x in independent.items()}==w
        positive=[symmetrizer_space(w,j) for j in (1,2,3)]
        negative=[symmetrizer_space(w,j,-1) for j in (1,2,3)]
        flags=[invertible_witness(space) is not None for space in positive]
        assert flags==[False,True,False] and all(invertible_witness(space) is None for space in negative)
        assert flatten(w,(0,2))==8
        fixtures.append({'H':h,'G':g,'rho':8,'positive_dimensions':[len(space) for space in positive],
                         'negative_dimensions':[len(space) for space in negative],'positive_invertible':flags})
    result={'result':'PASS','arithmetic':'exact polynomial identities in all12 parameters, plus rational fixtures',
            'swap02_identity_all_parameters':True,'swap13_identity_all_parameters':True,'fixtures':fixtures,
            'limits':['AS regularity and the rank formula use the existing signed-Ore theorem.','This does not cover every single-positive-witness AS tensor.']}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
