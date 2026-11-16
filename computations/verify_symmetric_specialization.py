#!/usr/bin/env python3
"""Exact regular symmetric CI specialization; individual orbit not closed."""
import argparse
import json
from pathlib import Path
from itertools import product,combinations
from fractions import Fraction as F
from math import factorial
import sympy as s
from frobenius_tensor import top_tensor,flatten
from tensor_symmetrization import symmetrize


def index(word):
    result=0
    for a in word:result=4*result+a
    return result


def relations(parameter):
    result=[]
    for a,b in combinations(range(4),2):
        relation={(b,a):s.Integer(1),(a,b):s.Integer(1)}
        if (a,b)==(0,1):relation[3,3]=-parameter**2
        result.append(relation)
    return result


def matrix(relations):
    return s.Matrix(16,6,lambda i,j:relations[j].get(divmod(i,4),0))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    t=s.symbols('t',nonzero=True);z=s.symbols('z0:4')
    tensor={}
    for word in product(range(4),repeat=4):
        if len(set(word))==4:tensor[index(word)]=s.Integer(1)
        elif sorted(word)==[2,3,3,3]:tensor[index(word)]=-t**2
    form=s.expand(sum(value*s.prod(z[a] for a in ((i//64)%4,(i//16)%4,(i//4)%4,i%4))
                      for i,value in tensor.items())/factorial(4))
    assert form==z[0]*z[1]*z[2]*z[3]-t**2*z[2]*z[3]**3/6
    change=s.diag(1/t,1/t,t**2,1)
    assert change.det()==1
    source=matrix(relations(s.Integer(1)));target=matrix(relations(t))
    scale=s.diag(*(change[a,a]*change[b,b] for a,b in combinations(range(4),2)))
    assert (s.kronecker_product(change,change)*source-target*scale).applyfunc(s.cancel)==s.zeros(16,6)
    for encoded,value in tensor.items():
        word=((encoded//64)%4,(encoded//16)%4,(encoded//4)%4,encoded%4)
        old=value.subs(t,1)
        assert s.cancel(old*s.prod(change[a,a] for a in word)-value)==0
    quadratic=z[0]*z[1]-t**2*z[3]**2/6
    assert s.hessian(quadratic,(z[0],z[1],z[3])).det()==t**2/3
    for value in (0,1,2):
        rel=[{pair:F(coefficient.subs(t,value)) for pair,coefficient in r.items() if coefficient.subs(t,value)!=0}
             for r in relations(t)]
        actual=top_tensor([rel]*4);normal=actual[27]
        actual={i:x/normal for i,x in actual.items()}
        expected={i:F(x.subs(t,value)) for i,x in tensor.items() if x.subs(t,value)!=0}
        assert actual==expected
        assert flatten(actual,(0,2))==6 and symmetrize(actual)[0] is not None
    certificate={'result':'PASS','arithmetic':'exact rational polynomial identities',
        'dual_quadrics':['e0^2','e1^2','e2^2','e3^2+t^2*e0*e1'],
        'associated_form':str(form),'SL4_gauge':'diag(t^-1,t^-1,t^2,1)',
        'graded_generic_isomorphism':'PASS','trace_tensor_generic_gauge':'PASS',
        'irreducible_quadratic_Hessian_det':'t^2/3','ranks_at_t_0_1_2':[6,6,6],
        'nonisomorphism_reason':'generic form has irreducible-factor degrees (1,1,2), special form (1,1,1,1); symmetric congruence criterion applies',
        'limits':['Individual orbit closure only, not failure of closedness of the entire symmetric chart image.',
                  'All fibers remain in the certified symmetric AS chart.']}
    if args.output:args.output.write_text(json.dumps(certificate,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(certificate,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
