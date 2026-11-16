"""Partially alternating rank-six normal-form tensors and Pfaffian matrices."""
from itertools import permutations
import sympy as s


def tensor(kind,canonical=True):
    w={}
    if kind=='volume':
        for word in permutations(range(4)):
            inversions=sum(word[i]>word[j] for i in range(4) for j in range(i+1,4))
            a,b,c,d=word;w[64*a+16*b+4*c+d]=(-1)**inversions
        return w
    assert kind in ('rank4','rank5')
    def add(pair,c,d):
        for a,b,value in [(pair[0],pair[1],1),(pair[1],pair[0],-1)]:
            encoded=64*a+16*c+4*b+d if canonical else 64*a+16*b+4*c+d
            w[encoded]=w.get(encoded,0)+value
    first=[(0,1),(0,2),(1,2)]
    second=first if kind=='rank4' else [(0,1),(0,3),(1,3)]
    for j,pair in enumerate(first,1):add(pair,0,j)
    for j,pair in enumerate(second,1):add(pair,j,0)
    if kind=='rank4':add((0,3),0,0)
    return {i:x for i,x in w.items() if x}


def pfaffian_matrices():
    a,b,c,d,e,f=s.symbols('z01 z02 z03 z12 z13 z23')
    X=s.Matrix([[0,a,b,c],[-a,0,d,e],[-b,-d,0,f],[-c,-e,-f,0]])
    Y=s.Matrix([[0,-f,e,-d],[f,0,-c,b],[-e,c,0,-a],[d,-b,a,0]])
    return (a,b,c,d,e,f),X,Y,a*f-b*e+c*d
