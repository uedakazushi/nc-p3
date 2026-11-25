"""Tangent bound at the unique quantum boundary and Laurent central sections for its Ore cells."""

from itertools import combinations, permutations, product

from flint import nmod_mat
from sympy.polys.domains import QQ_I

from verify_quantum_component_intersection import K, commutators

ONE, I = QQ_I.one, QQ_I(0,1)
ZERO = (0,0,0,0)
Y = (0,0,0,1)
WORDS = list(product(K,repeat=4))


def tangent():
    prime,root = 97,22
    assert root*root % prime == prime-1
    q = commutators((1,3,1))
    values = {}
    for word in permutations(K):
        inv = [(word[i],word[j]) for i,j in combinations(K,2) if word[i]>word[j]]
        values[word] = (-1)**len(inv)*pow(root,sum(q[i][j] for i,j in inv)%4,prime)
    rows = []
    for side in (0,1):
        def order(i,j,k,l):
            return (k,j,i,l) if side else (i,j,k,l)
        cut = nmod_mat([[values.get(order(i,j,k,l),0) for k,l in product(K,repeat=2)]
                        for i,j in product(K,repeat=2)],prime)
        left,a = cut.transpose().nullspace()
        right,b = cut.nullspace()
        assert a == b == 10
        for u,v in product(range(a),range(b)):
            rows.append([left[4*i+j,u]*right[4*k+l,v]
                         for word in WORDS for i,j,k,l in [order(*word)]])
    assert nmod_mat(rows,prime).rank()==177


def add(left,right,scale=ONE):
    out = dict(left)
    for key,val in right.items():
        out[key] = out.get(key,0)+scale*val
        if not out[key]:
            del out[key]
    return out


def multiply(left,right):
    out = {}
    for x,c in left.items():
        for y,d in right.items():
            key = tuple(i+j for i,j in zip(x,y,strict=True))
            sign = (-x[1]*y[0]+x[2]*y[0]-x[2]*y[1])%4
            out[key] = out.get(key,0)+c*d*I**sign
    return {k:z for k,z in out.items() if z}


def correction(p,a,s):
    b,c = (a+1)%3,(a+2)%3
    exp = [0]*4
    exp[a],exp[b],exp[c] = -1,s,2-s
    return {tuple(exp):(3*a+s+1)*I**(p*(s+1)%4)}


def generators(p,axes):
    out = {Y:ONE}
    for a in axes:
        for s in range(3):
            out=add(out,correction(p,a,s))
    return out


def power_section(p,j,axes):
    out = {ZERO:ONE}
    for n in range(p,p+j):
        out=multiply(out,generators(n,axes))
    return out


def polynomial_in_original(z,p):
    # Normal order with the ternary polynomial first and Ore letters on the right.
    remainder = dict(z)
    for j in range(4,-1,-1):
        leading = {(*e[:3],0):v for e,v in remainder.items() if e[3]==j}
        assert all(min(e)>=0 for e in leading)
        term=multiply(leading,power_section(p+4-j,j,range(3)))
        remainder=add(remainder,term,-ONE)
    assert not remainder


def verify_poles():
    y4 = {(0,0,0,4):ONE}
    glued = dict(y4)
    for a in range(3):
        z=power_section(0,4,[a])
        assert all(z==power_section(p,4,[a]) for p in K)
        assert all((e[0]-e[1])%4==(e[0]-e[2])%4==0 for e in z)
        glued=add(glued,add(z,y4,-ONE))
        for e in z:
            assert all(e[b]>=e[3] for b in range(3) if b!=a and e[3]<4)
    for p in K:
        polynomial_in_original(glued,p)
        polynomial_in_original({(1,1,1,1):ONE},p)


if __name__=="__main__":
    tangent()
    verify_poles()
    print("PASS: tangent dimension at most 79; all three coordinate-pole corrections glue")
