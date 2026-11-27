"""Recover the four intrinsic coordinate frames on the entire third rectangle."""
from itertools import combinations, permutations, product
from math import prod
import os
import subprocess

from flint import fmpz_mpoly_ctx
from verify_parallel_third_rectangle import K, WORDS, s, third_rectangle


def polynomials():
    v,h,b = s.symbols('v h b')
    ring = fmpz_mpoly_ctx.get(('v','h','b'))
    zero = ring.constant(0)

    def convert(f):
        return ring.from_dict({k:int(n) for k,n in s.Poly(f,v,h,b).terms()})

    units = convert(v*h*(h-1)*b*(b-1))

    def primitive(f):
        while f:
            g = f.gcd(units)
            if g.is_constant():
                return f//int(f.content())
            f = f//g
        return zero

    w = {k:convert(s.cancel(value*v*b*(h-1)*(b-1)))
         for k,value in third_rectangle(v,h/v,b).items()}

    def sign(seq):
        return (-1)**sum(i>j for i,j in combinations(seq,2))

    g = {}
    for key in WORDS:
        complements = [tuple(x for x in K if x!=i) for i in key]
        value = zero
        for tail in product(*(permutations(seq) for seq in complements[1:])):
            slots = [complements[0],*tail]
            x,y,z = zip(*slots)
            term = w.get(x,zero)*w.get(y,zero)*w.get(z,zero)
            if term:
                value += term*prod(sign((i,)+seq) for i,seq in zip(key,slots))
        g[key] = value
    for tensor in (w,g):
        coeffs = [a for a in tensor.values() if a]
        common = coeffs[0]
        for a in coeffs[1:]:
            common = common.gcd(a)
        assert primitive(common).is_constant()
        for key in tensor:
            tensor[key] = tensor[key]//common

    def second_moment(p,q):
        rest = tuple(j for j in K if j not in (p,q))
        def key(color,i,k):
            word = [0]*4
            for pos,value in zip((p,q)+rest,(i,i^color,k,k^color)):
                word[pos] = value
            return tuple(word)
        diagonal = [zero]*4
        for color in K:
            block = [[sum((w.get(key(color,i,k),zero)*g[key(color,j,k)] for k in K),zero)
                      for j in K] for i in K]
            for i in K:
                diagonal[i] += sum((block[i][j]*block[j][i] for j in K),zero)
        return diagonal

    weights = {k:value*g[k] for k,value in w.items()}
    def projected(p,chosen,letters):
        return [sum((a for k,a in weights.items() if k[p]==i
                     and tuple(k[j] for j in chosen)==letters),zero) for i in K]

    for p in (1,3,0,2):
        if p==1:
            entries = [second_moment(p,q) for q in (0,2,3)]
            source = second_moment(0,1)
            entries.append([sum((a*source[k[0]] for k,a in weights.items() if k[p]==i),zero)
                            for i in K])
        elif p==3:
            entries = [projected(p,(1,),(j,)) for j in K]
            entries.extend(second_moment(p,q) for q in (0,1))
        else:
            entries = [projected(p,(1,3),jl) for jl in product(K,repeat=2)]
        for i,j in combinations(K,2):
            yield (p,i,j),[str(primitive(row[i]-row[j])) for row in entries]


# Replacing a polynomial by its quotient by an invertible factor preserves its
# zero locus on the parameter torus. Gröbner bases are over the rationals.
SINGULAR = """
ring r=0,(v,h,b),dp;
poly u=v*h*(h-1)*b*(b-1);
proc stripUnits(poly f) {
    poly g=gcd(f,u);
    while (f!=0 && deg(g)>0) { f=f/g; g=gcd(f,u); }
    return(f);
}
proc removeUnits(ideal J) {
    ideal G=std(J); ideal S=G; int i; int changed=1;
    while(changed) {
        changed=0;
        for(i=1;i<=size(G);i++) {
            S[i]=stripUnits(G[i]);
            if(deg(S[i])<deg(G[i])) { changed=1; }
        }
        if(changed) { G=std(S); S=G; }
    }
    return(G);
}
proc checkPair(ideal P) {
    ideal J=std(ideal(0)); poly f; int i;
    for(i=1;i<=ncols(P);i++) {
        f=stripUnits(reduce(P[i],J));
        if(f!=0) { J=removeUnits(J+ideal(f)); attrib(J,"isSB",1); }
        if(J[1]==1) { return(1); }
    }
    return(0);
}
"""


def check():
    pairs = list(polynomials())
    commands = [f'print(checkPair(ideal({",".join(rows)})));' for _,rows in pairs]
    result = subprocess.run([os.environ.get('SINGULAR_BIN','Singular'),'-q'],
                            input=SINGULAR+'\n'.join(commands)+'\nquit;\n',
                            text=True,capture_output=True,check=True)
    output = [line.strip() for line in result.stdout.splitlines() if not line.startswith('//')]
    assert output==['1']*24,(pairs,result.stdout,result.stderr)
    print('PASS: all 24 pairs of letters separated over the rational parameter ring')


if __name__=='__main__':
    check()
