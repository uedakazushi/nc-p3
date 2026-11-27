"""Projective stabilizers of the fully classified sparse Klein rectangles."""
from itertools import product

from verify_klein_three_lines import character,smith
from verify_parallel_invertible_rectangles import K,WORDS,bit,candidates,s
from verify_parallel_septic_curve import window
from verify_third_rectangle_orbits import support_maps


def check():
    for i in (88,89,90,142,143,145):
        words = [k for k in WORDS if candidates(1)[i]&bit(k)]
        weights = s.Matrix([[int(k[p]==j) for p,j in product(K,repeat=2)] for k in words])
        diag = smith(weights)[0]
        assert [n for n in diag.diagonal() if n]==[1]*11+[2,2]
    # On the double boundary, count all successful permutations; the orbit
    # computation alone only needs the existence of one such permutation.
    y = s.Symbol('y')
    modulus = y*y-y+1
    def value(f):
        n,d = s.fraction(s.cancel(f))
        return s.rem(n*s.invert(d,modulus,y),modulus,y)
    for t in (1,-1):
        w = {k:value(c) for k,c in window(y,t).items() if value(c)!=0}
        words = sorted(w)
        weights = s.Matrix([[int(k[p]==j) for p,j in product(K,repeat=2)] for k in words])
        diag,left,_ = smith(weights)
        rows = left.tolist()[sum(bool(n) for n in diag.diagonal()):]
        good = []
        for pi in support_maps(w,w):
            ratios = [value(w[k]/w[tuple(p[j] for p,j in zip(pi,k))]) for k in words]
            if all(value(character(row,ratios)-1)==0 for row in rows):
                good.append(pi)
        assert good==[(tuple(K),)*4]
    print('PASS: exact diagonal stabilizers and the two period-one boundary stabilizers')


if __name__=='__main__':
    check()
