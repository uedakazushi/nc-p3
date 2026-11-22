"""Check the two positive exchanges and the cocycle relating the four exceptional sheets."""
from itertools import product
import sympy as sp
from verify_klein_elliptic_completion import gram_twist


def check():
    a,b,c,d,z=sp.symbols('a b c d z',nonzero=True)
    s=sp.Matrix([[a,b],[c,d]])
    t=sp.Matrix([[d,c/z],[z*b,a]])
    w,_,_=gram_twist(s,t,z,1/z)
    r=[1,1,z,1]
    v=[1,1/z,1,1]
    for i,j,k,l in product(range(4),repeat=4):
        x=w[4*i+j,4*k+l]
        assert sp.cancel(r[i]*x-r[k^1]*w[4*(k^1)+j,4*(i^1)+l])==0
        assert sp.cancel(v[j]*x-v[l]*w[4*i+l,4*k+j])==0
        if i^j^k^l:
            continue
        word=[i,j,k,l]
        twist=sum((x>>1)*(y&1) for n,x in enumerate(word) for y in word[n+1:])
        gauge=(j&1)+(k&1)+(j>>1)*(j&1)+(l>>1)*(l&1)
        target=((j^k)&1) if (i^j)<2 else 0
        assert (twist+gauge-target)%2==0
    print('Two symbolic opposite exchanges and the four-letter cocycle sign identity.')


if __name__=='__main__':
    check()
