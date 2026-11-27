"""Dense transverse (2,2,1,1) cuts: minimal primes, graded ranks and a symmetric twist."""
from itertools import combinations, permutations, product
from pathlib import Path
from subprocess import check_output
from tempfile import TemporaryDirectory

import sympy as s

from verify_klein_mixed_positive import fourth_character, matrix
from verify_klein_three_lines import K, WORDS, thin_system, twist

X = s.symbols('x:16')
HOLES = [w for w in WORDS if (w[0]^w[1])<2 and (w[1]^w[2]) in (0,2)]


def components():
    x = X
    cubic = -x[3]*x[5]*x[15]+x[1]*x[7]*x[15]-x[1]+x[3]+x[5]-x[7]
    indices = ([1,1,3,3,5,5,7,7]+[15]*8,
               [3]*4+[7]*4+[11]*4+[15]*4,
               [15,15,15,12,12,15,12,12,12,12,15,12,12,15,15,15])
    sub = [dict(zip(x[4:],[s.Integer(k)]*12)) for k in (-1,1)]
    sub += [{t:x[j] for t,j in zip(x,row) if t!=x[j]} for row in indices]
    primes = [[var-val for var,val in part.items()] for part in sub]
    primes[2].append(cubic)
    sub[2][x[1]]=s.solve(cubic,x[1])[0]
    return primes,sub


def index_maps():
    for rev,rot,t,h,k in product(range(2),K,K,(0,1),(0,2)):
        offset = (t,t^h,t^h^k,t^k)
        def perm(w):
            z = w[::-1] if rev else w
            z = z[rot:]+z[:rot]
            if rot%2:
                z = tuple((0,2,1,3)[i] for i in z)
            return tuple(i^j for i,j in zip(z,offset))
        yield {x:X[HOLES.index(perm(w))] for x,w in zip(X,HOLES)}


def decomposition(primes,w):
    def key(gens):
        return tuple(sorted(map(str,s.groebner(gens,*X,order='grevlex'))))
    seen = set()
    for prime,size in zip(primes,(16,16,8,4,4)):
        orbit = {key([f.xreplace(mapping) for f in prime]) for mapping in index_maps()}
        assert len(orbit)==size and not seen&orbit
        seen |= orbit
    gens = set()
    for p,group in enumerate(((0,1),(0,2))):
        for b in group:
            m = matrix(w,p,b)
            gens.update(m[I,J].det() for I,J in product(combinations(K,3),repeat=2))
    source = ('LIB "primdec.lib";\nring r=0,('+','.join(map(str,X))+'),dp;\nideal I='+
              ','.join(str(f).replace('**','^') for f in gens)+';\nlist L=minAssGTZ(I);\n'+
              'int j; for(j=1;j<=size(L);j++){print("PRIME"); print(L[j]);}; quit;')
    with TemporaryDirectory() as tmp:
        path = Path(tmp)/'decompose.sing'
        path.write_text(source)
        output = check_output(['Singular','-q',str(path)],text=True)
    actual = set()
    for chunk in output.split('PRIME\n')[1:]:
        actual.add(key([s.sympify(line.rstrip(',').replace('^','**')) for line in chunk.splitlines() if line]))
    assert actual==seen


def window(u,v):
    w = dict.fromkeys(WORDS,s.Integer(1))
    for word in HOLES:
        i,j,k,_ = word
        sign = (i&1)^(i>>1)^((i^j)*((j^k)>>1))
        w[word] = v if sign else u
    return w


def symmetric_model(w,u,v):
    sign = (1,1,1,-1)
    W = {}
    for t in WORDS:
        f = sum((-1)**sum((i//2)*j for i,j in zip(t,bits))*
                w.get(tuple((i%2)^(3*j) for i,j in zip(t,bits)),0) for bits in product(range(2),repeat=4))
        a,b,c,d = t
        W[a,b,c^1,d] = s.expand(twist(t)*sign[b]*sign[d]*f/2)
    for key,value in W.items():
        assert all(W[t]==value for t in permutations(key))
    m = [s.Matrix([[W.get((i,i^b,k,k^b^1),0) for k in K] for i in K]) for b in K]
    assert m[0][:,0]-m[0][:,2]==s.Matrix([4,4,-4,-4])
    assert m[0][:,0]+m[0][:,2]==2*s.Matrix([u+1,v+1,u+1,v+1])
    assert m[0][:,0]==m[0][:,1] and m[0][:,2]==m[0][:,3]
    assert m[2]==s.Matrix([u-1,v-1,u-1,v-1])*s.ones(1,4)
    assert m[3]==s.ones(4,1)*s.Matrix([[u-1,v-1,u-1,v-1]])
    assert m[1]==m[0].T
    # Check that the other Fourier coefficients vanish before applying the graded twist.
    from operator import xor
    from functools import reduce
    for t in product(K,repeat=4):
        if reduce(xor,t):
            val=sum((-1)**sum((i//2)*j for i,j in zip(t,bits))*
                    w.get(tuple((i%2)^(3*j) for i,j in zip(t,bits)),0) for bits in product(range(2),repeat=4))
            assert s.expand(val)==0


def main():
    thin,ds,_,gauge = thin_system([(2,3),(1,3)])
    assert len(thin)==48 and ds==[1]*35 and gauge.rank()==13
    w = dict.fromkeys(WORDS,s.Integer(1))
    w.update(zip(HOLES,X))
    primes,charts = components()
    decomposition(primes,w)
    for sub,expected in zip(charts,(11,11,11,13)):
        candidate = {key:s.factor(val.xreplace(sub).subs(sub,simultaneous=True)) for key,val in w.items()}
        assert fourth_character(candidate,0,0)==expected
    u,v = X[15],X[12]
    candidate = window(u,v)
    assert {key:val.xreplace(charts[-1]) for key,val in w.items()}==candidate
    symmetric_model(candidate,u,v)
    for word,val in candidate.items():
        rot = word[1:]+word[:1]
        shifted = tuple((0,2,1,3)[j]^k for j,k in zip(rot,(0,1,1,0)))
        assert candidate[shifted]==val
    for b in K:
        m = s.Matrix([[candidate[i,j,i^b,j^b] for j in K] for i in K])
        assert s.expand(m.det()-(u-1)*(v-1)*(u*v+u+v-3))==0
    print('PASS: all 48 components, four obstructions, the full regular surface and the one-shift action')


if __name__=='__main__':
    main()
