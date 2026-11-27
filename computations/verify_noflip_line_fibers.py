"""Smith ranks, Schur obstructions, and period strata of the unreflected Bell line."""
from itertools import product
from subprocess import run

import sympy as sp
from flint import nmod_mat
from sympy.polys.matrices import DomainMatrix

from verify_isolated_split_component import WORDS
from verify_reflected_collision_fibers import CharacterEquations
from verify_split_spectral_line import I, line_tensor

S = sp.Symbol('s', nonzero=True)
PRIME = 257


def linear(eq, words):
    rows = []
    for i,(left,right,_) in enumerate(eq.cuts):
        for a,b in product(left.to_list(),right.transpose().to_list()):
            rows.append([a[4*q[i]+q[(i+1)%4]]*b[4*q[(i+2)%4]+q[(i+3)%4]] for q in words])
    return eq.matrix(rows,len(words))


def ranks_and_obstructions():
    field = sp.QQ.frac_field(S)
    eq = CharacterEquations(line_tensor(S,(I,)*4),field)
    script = ['LIB "jacobson.lib";ring r=0,s,dp;',r"""
proc verifyRows(matrix a,int expected){
 matrix h=smith(transpose(a));int z,j;poly f;
 for(j=1;j<=ncols(a);j++){
    if(j>nrows(h) || j>ncols(h)){z++;}
    else{f=h[j,j];if(f==0){z++;}else{
        while(f%s==0){f=f/s;};if(deg(f)>0){print("NONUNIT");};
    };}
 }
 if(z!=expected){print("WRONG");};
}
"""]
    matrices = {}
    for n,predicted in ((3,20),(4,3)):
        words = [q for q in WORDS if sum(a<2 for a in q)==n]
        m = linear(eq,words)
        matrices[n] = words,m
        rows = []
        for row in m.to_list():
            if not any(row):
                continue
            vs = [field.to_sympy(a) for a in row]
            den = sp.lcm([sp.denom(a) for a in vs])
            assert all(f==S for f,_ in sp.factor_list(den)[1])
            rows.append([sp.factor(a*den) for a in vs])
        data = ','.join(str(a).replace('**','^') for row in rows for a in row)
        script.append(f'matrix a[{len(rows)}][{len(words)}]={data};verifyRows(a,{predicted});kill a;')
    out = run(['Singular','-q'],input='\n'.join(script+['quit;']),text=True,capture_output=True,check=True).stdout
    assert not any(a in out for a in ('?', 'WRONG', 'NONUNIT')), out
    words,m = matrices[3]
    orbit = eq.matrix([[eq.w.get(q[:j]+(t,)+q[j+1:],field.zero) if q[j]==a else field.zero
                        for j,a,t in product(range(4),range(2),range(2,4))] for q in words],16)
    rows = orbit.transpose().rref()[1]
    det = field.to_sympy(orbit.extract(rows,range(16)).det())
    num,den = sp.fraction(sp.factor(det))
    assert all(sp.denom(sp.cancel(S*(S**2-1)/p)).is_number
               for a in (num,den) for p,_ in sp.factor_list(a)[1])
    k = m.nullspace(divide_last=True).transpose()
    frame = orbit.hstack(k)
    ids = [j for j in frame.rref()[1] if j>=16]
    assert len(ids)==4
    cols = frame.extract(range(len(words)),ids).to_list()
    directions = [[eq.cut(dict(zip(words,[r[j] for r in cols])),i) for i in (0,1)] for j in range(4)]
    four = matrices[4][1]
    extra = four
    for j in range(4):
        for k in range(j,4):
            rhs = []
            for i,(left,right,g) in enumerate(eq.cuts):
                a,b = directions[j][i],directions[k][i]
                rhs.extend([v] for r in (left*(a*g*b+b*g*a)*right).to_list() for v in r)
            extra = extra.hstack(eq.matrix(rhs,1))
    assert four.rank()==extra.rank()==13


def matrix(rows):
    return nmod_mat(rows,PRIME)


def select(m,rows,cols):
    if not len(rows):
        return nmod_mat(0,len(cols),PRIME)
    return matrix([[int(m[i,j]) for j in cols] for i in rows])


def join(a,b):
    return matrix([[int(a[i,j]) for j in range(a.ncols())]+[int(b[i,j]) for j in range(b.ncols())]
                   for i in range(a.nrows())])


def identity(n):
    return matrix([[int(i==j) for j in range(n)] for i in range(n)])


def pivots(m):
    r,n = m.rref()
    return [next(j for j in range(m.ncols()) if r[i,j]) for i in range(n)]


def kernel(m):
    k,n = m.nullspace()
    return select(k,range(m.ncols()),range(n))


def representations(s,root,exceptional):
    raw = line_tensor(sp.Integer(s),(I,)*4)
    w = {q:int(sp.numer(a))*pow(int(sp.denom(a)),-1,PRIME)%PRIME for q,a in raw.items()}
    ks = []
    for i in range(4):
        m = nmod_mat(16,16,[w.get(q[-i:]+q[:-i] if i else q,0) for q in WORDS],PRIME)
        ks.append(kernel(m.transpose()))
    sectors = []
    for n in (3,4):
        words = [q for q in WORDS if sum(a<2 for a in q)==n]
        pos = {q:j for j,q in enumerate(words)}
        rows = []
        for i in (0,1):
            for x,y in product(range(10),repeat=2):
                row = [int(ks[i][4*q[i]+q[(i+1)%4],x]*ks[i+2][4*q[(i+2)%4]+q[(i+3)%4],y]) for q in words]
                if any(row):
                    rows.append(row)
        frame = kernel(matrix(rows))
        d = 0
        if n==3:
            g = matrix([[w.get(q[:i]+(t,)+q[i+1:],0) if q[i]==a else 0
                         for i,a,t in product(range(4),range(2),range(2,4))] for q in words])
            frame = join(g,frame)
            frame = select(frame,range(len(words)),pivots(frame))
            d = 16
        rs = pivots(frame.transpose())
        lift = select(frame,rs,range(frame.ncols())).inv()

        def action(fn):
            for q,a in w.items():
                target,c = fn(q,2)
                assert a*c%PRIME==w[target]
            image = nmod_mat(len(words),frame.ncols()-d,PRIME)
            for j,q in enumerate(words):
                target,c = fn(q,n)
                for b in range(d,frame.ncols()):
                    image[pos[target],b-d] = int(frame[j,b])*c%PRIME
            r = lift*select(image,rs,range(image.ncols()))
            assert frame*r==image
            return select(r,range(d,r.nrows()),range(r.ncols()))
        x = action(lambda q,n:(tuple(a^1 for a in q),1))
        a = action(lambda q,n:(q,(-1)**q.count(0)))
        b = action(lambda q,n:(q,(-1)**q.count(2)))
        def shift(q,n):
            r = q[1:]+q[:1]
            factor = root**(n-2)*pow(s,int(q[0]>=2)-int(q[2]<2),PRIME)
            return tuple(x^int(x>=2) for x in r),factor%PRIME
        sh = action(shift)
        assert sh**4==identity(sh.nrows())
        group = [x**i*a**j*b**k for i,j,k in product(range(2),repeat=3)]
        if exceptional:
            ph = [-s,1,-1,-s, s,1,s,1, 1,s,-s,1, s,1,s,1]
            e = action(lambda q,n:(tuple(x^int(x>=2) for x in q),
                                   int(sp.prod(ph[4*i+x] for i,x in enumerate(q)))%PRIME))
            group += [e*h for h in group]
        sectors.append((sh,group))
    return sectors


def components(s,root,exceptional):
    (sb,hb),(sc,hc) = representations(s,root,exceptional)
    group = list(zip(hb,hc))
    units = sorted({pow(3,4*j,PRIME) for j in range(64)})
    assert len(units)==64
    answer = []
    for m in (1,2):
        charts = []
        for gb,gc in group:
            mb,mc = gb*sb**m,gc*sc**m
            eigb = {a:kernel(mb-a*identity(4)) for a in units}
            eigc = {a:kernel(mc-a*identity(3)) for a in units}
            assert sum(v.ncols() for v in eigb.values())==4
            assert sum(v.ncols() for v in eigc.values())==3
            for a in units:
                if eigb[a].ncols():
                    charts.append((eigb[a],eigc[a*a%PRIME]))
                if eigc[a].ncols():
                    charts.append((nmod_mat(4,0,PRIME),eigc[a]))
        charts.sort(key=lambda z:-(z[0].ncols()+z[1].ncols()))
        reps = []
        def inside(a,b):
            return join(a,b).rank()==b.ncols()
        for b,c in charts:
            if not any(inside(b,gb*v) and inside(c,gc*w) for v,w in reps for gb,gc in group):
                reps.append((b,c))
        answer.append([b.ncols()+c.ncols()-1 for b,c in reps])
    return answer


if __name__=='__main__':
    ranks_and_obstructions()
    assert components(4,2,False)==[[4,4,3,3,3,3,2,2],[6]]
    assert components(16,4,True)==[[4,4,3,3,2],[6]]
    print('All no-flip line lifts, stabilizer comparisons and period components verified.')
