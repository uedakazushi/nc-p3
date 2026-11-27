"""All normalized unmarked vertical 112 orbits, including the Fourier frames on a=y."""
from itertools import permutations,product

from sympy.matrices.normalforms import hermite_normal_form
from verify_parallel_invertible_rectangles import K,T as THIN,X,bit,candidates,completed,s,smith,template,top
from verify_rectangle_partition_frames import vertical_parameters
from verify_third_rectangle_orbits import support_maps


def weight(k):
    return [int(k[p]==j) for p,j in product(K,repeat=2)]


def vertical_actions():
    a,y,d = s.symbols('a y d')
    base = template(112)
    holes = [k for k,f in base.items() if f in X]
    thin = [k for k in base if k not in holes]
    def exponent(f):
        return s.factor(THIN*s.diff(f,THIN)/f)
    gauge = s.Matrix([weight(k)+[exponent(base[k])] for k in thin])
    supports = {i:[k for k in base if candidates(1)[i]&bit(k)] for i in (118,124,129,134,139,142,145,147)}
    def maps(r,src,dst):
        return list(support_maps([k[r:]+k[:r] for k in src],dst))
    full = [maps(r,supports[147],supports[147]) for r in K]
    expected = {118:(2,0),124:(2,0),129:(1,1),134:(1,0),139:(1,0),142:(2,2),145:(2,0),147:(2,2)}
    for i,keys in supports.items():
        D = smith(s.Matrix(list(map(weight,keys))))[0]
        assert [n for n in D.diagonal() if n]==[1]*11+[2,2]
        for j,other in supports.items():
            if len(keys)!=len(other):
                continue
            for r in K:
                perms = maps(r,keys,other)
                assert len(perms)==(expected[i][r%2] if i==j else 0)
                assert set(perms)<=set(full[r])
    # The first half-turn permutation fixes every vertical window after rescaling.
    assert full[2][0]==(tuple(K),(1,0,3,2),tuple(K),(1,0,3,2))
    for t in (s.S.One,-s.S.One):
        raw = completed(112,t,1,1,t*a,1,y,t*d)
        for r,mp in enumerate(full):
            coordinates = []
            for pi in mp:
                def move(w):
                    return {tuple(p[i] for p,i in zip(pi,k[r:]+k[:r])):v for k,v in w.items()}
                z,zz = move(raw),move(base)
                vec = s.Matrix([exponent(zz[k]) for k in thin])
                sol = next(iter(s.linsolve((gauge,vec))))
                subs = dict.fromkeys(set().union(*(v.free_symbols for v in sol)),0)
                es = [v.xreplace(subs) for v in sol]
                tt = t**es[-1]
                def val(k):
                    return z[k]*t**-sum(es[4*p+i] for p,i in enumerate(k))
                assert tt==t and all(val(k)==base[k].subs(THIN,tt) for k in thin)
                normalized = [s.factor(val(k)/val(holes[4])) for k in holes[:8]]
                A,Y,D = normalized[2]/t,normalized[5],normalized[6]/t
                assert all(s.factor(f-g)==0 for f,g in zip(normalized,top(112,t,1,1,t*A,1,Y,t*D)))
                coordinates.append((A,Y,D))
            expected_coords = ((a,y,d),(y,a,d)) if r%2==0 else ((1-a,1-y,d+1-a-y),(1-y,1-a,d+1-a-y))
            assert set(coordinates)==set(expected_coords)
            if r==2:
                assert coordinates[0]==(a,y,d)


def fourier(raw):
    w = raw
    for p in K:
        u = 0 if p%2==0 else 2
        keys = {k[:p]+(j,)+k[p+1:] for k in w for j in K}
        def value(k):
            if k[p] not in (u,u+1):
                return w.get(k,0)
            return s.factor(w.get(k[:p]+(u,)+k[p+1:],0)+(1 if k[p]==u else -1)*w.get(k[:p]+(u+1,)+k[p+1:],0))
        w = {k:z for k in keys if (z:=value(k))}
    return w


def fourier_orbits():
    tensors = {}
    for t,i in product((1,-1),(118,124,142,145,147)):
        a,_,d = vertical_parameters(i,True)
        raw = completed(112,t,1,1,t*a,1,a,t*d)
        # The singleton projections give the same two rank-one points inside
        # each doublet. Check both adjacent pairs, all four singleton choices.
        for p in (0,2):
            matrices = []
            for k,l in product((2,3),(0,1)):
                matrices.append(s.Matrix([[raw.get((i,j,k,l) if p==0 else (k,l,i,j),0)
                                          for j in (2,3)] for i in (0,1)]))
            H = s.Matrix([[1,1],[1,-1]])
            rotated = [H*M*H for M in matrices]
            assert all(z[0,0]==z[1,1]==0 for z in rotated) or all(z[0,1]==z[1,0]==0 for z in rotated)
            assert s.Matrix([list(z) for z in rotated]).rank()==2
        w = fourier(raw)
        tensors[t,i] = w
        keys = sorted(w)
        # Every factor of a nonzero entry is a unit of the specified stratum.
        coeff = []
        for k in keys:
            n,den = s.fraction(w[k])
            ex = {}
            for sign,poly in ((1,n),(-1,den)):
                for f,e in s.factor_list(poly)[1]:
                    ex[f]=ex.get(f,0)+sign*e
            coeff.append(ex)
        variables = sorted(set().union(*(set(row) for row in coeff)),key=str)
        G = s.Matrix(list(map(weight,keys)))
        diag,left,_ = smith(G)
        ds = [n for n in diag.diagonal() if n]
        assert ds==([1]*11+[2,2] if t==1 else [1]*12)
        if variables:
            E = s.Matrix([[row.get(f,0) for f in variables] for row in coeff])
            lattice = left[len(ds):,:]*E
            assert hermite_normal_form(lattice.T)==s.eye(len(variables))
    per = [p for p in permutations(K) if set(p[:2])=={0,1}]
    for (t,i),w in tensors.items():
        for (s0,j),target in tensors.items():
            if len(w)!=len(target):
                continue
            for r in (0,1):
                src = [k[r:]+k[:r] for k in w]
                maps = []
                for perms in product(per,repeat=4):
                    pi = tuple(tuple(q^2 for q in p) for p in perms) if r else perms
                    if {tuple(p[k] for p,k in zip(pi,word)) for word in src}==set(target):
                        maps.append(pi)
                assert len(maps)==(2 if t==s0 and i==j and (r==0 or i in (142,147)) else 0)
                if r==0 and maps:
                    assert maps[0]==(tuple(K),)*4
                    assert maps[1]==((1,0,2,3),(0,1,3,2))*2


if __name__=='__main__':
    vertical_actions()
    fourier_orbits()
    print('PASS: every vertical 112 support, stabilizer weight lattice, Fourier frame and rotation')
