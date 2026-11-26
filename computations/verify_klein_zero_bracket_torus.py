"""Walsh cofactors, opposite pencils, cocycle signs and graded foldings at zero bracket."""
from itertools import combinations, product

import sympy as sp
from verify_klein_mixed_supports import character, smith, tensor


def coefficients(r,u,h):
    solution = {}
    for i in range(3):
        j,k = (i+1)%3,(i+2)%3
        values = [(u,1,u*r[i],r[i]),(1/r[j],h/r[k],-r[k]*r[i],-r[i]*r[j]/h)]
        for side,values_at_axis in enumerate(values):
            for slot,v in enumerate(values_at_axis):
                solution[side,i,slot] = sp.sympify(v)
    return solution


def walsh(w,x,y):
    rows,values = [],[]
    keys = []
    for i,b,e in product(range(4),range(1,4),range(1,4)):
        key = (i,i^b,i^b^e,i^e)
        keys.append(key)
        row = [0]*25
        for pos,label in enumerate(key):
            row[4*pos+label] = 1
        row[16+3*(b-1)+e-1] = -1
        rows.append(row)
        values.append((-1)**((i>>1)*((b^e)&1))/w[key])
    assert set(keys)==set(w)
    matrix = sp.Matrix(rows)
    d,l,_ = smith(matrix)
    rank = sum(v!=0 for v in d.diagonal())
    assert all(character(row,values)==1 for row in l.tolist()[rank:])
    wanted = sp.zeros(1,25)
    for idx,val in ((16,2),(20,2),(17,-2),(19,-2)):
        wanted[idx] = val
    certificate = sp.zeros(1,36)
    for idx,val in zip((9,10,12,13,18,19,21,22),(-1,1,1,-1,-1,1,1,-1)):
        certificate[idx] = val
    assert certificate*matrix==wanted
    assert sp.factor(character(certificate,values)-1-(x*y+x+y)/y**2)==0
    for j,k,r in ((2,3,x+y),(3,1,1+y),(1,2,1+x)):
        assert w[j,k,j,k]==w[k,j,k,j] and sp.factor(w[j,k,j,k]**2-r**2)==0


def pencils(w,x,y):
    t = [1]+[-1/(x*y)]*3
    v = [1,1,1/x,1/y]
    def entry(a,b,c,d):
        return sp.factor(t[c]*v[d]*w.get((a,b,c,d),0))
    for a,b,c,d in product(range(4),repeat=4):
        assert entry(a,b,c,d)==entry(c,b,a,d)==entry(a,d,c,b)
    def block(color):
        group = [(a,c) for a in range(4) for c in range(a,4) if a^c==color]
        return sp.Matrix([[entry(a,b,c,d) for b,d in group] for a,c in group])
    for i,target in enumerate((1/(x*y)**2,1/(x*y**2),1/(x**2*y)),1):
        assert sp.factor(block(i).det()-target)==0
    m = block(0)
    for a in (m,m.T):
        ns = sp.Matrix.hstack(*a.nullspace())
        assert a.rank()==2 and sp.factor(a[:2,:2].det())!=0
        assert ns[2:,:]==sp.eye(2)
        pluecker = [sp.factor(ns[[i,j],:].det()) for i,j in combinations(range(4),2)]
        factors = {z for f in pluecker for z,_ in sp.factor_list(sp.fraction(f)[0])[1]}
        assert factors-{x,y}=={x-1,y-1,x-y}
    # Full-support G and Q nets generate the coordinate projections from squares.
    # The det of y -> Q_y a and x -> G_x^t a is a nonzero monomial in all three coordinates.
    a = sp.symbols('a0:3')
    for side,expected in ((0,2),(1,2/(x*y))):
        columns = []
        for i in range(3):
            j,k = (i+1)%3,(i+2)%3
            col = [0]*3
            if side==0:
                col[j],col[k] = a[k],a[j]
            else:
                r = [1,x,y]
                col[j],col[k] = a[k]/r[k],a[j]/r[j]
            columns.append(col)
        assert sp.factor(sp.Matrix(columns).det()-expected*sp.prod(a))==0


def twist(w,x,y):
    def bichar(i,j):
        return (-1)**((i>>1)*(j&1))
    twisted = {k:v*bichar(k[2],k[2])*bichar(k[3],k[3])*sp.prod(bichar(a,b) for a,b in combinations(k,2))
               for k,v in w.items()}
    qr,gr = 1,1
    for i,r in enumerate((1,x,y),1):
        j,k = i%3+1,(i+1)%3+1
        assert sp.factor(twisted[k,i,j,0]*twisted[j,0,k,i]-r)==0
        qr *= twisted[j,i,k,0]/twisted[k,i,j,0]
        gr *= twisted[i,j,0,k]/twisted[i,k,0,j]
    assert sp.factor(qr)==sp.factor(gr)==1


def folding(w,sign):
    rows,values = [],[]
    for p,i,direction,swap in product(range(4),range(3),range(2),range(2)):
        nxt = (p+1)%4
        j,k = (i+1+swap)%3,(i+2-swap)%3
        key = (j+1,k+1)+((0,i+1) if not direction else (i+1,0))
        key = key[-p:]+key[:-p] if p else key
        row = [0]*18
        row[3*p+j] += 1
        row[3*nxt+k] += 1
        row[3*(p if direction else nxt)+i] -= 1
        row[12+2*i+(swap^direction)] -= 1
        rows.append(row)
        values.append(sign**direction/w[key])
    m = sp.Matrix(rows)
    d,l,_ = smith(m)
    rank = sum(v!=0 for v in d.diagonal())
    assert all(character(row,values)==1 for row in l.tolist()[rank:])


def main():
    x,y = sp.symbols('x y',nonzero=True)
    r = [sp.Integer(1),x,y]
    for u,h in product(map(sp.Integer,(1,-1)),repeat=2):
        w = tensor(coefficients(r,u,h))
        if u!=h:
            walsh(w,x,y)
        else:
            folding(w,u)
            if u==1:
                pencils(w,x,y)
            else:
                twist(w,x,y)
    print('PASS: both Walsh cofactors, two opposite pencils, both graded foldings and the Klein cocycle')


if __name__=='__main__':
    main()
