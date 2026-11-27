"""Certify the rank-two quadrics and exceptional identifications of the symmetric three-line sheet."""
from itertools import combinations, permutations, product

import sympy as s

K = range(4)


def main():
    a,b,c,d,t,x,y,z = s.symbols('a b c d t x y z')
    p = (a,b,c,d)
    Q = s.Matrix([[t*a,x,y,z],[x,t*b,-z,-y],[y,-z,t*c,-x],[z,-y,-x,t*d]])
    H = s.Matrix([[(-1)**((i&j).bit_count()) for j in K] for i in K])/2
    E = H*s.diag(-1,1,1,1)
    assert E*E.T==s.eye(4)
    thin = E*Q.subs(t,0)*E.T
    assert thin==s.diag(*thin.diagonal()) and thin.trace()==0
    r = (a+b-c-d,a+c-b-d,a+d-b-c)
    # The diagonalized fourth quadric is a Klein-circulant matrix.
    v = (sum(p),r[1],r[0],r[2])
    even = E*s.diag(*p)*E.T
    assert (4*even-s.Matrix([[v[i^j] for j in K] for i in K])).applyfunc(s.expand)==s.zeros(4)
    adj = Q.subs(t,1).adjugate()
    X = (x,y,z)
    for j in range(3):
        k,l = [i for i in range(3) if i!=j]
        m,n = [i+1 for i in range(3) if i!=j]
        f = 4*(adj[0,j+1]+adj[m,n])
        rhs = (sum(p)*r[j]-r[k]*r[l])*X[j]+4*r[j]*X[k]*X[l]
        assert s.expand(f-rhs)==0
    assert s.factor(Q.subs(t,0).det())==(x-y-z)*(x-y+z)*(x+y-z)*(x+y+z)
    # In the rank-one-frame case the 10 nodes are the differences of five basis vectors.
    nodes = [s.eye(5)[:,i]-s.eye(5)[:,j] for i,j in combinations(range(5),2)]
    planes = {ids for ids in combinations(range(10),6) if s.Matrix.hstack(*(nodes[i] for i in ids)).rank()==3}
    assert len(planes)==5
    # The action on the three nonconstant Fourier coefficients is exactly the even signed permutations.
    samples = (1,2,4,8)
    def fourier(values):
        return tuple(sum((-1)**((i&j).bit_count())*values[i] for i in K) for j in (1,2,3))
    v = fourier(samples)
    left = {fourier(w) for w in permutations(samples)}
    right = {tuple(j*k for j,k in zip(signs,pi)) for pi in permutations(v)
             for signs in product((1,-1),repeat=3) if s.prod(signs)==1}
    assert left==right and len(left)==24
    # In the two-pair case compute the graph on the two isotropic lines of each binary block.
    r = s.Symbol('r',nonzero=True)
    d0,d1,d2,d3 = s.symbols('d:4')
    time = (d0+d1+d2+d3)/4
    U = s.Matrix([d0+d1+2*r*time,d0+d1-2*r*time])
    V = s.Matrix([d2+d3+2*r*time,d2+d3-2*r*time])
    graph = s.Matrix([[1,-1-r],[1-r,-1]])/r
    assert (graph*U-V).applyfunc(s.expand)==s.zeros(2,1)
    assert s.factor(graph.det())==-1
    assert s.factor(graph[0,0]*graph[1,1]/(graph[0,1]*graph[1,0])-1/(1-r**2))==0
    print('PASS: diagonal net, all paired cofactors, five exceptional planes, cross-ratio involution')


if __name__=='__main__':
    main()
