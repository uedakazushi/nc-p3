"""Symbolic charts for the nodal and tacnodal split binary determinants."""
from sympy import Matrix, cancel, expand, symbols

lam, m = symbols("lambda m", nonzero=True)
# Normalize a=1; solve c+b*lam=-m, b+c*lam=-lam/m.
b = (lam/m - lam*m)/(lam**2-1)
c = -m-b*lam
n = lam/m
u,v,z,t = symbols("u v z t")
F = u*v + z*t + b*(u*z+v*t) + c*(u*t+v*z)
assert cancel(F.subs(v,lam*u)-(z-m*u)*(t-n*u)) == 0
assert cancel(F.subs({t:n*u,z:u/m})) == 0
assert cancel(F.subs({t:u/n,z:m*u})) == 0

A,B,C = symbols("A B C", nonzero=True)
theta = Matrix([[A,B],[C,-lam*A]])
scale = Matrix.diag(m,1)
eta = scale.inv()*theta*scale
assert eta == Matrix([[A,B/m],[m*C,-lam*A]])
inverse_graph = theta*Matrix.diag(1/lam,1)*theta
G = expand(u*(inverse_graph[1,0]*v+inverse_graph[1,1])-inverse_graph[0,0]*v-inverse_graph[0,1])
H = m*C*u*v - n*A*u - m*A*v - B/m
# For B,C nonzero, proportionality forces m^2=1.
assert cancel(H.coeff(u*v)*G.subs({u:0,v:0}) - G.coeff(u*v)*H.subs({u:0,v:0}) - A*B*C*(lam-1/lam)*(m-1/m)) == 0
# Upper triangular case: u/v comparison gives m^2=-1/lam^2,
# and the v/constant comparison gives m^2*(1-lam^2)=1.
g,h = expand(G.subs(C,0)), expand(H.subs(C,0))
assert cancel(h.coeff(u)*g.coeff(v)-g.coeff(u)*h.coeff(v) - A**3*(1/m+lam**2*m)) == 0
assert cancel(h.coeff(v)*g.subs({u:0,v:0})-g.coeff(v)*h.subs({u:0,v:0}) - A**2*B*(m*(1-lam**2)-1/m)/lam) == 0

e,f,k = symbols("e f k", nonzero=True)
Ft = (u-v)*(t-v)+(z-v)*(-f*(u-v)+e*(t-v)-e*f)/(f-e)
assert cancel(Ft.subs(v,u) - e*(z-u)*(t-u-f)/(f-e)) == 0
assert cancel(Ft.subs(v,u+e) - e*(z-u-f)*(t-u-e)/(f-e)) == 0
assert cancel(Ft.subs({z:u+f,t:u+f})) == 0
Theta = Matrix([[0,k],[1,0]])
shift = Matrix([[1,f],[0,1]])
Eta = shift*Theta*shift.inv()
assert Eta == Matrix([[f,k-f**2],[1,-f]])
residual = Theta*Matrix([[1,e],[0,1]])*Theta
assert residual == Matrix([[k,0],[e,k]])
assert (u*(u+f)-f*u-f*(u+f)-(k-f**2)).subs(u,0) == -k
print("PASS: two nodal transfers, mixed traces and incompatibility, additive interpolation and tacnode residual")
