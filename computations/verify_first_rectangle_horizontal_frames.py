"""All regular-candidate horizontal 76 frames by two-plane reconstruction."""
from itertools import combinations,product
import os
import subprocess

from flint import fmpz_mpoly_ctx
from verify_parallel_invertible_rectangles import K,completed,s
from verify_rectangle_partition_frames import stable_partitions
from verify_third_rectangle_frames import SINGULAR,frame_operators


def check():
    b,t,x = s.symbols('b t x')
    c = (t*t-b)/(1-b)
    raw = completed(76,t,b,c,x,1,1,c*(1-x)/t**2)
    multiplier = s.lcm([s.fraction(s.cancel(z))[1] for z in raw.values()])
    units = s.prod({f for value in list(raw.values())+[b+c,b-1]
                    for poly in s.fraction(s.factor(value)) for f,_ in s.factor_list(poly)[1]})
    target = tuple(((0,h),(1,1^h)) for h in (2,3,2,3))
    whole = ((tuple(K),),)*4
    cases = [P for P in product(*((T,(tuple(K),)) for T in target)) if P not in (target,whole)]
    stable_partitions(raw,(b,t,x),multiplier,units,cases)

    primitive,moment,_,_,_ = frame_operators(raw,(b,t,x),multiplier,units)
    ring = fmpz_mpoly_ctx.get(('b','t'))
    coeffs = []
    for p,q in product(K,repeat=2):
        if p!=q:
            row = moment(p,q)
            f = primitive(row[0]-row[1])
            assert f.degrees()[2]<=2
            coeffs.append([ring.from_dict({(i,j):z for (i,j,h),z in f.to_dict().items() if h==n}) for n in range(3)])
    bb,tt = ring.gens()
    base_unit = bb*tt*(bb-1)*(tt-bb)*(tt+bb)*(tt*tt-bb)
    minors = set()
    for (a,b1,c1),(d,e,f),(g,h,i) in combinations(coeffs,3):
        m = a*(e*i-f*h)-b1*(d*i-f*g)+c1*(d*h-e*g)
        if m:
            while not (v:=m.gcd(base_unit)).is_constant():
                m = m//v
            m = m//int(m.content())
            minors.add(str(m if m.leading_coefficient()>0 else -m))
    program = 'LIB "elim.lib";\n'+SINGULAR.replace('v,h,b','b,t').replace('v*h*(h-1)*b*(b-1)',str(base_unit))
    program += 'ideal H=std(pairIdeal(ideal('+','.join(sorted(minors,key=lambda z:(len(z),z)))+')));\n'
    program += 'H=std(sat(H,ideal(u))[1]);print(H);quit;'
    res = subprocess.run([os.environ.get('SINGULAR_BIN','Singular'),'-q'],input=program,
                         text=True,capture_output=True,check=True)
    assert [z for z in res.stdout.splitlines() if z.strip() and not z.startswith('//')]==['1']
    for ids in ((0,0,1,0),(1,1,0,1)):
        P = [(i,i^h) for i,h in zip(ids,(2,3,2,3))]
        M = s.Matrix([[raw.get((i,j,k,l),0) for j,l in product(P[1],P[3])] for i,k in product(P[0],P[2])])
        assert M==s.Matrix([[0,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,0]])
    print('PASS: all partial plane partitions and full-rank quadratic coefficient minors on the horizontal open')


if __name__=='__main__':
    check()
