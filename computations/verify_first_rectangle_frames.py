"""Recover the four intrinsic letter frames for vertical 76 windows and on a horizontal open."""
from itertools import combinations,product
import os
import subprocess

from verify_parallel_invertible_rectangles import K, WORDS, completed, matrix, s
from verify_third_rectangle_frames import SINGULAR,frame_differences
from verify_transverse_rank_twelve_isomorphism import alternating_gradient


def pair_planes(raw):
    w = {k:raw.get(k,0) for k in WORDS}
    g = alternating_gradient(w)
    for p in K:
        rows = []
        for q in K:
            if p==q:
                continue
            slots = (p,q,*sorted(set(K)-{p,q}))
            W,G = [{tuple(k[j] for j in slots):v for k,v in z.items()} for z in (w,g)]
            blocks = [matrix(W,0,h)*matrix(G,0,h).T for h in K]
            rows.append([s.factor(sum((M*M)[i,i] for M in blocks)) for i in K])
        shift = (2,3,2,3)[p]
        assert all(row[i]==row[i^shift] for row in rows for i in K)
        g0 = s.gcd_list([s.fraction(s.factor(row[0]-row[1]))[0] for row in rows])
        assert g0!=0 and g0.is_number
    # These two complementary plane quadruples use only the thin entries.
    # A two-dimensional image or row space has just two rank-one lines.
    for ids in ((0,0,1,0),(1,1,0,1)):
        pairs = [(i,i^u) for i,u in zip(ids,(2,3,2,3))]
        block = s.Matrix([[w.get((i,j,k,l),0) for j,l in product(pairs[1],pairs[3])]
                          for i,k in product(pairs[0],pairs[2])])
        rows = [i for i in K if any(block.row(i))]
        cols = [i for i in K if any(block.col(i))]
        assert len(rows)==len(cols)==2 and rows[0]^rows[1]==cols[0]^cols[1]==3
        det = s.factor(block.extract(rows,cols).det())
        assert det!=0 and det.is_number


def check():
    a,d = s.symbols('a d')
    unit = a*d*(a-1)*(d-1)*(a-d)
    for t in (1,-1):
        raw = completed(76,t,1,1,a,1,1,d)
        pairs = list(frame_differences(raw,(a,d),a-d,unit))
        program = SINGULAR.replace('v,h,b','a,d').replace('v*h*(h-1)*b*(b-1)',str(unit).replace('**','^'))
        for _,row in pairs:
            program += f'print(pairIdeal(ideal({",".join(row)})));\n'
        result = subprocess.run([os.getenv('SINGULAR_BIN','Singular'),'-q'],input=program+'quit;\n',
                                text=True,capture_output=True,check=True)
        output = [v.strip() for v in result.stdout.splitlines() if v.strip() and not v.startswith('//')]
        assert output==['a+d-1' if p==1 and (i,j) in ((0,3),(1,2)) else '1' for (p,i,j),_ in pairs]
        pair_planes(completed(76,t,1,1,a,1,1,1-a))
        edge = completed(76,t,1,1,1,1,1,d)
        for _,vals in frame_differences(edge,(d,),1-d,d*(d-1)):
            gcd = s.gcd_list([s.sympify(v.replace('^','**')) for v in vals])
            assert gcd!=0 and gcd.is_number
    pair_planes(completed(76,2,3,s.Rational(-1,2),7,1,1,s.Rational(3,4)))


if __name__=='__main__':
    check()
    print('PASS: four unmarked frames throughout both vertical 76 charts and at a horizontal generic witness')
