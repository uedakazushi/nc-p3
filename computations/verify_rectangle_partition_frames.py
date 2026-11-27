"""Intrinsic coordinate subspaces from stable partitions under tensor contractions."""
from functools import cache
from itertools import product
import os
import subprocess

from verify_parallel_invertible_rectangles import K,completed,s
from verify_third_rectangle_frames import SINGULAR,frame_operators


def partitions(elements):
    if not elements:
        yield ()
        return
    first = elements[0]
    for rest in partitions(elements[1:]):
        yield ((first,),)+rest
        for i,block in enumerate(rest):
            yield rest[:i]+((first,)+block,)+rest[i+1:]


def stable_partitions(raw,variables,multiplier,units,possibilities):
    primitive,moment,_,weights,zero = frame_operators(raw,variables,multiplier,units)
    @cache
    def difference(p,i,j,blocks):
        other = [r for r in K if r!=p]
        return primitive(sum((f*((k[p]==i)-(k[p]==j)) for k,f in weights.items()
                              if all(k[r] in B for r,B in zip(other,blocks))),zero))
    @cache
    def moment_difference(p,i,j,q):
        row = moment(p,q)
        return primitive(row[i]-row[j])
    def insert(rows,f):
        if f:
            rows.add(str(f if f.leading_coefficient()>0 else -f))
    polynomials = {}
    commands = []
    for parts in possibilities:
        rows = set()
        for p,P in enumerate(parts):
            other = [r for r in K if r!=p]
            for block in P:
                for i,j in zip(block,block[1:]):
                    for blocks in product(*(parts[q] for q in other)):
                        insert(rows,difference(p,i,j,blocks))
        if '1' in rows:
            continue
        for p,P in enumerate(parts):
            for block in P:
                for i,j in zip(block,block[1:]):
                    for q in K:
                        if q!=p:
                            insert(rows,moment_difference(p,i,j,q))
        if '1' in rows:
            continue
        names = [polynomials.setdefault(f,'f'+str(len(polynomials))) for f in sorted(rows,key=lambda z:(len(z),z))]
        commands.append('print(localUnit(ideal('+(','.join(names) or '0')+')));')
    program = SINGULAR.replace('v,h,b',','.join(map(str,variables))).replace('v*h*(h-1)*b*(b-1)',str(units).replace('**','^'))
    program = 'LIB "elim.lib";\n'+program+'''
proc localUnit(ideal I) {
    ideal H=std(pairIdeal(I));
    if(H[1]==1) {return(1);}
    poly w=reduce(u,H);
    int i;
    for(i=0;i<10 && w!=0;i++) {w=reduce(w*w,H);}
    if(w==0) {return(1);}
    H=std(sat(H,ideal(u))[1]);
    return(int(H[1]==1));
}
'''
    program += '\n'.join(f'poly {n}={f};' for f,n in polynomials.items())+'\n'+'\n'.join(commands)+'\nquit;'
    result = subprocess.run([os.environ.get('SINGULAR_BIN','Singular'),'-q'],input=program,text=True,capture_output=True,check=True)
    out = [line.strip() for line in result.stdout.splitlines() if line.strip() and not line.startswith('//')]
    assert out==['1']*len(commands),(out,result.stderr)
    return len(commands)


def vertical_parameters(index,equal=False):
    a,y,d = s.symbols('a y d')
    if equal:
        return {118:(1,1,0),124:(1,1,d),142:(s.S.Half,s.S.Half,0),145:(a,a,2*a-1),147:(a,a,d)}[index]
    return {129:(1,0,d),134:(1,y,0),139:(1,y,d),142:(a,1-a,0),145:(a,y,a+y-1),147:(a,y,d)}[index]


def vertical_check():
    for t in (1,-1):
        for equal,indices in ((True,(118,124,142,145,147)),(False,(129,134,139,142,145,147))):
            target = ((((0,1),(2,),(3,)),((0,),(1,),(2,3))) if equal else (((0,),(1,),(2,),(3,)),)*2)
            choices = [[P for P in partitions(K) if all(any(set(B)<=set(C) for C in P) for B in T)] for T in target]
            cases = [(A,B,A,tuple(tuple(i^1 for i in block) for block in B)) for A,B in product(*choices)
                     if (len(A),len(B))!=tuple(map(len,target))]
            for idx in indices:
                a,y,d = map(s.sympify,vertical_parameters(idx,equal))
                w = completed(112,t,1,1,t*a,1,y,t*d)
                delta = a*y-d
                variables = tuple(sorted(a.free_symbols|y.free_symbols|d.free_symbols,key=str)) or (s.Symbol('v'),)
                units = s.prod({f for value in list(w.values())+[delta,a-y]
                                for poly in s.fraction(s.factor(value)) for f,_ in s.factor_list(poly)[1]})
                # Clear rational constants too, e.g. a=y=1/2 in stratum 142.
                scale = s.ilcm(1,*[s.Poly(s.cancel(delta*z),*variables).clear_denoms()[0] for z in w.values()])
                stable_partitions(w,variables,scale*delta,units,cases)
        print('PASS: intrinsic partitions for all vertical 112 strata at sign',t,flush=True)


if __name__=='__main__':
    vertical_check()
