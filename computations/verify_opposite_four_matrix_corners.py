"""Integer dimension bounds for universal extensions at two vertices of a four-cycle."""

from itertools import combinations, product


def surviving_extensions(d, marked, *, separated=False):
    h = [2*d[i] - 4*d[i-1] - 4*d[(i+1)%4] + 6*d[(i+2)%4] for i in range(4)]
    euler = sum(x*y for x,y in zip(d, h))
    strong = separated or any(2*d[t] > 4 for t in marked)
    if euler < 0:
        return []
    options = []
    for i in range(4):
        left, right = 3*(d[i-1]-1)-d[i], 3*(d[(i+1)%4]-1)-d[i]
        options.append([(a,b) for a,b in product(range(left+1),range(right+1))
                        if a+b+h[i] >= 0])
    survivors = []
    for profile in product(*options):
        a,b = zip(*profile)
        if sum(a[t]*b[t] for t in range(4) if strong or t in marked) > euler:
            continue
        for i in range(4):
            j,k,f = (i+1)%4,(i+2)%4,d[(i+3)%4]
            A,B = a[i],b[k]
            p,l = d[i]+A,d[k]+B
            threshold = 4 + max(0,3*(f-3))
            rank_a = max(0,p-2*(f-2)) if l>threshold else 0
            rank_b = max(0,l-2*(f-2)) if p>threshold else 0
            defects = h[j]+a[j]+b[j]
            u = a[j]+max(0,4*A-defects)
            w = b[j]+max(0,4*B-defects)
            via_i = u if strong or j in marked else u-a[j] if i in marked else 0
            via_k = w if strong or j in marked else w-b[j] if k in marked else 0
            loss_i = 0 if strong or i in marked else A*b[i]
            loss_k = 0 if strong or k in marked else B*a[k]
            if via_i*b[j] > euler+A*(h[i]+a[i])+loss_i:
                break
            if a[j]*via_k > euler+B*(h[k]+b[k])+loss_k:
                break
            lower = u*w if strong or j in marked else u*w-(u-via_i)*(w-via_k)
            upper = (euler+A*(h[i]+a[i])+B*(h[k]+b[k])+6*A*B
                     -max(A*rank_a,B*rank_b)+loss_i+loss_k)
            if lower>upper:
                break
        else:
            survivors.append(profile)
    return survivors


def main():
    for e_size in (2,3):
        remaining = set()
        for x in range(2,3*(e_size-1)+1):
            for z in range(x,3*(e_size-1)+1):
                for y in range(2,3*(x-1)+1):
                    d=(e_size,x,y,z)
                    if surviving_extensions(d,{0}):
                        remaining.add(d)
        expected={(3,3,3,3)} if e_size==3 else {
            (2,2,2,2),(2,2,2,3),(2,2,3,2),(2,3,2,3),(2,3,3,3)}
        assert remaining==expected, remaining
    for d in ((2,2,2,3),(2,3,2,3),(2,3,3,3)):
        for marked in combinations(range(4),2):
            assert not surviving_extensions(d,set(marked)), (d,marked)
    print("Single four-matrix corner profiles and all six pairs of marked corners verified.")


if __name__=="__main__":
    main()
