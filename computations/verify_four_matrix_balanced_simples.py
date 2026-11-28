"""Euler and obstruction dimensions in the nonsimple endpoint calculation."""


def euler(u, v):
    return sum(x*(2*v[i]-4*v[i-1]+6*v[(i+2)%4]-4*v[(i+1)%4]) for i,x in enumerate(u))


def main():
    s = (2,2,2,3)
    u, v = (6,2,2,3), (2,2,6,3)
    assert euler(s,s) == 2 and euler(u,v) == 66
    assert [(a,b) for a in range(5) for b in range(2)
            if a+b>=4 and a*b<=2] == [(4,0)]
    assert [(a,b) for a in range(2) for b in range(2)
            if (10-b)*b<=2 and a*(10-a)<=2] == [(0,0)]
    assert (4*4-6)**2 == 100 > euler(u,v)+2*(4*4) == 98
    print("Universal head and socle extensions give incompatible second Ext bounds 100 and 98.")


if __name__ == "__main__":
    main()
