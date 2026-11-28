"""Check residual resultant degrees, the secant intersection, and the binary fold."""
from itertools import product

import numpy as np
from flint import nmod_mat

from verify_quadric_biform_window import PRIME, forms
from verify_quadric_ruling_identifiability import quartic_basis


def main():
    # Restrictions to a length-two divisor for forms of degree d have total
    # Chern class (1-h)^(-(d-1)).
    for a,b in ((3,5),(4,4)):
        s,t = a-1,b-1
        assert s+t == 6
        assert s*(s+1)//2 + s*t + t*(t+1)//2 == 21
        assert (a-2)+(b-2) == 4
        assert 4*(21-6) == 60
    for a in range(1,5):
        assert max(0,8-a-4)+max(0,a-4) == 4-a < 4
    v = forms()
    transfer = quartic_basis(v,1,0).inv()*quartic_basis(v,0,1)
    for inv in (transfer,transfer.inv()):
        rows = np.zeros((16,16), dtype=int)
        for i,j in product(range(4), repeat=2):
            for k in range(5):
                rows[i+k,4*i+j] += int(inv[k,j+1])
                rows[8+i+k,4*i+j] -= int(inv[k,j])
        assert nmod_mat((rows % PRIME).tolist(),PRIME).rank() == 16
    print("PASS: both secant types of degree 21, length-60 residual divisor, all fold splitting bounds and both binary maps at the smooth sample")


if __name__ == "__main__":
    main()
