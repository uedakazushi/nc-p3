"""All regular fibers of the first rectangle on u=1,z=-1 and its rotations and twists."""

from flint import fmpz_poly
from flint_rational_function import rational_function_field
from verify_first_rectangle_edge_center import edge_rules
from verify_parallel_curve_residues import residue_test
from verify_parallel_invertible_rectangles import K,s
from verify_rectangular_surface_cut_charts import chart,relation_cut_data


def check():
    v = s.Symbol('v')
    field = rational_function_field(('v',))
    raw = edge_rules(v,-1,field.from_sympy,(0,2,2,0))
    data = relation_cut_data(raw,field,field.from_sympy(v*(v*v-1)).n)
    factors=[]
    for p in K:
        polys = [s.Poly(s.sympify(f.replace('^','**')),v) for f in chart(p,data)]
        factors.append({tuple(map(int,reversed(f.all_coeffs()))) for f in polys})
    assert list(map(len,factors))==[196,183,185,175]
    def reverse(coeff):
        return tuple(z if coeff[0]>0 else -z for z in reversed(coeff))
    bad = [factors[p]&{reverse(f) for f in factors[p+2]} for p in (0,1)]
    assert [(len(B),sum(len(f)-1 for f in B)) for B in bad]==[(14,31),(19,38)]
    combined = set.union(*bad)
    assert len(combined)==26 and max(len(f)-1 for f in combined)==4
    def modular(value,prime):
        def convert(z):
            n,d = s.fraction(z)
            return int(n)*pow(int(d),-1,prime)%prime
        return edge_rules(s.Integer(value),-1,convert,(0,0,0,0))
    for coeff in sorted(combined):
        f = fmpz_poly(list(coeff))
        for P in s.primerange(1009,20000):
            if not int(f.leading_coefficient())%P:
                continue
            roots = [i for i in range(2,P-1) if int(f(i))%P==0]
            if not roots:
                continue
            try:
                residue_test(roots[0],P,modular)
            except (ValueError,AssertionError):
                continue
            break
        else:
            raise AssertionError(('no successful finite-place certificate',coeff))
    print('PASS: every generic cut and all 26 exceptional number fields on the rank-ten curve')


if __name__=='__main__':
    check()
