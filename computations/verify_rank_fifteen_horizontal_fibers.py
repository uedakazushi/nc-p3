"""The hyperelliptic models and degenerations of the nonzero horizontal fibers."""

import sympy as sp

x, t, c = sp.symbols("x t c")
HALF = sp.Rational(1, 2)


def verify():
    cubic = x**3-27*c
    quadratic = (1-x)**2-4*c
    assert sp.expand(sp.resultant(cubic, quadratic, x)-(c-1)**2*(1-64*c)) == 0
    x_at_t = 1-(t+4*c/t)*HALF
    q_at_t = (t-4*c/t)*HALF
    assert sp.cancel(1-x_at_t+q_at_t-t) == 0
    assert sp.cancel(q_at_t**2-quadratic.subs(x, x_at_t)) == 0
    sextic = (-t**2+2*t-4*c)**3-216*c*t**3
    branch = 2*t*sextic
    assert sp.cancel((4*t**2)**2*cubic.subs(x, x_at_t)-branch) == 0
    assert sp.cancel(branch.subs(t, 4*c/t)-256*c**4*branch/t**8) == 0
    discriminant = sp.factor(sp.discriminant(sextic, t))
    assert sp.cancel(discriminant/(c**10*(c-1)**2*(64*c-1))).is_Integer
    assert discriminant != 0
    for level, common_x, double_t in ((1, 3, -2), (sp.Rational(1, 64), sp.Rational(3, 4), sp.Rational(1, 4))):
        assert sp.gcd(cubic.subs(c, level), quadratic.subs(c, level)) == x-common_x
        specialized = sp.Poly(branch.subs(c, level), t)
        assert specialized.gcd(specialized.diff()).monic().as_expr() == t-double_t
        normalized = specialized.exquo(sp.Poly((t-double_t)**2, t))
        assert normalized.degree() == 5 and normalized.discriminant() != 0
        shift = sp.Poly((cubic*quadratic).subs(c, level), x)
        quotient = shift.exquo(sp.Poly((x-common_x)**2, x))
        assert quotient.degree() == 3 and quotient.discriminant() != 0
    print("Horizontal genera 3 and 2; both ordinary nodes and shift-quotient genera verified.")


if __name__ == "__main__":
    verify()
