from sage.rings.rational_field import QQ
from sage.rings.qqbar import QQbar, number_field_elements_from_algebraics
from sage.arith.functions import lcm
from sage.arith.misc import valuation
from ore_algebra import DifferentialOperators
from ore_algebra.analytic.differential_operator import DifferentialOperator
from sage.rings.polynomial.polynomial_ring import PolynomialRing_field
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing



def with_polynomial_coefficients(dop):

    """
    Return a polynomial q and a differential operator L with polynomial
    coefficients such that q*dop = L.
    """

    if not isinstance(dop.base_ring(), PolynomialRing_field):
        q = lcm([c.denominator() for c in dop])
        L = dop.parent().change_ring(dop.base_ring().base())(q*dop)
        return q, L
    else:
        return dop.base_ring().one(), dop



def is_fuchsian(dop):

    r"""
    Return True if "dop" is fuchian, False if not.

    Use Fuch's criterion: Let L = a_n*Dz^n + ... + a_0 be an operator whose
    coefficients are polynomial and p be a finite singular point of L. Then p is
    regular iff no (z-p)^{n-k}*a_k/a_n admits p as pole.

    """

    z, n = dop.base_ring().gen(), dop.order()

    dop = with_polynomial_coefficients(dop)[1]

    coeffs = dop.coefficients()
    fac = coeffs.pop().factor()
    for (f, m) in fac:
        for k, ak in enumerate(coeffs):
            mk = valuation(ak, f)
            if mk - m < k - n: return False

    dop = dop.annihilator_of_composition(1/z)
    for k, frac in enumerate(dop.monic().coefficients()[:-1]):
        d = (z**(n-k)*frac).denominator()
        if d(0)==0: return False

    return True


def _E(dop, S=False):

    """
    Compute the largest modulus of the local exponents of "dop" at infinity and
    at its finite non-apparent singular points. See [BRS-19-EDB].
    With low probability, the result can be bigger than expected.
    Assumption: "dop" has only polynomial coefficients.
    """

    if not isinstance(dop.base_ring(), PolynomialRing_field):
        raise TypeError("The operator must to have polynomial coefficients.")

    z = dop.base_ring().gen()
    pol = dop.indicial_polynomial(1/z)
    E = max(abs(r) for r in pol.roots(QQbar, multiplicities=False)).ceil()

    sing = DifferentialOperator(dop.desingularize())._singularities(QQbar)
    singQQ = [QQ(x) for x in sing if x in QQ]
    singQQbar = [x for x in sing if x not in QQ]

    for s in singQQ:
        pol = dop.indicial_polynomial(z - s)
        m = max(abs(r) for r in pol.roots(QQbar, multiplicities=False)).ceil()
        E = max(E, m)

    for s in singQQbar:
        K, Ks, mor = number_field_elements_from_algebraics(s)
        Dops, z, Dz = DifferentialOperators(K, 'z')
        pol = Dops(dop).indicial_polynomial(z - Ks) # improvement: same computation for all conjugates
        pol = pol.parent().change_ring(QQbar)([mor(c) for c in pol])
        m = max(abs(r) for r in pol.roots(QQbar, multiplicities=False)).ceil()
        E = max(E, m)

    if not S:
        return E
    else:
        return E, len(sing)


def explicit_degree_bound(dop, fuchsian=None):
    """ Assumptions: Ground field is QQ, "dop" has only polynomial coefficients
    and is Fuchsian. """

    if not isinstance(L.base_ring(), PolynomialRing_field):
        raise TypeError("The operator must to have polynomial coefficients.")

    if fuchsian==None:
        fuchsian = is_fuchsian(dop)

    E, S = _E(dop, S=True)

    if fuchsian:
        db = lambda r: (r**2)*(S + 1)*E + r*S + (1/2)*(r**2)*(r - 1)*(S - 1)

    return db
