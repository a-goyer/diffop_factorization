from sage.rings.real_mpfr import RealField
from sage.rings.rational_field import QQ
from sage.rings.qqbar import QQbar
from sage.matrix.constructor import matrix
from sage.matrix.matrix_dense import Matrix_dense
from sage.modules.free_module_element import FreeModuleElement_generic_dense
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.functions.other import floor
from sage.functions.log import log



def hp_approx (F, d):

    r"""
    Return the Hermite-Padé approximants of F at order d.

    Let F = [f1, ..., fm]. This function returns a list of polynomials P =
    [p1, ..., pm] such that:
    - max(deg(p1), ..., deg(pm)) is minimal,
    - p1*f1 + ... + pm*fm = O(x^d).

    Note that this function calls some methods of the Library of Polynomial
    Matrices, see https://github.com/vneiger/pml to install it.

    INPUT:
     - "F" - a list of polynomials or series

    OUTPUT:
     - "P" - a list of polynomials

    EXAMPLES::

        sage: from diffop_factorization.guessing_tools import hp_approx
        sage: f = taylor(log(1+x), x, 0, 8).series(x).truncate().polynomial(QQ); f
        -1/8*x^8 + 1/7*x^7 - 1/6*x^6 + 1/5*x^5 - 1/4*x^4 + 1/3*x^3 - 1/2*x^2 + x
        sage: F = [f, f.derivative(), f.derivative().derivative()]
        sage: P = hp_approx(F, 5); P
        (0, 1, x + 1)
        sage: from ore_algebra import OreAlgebra
        sage: Pol.<x> = QQ[]; OA.<Dx> = OreAlgebra(Pol)
        sage: diffop = OA(list(P)); diffop
        (x + 1)*Dx^2 + Dx
        sage: diffop(log(1+x))
        0

    """

    try:
        F = [f.truncate() for f in F]
    except: pass

    mat = matrix(len(F), 1, F)
    basis = mat.minimal_approximant_basis(d)
    rdeg = basis.row_degrees()
    i = min(range(len(rdeg)), key = lambda i: rdeg[i])
    P = basis[i]

    return P



def guess_rational(x, p=None):

    r"""
    Guess rational coefficients for a vector or a matrix or a polynomial or a
    list or just a complex number.

    Note: this function is designed for ComplexOptimisticField as base ring.

    INPUT:
     - 'x' - object with approximate coefficients

    OUTPUT:
     - 'r' - object with rational coefficients

    EXAMPLES::

        sage: from diffop_factorization.complex_optimistic_field import ComplexOptimisticField
        sage: from diffop_factorization.guessing import guess_rational
        sage: C = ComplexOptimisticField(30, 2^-10)
        sage: a = 1/3 - C(1+I)*C(2^-20)
        sage: P.<x> = C[]; pol = (1/a)*x + a; pol
        ([3.0000086 +/- 2.86e-8] + [8.5831180e-6 +/- 7.79e-14]*I)*x + [0.333332379 +/- 8.15e-10] - [9.53674316e-7 +/- 4.07e-16]*I
        sage: guess_rational(pol)
        3*x + 1/3

    """

    if isinstance(x, list) :
        r = [guess_rational(c, p=p) for c in x]
        return r

    if isinstance(x, FreeModuleElement_generic_dense) or \
    isinstance(x, Matrix_dense) or isinstance(x, Polynomial):
        r = x.parent().change_ring(QQ)(guess_rational(x.list(), p=p))
        return r

    if p is None:
        eps = x.parent().eps
        p = floor(-log(eps, 2))
    else:
        eps = RealField(30).one() >> p

    if not x.imag().above_abs().mid()<eps:
        raise ValueError('This number does not seem a rational number.')

    x = x.real().mid()
    r = x.nearby_rational(max_error=x.parent()(eps))

    return r



def guess_algebraic(x, p = 30, d = 3):

    r"""
    Return the simplest algebraic number of degree at most d equals to x to the
    accuracy p.

    INPUT:
     - 'x' - a complex number or a list of complex numbers

    OUTPUT:
     - 'a' - an algebraic number or a list of algebraic numbers

    EXAMPLES::

        sage: from diffop_factorization.guessing_tools import guess_algebraic
        sage: a = CC(sqrt(2))
        sage: guess_algebraic(a)
        1.414213562373095?
        sage: _.minpoly()
        x^2 - 2

    """

    if isinstance(x, list) :
        a = [guess_algebraic(c) for c in x]
        return a

    minpol = algdep(x, degree = d, known_bits = p)
    roots = minpol.roots(QQbar, multiplicities=False)
    l = [r-x for r in roots]
    i = min(range(len(l)), key = lambda i: abs(l[i])); i
    a = roots[i]

    return a
