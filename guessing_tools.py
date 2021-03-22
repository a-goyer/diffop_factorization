from sage.all import *

def hermite_pade_approximants (F, d):

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

        sage: from diffop_factorization.guessing_tools import hermite_pade_approximants
        sage: f = taylor(log(1+x), x, 0, 8).series(x).truncate().polynomial(QQ); f
        -1/8*x^8 + 1/7*x^7 - 1/6*x^6 + 1/5*x^5 - 1/4*x^4 + 1/3*x^3 - 1/2*x^2 + x
        sage: F = [f, f.derivative(), f.derivative().derivative()]
        sage: P = hermite_pade_approximants(F, 5); P
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


def guess_rational(x, p = 30):

    r"""
    Return the simplest rational number equals to x to the accuracy p.

    INPUT:
     - 'x' - a complex number or a list of complex numbers

    OUTPUT:
     - 'r' - a rational number or a list of rational numbers

    EXAMPLES::

        sage: from diffop_factorization.guessing_tools import guess_rational
        sage: a = CC(sqrt(2))^2; a == 2
        False
        sage: guess_rational(a)
        2
        sage: b = a + 2^(-40); b
        2.00000000000091
        sage: guess_rational(b)
        2
        sage: c = a + 2^(-20); c
        2.00000095367432
        sage: guess_rational(c)
        2095107/1047553
        sage: guess_rational(c, 20)
        2

    """

    if isinstance(x, list) :
        r = [guess_rational(c) for c in x]
        return r

    im = x.imag_part()
    eps = 2.**(-p)
    if not im <= eps:
        raise ValueError('This number does not seem a rational number.')

    x = x.real_part()
    r = x.nearby_rational(max_error = eps)

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
