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

    m = len(F)
    K = F[0].base_ring()

    Pol, x = PolynomialRing(K, 'x').objgen()
    F = [Pol(f) for f in F]

    mat = matrix(m, 1, F)
    basis = mat.minimal_approximant_basis(d)
    rdeg = basis.row_degrees()
    i = min(range(len(rdeg)), key = lambda i: rdeg[i])
    P = basis[i]

    return P
