

def _clean(pol):

    l = pol.coefficients()
    while len(l)>0 and l[-1].contains_zero(): l.pop()
    cpol = pol.parent()(l)

    return cpol



def GCD(a, b):

    r"""
    Return a *non-rigorous* gcd of the polynomials "a" and "b".

    Note: this function is designed for BallField as base ring.

    Some words about the correction of this function:
    Let a· and b· be fixed. If a and b are precise enough, GCD(a, b) contains
    the gcd of a and b.

    INPUT:

     -- "a" -- polynomial
     -- "b" -- polynomial

    OUTPUT:

     -- "a" -- polynomial


    EXAMPLE::

        sage: from diffop_factorization.polynomials import GCD
        sage: P.<x> = CBF[]; a = CBF(pi)
        sage: p, q = (x-1)*(x-2)**2, (x-2)*(x-3)**2
        sage: p, q = p(x*a), q(x*a)
        sage: d = GCD(p, q); d(x/a).monic()
        ([1.0000000000 +/- 1.34e-12])*x + [-2.00000000000 +/- 1.94e-12]

    """

    a, b = _clean(a), _clean(b)
    if a==0: return b
    if b==0: return a
    if a.degree() < b.degree(): return GCD(b, a)

    while b != 0:
        a, b = b, a.quo_rem(b)[1]
        b = _clean(b)

    return a



def XGCD(a, b):

    r"""
    Return the *non-rigorous* monic gcd of the polynomials "a" and "b" and the
    coefficients in the Bezout identity.

    Note: this function is designed for BallField as base ring.

    Some words about the correction of this function:
    Let a· and b· be fixed.

    INPUT:

     -- "a" -- polynomial
     -- "b" -- polynomial

    OUTPUT:

     -- "d" -- polynomial
     -- "u" -- polynomial
     -- "v" -- polynomial


    EXAMPLE::

        sage: from diffop_factorization.polynomials import XGCD, _clean
        sage: P.<x> = CBF[]; a = CBF(pi)
        sage: p, q = (x-1)*(x-2)**2, (x-2)*(x-3)**2
        sage: p, q = p(x*a), q(x*a)
        sage: d, u, v = XGCD(p, q); d, 2/a
        (([1.0000000000 +/- 1.33e-12])*x + [-0.63661977237 +/- 3.04e-12],
         [0.636619772367581 +/- 4.28e-16])
        sage: _clean(u*p + v*q)
        ([1.000000000 +/- 3.39e-11])*x + [-0.63661977237 +/- 6.95e-12]

    """

    P = a.parent()

    a, b = _clean(a), _clean(b)
    if a==0: return b, P.zero(), P.one()
    if b==0: return a, P.one(), P.zero()
    if a.degree() < b.degree():
        d, v, u = XGCD(b, a)
        return d, u, v

    r0, u0, v0, r1, u1, v1 = a, P.one(), P.zero(), b, P.zero(), P.one()
    while r1!=0:
        r0, (q, r1) = r1, r0.quo_rem(r1)
        u0, v0, u1, v1 = u1, v1, u0 - q*u1, v0 - q*v1
        r1 = _clean(r1)

    lc = r0.leading_coefficient()
    d, u, v = r0.monic(), _clean(u0/lc), _clean(v0/lc)

    return d, u, v



def radical(pol):

    r"""
    Return the *non-rigorous* radical of the polynomial "pol".

    Note: this function is designed for BallField as base ring.

    Some words about the correction of this function:
    Let pol· be fixed. If pol is precise enough, radical(pol)) contains the
    radical of pol·.


    INPUT:

     -- "pol" -- polynomial


    OUTPUT:

     -- "rad" -- polynomial


    EXAMPLE::

        sage: from diffop_factorization.polynomials import radical
        sage: P.<x> = CBF[]; a = CBF(pi)
        sage: p = (x-1)*(x-2)**2
        sage: p = p(x*a).monic(); p
        ([1.0000000000000 +/- 1.57e-15])*x^3 + ([-1.59154943091895 +/- 5.18e-15])*x^2 + ([0.81056946913870 +/- 2.96e-15])*x + [-0.129006137732798 +/- 1.12e-16]
        sage: rad = radical(p); rad
        ([-44.4132198049 +/- 2.39e-11])*x^2 + ([42.4115008235 +/- 7.69e-11])*x + [-9.0000000000 +/- 3.38e-11]
        sage: rad.roots(multiplicities=False)
        [[0.3183098862 +/- 2.31e-11] + [+/- 6.86e-12]*I,
         [0.6366197724 +/- 4.20e-11] + [+/- 9.55e-12]*I]
        sage: [1/a, 2/a]
        [[0.318309886183791 +/- 4.43e-16], [0.636619772367581 +/- 4.28e-16]]

    """

    d = GCD(pol, pol.derivative())
    rad = _clean(pol.quo_rem(d)[0])

    return rad
