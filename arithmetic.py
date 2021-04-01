

def _clean_pol(pol):
    l = pol.coefficients()
    while len(l)>0 and l[-1].contains_zero(): l.pop()
    return pol.parent()(l)

def GCD(a, b):

    """
    Return a *non-rigorous* gcd of the polynomials a and b.

    Note: this function is designed for BallField as base ring.

    Some words about the correction of this function:
    Let a· and b· be fixed. If a and b are precise enough, GCD(a, b) contains
    the gcd of a and b.

    INPUT:

    OUTPUT:
    """

    a, b = _clean_pol(a), _clean_pol(b)
    if a==0: return b
    if b==0: return a
    if a.degree() < b.degree(): return GCD(b, a)

    while b != 0:
        a, b = b, a.quo_rem(b)[1]
        a, b = _clean_pol(a), _clean_pol(b)

    return a
