from .precision_error import PrecisionError
from .polynomials import radical, XGCD
from .reduction import ker
from sage.matrix.special import identity_matrix
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.misc_c import prod

def eigenvalues(mat, multiplicities=False):

    r"""
    Return the eigenvalues of "mat".

    Note: this function is designed for ComplexBallField as base ring.

    Some words about the correction of this algorithm:


    INPUT:

     -- "mat"            -- n×n matrix
     -- "multiplicities" -- boolean


    OUTPUT:

     -- "s" -- list of complex numbers

    If 'multiplicities=True' is specified, "s" is a list of couples (ev, m) with
    ev a complex number and m a positive integer.


    EXAMPLES:
    """

    p = mat.charpoly()
    try:
        s = radical(p).roots(multiplicities=False)
    except ValueError:
        raise PrecisionError("Cannot compute eigenvalues.")

    if not multiplicities: return s

    n = mat.nrows()
    derivatives = [p]
    for i in range(n):
        p = p.derivative()
        derivatives.append(p)

    for j, ev in enumerate(s):
        m = 1
        evaluations = [pol(ev) for pol in derivatives]
        while evaluations[m].contains_zero():
            m = m + 1
        s[j] = (ev, m)

    if sum(m for _, m in s)<n:
        raise PrecisionError("Cannot compute multiplicities.")

    return s



def gen_eigenspaces(mat, projections=False):

    r"""
    Return the generalized eigenspaces of "mat".

    Note: this function is designed for ComplexBallField as base ring.

    Some words about the correction of this algorithm:
    Let GenEigSpaces be the output for gen_eigenspaces(mat, projections=True).
    For any mat· in mat, there is a selection [space1·, ..., spacek·] in
    GenEigSpaces=[space1, ..., spacek] where each spacei· is a selection
    {'eigenvalue' : li·, 'multiplicity' : mi, 'eigenvectors' : bi·,
    'projection' : pi·} in spacei, such that:
    i)
    ii) the li· are pairwise disjoints,
    iii) the sum of the mi is equal to the dimension,
    iv) for each i, bi· is a list od mi linearly independent vectors.
    Reversely, let mat· be fixed. If mat is precise enough, no PrecisionError
    is raised and the spaces in gen_eigspaces(mat) correspond to the generalized
    eigenspaces of mat one-to-one.

    INPUT:

     -- "mat"         -- n×n matrix
     -- "projections" -- boolean


    OUTPUT:

     -- "GenEigSpaces" -- list of dictionary

    Each dictionary of "GenEigSpaces" represents a generalized eigenspace of
    "mat", whose keys are the following strings:
     - 'eigenvalue'   : complex number
     - 'multiplicity' : integer
     - 'basis' : list of vectors
     - 'projection'   : polynomial (if 'projections=True' is specified).


    EXAMPLES:
    """

    n, C = mat.nrows(), mat.base_ring()
    I = identity_matrix(C, n)
    Pol, x = PolynomialRing(C, 'x').objgen()

    s = eigenvalues(mat, multiplicities=True)

    if projections:
        k, l = len(s), len(s)//2
        P = [(x-ev)*m for ev, m in s]
        Q = [Pol(prod(P[j] for j in range(k) if j != i)) for i in range(k)]
        d, u1, u2 = XGCD(sum(Q[:l]), sum(Q[l:]))
        proj = [u*q for u, q in zip([u1]*l + [u2]*(k-l), Q)]

    GenEigSpaces = []
    for i, (ev, m) in enumerate(s):
        b = ker((mat - ev*I)**m)
        if len(b)<m:
            raise PrecisionError("Cannot compute a basis of this generalized eigenspace. ")
        space = {'eigenvalue' : ev, 'multiplicity' : m, 'basis' : b}
        if projections:
            space['projection'] = proj[i]
        GenEigSpaces.append(space)

    return GenEigSpaces
