from .precision_error import PrecisionError
from .polynomials import radical, XGCD



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
    for i in range(n-1):
        p = p.derivative()
        derivatives.append(p)

    for j, ev in enumerate(s):
        m = 1
        evaluations = [pol(ev) for pol in derivatives]
        while evaluations[m-1] == 0:
            m = m + 1
        s[j] = (ev, m)

    if n > sum(m for _, m in s):
        raise PrecisionError("Cannot compute multiplicities.")

    return s



def gen_eigspaces(mat, projections=False):

    r"""
    Return the generalized eigenspaces of "mat".

    Note: this function is designed for ComplexBallField as base ring.

    Some words about the correction of this algorithm:
    Let GenEigSpaces be the output for gen_eigspaces(mat, projections=True). For
    any mat· in mat, there is a selection [space1·, ..., spacek·] in
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
     - 'eigenvectors' : list of vectors
     - 'projection'   : polynomial (if 'projections=True' is specified).


    EXAMPLES:
    """

    n, C = mat.nrows(), mat.base_ring()

    s = eigenvalues(mat, multiplicities=True)

    GenEigSpaces = []





    """
    Pol, x = PolynomialRing(C, 'x').objgen()

    # call 'gen_eigspaces' function
    GenEigSpaces = gen_eigspaces(mat)
    s = len(GenEigSpaces)
    if s == 1: return [tuple(list(GenEigSpaces[0])+[Pol(1)])]

    # compute cofactors
    P = [(x-eig)**mult for eig, mult, _ in GenEigSpaces]
    Q = [prod(P[j] for j in range(s) if j != i) for i in range(s)]
    r = s//2
    D, U1, U2 = xgcd(sum(Q[:r]), sum(Q[r:]))
    U = [U1]*r + [U2]*(s-r)

    dec = [space + (u*q,) for space, u, q in zip(GenEigSpaces, U, Q)]






    # computation of a basis of the corresponding generalized eigenspace
    b = ((mat - ev*identity_matrix(C, nrows))**m).right_kernel().basis()
    if len(b) != m:
        raise PrecisionError("Cannot compute a basis of this generalized eigenspace.")
    b = [V(u) for u in b]

    GenEigSpaces.append((ev, m, b))"""

    return GenEigSpaces
