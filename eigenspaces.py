from .precision_error import PrecisionError
from .polynomials import radical, XGCD
from .reduction import ker
from sage.matrix.special import identity_matrix
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.misc_c import prod

def eigenvalues(mat, *, multiplicities=False):

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


    EXAMPLE::

        sage: from diffop_factorization.eigenspaces import eigenvalues
        sage: r2, r3 = CBF(sqrt(2)), CBF(sqrt(3))
        sage: mat = jordan_block(r2, 2).block_sum(matrix([r3]))
        sage: ran = MatrixSpace(CC, 3).random_element().change_ring(CBF)
        sage: mat = ~ran * mat * ran
        sage: eigenvalues(mat, multiplicities=True)
        [([1.4142136 +/- 5.47e-8] + [+/- 1.71e-8]*I, 2),
         ([1.7320508 +/- 2.65e-8] + [+/- 1.89e-8]*I, 1)]

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



def gen_eigenspaces(mat, *, projections=False):

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
     - 'basis'        : list of vectors
     - 'projection'   : polynomial (if 'projections=True' is specified).


    EXAMPLES:

    A generic example ::

        sage: from diffop_factorization.eigenspaces import gen_eigenspaces
        sage: mat = MatrixSpace(CC, 3).random_element().change_ring(CBF)
        sage: gen_eigenspaces(mat)
        [{'eigenvalue': [0.008731294334 +/- 5.79e-13] + [-0.457554487127 +/- 3.74e-13]*I,
          'multiplicity': 1,
          'basis': [([-0.081853906993 +/- 6.52e-13] + [-0.087179008676 +/- 7.81e-13]*I, 1.000000000000000, [-0.342339920094 +/- 3.62e-13] + [0.098986705031 +/- 3.30e-13]*I)]},
         {'eigenvalue': [0.339079560015 +/- 4.32e-13] + [0.888549954437 +/- 3.85e-13]*I,
          'multiplicity': 1,
          'basis': [([-0.409699282668 +/- 6.41e-13] + [0.247473920868 +/- 5.86e-13]*I, [-0.066871473678 +/- 7.28e-13] + [0.685673697712 +/- 7.82e-13]*I, 1.000000000000000)]},
         {'eigenvalue': [0.606834685948 +/- 7.61e-13] + [-0.992982372445 +/- 8.72e-13]*I,
          'multiplicity': 1,
          'basis': [([-0.63246660126 +/- 2.18e-12] + [-0.05887072347 +/- 3.44e-12]*I, 1.000000000000000, [-0.59541382749 +/- 5.20e-12] + [0.10852249946 +/- 5.37e-12]*I)]}]

    An example with a multiple eigenvalue ::

        sage: from diffop_factorization.eigenspaces import gen_eigenspaces
        sage: r2, r3 = CBF(sqrt(2)), CBF(sqrt(3))
        sage: mat = jordan_block(r2, 2).block_sum(matrix([r3]))
        sage: ran = MatrixSpace(CC, 3).random_element().change_ring(CBF)
        sage: mat = ~ran * mat * ran
        sage: GenEigSpaces = gen_eigenspaces(mat, projections=True)
        sage: [(space['eigenvalue'], space['multiplicity']) for space in GenEigSpaces]
        [([1.41421356 +/- 7.84e-9] + [+/- 5.47e-9]*I, 2),
         ([1.73205081 +/- 8.48e-9] + [+/- 6.05e-9]*I, 1)]
        sage: ev, vec = GenEigSpaces[1]['eigenvalue'], GenEigSpaces[1]['basis'][0]
        sage: (mat - ev*identity_matrix(CBF, 3))*vec
        ([+/- 1.37e-6] + [+/- 1.37e-6]*I, [+/- 4.33e-7] + [+/- 4.33e-7]*I, [+/- 1.17e-6] + [+/- 1.17e-6]*I)
        sage: T = matrix(GenEigSpaces[0]['basis'] + GenEigSpaces[1]['basis']).transpose()
        sage: pol = GenEigSpaces[0]['projection']
        sage: P = ~T * pol(mat) * T; P
        [  [1.00 +/- 1.19e-4] + [+/- 1.19e-4]*I        [+/- 8.62e-5] + [+/- 8.62e-5]*I        [+/- 1.53e-4] + [+/- 1.53e-4]*I]
        [       [+/- 3.97e-5] + [+/- 3.97e-5]*I [1.0000 +/- 4.89e-5] + [+/- 4.89e-5]*I        [+/- 5.52e-5] + [+/- 5.52e-5]*I]
        [       [+/- 9.55e-5] + [+/- 9.55e-5]*I        [+/- 7.19e-5] + [+/- 7.19e-5]*I        [+/- 1.25e-4] + [+/- 1.25e-4]*I]

    """

    n, C = mat.nrows(), mat.base_ring()
    I = identity_matrix(C, n)
    Pol, x = PolynomialRing(C, 'x').objgen()

    s = eigenvalues(mat, multiplicities=True)

    if projections:
        k, l = len(s), len(s)//2
        P = [(x-ev)**m for ev, m in s]
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
