from sage.all import *
from sage.rings.complex_field import ComplexField_class



# This module provides some functions for linear algebra with ComplexField with
# an approximate zero test as base field.
# For more information about this field, see:
# https://git.sagemath.org/sage.git/commit/?h=u/mmezzarobba/CC_eps&id=442837534e627be3047f3faf1a11501b614b15ae
# Here an example of use:
#    sage: from sage.rings.complex_mpfr import ComplexField_class
#    sage: my_CC = ComplexField_class(eps=1e-10)
#    sage: a = my_CC(1)
#    sage: a == a + 1e-12
#    True





def gen_eigspaces (mat) :

    """
    Return the generalized eigenspaces (eigval, mult, basis) of "mat".


    This function uses the square-free decomposition of the characteristic
    polynomial of "mat" to compute the eigenvalues (and corresponding
    multiplicities) with a good number of correct decimals. No guarantees are
    made about the accuracy of the output even though it is good in practice.

    Note that the bases of the generalized eigenspaces are not computed such
    that the matrix "mat" written with respect to these bases is triangular.


    INPUT :

     -- "mat" -- an n×n matrix

    OUTPUT:

     -- "GenEigSpaces" -- a list of triplets ("ev", "m", "b") with:
                           - "ev" - a complex number
                           - "m"  - a positive integer
                           - "b"  - a list of "m" vectors


    EXAMPLE:

        sage: C = ComplexField_class(30, eps = 2**(-15))
        sage: r2, r3 = C(sqrt(2)), C(sqrt(3))
        sage: mat = matrix(C, [[r2, 1, 0], [0, r2, 0], [0, 0, r3]]); mat
        [ 1.4142136  1.0000000 0.00000000]
        [0.00000000  1.4142136 0.00000000]
        [0.00000000 0.00000000  1.7320508]
        sage: T = matrix(C, 3, [CC.random_element() for i in range(9)])
        sage: mat = ~T * mat * T; mat
        [    2.0104457 - 1.3085006*I    -1.4221239 + 1.2637281*I    1.3569532 + 0.86870050*I]
        [-0.026154678 + 0.65325782*I    2.3127517 - 0.75237064*I  -0.73683954 - 0.30247458*I]
        [   1.9239033 + 0.57314357*I    -1.9719328 - 2.4405489*I    0.23728051 + 2.0608712*I]
        sage: GenEigSpaces = gen_eigspaces(mat); GenEigSpaces
        [(1.4142147 - 1.8442147e-8*I,
          2,
          [(1.0000000 - 1.1641532e-10*I, -7.1150768e-20 + 3.4924597e-10*I, 2.0671185 - 3.6167470*I),
           (0.00000000, 1.0000000 + 9.3132257e-10*I, 6.1766047 - 9.5560085*I)]),
         (1.7320508 + 8.0330744e-8*I,
          1,
          [(1.0000000 - 9.3132257e-10*I, -0.61437803 + 0.58751076*I, 0.76313802 + 1.6636404*I)])]
        sage: GenEigSpaces[0][0] == r2
        True


    """


    nrows, ncols = mat.dimensions()

    C = mat.base_ring()
    V = VectorSpace(C, nrows)

    # computation of eigenvalues
    charpol = mat.charpoly()
    p = charpol.derivative()
    eigvals = (charpol//gcd(charpol, p)).roots(multiplicities=False)

    # the list of derivatives of the characteristic polynomial to detect multiplicities
    derivatives_of_charpol = [p]
    for i in range(nrows-1) :
        p = p.derivative()
        derivatives_of_charpol.append(p)

    GenEigSpaces = []
    for ev in eigvals :
        # computation of the multiplicity of ev
        m = 1
        evaluations = [pol(ev) for pol in derivatives_of_charpol]
        while evaluations[m-1] == 0: m = m + 1

        # computation of a basis of the corresponding generalized eigenspace
        b = ((mat - ev*identity_matrix(C, nrows))**m).right_kernel().basis()
        if len(b) != m:
            raise FloatingPointError("Cannot compute a basis of the generalized eigenspace.")
        b = [V(u) for u in b]

        GenEigSpaces.append((ev, m, b))

    return GenEigSpaces





def gen_eigspaces_with_proj (mat) :

    """
    Return the generalized eigenspaces with projections as polynomials.


    This function returns a list of 4-tuples of the form (eigenvalue,
    multiplicity, basis, projection) which represent generalized eigenspaces
    one-to-one. The polynomial "projection" evaluated at the matrix "mat" is
    the projection onto the corresponding generalized eigenspace along the
    others.

    No guarantees are made about the accuracy of the output.


    INPUT :
     -- "mat" -- an n×n matrix

    OUTPUT:
     -- "dec" -- a list of 4-tuples ("ev", "m", "b", "P") with:
                  - "ev" - a complex number
                  - "m"  - a positive integer
                  - "b"  - a list of "m" vectors
                  - "P"  - a polynomial with complex coefficients


    EXAMPLE:

        sage: C = ComplexField_class(30, eps = 2**(-15))
        sage: r2, r3 = C(sqrt(2)), C(sqrt(3))
        sage: mat = matrix(C, [[r2, 1, 0], [0, r2, 0], [0, 0, r3]]); mat
        [ 1.4142136  1.0000000 0.00000000]
        [0.00000000  1.4142136 0.00000000]
        [0.00000000 0.00000000  1.7320508]
        sage: T = matrix(C, 3, [CC.random_element() for i in range(9)])
        sage: mat = ~T * mat * T; mat
        [   1.6663521 - 0.44523535*I  0.72289035 - 0.054920460*I  0.44581260 + 0.074757800*I]
        [ 0.015141091 + 0.23939887*I   1.3881664 + 0.044850123*I -0.23020840 + 0.042792508*I]
        [ 0.37761244 + 0.070053958*I  -0.15142170 + 0.64512325*I    1.5059595 + 0.40038523*I]
        sage: dec = gen_eigspaces_with_proj(mat); dec
        [(1.4142140 + 1.2668544e-7*I,
          2,
          [(1.0000000, 0.00000000, 0.014412320 + 1.0637495*I),
           (0.00000000, 1.0000000, 0.58391053 + 0.43399576*I)],
          (-9.8990061 - 7.7844543e-6*I)*x^2 + (27.998626 + 0.000024525888*I)*x - 18.798024 - 0.000019115936*I),
         (1.7320508 + 1.7139994e-9*I,
          1,
          [(1.0000000, 0.29724229 + 0.87487227*I, -0.49280685 - 0.30065441*I)],
          (9.8990061 + 7.7844543e-6*I)*x^2 + (-27.998626 - 0.000024525888*I)*x + 19.798024 + 0.000019115936*I)]
        sage: adapted_basis, proj = dec[0][2] + dec[1][2], dec[0][3]
        sage: P = matrix(C, adapted_basis).transpose()
        sage: proj(P**(-1) * mat * P) == matrix(C, [[1, 0, 0], [0, 1, 0], [0, 0, 0]])
        True


    """


    C = mat.base_ring()
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

    return dec





def approx_echelon_form (mat, *, rank=False, take_over=False) :

    """
    Return the row echelon form of "mat" and the transition matrix.


    The row echelon form is 'almost reduced': the pivots are set to 1 and the
    coefficients below the pivots are set to 0 but the coefficients above the
    pivots may be nonzero.

    If the keyword option "rank=True" is specified, this functions returns
    the rank of "mat" (computed as the number of used pivots).

    If the keyword option "take_over=True" is specified, this function begins
    with checking if the columns has already been reduced. For example, if
    "mat" is already in almost reduced row echelon form, no transformation is
    made (just some checkings). Moreover, it stores and returns the indices of
    the rows which have been modified.

    Denoting by "mat_REF" the reduced row echelon form of "mat" and by "T" the
    transition matrix, we have: "T" × "mat" = "E".

    No guarantees are made about the accuracy of the output.


    INPUT:

     -- "mat"       -- an n×p matrix
     -- "rank"      -- a boolean (optional, default: 'False')
     -- "take_over" -- a boolean (optional, default: 'False')

    OUTPUT:

     -- "mat_REF"       -- an n×p matrix
     -- "T"             -- an n×n matrix
     -- "r"             -- an integer (only if "rank=True" is specified)
     -- "modified_rows" -- a list of integers (only if "take_over=True" is specified)


    EXAMPLE:

        sage: C = ComplexField_class(20, eps = 2**(-10))
        sage: M = MatrixSpace(C, 3)
        sage: mat = M.random_element(); mat
        [ -0.53347 - 0.58702*I   0.10785 - 0.96536*I -0.051985 + 0.13611*I]
        [  0.44753 + 0.22854*I   0.78818 + 0.69078*I  0.097000 + 0.92620*I]
        [ -0.69286 - 0.38037*I  -0.71133 - 0.48489*I -0.030800 - 0.34716*I]
        sage: mat_REF, T = approx_echelon_form(mat)
        sage: T * mat == mat_REF
        True
        sage: all(mat_REF[i,j] == M.one()[i,j] for j in range(3) for i in range(3) if j <= i)
        True

    """

    nrows, ncols = mat.dimensions()
    K = (mat[0,0]).parent()

    T = identity_matrix(K, nrows)

    positions_of_pivots = []

    if take_over: modified_rows = []

    r = 0 # number of used pivots
    j = 0 # index of current column
    while j < ncols and r < nrows:

        C = T * vector(mat[:,j])

        if take_over and mat[r,j] == K(1):
            # put zeros under the pivot
            l = r + 1
            while l < nrows and mat[l,j] == K(0): l = l + 1

            for i in range(l, nrows):
                T[i] = [T[i,k] - C[i]*T[r,k] for k in range(nrows)]
                if not i in modified_rows: modified_rows.append(i)

            positions_of_pivots.append(j)
            r = r + 1

        else:
            # find a good pivot in the current column
            maximum = RealField(K.precision())(0)
            for i in range(r, nrows) :
                coeff = C[i]
                if not coeff == 0 :
                    x = coeff.abs()
                    if x > maximum :
                        pivot_row = i
                        maximum = x

            if maximum > 0 :

                # make the operations only on T
                if r != pivot_row :
                    T[r], T[pivot_row] = T[pivot_row], T[r]
                    C[r], C[pivot_row] = C[pivot_row], C[r]
                T[r] = [T[r,k]/C[r] for k in range(nrows)]
                if take_over and not r in modified_rows: modified_rows.append(r)

                for i in range(r+1, nrows) :
                    T[i] = [T[i,k] - C[i]*T[r,k] for k in range(nrows)]
                    if take_over and not i in modified_rows: modified_rows.append(i)

                positions_of_pivots.append(j)
                r = r + 1

        j = j + 1

    # compute mat_REF from T
    mat_REF = T * mat
    for j, p in enumerate(positions_of_pivots):
        mat_REF[j,p] = K(1)
        for i in range(j+1, nrows) : mat_REF[i,p] = K(0)

    if det(T) == 0 : raise FloatingPointError("Cannot compute an invertible transition matrix.")

    if rank and take_over:
        return mat_REF, T, r, modified_rows
    elif rank:
        return mat_REF, T, r
    elif take_over:
        return mat_REF, T, modified_rows

    return mat_REF, T





def Inv_with_T (listM, v) :

    """
    Return the smallest subspace containing "v" and invariant under "listM".


    This function returns a basis of the smallest subspace which contains "v"
    and which is invariant under the action of all the matrices of "listM" and
    the matrices "transition_matrices" which pass from "v" to the vectors of
    this basis.

    Write "basis" = [v1, ..., vs] and "transition_matrices" = [T1, ..., Ts]. We
    have that the Tj are polynomials in the matrices of "listM" and
    vj = Tj * v.

    No guarantees are made about the accuracy of the output.


    INPUT:

     -- "listM"          -- a list of n×n matrices
     -- "v"              -- a vector of size n

    OUTPUT:

     -- "basis"           -- a list of vectors of size n
     -- "transition_matrices" -- a list of n×n matrices


    EXAMPLE:

        sage: C = ComplexField_class(20, eps = 2**(-10))
        sage: mat = matrix(C, [[1, 1, 0], [0, 1, 1], [0, 0, 1]])
        sage: u, v = vector(C, [1, 0, 0]), vector(C, [0, 1, 0])
        sage: Ran = matrix(C, 3, [CC.random_element() for i in range(9)])
        sage: u, v, mat = Ran * u, Ran * v, Ran * mat * (~Ran)
        sage: len(Inv_with_T([mat], u)[0]), len(Inv_with_T([mat], v)[0])
        (1, 2)
        sage: [[v0, v1], [T0, T1]] = Inv_with_T([mat], v)
        sage: (T0 * v - v0).norm(), (T1 * v - v1).norm()
        (6.8480e-7, 3.1720e-6)


    """

    nrows, ncols = listM[0].dimensions()
    n = len(v)
    C = v.parent().base()

    basis = [v]
    modified_vectors = [v]

    transition_matrices = [identity_matrix(C, n)]
    modified_transition_mat = [identity_matrix(C, n)]

    old_dim, dim = 0, 1
    while old_dim < dim and dim < n :

        # add the products "matrix × vector"
        basis.extend([M*u for M in listM for u in modified_vectors])
        transition_matrices.extend([M*T for M in listM for T in modified_transition_mat])

        # keep only a list of linearly independent vectors
        mat = matrix(C, basis)
        mat, T, r , modified_rows = approx_echelon_form(mat, rank=True, take_over=True)
        basis = list(mat[:r])

        modified_vectors = [basis[i] for i in range(r) if i in modified_rows]

        l = len(transition_matrices)
        transition_matrices = [sum(T[i,j]*transition_matrices[j] for j in range(l)) for i in range(r)]
        modified_transition_mat =[transition_matrices[i] for i in range(r) if i in modified_rows]

        old_dim, dim = dim, r

    return basis, transition_matrices





def unique_eigenvalue(mat):

    """
    Return the unique eigenvalue of "mat".


    Assumption: the matrix "mat" has only one eigenvalue.

    This function returns the unique eigenvalue of "mat" if "mat" or raises a
    FloatingPointError if "mat" seems have multiple eigenvalues.


    INPUT:

     -- "mat" -- a square matrix

    OUTPUT:

     -- "eigval" -- a complex number


    EXAMPLES:

    An example with a matrix which has only one eigenvalue.

        sage: C = ComplexField_class(20, eps = 2**(-10))
        sage: ev = C.random_element(); ev
        0.066484 + 0.98568*I
        sage: Ran = matrix(C, 2, [CC.random_element() for i in range(4)])
        sage: mat = ~Ran * matrix(C, [[ev, 1], [0, ev]]) * Ran; mat
        [-0.64527 + 0.42409*I -0.18501 - 0.53142*I]
        [  1.4534 + 0.14619*I   0.77824 + 1.5473*I]
        sage: unique_eigenvalue(mat)
        0.066485 + 0.98568*I

    An example with a matrix which has two eigenvalues.

        sage: C = ComplexField_class(20, 2^(-10))
        sage: ev1, ev2 = C.random_element(), C.random_element(); ev1, ev2
        (-0.42558 + 0.70784*I, -0.38403 - 0.49669*I)
        sage: Ran = matrix(C, 2, [CC.random_element() for i in range(4)])
        sage: mat = ~Ran * matrix(C, [[ev1, 1], [0, ev2]]) * Ran; mat
        [ -0.78315 - 1.0075*I   1.2594 + 0.94418*I]
        [0.042343 - 0.72041*I -0.026468 + 1.2187*I]
        sage: unique_eigenvalue(mat)
        ---------------------------------------------------------------------------
        FloatingPointError                        Traceback (most recent call last)
        <ipython-input-6-ac384bd88e57> in <module>
        ----> 1 unique_eigenvalue(mat)

        ~/sage-9.2/MyPrograms/invariant_subspace/approximate_arithmetic_module.py in unique_eigenvalue(mat)
            473     pol = mat.charpoly()
            474     pol = pol = pol//gcd(pol, pol.derivative())
        --> 475     if pol.degree() > 1 : raise FloatingPointError("The matrix seems have multiple eigenvalues.")
            476     eigval = -pol[0]
            477

        FloatingPointError: The matrix seems have multiple eigenvalues.


    """


    pol = mat.charpoly()
    pol = pol = pol//gcd(pol, pol.derivative())
    if pol.degree() > 1 : raise FloatingPointError("The matrix seems have multiple eigenvalues.")
    eigval = -pol[0]

    return eigval
