from sage.all import *



# This module provides some functions for linear algebra with ComplexBallField
# (a binding to the Arb library, see https://arblib.org/ for more information)
# as base field.





def ball_echelonize (mat, *, rank=False, take_over=False) :

    """
    Transform "mat" in row echelon form in place.


    The row echelon form is 'almost reduced': the pivots are set to 1 and the
    coefficients below the pivots are set to 0 but the coefficients above the
    pivots may be nonzero.

    If the keyword option "rank=True" is specified, this functions returns
    the rank of "mat".

    If the keyword option "take_over=True" is specified, this function begins
    with checking if the columns has already been reduced. For example, if
    "mat" is already in almost reduced row echelon form, no transformation is
    made (just some checkings). Moreover, it stores and returns the indices of
    the rows which have been modified.

    Some words about correction of this algorithm:
0123456789012345678901234567890123456789012345678901234567890123456789012345678
    The computed rank is the number of pivots used in this algorithm. It cannot
    exceed the rank of a point-matrix in "mat" but it can be strictly smaller.
    If the precision is large enough, "rank" is equal to the minimal rank of
    the point-matrices in "mat".
    Denote by "mat_REF" the matrix computed by this algorithm. Then "mat_REF"
    is in almost reduced row echelon form within the ball-arithmetic meaning.
    Moreover, whatever the point-matrix M in "mat", there are a point-matrix R
    in "mat_REF" and an invertible point-matrix T such that R = T * M. But R is
    not necessarily in almost reduced row echelon form.


    INPUT:

     -- "mat"       -- a matrix
     -- "rank"      -- a boolean (optional, default: 'False')
     -- "take_over" -- a boolean (optional, default: 'False')


    OUTPUT:

    The matrix "mat" is put into row echelon form. Nothing is returned unless:
     - the keyword option "rank=True" is specified, in which case "r" (a non
     negative integer) is returned or
     - the keyword option "take_over=True" is specified, in which case
     "modified_rows" (list of non negative integers) is returned or
     - the both options "rank=True" and "take_over=True" are specified, in
     which case the couple ("r", "modified_rows") is returned.


    EXAMPLES:

    An example with a non singular square matrix.

        sage: R = RealBallField(20)
        sage: mat =  matrix(R, 3, [RR.random_element() for i in range(9)]); mat
        [[-0.351698 +/- 1.30e-8]  [0.702078 +/- 1.33e-7]  [0.268266 +/- 4.78e-7]]
        [ [0.782758 +/- 1.63e-7] [-0.385399 +/- 3.08e-7] [-0.613671 +/- 3.50e-7]]
        [[0.0882850 +/- 1.02e-8] [-0.481625 +/- 9.13e-8]  [0.331743 +/- 2.37e-7]]
        sage: ball_echelonize(mat, take_over=True)
        [0, 1, 2]
        sage: mat
        [                1.00000 [-0.492361 +/- 8.85e-7]  [-0.78398 +/- 5.14e-6]]
        [                      0                 1.00000  [-0.01410 +/- 3.81e-6]]
        [                      0                       0                 1.00000]
        sage: ball_echelonize(mat, take_over=True)
        []
        sage: mat
        [                1.00000 [-0.492361 +/- 8.85e-7]  [-0.78398 +/- 5.14e-6]]
        [                      0                 1.00000  [-0.01410 +/- 3.81e-6]]
        [                      0                       0                 1.00000]


    An example with a singular square matrix.

        sage: R = RealBallField(20)
        sage: mat = random_matrix(QQ, 3, 3, algorithm='echelonizable', rank=2)
        sage: T = matrix(R, 3, [RR.random_element() for i in range(9)])
        sage: mat = ~T * mat.change_ring(R) * T; mat
        [ [-9.67 +/- 6.46e-3]  [51.74 +/- 6.42e-3] [-48.80 +/- 6.04e-3]]
        [  [0.40 +/- 4.17e-3]  [-4.04 +/- 1.45e-3]   [4.39 +/- 5.03e-3]]
        [ [-1.43 +/- 3.03e-3] [10.119 +/- 7.86e-4] [-10.29 +/- 5.41e-3]]
        sage: ball_echelonize(mat, rank=True), mat
        (
           [            1.00000 [-5.35 +/- 3.93e-3]  [5.05 +/- 1.91e-3]]
           [                  0             1.00000 [-1.25 +/- 6.20e-3]]
        2, [                  0                   0        [+/- 0.0288]]
        )


    Two examples with generic rectangular matrices.

        sage: R = RealBallField(10)
        sage: mat =  matrix(R, 3, 4, [RR.random_element() for i in range(12)])
        sage: ball_echelonize(mat); mat
        [                1.00 [-0.107 +/- 3.11e-4]  [0.111 +/- 4.51e-4] [-0.156 +/- 2.51e-4]]
        [                   0                 1.00   [2.4e+1 +/- 0.535]  [-2.0e+1 +/- 0.235]]
        [                   0                    0                 1.00    [-0.3 +/- 0.0510]]
        sage: mat =  matrix(R, 4, 3, [RR.random_element() for i in range(12)])
        sage: ball_echelonize(mat); mat
        [               1.00  [1.65 +/- 2.35e-3] [-0.53 +/- 2.23e-3]]
        [                  0                1.00 [-0.59 +/- 7.14e-3]]
        [                  0                   0                1.00]
        [                  0                   0                   0]


    """


    nrows, ncols = mat.dimensions()
    C = mat.base_ring()

    if take_over: modified_rows = []

    r = 0
    for j in range (ncols) :

        if r < nrows:

            if take_over and mat[r,j].is_one():
                # put zeros under the pivot
                for i in range(r+1, nrows):
                    if not mat[i,j].is_zero():
                        mat[i] = [mat[i,k] - mat[i,j]*mat[r,k] for k in range(ncols)]
                        mat[i,j] = C(0)
                        if not i in modified_rows: modified_rows.append(i)

                r = r + 1

            else:
                # find a good pivot in the j-th column
                pivot_row = max((i for i in range(r, nrows) if mat[i,j].is_nonzero()), key=lambda i: mat[i,j].below_abs(), default=None)

                if not pivot_row is None:
                    # move the selected row to the beginning
                    mat[r], mat[pivot_row] = mat[pivot_row], mat[r]
                    if take_over and not r in modified_rows: modified_rows.append(r)

                    # normalize the i-th row with pivot = 1
                    pivot = mat[r,j]
                    mat[r] = [mat[r,k]/pivot for k in range(ncols)]
                    mat[r,j] = C(1)

                    # put zeros under the pivot
                    for i in range(r+1, nrows) :
                        mat[i] = [mat[i,k] - mat[i,j]*mat[r,k] for k in range(ncols)]
                        mat[i,j] = C(0)
                        if take_over and not i in modified_rows: modified_rows.append(i)

                    r = r + 1

    if rank and take_over:
        return r, modified_rows
    elif rank:
        return r
    elif take_over:
        return modified_rows

    return





def Inv (listM, v) :

    """
    Return the smallest subspace containing "v" and invariant under "listM".


    This function computes a linearly independent sequence of ball-vectors
    "basis" such that the subspace generated by these vectors contains "v" and
    is invariant under the action of the matrices of "listM".

    Some words about correction:

    Whatever the point-vector vp in "v" and whatever the selection of point-
    matrices listMp in "listM" (one for each ball-matrix), there is a choice
    of point-vectors vp1, ..., vpr in "basis" (one for each ball-vector) such
    that Span(vp1, ..., vpr) is contained in the smallest invariant subspace
    containing vp and is invariant under the action of the matrices of listMp.
    In particular, if "basis" has lenght n, then whatever the point- vector
    vp in "v" and whatever the selection of point-matrices listMp in "listM",
    this proves that Inv_listMp(vp) spans the entire space.


    INPUT:

     -- "listM" -- a list of n×n matrices
     -- "v"     -- a vector of size n

    OUTPUT:

     -- "basis" -- a list of vectors of size n


    EXAMPLES:

    An example with just one matrix. ::

        sage: R = RealBallField(20)
        sage: mat = matrix(R, [[1, 1, 0], [0, 1, 1], [0, 0, 1]]); mat
        [1.00000 1.00000       0]
        [      0 1.00000 1.00000]
        [      0       0 1.00000]
        sage: u, v, w = vector(R, [1, 0, 0]), vector(R, [0, 1, 0]), vector(R, [0, 0, 1])
        sage: T = matrix(R, 3, [RR.random_element() for i in range(9)])
        sage: u, v, w, mat = T*u, T*v, T*w, T * mat * ~T
        sage: len(Inv([mat], u)), len(Inv([mat], v)), len(Inv([mat], w))
        (1, 2, 3)
        sage: u1 = Inv([mat], u); u1
        [(1.00000, [-0.59097 +/- 7.29e-6], [1.07235 +/- 5.26e-6])]
        sage: u/u[0]
        ([1.0000 +/- 6.01e-6], [-0.59097 +/- 8.06e-6], [1.07235 +/- 6.66e-6])

    An exemple with two matrices. ::

        sage: R = RealBallField(20)
        sage: mat1 = matrix(R, [[0, 0, 1], [0, 0, 0], [0, 0, 0]])
        sage: mat2 = matrix(R, [[0, 0, 0], [0, 0, 1], [0, 0, 0]])
        sage: u = vector(R, [0, 0, 1])
        sage: print(mat1, '\n'); print(mat2, '\n'); print(u)
        [      0       0 1.00000]
        [      0       0       0]
        [      0       0       0]

        [      0       0       0]
        [      0       0 1.00000]
        [      0       0       0]

        (0, 0, 1.00000)
        sage: T = matrix(R, 3, [RR.random_element() for i in range(9)])
        sage: u, mat1, mat2 =  T*u, T * mat1 * ~T, T * mat2 * ~T
        sage: len(Inv([mat1], u)), len(Inv([mat2], u))
        (2, 2)
        sage: len(Inv([mat1, mat2], u))
        3


    """


    n = len(v)

    basis = [v]
    modified_vectors = [v]

    old_dim, dim = 0, 1
    while old_dim < dim and dim < n :

        # add the products "matrix × vector"
        basis.extend([M*u for M in listM for u in modified_vectors])

        # only keep a list of linearly independent vectors
        mat = matrix(basis)
        r, modifified_rows = ball_echelonize(mat, rank=True, take_over=True)
        basis = list(mat[:r])

        modified_vectors = [basis[i] for i in range(r) if i in modifified_rows]

        old_dim, dim = dim, r

    return basis





def basis_of_algebra (list_of_matrices):

    """
    Return a basis of the algebra generated by "list_of_matrices".


    Computes a sequence of linearly independent ball-matrices "basis" such that
    the algebra generated by the matrices of "list_of_matrices" is spanned by
    the matrices of "basis".

    Some words about correction:
    Whatever the selection of point-matrices Mp1, ..., Mpr in
    "list_of_matrices" (one for each ball-matrix), there is a selection of
    point-matrices Bp1, ..., Bps in "basis" such that Vect(Bp1, ..., Bps) is
    include in Alg(Mp1, ..., Mpr).
    Moreover, any selection of point-matrices in "basis" is a sequence of
    linearly independent matrices. So if this function returns n×n ball
    matrices, it is certified that the algebra generated by  any selection of
    point-matrices Mp1, ..., Mpr in "list_of_matrices" is all the matrices.


    INPUT:

     -- "list_of_matrices" -- a list of n×n matrices

    OUTPUT:

     -- "basis" -- a list of n×n matrices


    EXAMPLES:

        sage: R = RealBallField()
        sage: dim = 5
        sage: mat = matrix(R, dim, [RR.random_element() for i in range(dim**2)])
        sage: mat1, mat2 = mat**2, mat**3
        sage: len(basis_of_algebra([mat1, mat2])) == dim
        True
        sage: mat1 = matrix(R, dim, [RR.random_element() for i in range(dim**2)])
        sage: mat2 = matrix(R, dim, [RR.random_element() for i in range(dim**2)])
        sage: len(basis_of_algebra([mat1, mat2])) == dim**2
        True


    """


    n = list_of_matrices[0].dimensions()[0]
    C = list_of_matrices[0].base_ring()

    basis = list_of_matrices.copy()
    modified_matrices = list_of_matrices.copy()

    old_dim, dim = 0, 1
    while old_dim < dim and dim < n*n :

        # add the products "matrix × matrix"
        basis.extend([M*N for M in basis for N in modified_matrices])

        # keep only a list of linearly independent matrices
        basis = [M.list() for M in basis]
        mat = matrix(basis)
        r, modifified_rows = ball_echelonize(mat, rank=True, take_over=True)
        basis = list(mat[:r])
        basis = [matrix(C, n, n, R) for R in basis]

        modified_matrices = [basis[i] for i in range(r) if i in modifified_rows]

        old_dim, dim = dim, r

    return basis
