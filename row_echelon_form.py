from .precision_error import PrecisionError


def REF(mat, transformation=False, struct=False, prec_struct={}):

    r"""
    Return a row echelon form of "mat".


    The computed row echelon form is 'almost reduced': the pivots are set to 1
    and the coefficients below the pivots are set to 0 but the coefficients
    above the pivots may be nonzero.

    Some words about the correction of this function:
    Let (R, T, s) be the output for REF(mat).
    1. For any matp in mat, there are Rp in R, Tp in T such that Rp = Tp*matp.
    2. For any Tp in T, Tp is invertible.
    3. The rank of any matp i mat is greater than or equal to the lenght of s.
    4. For each key j of s:
        4.a. R[i,j] = 1 (exactly) for i = s[j]['pivot'],
        4.b. R[i,j] = 0 (exactly) for each i in s[j]['zeros'].


    INPUT:

     -- "mat"            -- an n×p matrix of balls
     -- "transformation" -- a boolean (optional, default: False)
     -- "struct"         -- a boolean (optional, default: False)
     -- "prec_struct"    -- a dictionary (optional, default: {})


    OUTPUT:

     -- "R" -- an n×p matrix of balls
     -- "T" -- an n×n matrix of balls (if transformation=True is specified)
     -- "s" -- a dictionary (if struct=True is specified, more details below)

    The keys of the dictionary "s" are integers between 0 and p-1. For each such
    key j, s[j] is a dictionary with two keys 'pivot' and 'zeros': s[j]['pivot']
    is an integer between 0 and n-1 and s[j]['zeros'] is a list of integers
    between 0 and n-1.


    EXAMPLES::
    """

    n, p, C = mat.nrows(), mat.ncols(), mat.base_ring()
    T = identity_matrix(C, n)
    s = prec_struct

    for j in range(p):
        

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



    if transformation:
        if det(T).contains_zero():
            raise PrecisionError("Cannot compute an invertible matrix.")
        if struct:
            return R, T, s
        else:
            return R, T
    if struct:
        return R, s

    return R
