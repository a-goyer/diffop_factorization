from sage.all import *
from sage.rings.complex_mpfr import ComplexField_class
from .ball_arithmetic import *
from .approximate_arithmetic import *




# Here is the algorithm for computing an invariant subspace.
#
# Let [M1, ..., Mr] be a list of approximate n×n matrices with rigorous bounds.
# If the bounds are small enough, the 'invariant_subspace' function:
#  - either computes an approximate nontrivial subspace (without error bounds)
#  which is invariant (at the working precision) under the action of these
#  matrices,
#  - or certifies that, whatever the exact matrices [eM1, ..., eMr]
#  represented by [M1, ..., Mr], there is no nontrivial subspace invariant
#  under eM1, ..., eMr.
# Thus only the second case is rigorous. If the precision is not large enough,
# this function raises a PrecisionError.
#
# We use a simpler version of van der Hoeven's algorithm (see the pseudo-code
# in simple_version_for_invariant_subspace.pdf).



class Splitting ():

    r"""
    This class handles the current splitting.


    Let M = {M1, ..., Mr} be a list of n×n matrices with complex coefficients.
    An M-splitting is a decomposition of C^n as a direct sum E1 + ... + Ek such
    that each projection pj onto Ej along the others is polynomial in M.

    The attributes of a Splitting object are:
     - the list of matrices, viewed with respect to a basis B adapted to the
     current splitting,
     - the partition of the dimension (= the list of the dimensions of the Ej),
     - the basis B with respect to the canonical basis,
     - the projection matrices with respect to the basis B, computed from the
     input matrices.

    (In an exact world, the projections are just the matrices with only zero
    coefficients except for some diagonal coefficients which are equal to 1.)

    """



    def __init__ (self, list_of_matrices):

        mat = list_of_matrices[0] # what if [] ?
        self._n = mat.nrows()
        self._C = mat.base_ring()
        self._prec = self._C.precision()

        self._field_approx = ComplexField_class(self._prec, eps = 2**(-self._prec//2))

        self._matrices = list_of_matrices.copy()
        self._partition = [self._n]
        self._basis = list(identity_matrix(self._C, self._n))
        self._projections = [identity_matrix(self._C, self._n)]





    def refine (self, mat):

        """
        Replace the ambient splitting by a finer splitting adapted to "mat".


        An M-splitting C^n = F1 + .. + Fl is finer than C^n = E1 + ... + Ek if
        each Fi is equal to a sum of some Ej.

        An M-splitting F1 + .. + Fl is adpated to "mat" if the endomorphism
        induced by "mat" on each Fi have a single eigenvalue.

        If all the induced endomorphisms have only one eigenvalue, this method
        cannot strictly refine the splitting and 'False' is returned.

        This method uses the generalized eigenspaces of the induced
        endomorphisms to split each Ej.

        Assumption: "mat" is polynomial in self._matrices.


        INPUT:

         -- 'mat' -- a matrix

        OUTPUT:

         -- 'strict_refining' -- a boolean


        EXAMPLES:
        """

        Capprox = self._field_approx
        Pol, x = PolynomialRing(self._C, 'x').objgen()

        new_partition = []
        new_basis = [] # written with respect to self._basis
        new_projections = []


        # decompose each induced endomorphism
        strict_refining=False
        s = 0
        for j, mj in enumerate(self._partition) :

            matj = mat.submatrix(s, s, mj, mj)
            matj = matj.change_ring(Capprox)
            dec = gen_eigspaces_with_proj(matj)

            if len(dec) > 1 :
                strict_refining=True

                new_partition.extend([space[1] for space in dec])

                basisj = []
                for space in dec : basisj.extend(space[2])
                basisj = [list(v.change_ring(self._C)) for v in basisj]
                new_basis.extend([[self._C(0)]*s + v + [self._C(0)]*(self._n-s-mj) for v in basisj])

                polj = [space[3].change_ring(self._C) for space in dec]
                projj = self._projections[j]
                new_projections.extend([projj*p(mat)*projj for p in polj])

            else :
                new_partition.append(mj)
                new_basis.extend(list(identity_matrix(self._C, self._n))[s:(s+mj)])
                new_projections.append(self._projections[j])

            s = s + mj

        # update self replacing the old splitting by the new splitting
        if strict_refining:
            self._partition = new_partition
            self._projections = new_projections
            T = matrix(self._C, new_basis)
            self._basis = list(T * matrix(self._C, self._basis))

            # write list_of_matrices and projections with respect to the new basis
            T = T.transpose()
            for i in range(len(self._matrices)) :
                M = self._matrices[i]
                M = ~T * M * T
                self._matrices[i] = M
            for i in range(len(self._projections)) :
                proj = self._projections[i]
                proj = ~T * proj * T
                self._projections[i] = proj


        return strict_refining


    def invariant_subspace(self):

        """

        EXAMPLES:
        """


        Capprox = self._field_approx
        In = identity_matrix(self._C, self._n)

        ### Computation of a maximal splitting ###

        mat = sum(self._C(Capprox.random_element()) * M for M in self._matrices)
        self.refine(mat)



        ### Check the 1-dimensional subspaces ###
        s = 0
        for m in self._partition :
            if m == 1 :

                # compute the error bound
                mat = self._projections[s]
                mat[s,s] = mat[s,s] - 1
                epsilon = max(sum(mat[i,j].above_abs() for j in range(self._n)) for i in range(self._n)) # upper bound for the operator norm associated to Infinity-norm
                eps_ball = self._C(0).add_error(epsilon)
                v = In[s] + vector([eps_ball]*self._n)
                if all(coeff.contains_zero() for coeff in v):
                    raise FloatingPointError("The projections are not accurate enough.")

                Inv_v = Inv(self._matrices, v)
                if len(Inv_v) < self._n :
                    T = matrix(Capprox, self._basis).transpose()
                    Inv_v = [T * u.change_ring(Capprox) for u in Inv_v]
                    return Inv_v

            s = s + m

        if len(self._partition) == self._n : return None


        ### Find a 'good vector' (i.e. s.t. Inv(v) != C^n) in a d-dimensional subspace, d>1 ###

        approx_list_of_matrices = [mat.change_ring(Capprox) for mat in self._matrices]

        # compute an index j for which the j-th element of the partition is minimal among those > 1
        mj = max(self._partition)
        j = self._partition.index(mj)
        for i, mi in enumerate(self._partition) :
            if mi > 1 and mi < mj : j, mj = i, mi

        s = sum(self._partition[i] for i in range(j))
        I = range(s, s + mj)

        # compute K the space of commun eigenvectors
        K = VectorSpace(Capprox, mj)
        for mat in approx_list_of_matrices :
            matj = mat.matrix_from_rows_and_columns(I, I)
            lambdaj = unique_eigenvalue(matj)
            K = K.intersection((matj - lambdaj*identity_matrix(Capprox, mj)).right_kernel())

        while K.dimension() > 0:
            vj = K.random_element()
            v = vector(Capprox, [0]*s + list(vj) + [0]*(self._n-s-mj))
            Inv_v, T = Inv_with_T(approx_list_of_matrices, v)

            if len(Inv_v) < self._n :
                T = matrix(Capprox, self._basis).transpose()
                Inv_v = [T * u for u in Inv_v]
                return Inv_v

            matj = T[s].matrix_from_rows_and_columns(I, I)
            lambdaj = unique_eigenvalue(matj)
            K = K.intersection((matj - lambdaj*identity_matrix(Capprox, mj)).right_kernel())

        self._matrices = basis_of_algebra(self._matrices)
        if len(self._matrices) < self._n**2 :
            return self.invariant_subspace()
        else : return None





def invariant_subspace(L):

    """


    """

    split = Splitting(L)
    V = split.invariant_subspace()

    return V
