from .precision_error import PrecisionError
from .linear_algebra import REF, orbit, generated_algebra
from .linear_algebra import eigenvalues, gen_eigenspaces
from .useful_functions import customized_accuracy
from .complex_optimistic_field import ComplexOptimisticField

try:
    from sage.rings.complex_mpfr import ComplexField
except ModuleNotFoundError: # versions of sage older than 9.3
    from sage.rings.complex_field import ComplexField

from sage.rings.real_mpfr import RealField
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix

from sage.misc.misc_c import prod

class Splitting():

    r"""
    """

    def __init__(self, Mats):

        self.n = Mats[0].nrows()
        self.C = Mats[0].base_ring()
        self.I = identity_matrix(self.C, self.n)

        self.matrices = Mats.copy()
        self.partition = [self.n]
        self.basis = self.I          # column-wise
        self.projections = [self.I]



    def refine(self, mat):

        new_dec, s = [], 0
        for j, nj in enumerate(self.partition):
            new_dec.append(gen_eigenspaces(mat.submatrix(s, s, nj, nj), projections=True))
            s = s + nj

        self.partition = [space['multiplicity'] for bloc in new_dec for space in bloc]
        self.projections = [p*space['projection'](mat)*p for j, p in enumerate(self.projections) for space in new_dec[j]]

        T = matrix()
        for bloc in new_dec:
            basis = []
            for s in bloc:
                basis.extend(s['basis'])
            T = T.block_sum(matrix(basis).transpose())
        self.basis = T*self.basis
        try:
            invT = ~T
        except ZeroDivisionError:
            raise PrecisionError("Cannot compute the transition to the old basis from the new one.")
        self.matrices = [invT*M*T for M in self.matrices]
        self.projections = [invT*p*T for p in self.projections]

        return



    def check_lines(self):

        s = 0
        for j, nj in enumerate(self.partition):
            if nj==1:
                p = self.projections[j]
                p[s,s] = p[s,s] - self.C.one()
                err = max(sum(p[i,j].above_abs() for j in range(self.n)) for i in range(self.n))
                err = self.C.zero().add_error(err)
                vec = self.I[s] + vector([err]*self.n)
                V = orbit(self.matrices, vec)
                if len(V)==0:
                    raise PrecisionError('Projections are not precise enough.')
                if len(V)<self.n:
                    T = self.basis
                    V = [T*v for v in V]
                    return (True, V)
            s = s + nj

        return (False, None)



    def COF_version(self):

        prec1 = self.C.precision()
        prec2 = min(customized_accuracy(M) for M in self.matrices)
        prec2 = min(prec2, customized_accuracy(self.basis))

        if 2*prec2 < prec1: raise PrecisionError("Losing too much precision to continue.")

        COF = ComplexOptimisticField(prec1, eps = RealField(30).one()>>(3*prec1//8))

        Mats = [M.change_ring(COF) for M in self.matrices]
        b = self.basis.change_ring(COF)

        return COF, Mats, b



    def check_nolines(self, verbose=False):

        (COF, Mats, b), p = self.COF_version(), self.partition

        s=0
        for j, nj in enumerate(p):
            if nj>1:
                if verbose: print('Check in a subspace of dimension', nj)
                ind = range(s, s + nj)
                basis =  identity_matrix(COF, self.n)[s:s+nj]
                K = VectorSpace(COF, nj)
                for M in Mats:
                    mat = M.matrix_from_rows_and_columns(ind, ind)
                    K = intersect_eigenvectors(K, mat)
                while K.dimension()>0:
                    if verbose: print('dim K =', K.dimension())
                    vec0 = vector(COF, [0]*s + list(K.basis()[0]) + [0]*(self.n - s - nj))
                    V, T, p = orbit(Mats, vec0, transition=True, pivots=True)
                    if len(V)<self.n:
                        V = [(b*v).change_ring(self.C) for v in V]
                        return (True, V)
                    vec1 = basis[0]
                    if len(REF(matrix([vec0, vec1]), pivots=True)[1])==1: vec1 = basis[1]
                    lc = linear_combination(vec1, V, p)
                    M = sum(cj*T[j] for j, cj in enumerate(lc))
                    mat = M.matrix_from_rows_and_columns(ind, ind)
                    if len(eigenvalues(mat))>1: return ('new_matrix', M)
                    K = intersect_eigenvectors(K, mat)

            s = s + nj

        return (False, None)



def linear_combination(vec, Vecs, p):

    n = len(Vecs)
    p = {value:key for key, value in p.items()}
    lc = [0]*n
    for i in range(n-1, -1, -1):
        x = vec[p[i]]
        vec = vec - x*Vecs[i]
        lc[p[i]] = x

    return lc



def intersect_eigenvectors(K, mat):

    eigvals = eigenvalues(mat)
    if len(eigvals)>1:
        raise PrecisionError('This matrix seems have several eigenvalues.')
    K = K.intersection((mat-eigvals[0]*(mat.parent().one())).right_kernel())

    return K



def invariant_subspace(Mats, *, verbose=False):

    r"""
    Return either a nontrivial subspace invariant under the action of the
    matrices of "Mats" or None if there is none.

    Note: this function is designed for BallField as base ring.

    Note: only the output None is rigorous, in the following sense. If
    invariant_subspace(Mats) is None than for any [M1·, ..., Mr·] in
    Mats=[M1, ..., Mr], there is no nontrivial subpace invariant under the
    action of M1·, ..., Mr·.

    INPUT:

     -- "Mats" -- list of n×n matrices


    OUTPUT:

     -- "V" -- list of vectors of size n or None


    EXAMPLE::

        sage: from diffop_factorization import invariant_subspace
        sage: C = ComplexBallField(30)
        sage: mat1 = matrix(C, [[1, 1], [1, 1]])
        sage: mat2 = matrix(C, [[3, -1], [0 , 2]])
        sage: invariant_subspace([mat1, mat2])
        [([1.00 +/- 1.17e-3] + [+/- 1.17e-3]*I, [1.00 +/- 9.96e-4] + [+/- 9.96e-4]*I)]

    """

    if Mats==[]: raise TypeError("This function requires at least one matrix.")

    n, C = Mats[0].nrows(), Mats[0].base_ring()
    split = Splitting(Mats)
    mat = sum(C(ComplexField().random_element())*M for M in split.matrices)
    split.refine(mat)

    hope = True
    while hope:
        if verbose: print("The partition is currently " + str(split.partition) + ".")
        b, V = split.check_lines()
        if b: return V
        if verbose: print("Lines checked.")

        if len(split.partition)==n: return None
        if verbose: print("Need to check nolines.")

        b, x = split.check_nolines(verbose=verbose)
        if b=='new_matrix': split.refine(x)
        elif b: return x
        else: hope=False

    if verbose: print("Need to compute a basis of the algebra.")
    Mats = generated_algebra(Mats)
    if len(Mats)==n**2:
        return None
    else:
        if verbose: print("Restart with the basis of the algebra.")
        return invariant_subspace(Mats)
