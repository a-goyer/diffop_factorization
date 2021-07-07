import collections

from .precision_error import PrecisionError

from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RealField
from sage.rings.qqbar import QQbar
from .complex_optimistic_field import ComplexOptimisticField

from sage.rings.power_series_ring import PowerSeriesRing
from ore_algebra.analytic.differential_operator import PlainDifferentialOperator

from sage.modules.free_module_element import vector

from sage.misc.misc_c import prod
from sage.arith.functions import lcm
from sage.arith.misc import valuation
from .useful_functions import customized_accuracy, power_series_coerce, derivatives
from ore_algebra.analytic import accuracy
from ore_algebra.analytic.monodromy import _monodromy_matrices
from .splitting import invariant_subspace
from .guessing import hp_approximants, guess_exact_numbers


Radii = RealField(30)

MonoData = collections.namedtuple("MonoData", ["precision", "matrices", "points", "loss"])


class LinearDifferentialOperator(PlainDifferentialOperator):

    r"""
    A subclass of differential operators for internal use.
    Assumption: polynomial coefficients and 0 is an ordinary point.
    """

    def __init__(self, dop):

        super(LinearDifferentialOperator, self).__init__(dop)

        self.z = self.base_ring().gen()
        self.Dz = self.parent().gen()
        self.n = self.order()

        self.order_of_truncation = max(100, 2*self.degree() + self.n)
        self.algebraicity_degree = self.base_ring().base_ring().degree()
        self.precision = 100*self.algebraicity_degree

        self.monodromy_data = MonoData(0, [], None, 0)

        self.fuchsian_info = None


    def is_fuchsian(self):

        r"""
        Return True if "self" is fuchian, False otherwise.

        Fuch's criterion: p is a regular point of a_n*Dz^n + ... + a_0 (with a_i
        polynomial) iff no (z-p)^{n-k}*a_k/a_n admits p as pole.
        """

        coeffs = self.coefficients()
        fac = coeffs.pop().factor()
        for (f, m) in fac:
            for k, ak in enumerate(coeffs):
                mk = valuation(ak, f)
                if mk - m < k - self.n: return False

        dop = self.annihilator_of_composition(1/self.z)
        for k, frac in enumerate(dop.monic().coefficients()[:-1]):
            d = (self.z**(self.n - k)*frac).denominator()
            if d(0)==0: return False

        return True


    def monodromy(self, precision, verbose=False):

        r"""
        Compute a generating set of matrices for the monodromy group of "self"
        at 0, such that the (customized) precision of each coefficient is at
        least equal to "precision".
        """

        if verbose: print("Monodromy computation with wanted precision = " + str(precision) + ".")
        if self.monodromy_data.precision<precision:
            success, increment, loss = False, 50, self.monodromy_data.loss
            if self.monodromy_data.points==None:
                useful_singularities = self._singularities(QQbar, include_apparent=False)
            else:
                useful_singularities = self.monodromy_data.points
            while not success:
                try:
                    p = precision + loss + increment
                    if verbose: print("Try with precision = " + str(p) + ".")
                    it = _monodromy_matrices(self, 0, eps=Radii.one()>>p, sing=useful_singularities)
                    points, matrices = [], []
                    for pt, mat, is_scalar in it:
                        if not is_scalar: matrices.append(mat); points.append(pt)
                    output_precision = min(min([customized_accuracy(mat.list()) for mat in matrices], default=p), p)
                    if output_precision<precision:
                        if verbose: print("Insufficient precision, loss = " + str(p - output_precision) + ".")
                        increment = 50 if loss==0 else increment<<1
                    else: success=True
                    loss = max(loss, p - output_precision)
                except (ZeroDivisionError, accuracy.PrecisionError):
                    if verbose: print("Insufficient precision for computing monodromy.")
                    increment = increment<<1
            self.monodromy_data =  MonoData(output_precision, matrices, points, loss)


    def _symbolic_guessing(self):

        """
        Return a non-trivial right factor under the assumtion that the elements
        of the differential Galois group of "self" are homotheties.
        """

        T = self.order_of_truncation
        R = self.base_ring().base_ring()

        while True:

            S = PowerSeriesRing(R, default_prec=T + 1)
            basis = self.local_basis_expansions(QQ.zero(), T + 1) # computing only the first one?
            f = power_series_coerce(basis[0], S)
            pols = hp_approximants([f, f.derivative()], T)
            dop = self.parent()(pols)
            if self%dop==0: return dop
            T = T<<1


    def _guessing(self, v, m):

        """
        Return a non-trivial right factor thanks to an oracle that indicates a
        good linear combination of the solutions of "self" at 0, that is, a
        solution annihilating an operator of smaller order than "self".
        """

        T0, T, d0 = 25, self.order_of_truncation, 0
        p = customized_accuracy(v)
        if p<50: raise PrecisionError("Loosing too much precision to attempt the guessing part.")
        C = ComplexOptimisticField(p, eps=Radii.one()>>p//3)
        v = v.change_ring(C)

        while T0<=2*T:

            S = PowerSeriesRing(C, default_prec=T0 + m)
            basis = self.local_basis_expansions(QQ.zero(), T0 + m) # avoiding the re-computation of the first terms?
            basis = power_series_coerce(basis, S)
            f = v*vector(basis)
            pols = hp_approximants(derivatives(f, m), T0)
            p, d1 = customized_accuracy(pols), max(pol.degree() for pol in pols)
            if d1==d0:
                alg_deg = self.algebraicity_degree
                while 50*alg_deg<p:
                    try:
                        exact_pols = guess_exact_numbers(pols, alg_deg)
                        coeffs = [c for pol in exact_pols for c in pol.coefficients()]
                        selftilde = self.extend_scalars(*coeffs)[0] # not recursive, need an embedding + monodromy problem
                        dop = selftilde.parent()(exact_pols)
                        if selftilde%dop==0:
                            self, self.algebraicity_degree = selftilde, alg_deg
                            self.z, self.Dz = self.base_ring().gen(), self.parent().gen()
                            return dop
                    except PrecisionError: pass
                    alg_deg = alg_deg + 1; print('alg_deg=',alg_deg)
            d0, T0 = d1, T0<<1

        self.order_of_truncation = self.order_of_truncation<<1

        raise PrecisionError("Insufficient precision for the guessing part.")


    def right_Dfactor(self, verbose=False):

        r"""
        Return either a non-trivial right factor of "self" or the string
        'irreducible' if "self" is irreducible.
        """

        if self.precision > 20000: raise NotImplementedError

        if self.n<2: return 'irreducible'
        if verbose: print("Try to factorize an operator of order " + str(self.n) + ".")
        if self.fuchsian_info==None:
            self.fuchsian_info = self.is_fuchsian()
            if not self.fuchsian_info: print("WARNING: The operator is not fuchsian: termination is not guaranteed.")

        self.monodromy(self.precision, verbose=verbose)
        self.precision = self.monodromy_data.precision
        matrices = self.monodromy_data.matrices
        if verbose: print("Monodromy computed with precision = " + str(self.precision) + ".")

        if matrices==[]:
            if verbose: print("Any subspace is invariant --> symbolic guessing.")
            dop = self._symbolic_guessing()
        else:
            try:
                V = invariant_subspace(matrices, verbose=verbose)
                if V is None: return 'irreducible'
                if verbose: print("Find an invariant subspace of dimension " + str(len(V)) + " --> guessing.")
                dop = self._guessing(V[0], len(V))
            except PrecisionError:
                if verbose: print("Insufficient precision.")
                self.precision = self.precision<<1
                return self.right_Dfactor(verbose=verbose)

        return dop


def right_Dfactor(dop, verbose=False):

    r"""
    Return either a non-trivial right factor of "dop" or the string
    'irreducible' if "dop" is irreducible.
    """

    if dop.order()==1: return 'irreducible'
    success, rfactor = try_rational(dop)
    if success: return rfactor

    coeffs, z0, z = dop.monic().coefficients(), QQ.zero(), dop.base_ring().gen()
    while min(c.valuation(z - z0) for c in coeffs)<0: z0 = z0 + QQ.one()
    shifted_dop = dop.annihilator_of_composition(z + z0)

    LDO = LinearDifferentialOperator(shifted_dop)
    result = LDO.right_Dfactor(verbose=verbose)
    if result=='irreducible': return 'irreducible'
    result = result.annihilator_of_composition(z - z0)

    return result


def Dfactor(dop, verbose=False):

    r"""
    Return a list of irreductible operators [L1, L2, ..., Lr] such that L is
    equal to the composition L1.L2...Lr.
    """

    rfactor = right_Dfactor(dop, verbose=verbose)
    if rfactor=='irreducible': return [dop]
    lfactor = rfactor.parent()(dop)//rfactor

    return Dfactor(lfactor, verbose=verbose) + Dfactor(rfactor, verbose=verbose)


def try_rational(dop):
    for (f,) in dop.rational_solutions():
        d = f.gcd(f.derivative())
        rfactor = (f/d)*dop.parent().gen() - f.derivative()/d
        return True, rfactor
    return False, None
