from .precision_error import PrecisionError
from .useful_functions import customized_accuracy
from .invariant_subspace import InvSub
from .complex_optimistic_field import ComplexOptimisticField
from .guessing import hp_approx, guess_rational
from .diffops import is_fuchsian
from ore_algebra.analytic import monodromy_matrices
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.real_mpfr import RealField
from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import vector
from sage.rings.polynomial.polynomial_ring import PolynomialRing_field

Radii = RealField(30)



def right_dfactor(L, prec=100, T0=20, verbose=False, fuchsian=None):

    r"""
    Return a nontrivial right-hand factor of the linear differential operator L
    or None if there is none.
    """

    n = L.order()
    if n<2: return None

    if fuchsian==None:
        fuchsian = is_fuchsian(L)
        if not fuchsian:
            print("WARNING! The operator is not fuchsian: termination is not guaranteed.")

    if verbose: print("Try to factorize an operator of order", n)

    OA = L.parent()
    z, Dz = L.base_ring().gen(), OA.gen()

    rat_sol = L.rational_solutions()
    if len(rat_sol)!=0:
        f = rat_sol[0][0]
        R = f*Dz - f.derivative()
        return R

    coeffs = L.monic().coefficients()
    z0 = 0
    while min(c.valuation(z-z0) for c in coeffs)<0:
        z0 += 1
    K = OA(L.annihilator_of_composition(z + z0))

    if verbose: print("Monodromy computation with prec =", prec)

    p = prec
    success = False
    cpt = 0
    while not success:
        eps = Radii.one() >> p
        try:
            mono = monodromy_matrices(K, 0, eps=eps)
            output_prec = min([customized_accuracy(mat.list()) for mat in mono], default=p)
            if output_prec<prec//2:
                p = prec//2 + p - output_prec + (100<<cpt)
                cpt = cpt + 1
                if verbose: print("Need to increase precision :", p)
            else:
                if verbose: print("Monodromy computed")
                success = True
        except ZeroDivisionError:
            p = 2*p

    if len(mono)==0:
        raise NotImplementedError("Cannot continue. The operator is not Fuschian.")

    try:
        V = InvSub(mono)
    except PrecisionError:
        if verbose: print("Precision not good enough to detect an invariant subspace.")
        return right_dfactor(L, prec=2*prec, T0=T0, verbose=verbose, fuchsian=fuchsian)

    if V is None:
        return None
    d = len(V)

    if verbose: print("Found an invariant subspace of dimension", d)
    try:
        p = min([customized_accuracy(v.list()) for v in V], default=output_prec)
        if p<20: raise PrecisionError("Loosing too much precision to continue.")
        C = ComplexOptimisticField(p, eps = Radii.one()>>p//2)
        if isinstance(K.base_ring(), PolynomialRing_field):
            T = max(coeff.degree() for coeff in K) + K.order() + T0
        else:
            T = max(max(coeff.numerator().degree() for coeff in K), max(coeff.denominator().degree() for coeff in K)) + T0
        if verbose: print("T =", T)
        basis = K.local_basis_expansions(ZZ(0), ZZ(T+d))
        S = PowerSeriesRing(C, default_prec=T+d)
        for k, sol in enumerate(basis):
            sol2 = S.zero()
            for c, mon in sol:
                if c!=0:
                    sol2 += c*S.gen()**mon.n
            basis[k] = sol2
        f = V[0].change_ring(C)*vector(basis)
        df = [f]
        for k in range(d):
            f = f.derivative()
            df.append(f)
        P = [guess_rational(pol) for pol in hp_approx(df, T)]

    except (PrecisionError, ZeroDivisionError):
        if verbose: print("Precision not good enough to guess a candidate right factor.")
        return right_dfactor(L, prec=2*prec, T0=2*T0, verbose=verbose, fuchsian=fuchsian)

    R = OA(P)
    if K%R==0:
        R = OA(R.annihilator_of_composition(z - z0))
        return R
    else:
        if verbose: print("Candidate right factor not good.")
        return right_dfactor(L, prec=2*prec, T0=2*T0, verbose=verbose, fuchsian=fuchsian)



def dfactor(L, verbose=False, fuchsian=None):

    r"""
    Return a list of irreductibles operators [L1, L2, ..., Lr] such that L =
    L1.L2...Lr.
    """

    R = right_dfactor(L, verbose=verbose, fuchsian=fuchsian)
    if R is None:
        return [L]
    else:
        Q = L//R
        return dfactor(Q, verbose=verbose) + dfactor(R, verbose=verbose)
