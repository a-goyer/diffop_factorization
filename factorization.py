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


def _mono_computation(L, prec, loss, verbose):
    k, success = 50, False
    while not success:
        try:
            p = 2*prec + loss + k
            mono = monodromy_matrices(L, 0, Radii.one()>>p)
            output_prec = min([customized_accuracy(mat.list()) for mat in mono], default=p)
            loss = max(loss, p - output_prec)
            if output_prec<2*prec:
                k = k<<1
            else:
                success = True
                prec = output_prec
        except ZeroDivisionError:
            k = k<<1
    return mono, prec, loss



def _T(L):
    if isinstance(L.base_ring(), PolynomialRing_field):
        deg = max(coeff.degree() for coeff in L)
    else:
        deg = max(max(coeff.numerator().degree() for coeff in L), max(coeff.denominator().degree() for coeff in L))
    T = max(100, 2*deg + L.order())
    return T



def _guess(L, V, C, p, T):
    """ return 'Fail' if no PrecisionError is raised and no result is found. """

    current_p = min([customized_accuracy(v.list()) for v in V], default=p)
    if current_p<20: raise PrecisionError("Loosing too much precision to attempt the guessing part.")

    T0, D0, d = 25, 0, len(V)

    while T0<=T:
        basis = L.local_basis_expansions(ZZ(0), ZZ(T0+d))
        S = PowerSeriesRing(C, default_prec=T0+d)
        for k, sol in enumerate(basis):
            sol2 = S.zero()
            for c, mon in sol:
                if c!=0:
                    sol2 += c*S.gen()**mon.n
            basis[k] = sol2
        f = V[0]*vector(basis)
        df = [f]
        for k in range(d):
            f = f.derivative()
            df.append(f)
        P = hp_approx(df, T0)
        D1 = max(pol.degree() for pol in P)
        if D0==D1: # improvement: increase algebraicity degree.
            P = [guess_rational(pol) for pol in P]
            L2 = L.parent()(P)
            if L%L2==0:
                return L2
            else:
                print("WARNING: problem seems to come from degree of algebraicity (not implemented yet).")
                raise PrecisionError("Candidate right factor not good.")
        else:
            D0 = D1
            T0 = 2*T0
    print('cc')
    return 'Fail'



def right_dfactor(L, prec=50, loss=0, T=None, fuchsian=None, verbose=False):

    n = L.order()
    if n<2: return None

    if fuchsian==None:
        if verbose: print("Try to factorize an operator of order " + str(n) + ".")
        fuchsian = is_fuchsian(L)
        if not fuchsian:
            print("WARNING! The operator is not fuchsian: termination is not guaranteed.")

    z, Dz = L.base_ring().gen(), L.parent().gen()

    rat_sol = L.rational_solutions()
    if len(rat_sol)!=0:
        f = rat_sol[0][0]
        L2 = f*Dz - f.derivative()
        return L2

    coeffs, z0 = L.monic().coefficients(), 0
    while min(c.valuation(z - z0) for c in coeffs)<0: z0 += 1
    K = L.annihilator_of_composition(z + z0)

    if verbose: print("Monodromy computation with wanted prec = " + str(2*prec) + ".")
    mono, prec, loss = _mono_computation(K, prec, loss, verbose=verbose)
    if len(mono)==0:
        raise NotImplementedError("Cannot continue. The operator is not Fuschian.")

    if verbose: print("Monodromy computed with prec = " + str(prec) + ".")
    C = ComplexOptimisticField(prec, eps=Radii.one()>>(prec//4))
    #mono = [mat.change_ring(C) for mat in mono]

    try:
        V = InvSub(mono)
    except PrecisionError:
        if verbose: print("Precision is not good enough for invariant subspace computation.")
        return right_dfactor(L, prec, loss, T, fuchsian, verbose)

    if V==None: return None
    if verbose: print("Found an invariant subspace of dimension "+str(len(V))+".")

    V = [v.change_ring(C) for v in V]

    if T==None:
        T = _T(K)
        if verbose: print("Order of truncation = " + str(T) + ".")
    try:
        L2 = _guess(K, V, C, prec//4, T)
    except PrecisionError:
        if verbose: print("Precision not good enough for the guessing part.")
        return right_dfactor(L, prec, loss, T, fuchsian, verbose)

    if L2=='Fail':
        if verbose: print("New maximal order of truncation = " + str(2*T) + ".")
        return right_dfactor(L, prec, loss, 2*T, fuchsian, verbose)
    else:
        L2 = L2.annihilator_of_composition(z - z0)
        return L2


def dfactor(L, verbose=False):

    r"""
    Return a list of irreductible operators [L1, L2, ..., Lr] such that L =
    L1.L2...Lr.
    """

    R = right_dfactor(L, verbose=verbose)
    if R is None:
        return [L]
    else:
        Q = L//R
        return dfactor(Q, verbose=verbose) + dfactor(R, verbose=verbose)
