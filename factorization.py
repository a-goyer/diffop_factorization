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



def right_factor(L, prec=100, T0=5, verbose=False):

    r"""
    Return a nontrivial right-hand factor of the linear differential operator L
    or None if there is none.
    """

    n = L.order()
    if n<2: return None

    fuchsian = is_fuchsian(L)
    if not fuchsian:
        print("WARNING: Operator not fuchsian! Termination is not guaranteed.")

    if verbose: print("Try to factorize an operator of order", n)

    OA = L.parent()
    z, Dz = L.base_ring().gen(), OA.gen()

    rat_sol = L.rational_solutions()
    if len(rat_sol)!=0:
        f = rat_sol[0][0]
        R = f*Dz - f.derivative()
        return R

    z0 = 0
    while z0 in L.singularities():
        z0 += 1
    K = OA(L.annihilator_of_composition(z + z0))

    if verbose: print("Monodromy computation with prec =", prec)

    p = prec
    success = False
    while not success:
        eps = Radii.one() >> p
        try:
            mono = monodromy_matrices(K, 0, eps=eps)
            output_prec = min([customized_accuracy(mat.list()) for mat in mono], default=prec)
            if output_prec<0:
                p = p - output_prec + 10
                if verbose: print("Need to increase precision :", p)
            else:
                if verbose: print("Monodromy computed")
                success = True
        except ZeroDivisionError:
            p = 2*p

    if len(mono)==0:
        raise NotImplementedError("The operator is not Fuschian.")

    try:
        V = InvSub(mono)
    except PrecisionError:
        if verbose: print("Precision not good enough to detect an invariant subspace.")
        return right_factor(L, prec=2*prec, verbose=verbose)

    if V is None:
        return None
    d = len(V)

    if verbose: print("Found an invariant subspace of dimension", d)

    try:
        C = ComplexOptimisticField(prec, eps = Radii.one() >> prec//2)
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
        return right_factor(L, prec=2*prec, T0=2*T0, verbose=verbose)

    R = OA(P)
    if K%R==0:
        R = OA(R.annihilator_of_composition(z - z0))
        return R
    else:
        return right_factor(L, prec=2*prec, T0=2*T0, verbose=verbose)



def factors(L, verbose=False):

    R = right_factor(L, verbose=verbose)
    if R is None:
        return [L]
    else:
        Q = L//R
        return factors(Q, verbose=verbose) + factors(R, verbose=verbose)



#Free.<z,Dz> = FreeAlgebra(QQ)
#Dop.<z,Dz> = Free.g_algebra(relations={Dz*z: z*Dz+1})
#dop = (z^2*Dz+3)*((z-3)*Dz+4*z^5)
