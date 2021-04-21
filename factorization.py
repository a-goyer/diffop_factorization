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

Radii = RealField(30)



def right_factor(L, prec=100, T0=5, verbose=False):

    r"""
    Return a nontrivial right-hand factor of the linear differential operator L
    or None if there is none.
    """

    if not is_fuchsian(L):
        print("WARNING: Operator not fuchsian! Termination is not guaranteed.")

    if L.order()==1: return None

    if verbose: print("Operator =", L)

    OA = L.parent()
    z, Dz = L.base_ring().gen(), OA.gen()

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
            if 4*min(customized_accuracy(mat.list()) for mat in mono)<3*prec:
                p = 2*p
                if verbose: print("Need to increase precision :", p)
            else:
                success = True
        except ZeroDivisionError:
            p = 2*p

    try:
        V = InvSub(mono)
    except PrecisionError:
        return right_factor(L, prec=2*prec, verbose=verbose)

    if V is None:
        return None
    d = len(V)

    if verbose: print("Found an invariant subspace of dimension", d)

    C = ComplexOptimisticField(prec, eps = Radii.one() >> prec//2)
    T = max(pol.degree() for pol in K) + K.order() + T0
    if verbose: print("T =", T)
    S = PowerSeriesRing(C, 'z', T+d)
    basis = [S(str(v)) for v in K.local_basis_expansions(ZZ(0), ZZ(T+d))]
    f = V[0].change_ring(C)*vector(basis)
    df = [f]
    for k in range(d):
        f = f.derivative()
        df.append(f)
    try:
        P = [guess_rational(pol) for pol in hp_approx(df, T)]
    except PrecisionError:
        return right_factor(L, prec=2*prec, T0=2*T0, verbose=verbose)

    R = OA(P)
    if K%R==0:
        R = OA(R.annihilator_of_composition(z - z0))
        return R
    else:
        return right_factor(L, prec=2*prec, T0=2*T0, verbose=verbose)


def factor(L, verbose=False):

    R = right_factor(L, verbose=verbose)
    if R is None:
        return [L]
    else:
        Q = L//R
        return factor(Q, verbose=verbose) + factor(R, verbose=verbose)


#Free.<z,Dz> = FreeAlgebra(QQ)
#Dop.<z,Dz> = Free.g_algebra(relations={Dz*z: z*Dz+1})
#dop = (z^2*Dz+3)*((z-3)*Dz+4*z^5)
