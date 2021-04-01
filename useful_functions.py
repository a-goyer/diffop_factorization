from sage.all import *
from sage.matrix.matrix_complex_ball_dense import Matrix_complex_ball_dense
from sage.modules.free_module_element import FreeModuleElement_generic_dense


def overlaps(a, b):
    """
    Function overlaps designed for vectors and matrices.
    """
    l = len(a.list())
    return all(a.list()[i].overlaps(b.list()[i]) for i in range(l))

def customized_accuracy(x):

    """
    Return either the absolute accuracy of x if x contains 0 or the relative
    accuracy of x if x does not contains 0.

    Note that works also if x is a list or a vector or a matrix (minimal
    accuracy of the coefficients).

    INPUT:
     - 'x' - a complex ball

    OUTPUT:
     - 'acc' - a nonnegative integer
    """

    if isinstance(x, Matrix_complex_ball_dense) or isinstance(x, FreeModuleElement_generic_dense):
        x = x.list()
        acc = min(customized_accuracy(c) for c in x)
        return acc

    if isinstance(x, list):
        acc = min(customized_accuracy(c) for c in x)
        return acc

    if x.contains_zero() :
        acc = -log(x.rad(), 2)
        if not acc == +infinity:
            acc = int(acc)
    else:
        acc = x.accuracy()

    return acc
