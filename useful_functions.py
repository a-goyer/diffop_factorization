from sage.all import *


def customized_accuracy(x):

    """
    Return either the absolute accuracy of x if x contains 0 or the relative
    accuracy of x if x does not contains 0.

    Note that works also if x is a list (minimal accuracy of the coefficients).

    INPUT:
     - 'x' - a complex ball or a list of complex balls

    OUTPUT:
     - 'acc' - a nonnegative integer
    """

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
