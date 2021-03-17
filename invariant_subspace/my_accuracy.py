from sage.all import *

def my_acc(z):

    r"""
    Return z.accuracy() if z does not contain 0, int(-log(z.rad(), 2) otherwise.

    Examples::
        sage: a = CBF(2).add_error(0.1); a
        [2e+0 +/- 0.101] + [+/- 0.101]*I
        sage: my_acc(a)
        4
        sage: a = a-2; my_acc(a)
        2

    """

    if z.contains_zero():
        acc = int(-log(z.rad(), 2))
    else:
        acc = z.accuracy()

    return acc

a = 2
print("coucou!")
