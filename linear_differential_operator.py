from sage.arith.all import lcm
from sage.rings.qqbar import QQbar
from sage.rings.rational_field import QQ
from ore_algebra.ore_operator_1_1 import UnivariateDifferentialOperatorOverUnivariateRing
from ore_algebra.analytic.differential_operator import DifferentialOperator

class LinDiffOp(UnivariateDifferentialOperatorOverUnivariateRing):
    r"""
    A subclass of differential operators for internal use.
    """

    def __init__(self, dop):
        if not dop:
            raise ValueError("operator must be nonzero")
        if not dop.parent().is_D():
            raise ValueError("expected an operator in K(x)[D]")
        _, _, _, dop = dop.numerator()._normalize_base_ring()
        den = lcm(c.denominator() for c in dop)
        dop *= den
        super(LinDiffOp, self).__init__(
                dop.parent(), dop)

        self.singularities = DifferentialOperator(dop)._singularities(QQbar)
        z = QQ.zero()
        while z in self.singularities: z = z + QQ.one()
        self.base_point = z
