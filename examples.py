# coding: utf-8
r"""
This file contains some examples: the one raised by Salvy, ...
"""

from sage.rings.rational_field import QQ
from ore_algebra import DifferentialOperators

diffops, z, Dz = DifferentialOperators(QQ, 'z')



# DEtools[DFactor] (Maple, diffop package) fails with the following operator
# (in 2020, Dec).
salvy = (z**2*Dz + 3) * ((z - 3)*Dz + 4*z**5)
