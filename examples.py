# coding: utf-8
r"""
More examples can be found in ore_algebra package.
"""

from sage.rings.rational_field import QQ
from ore_algebra import DifferentialOperators

Diffops, z, Dz = DifferentialOperators(QQ, 'z')



# from http://koutschan.de/data/fcc1/. (fcc4, fcc5, fcc6 are available in
# ore_algebra's examples).
fcc3 = 2*(-1+z)*z**2*(3+z)**2*Dz**3+3*z*(3+z)*(-6+5*z+5*z**2)*Dz**2+6*(-3+3*z+12*z**2+4*z**3)*Dz+6*z*(2+z)


# DEtools[DFactor] (of Maple, diffop package) fails with the following operator
# (example provided by Bruno Salvy). We suspect that the large exponent (=-972)
# at point 3 is the cause. Update: 2020, Dec
DiffopWithLargeExponent = (z**2*Dz + 3)*((z - 3)*Dz + 4*z**5)

# The only right factor of the following operator has a degree k (a parameter)
# while the degree of the operator is 2. For more details, see the article
# "Explicit degree bounds for right factors of linear differential operators" by
# Alin Bostan, Tanguy Rivoal and Bruno Salvy (2020).
DiffopWithLargeRightFactor = lambda k: z*Dz**2 + (2-z)*Dz + k
