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

# Fuchsian operators with linear factors but without rational solution.
# The first one is the annihilator of sqrt(1+z) and sqrt(1+2z).
# Thanks to Emre Sertoz for the second one (arosing from actual computations).
readme_ex = (4*z**2 + 6*z + 2)*Dz**2 + (4*z + 3)*Dz - 1
sertoz_ex = (-128*z**2 + 8)*(z*Dz)**3 + (-256*z**2-24)*(z*Dz)**2 + (32*z**2 + 10)*(z*Dz)+ 64*z**2

# DEtools[DFactor] (of Maple, diffop package) fails with the following operator
# (thanks to Bruno Salvy for reporting it). We suspect that the large exponent
# (=-972) at point 3 is the cause. !Not Fuchsian! (Update: 2020, Dec)
salvy_ex = (z**2*Dz + 3)*((z - 3)*Dz + 4*z**5)

# The only right factor of the following operator has a degree k (a parameter)
# while the degree of the full operator is 2. For more details, see the article
# "Explicit degree bounds for right factors of linear differential operators" by
# Alin Bostan, Tanguy Rivoal and Bruno Salvy (2020). !Not Fuchsian!
bostan_ex = lambda k: z*Dz**2 + (2-z)*Dz + k

# This example is from van Hoeij's phd thesis (section 3.1). It seems that its
# only right factor has a degree n^2. Not Fuchsian!
vanhoeij_ex = lambda n: Dz**2 - (1/n)*Dz + n/z
