from sage.all import *
from ore_algebra import OreAlgebra

P, z = PolynomialRing(QQ, 'z').objgen()
OA, Dz = OreAlgebra(P, 'Dz').objgen()


def gens():
    """
    Return (z, Dz).
    """
    return P.gen(), OA.gen()


### Salvy's example (>9h in Maple) ###
salvy_example = {}
salvy_example["left_factor"] = z**2*Dz + 3
salvy_example["right_factor"] = (z - 3)*Dz + 4*z**5
salvy_example["operator"] = salvy_example.get("left_factor") * salvy_example.get("right_factor")
