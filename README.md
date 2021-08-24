# Package for the factorization of linear differential operators whose coefficients are rational functions.

Distributed under the terms of the GNU General Public License

## Example

```
sage: from ore_algebra import DifferentialOperators
sage: from factorization import factor
sage: Diffops, z, Dz = DifferentialOperators(QQ, 'z')
sage: L = (4*z^2 + 6*z + 2)*Dz^2 + (4*z + 3)*Dz - 1
sage: factor(L)
[(4*z + 4)*Dz + 2, (z + 1/2)*Dz - 1/2]
sage: prod(_) == L
True

```

## Installation

__Requirements:__
- SageMath ([https://www.sagemath.org/](https://www.sagemath.org/))
- ore_algebra ([https://github.com/mkauers/ore_algebra/](https://github.com/mkauers/ore_algebra/))

__How to import this package?__
- download this repository in your pc
- retain the path of the PARENT folder
- execute the following lines in a Sage interpreter (replacing /path/to/parent/folder by the corresponding path)
```
sage: import sys
sage: sys.path.append('/path/to/parent/folder')
```
