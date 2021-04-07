<p align="center">  <img src="https://histcultcine.hypotheses.org/files/2020/01/work-in-progress-wip.jpg" width="130" height="130" /> </p align="center">

# Package for the factorization of linear differential operators whose coefficients are polynomials.

## Examples (automation to come)

#### Find a factorization

```                                                                                               
sage: from ore_algebra import *
sage: from ore_algebra.analytic import monodromy_matrices
sage: from diffop_factorization import InvSub
sage: from diffop_factorization.guessing import guess_rational, hp_approx
sage: Diffops, z, Dz = DifferentialOperators(QQ, var='z')
sage: L1, L2 = (z + 1)*Dz + 1, Dz - 3
sage: L = L1*L2; L
(z + 1)*Dz^2 + (-3*z - 2)*Dz - 3
sage: mono = monodromy_matrices(L, 0)
sage: V = InvSub(mono); V
[(1.0000000000000000, [3.000000000000000 +/- 6.92e-16] + [+/- 1.92e-16]*I)]
sage: v = guess_rational(V[0], p=40); v
(1, 3)
sage: S = PowerSeriesRing(QQ, 'z', 15)
sage: basis = [S(str(x)) for x in L.local_basis_expansions(0, 15)]
sage: f = v*vector(basis)
sage: P = hp_approx([f, f.derivative()], 10); P
(-3, 1)
sage: K2 = Diffops(P.list()); K1 = L//K2; K1, K2
((z + 1)*Dz + 1, Dz - 3)
```

#### Prove the irreducibility (the following example takes 1 minute in Maple)
```
sage: from ore_algebra import *                                                    
sage: from ore_algebra.examples.fcc import dop5                                          
sage: from ore_algebra.analytic import monodromy_matrices                          
sage: from diffop_factorization import InvSub                     
sage: dop5.singularities()                                                         
{-5, 0, 1, 5}
sage: %time mono = monodromy_matrices(dop5, -1, eps=1e-200)                      
CPU times: user 6.4 s, sys: 0 ns, total: 6.4 s
Wall time: 6.4 s
sage: InvSub(mono) is None                                                         
True
```

## Help to import this package

- clone diffop_factorization in your machine
- retain the path of the PARENT folder that contains this package
- execute the following lines replacing /path/to/parent/folder by the corresponding path
```
sage: import sys
sage: sys.path.append('/path/to/parent/folder')
```

_Example:_
- put diffop_factorization in your "Downloads" folder
- execute the following lines replacing 'yourname' by your name
```
sage: import sys
sage: sys.path.append('/home/yourname/Downloads')
```
