# tempered_diffusion
This code was created with MATLAB_R2017A in collaboration with [Jim Kelly](https://msu.edu/~kellyja8/) to solve tempered fractional diffusion problems in 1D with reflecting boundary conditions.
Both explicit and implicit Euler methods are implemented for both normalized and centered normalized tempered fractional diffusion equations.
The implementation includes functionality for left- and right-sided tempered Riemann-Liouville derivatives and a convex combination of the two.
To run the code, set the desired parameter values in simulate_tempered.m and run.
Note: tempered_stable.m relies on the function stblpdf.m, written by Mark Veillette and downloaded as part of the [stbl library](https://github.com/markveillette/stbl).
Last update: November 4, 2019 Anna Lischke
