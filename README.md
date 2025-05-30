# HeisenbergXX
### About

Repository containing optimised codes to study magnetic quantum systems. Code written to perform calculations for Physics Bachelor's Thesis.

Built with Fortran 90 and Python 3.13.

N.B.: Requires Python installation, access to a Python package manager, Fortran compiler (`gfortran` is assumed) as well as the `BLAS` & `LAPACK` numerical libraries.

### Structure

This repository contains folders with data files &ndash; namely `correl_data/`, `spin_chain_data/`  &ndash; which record the output of numerical calculations. The problem parameters are set in `spin_chain.dat`.

All `.f90` files do the heavy number-crunching, while `utils.py` defines some convenient Python wrappers. The `figures.ipynb` notebook produces many figures, saved in the `figures/` folder, most of which are included in the actual thesis. The `examples.ipynb` notebook explains how to use this repository's code.

```
HeisenbergXX
.
|-- src/
|  |-- correl_data/
|  |  |-- correl.dat
|  |  |-- eigvals.dat
|  |-- spin_chain_data/
|  |  |-- eigenvectors.dat
|  |  |-- poly_data.dat
|  |  |-- polynomials.dat
|  |-- correl.f90
|  |-- spin_chain.dat
|  |-- spin_chain.f90
|  |-- utils.py
|-- figures/
|-- .gitignore
|-- examples.ipynb
|-- figures.ipynb
|-- README.md
```