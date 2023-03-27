
# HW4 README
`Driver_LinAl.f90`:
Purpose: A driver file that uses the module file `LinAl.f90` for its subroutines.
Run Command: Type 'make' to command window to create executable 'LinAl.exe'. Type './LinAl.exe' to run executable.

`LinAl.f90`:
Purpose: Stores (relevant) subroutines:
- ChoDecomp: Performs Cholesky Decomposition
- ChoBacksub: Performs back-substitution after Cholesky Decomposition
- vanderMat: Creates a Vandermonde matrix given a vector X
- normalEq: Forms the normal equation
- HouseQR: Performs a Householder based QR Decomposition
- Lines 251-413 are the new lines of code for this homework
- NOTE: gaussBacksub was used for QR back-substitution since R is an upper triangular matrix and Rx=Q^T*b

`plotter.py`:
Purpose: Plots the fitted curve to the given data points in 'atkinson.dat' for degree 3 and 5 polynomials using Cholesky and QR decomposition.
Run Command: Type 'python3 plotter.py' to command window to run script.
