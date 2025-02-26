# HW5 README
`Driver_LinAl.f90`:
Purpose: A driver file that uses the module file `LinAl.f90` for its subroutines.
Run Command: Type 'make' to command window to create executable 'LinAl.exe'. Type './LinAl.exe' to run executable.

`LinAl.f90`:
Purpose: Stores (relevant) subroutines:
- hessenDecomp: Performs Hessenberg Decomposition
- QRNoShift: Performs QR Decomposition without shifts to calculate the eigenvectors/values of matrix A
- QRShift: Performs QR Decomposition with shifts to calculate the eigenvectors/values of matrix A
- inverseIter: Performs an Inverse Iteration algorithm to calculate an eigenvector of matrix A
- Lines 435-576 are the new lines of code for this homework
