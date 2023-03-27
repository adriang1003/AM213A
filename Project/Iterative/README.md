# Final Project: Iterative README
`Driver_LinAl.f90`:
Purpose: A driver file that uses the module file `LinAl.f90` for its subroutines.
Run Command: Type 'make' to command window to create executable 'LinAl.exe'. Type './LinAl.exe' to run executable.
**NOTE** -- 'LinAl.exe' takes in 4 inputs from the user:
- 'D = ' integer
- 'Type 1, 2, or 3 to solve Ax = b using a Gauss-Jacobi (1), Gauss-Seidel (2), or Conjugate Gradient algorithm (3):' integer
- 'Do you want to run the a_ii = i case for problem 1? (yes (1) or no (0)):' integer
- 'Do you want to run the a_ii = i case for problem 2? (yes (1) or no (0)):' integer
So, to run 'LinAl.exe', the user needs to supply those inputs

`LinAl.f90`:
Purpose: Stores (relevant) subroutines:
- gaussJacobi: Solves the equation Ax = b using the Gauss-Jacobi algorithm
- gaussSeidel: Solves the equation Ax = b using the Gauss-Seidel algorithm
- CGNoCond: Performs the smart conjugate gradient algorithm without pre-conditioning to solve the problem Ax = b
- CGPreCond: Performs the smart conjugate gradient algorithm with pre-conditioning to solve the problem Ax = b
- Lines 576-836 are the new lines of code for this part of the project

`plotter.py`:
Purpose: Plots the Gauss-Jacobi and Gauss-Seidel algorithms for D = [10,100,1000].
Run Command: Type 'python3 plotter.py' to command window, see 'ConvergenceXXXX.jpg' for the resulting images
**NOTE** -- To run 'plotter.py', you would need to run 'LinAl.exe' for D = [10,100,1000] cases for both the Gauss-Jacobi and Gauss-Seidel algorithms.
