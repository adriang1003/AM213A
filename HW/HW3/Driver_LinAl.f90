! File: Driver_LinAL.f90
! Author: Adrian Garcia
! Date: 02/10/2023
! DISCLAIMER: The original template code was given by AM213A Professor Dongwook Lee
! Purpose: Driver file for Homework 3 Part 1: Code.
! Run Command: Type 'make' to command window to create executable 'LinAl.exe'. Type './LinAl.exe' to run executable.
Program Driver_LinAl
  use LinAl
  implicit none
  character(len = 100) :: myFileNameA,myFileNameB,myFileNameC,myFileNameD
  integer :: i,j
  logical :: singular
  real(kind = kind(0.d0)) :: traceA,enormX
  myFileNameA = 'Amat.dat'
  myFileNameB = 'Bmat.dat'
  open(10,file = myFileNameA)
  read(10,*) msize,nsize
  close(10)
  allocate(mat(msize,nsize),A(msize,nsize))
  ! Initialize matrices as 0
  mat = 0.0
  A = 0.0
  call readMat(myFileNameA)
  A = mat
  Am = msize
  An = nsize
  deallocate(mat)
  !*********************Part 1: Problem 1********************
  write(*,*) '*********************Part 1: Problem 1********************'
  ! Print matrix A
  write(*,*) 'Matrix A ='
  ! Function call for printing matrix
  call printMat(A,Am,An)
  write(*,*) '**********************************************************'
  ! Function call for Trace
  call traceMat(A,Am,traceA)
  ! Print Trace
  write(*,*) 'Trace(A) =',traceA
  write(*,*) '**********************************************************'
  do i = 1,An
    ! Function call for Euclidean norm
    call enormVector(A(:,i),Am,enormX)
    ! Print Euclidean norm (for all columns)
    write(*,*) 'Column : ||x||_2 =',i,':',enormX
  end do
  write(*,*) '**********************************************************'
  !*********************Part 1: Problem 2********************
  write(*,*) '*********************Part 1: Problem 2********************'
  open(10,file = myFileNameB)
  read(10,*) msize,nsize
  close(10)
  allocate(mat(msize,nsize),B(msize,nsize),X(msize,nsize),E(msize,nsize))
  ! Initialize matrices as 0
  mat = 0.0
  B = 0.0
  X = 0.0
  E = 0.0
  call readMat(myFileNameB)
  B = mat
  Bm = msize
  Bn = nsize
  Xm = Bm
  Xn = Bn
  Em = Xm
  En = Xn
  deallocate(mat)
  ! Print matrix A and B (before Gaussian Elimination)
  write(*,*) 'Matrix A ='
  call printMat(A,Am,An)
  write(*,*) '**********************************************************'
  write(*,*) 'Matrix B ='
  call printMat(B,Bm,Bn)
  write(*,*) '********Gaussian Elimination with partial pivoting********'
  ! Function call for Gaussian Elimination with partial pivoting
  call gaussElim(A,Am,An,B,Bm,Bn,singular)
  ! Print matrix A (modified: upper triangular) and matrix B (modified)
  write(*,*) 'Matrix A ='
  call printMat(A,Am,An)
  write(*,*) '**********************************************************'
  write(*,*) 'Matrix B ='
  call printMat(B,Bm,Bn)
  write(*,*) '*******Back substitution after Gaussian Elimination*******'
  ! Function call for back substitution after Gaussian Elimination
  call gaussBacksub(A,B,X)
  ! Print solution matrix X
  write(*,*) 'Matrix X ='
  call printMat(X,Xm,Xn)
  ! Calculate error matrix
  E = matmul(A,X) - B
  write(*,*) '**********************************************************'
  ! Print error matrix E
  write(*,*) 'Error Matrix E ='
  call printMat(E,Em,En)
  write(*,*) '**********************************************************'
  do i = 1,En
    ! Function call for Euclidean norm
    call enormVector(E(:,i),Em,enormX)
    ! Print Euclidean norm (for all columns)
    write(*,*) 'Column : ||x||_2 =',i,':',enormX
  end do
  write(*,*) '**********************************************************'
  deallocate(A,B,X,E)
  !*********************Part 1: Problem 3********************
  write(*,*) '*********************Part 1: Problem 3********************'
  open(10,file = myFileNameA)
  read(10,*) msize,nsize
  close(10)
  allocate(mat(msize,nsize),A(msize,nsize),ALU(msize,nsize),L(msize,nsize),U(msize,nsize),s(msize))
  ! Initialize matrices as 0
  mat = 0.0
  A = 0.0
  ALU = 0.0
  L = 0.0
  U = 0.0
  s = 0
  call readMat(myFileNameA)
  A = mat
  ALU = mat
  Am = msize
  An = nsize
  deallocate(mat)
  open(10,file = myFileNameB)
  read(10,*) msize,nsize
  close(10)
  allocate(mat(msize,nsize),B(msize,nsize),X(msize,nsize),E(msize,nsize))
  ! Initialize matrices as 0
  mat = 0.0
  B = 0.0
  X = 0.0
  E = 0.0
  call readMat(myFileNameB)
  B = mat
  Bm = msize
  Bn = nsize
  Xm = Bm
  Xn = Bn
  Em = Xm
  En = Xn
  deallocate(mat)
  ! Print matrix A (before LU Decomposition)
  write(*,*) 'Matrix A ='
  call printMat(ALU,Am,An)
  write(*,*) '**********LU Decomposition with partial pivoting**********'
  ! Function call for LU Decomposition with partial pivoting
  call LUDecomp(ALU,Am,s,singular)
  ! Print matrix A (modified: LU Decomposed)
  write(*,*) 'Matrix A ='
  call printMat(ALU,Am,An)
  ! Find L and U
  do j = 1,Am
    do i = 1,Am
      if (i .eq. j) then
        U(i,j) = ALU(i,j)
        L(i,j) = 1
      else if (i < j) then
        U(i,j) = ALU(i,j)
      else
        L(i,j) = ALU(i,j)
      endif
    enddo
  enddo
  write(*,*) '**********************************************************'
  write(*,*) "Matrix L ="
  call printMat(L, Am, An)
  write(*,*) '**********************************************************'
  write(*,*) "Matrix U ="
  call printMat(U, Am, An) 
  write(*,*) '*********Back substitution after LU Decomposition*********'
  ! Function call for back substitution after LU Decomposition
  call LUBacksub(ALU,Am,B,s,X,singular)
  ! Print solution matrix X
  write(*,*) 'Matrix X ='
  call printMat(X,Xm,Xn)
  ! Calculate error matrix
  enormX = 0.0
  E = matmul(A,X) - B
  write(*,*) '**********************************************************'
  ! Print error matrix E
  write(*,*) 'Error Matrix E ='
  call printMat(E,Em,En)
  write(*,*) '**********************************************************'
  do i = 1,En
    ! Function call for Euclidean norm
    call enormVector(E(:,i),Em,enormX)
    ! Print Euclidean norm (for all columns)
    write(*,*) 'Column : ||x||_2 =',i,':',enormX
  end do
  write(*,*) '**********************************************************'
  deallocate(A,ALU,B,X,E)
  !*********************Part 1: Problem 3********************
  write(*,*) '*********************Part 1: Problem 3********************'
  myFileNameC = 'Cmat.dat'
  myFileNameD = 'Dmat.dat'
  open(10,file = myFileNameC)
  read(10,*) msize,nsize
  close(10)
  allocate(mat(msize,nsize),A(msize,nsize))
  ! Initialize matrices as 0
  mat = 0.0
  A = 0.0
  call readMat(myFileNameC)
  A = mat
  Am = msize
  An = nsize
  deallocate(mat)
  open(10,file = myFileNameD)
  read(10,*) msize,nsize
  close(10)
  allocate(mat(msize,nsize),B(msize,nsize),X(msize,nsize))
  ! Initialize matrices as 0
  mat = 0.0
  B = 0.0
  X = 0.0
  call readMat(myFileNameD)
  B = mat
  Bm = msize
  Bn = nsize
  Xm = Bm
  Xn = Bn
  deallocate(mat)
  ! Print matrix A and B (before Gaussian Elimination)
  write(*,*) 'Matrix A ='
  call printMat(A,Am,An)
  write(*,*) '**********************************************************'
  write(*,*) 'Matrix B ='
  call printMat(B,Bm,Bn)
  write(*,*) '********Gaussian Elimination with partial pivoting********'
  ! Function call for Gaussian Elimination with partial pivoting
  call gaussElim(A,Am,An,B,Bm,Bn,singular)
  ! Print matrix A (modified: upper triangular) and matrix B (modified)
  write(*,*) 'Matrix A ='
  call printMat(A,Am,An)
  write(*,*) '**********************************************************'
  write(*,*) 'Matrix B ='
  call printMat(B,Bm,Bn)
  write(*,*) '*******Back substitution after Gaussian Elimination*******'
  ! Function call for back substitution after Gaussian Elimination
  call gaussBacksub(A,B,X)
  ! Print solution matrix X
  write(*,*) 'Matrix X ='
  call printMat(X,Xm,Xn)
  open(10,file = 'myplot.dat')
  do i = 1,size(X,dim=1)
    write(10,*) X(i,:)
  enddo
  close(10)
  deallocate(A,B,X)
End Program Driver_LinAl