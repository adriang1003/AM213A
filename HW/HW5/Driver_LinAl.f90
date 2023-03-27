! File: Driver_LinAL.f90
! Author: Adrian Garcia
! Date: 03/10/2023
! Purpose: Driver file for Homework 5 Part 1: Code
Program Driver_LinAl
  use LinAl
  implicit none
  character(len = 100) :: myFileName
  integer :: i
  real(kind = kind(0.d0)) :: j
  real(kind = kind(0.d0)),dimension(4) :: mu
  !*********************Part 1: Problem 1********************
  myFileName = 'Amat.dat'
  open(10,file = myFileName)
  read(10,*) msize,nsize
  close(10)
  allocate(mat(msize,nsize),A(msize,nsize))
  ! Initialize matrices as 0
  mat = 0.0
  A = 0.0
  call readMat(myFileName)
  A = mat
  write(*,*) '*********************Part 1: Problem 1********************'
  ! Print matrix A
  write(*,*) 'Matrix A ='
  ! Function call for printing matrix
  call printMat(A,msize,nsize)
  write(*,*) '***************After Hessenberg Decomposition*************'
  ! Function call for Hessenberg Decomposition
  call hessenDecomp(A,msize,nsize)
  ! Print matrix A (modified: Hessenberg form)
  write(*,*) 'Matrix A ='
  call printMat(A,msize,nsize)
  write(*,*) '**********************************************************'
  deallocate(mat,A)
  !*********************Part 1: Problem 2********************
  myFileName = 'Bmat.dat'
  open(10,file = myFileName)
  read(10,*) msize,nsize
  close(10)
  allocate(mat(msize,nsize),A(msize,nsize),V(msize,nsize))
  ! Initialize matrices as 0
  mat = 0.0
  A = 0.0
  V = 0.0
  call readMat(myFileName)
  A = mat
  Am = msize
  An = nsize
  Vm = Am
  Vn = An
  write(*,*) '*********************Part 1: Problem 2********************'
  ! Print matrix A
  write(*,*) 'Matrix A ='
  ! Function call for printing matrix
  call printMat(A,Am,An)
  write(*,*) '****************After QR algorithm w/o shift**************'
  ! Function call for QR algorithm w/o shift
  call QRNoShift(A,Am,An,V)
  ! Print eigenvectors (matrix V)
  write(*,*) 'eigenvectors ='
  call printMat(V,Vm,Vn)
  ! Print eigenvalues (diagonal entries of matrix A)
  write(*,*) 'eigenvalues ='
  do i = 1,An
    write(*,*) A(i,i)
  enddo
  write(*,*) '*****************After QR algorithm w/ shift**************'
  A = mat
  V = 0.0
  ! Function call for QR algorithm w/ shift
  call QRShift(A,Am,An,V)
  ! Print eigenvectors (matrix V)
  write(*,*) 'eigenvectors ='
  call printMat(V,Vm,Vn)
  ! Print eigenvalues (diagonal entries of matrix A)
  write(*,*) 'eigenvalues ='
  do i = 1,An
    write(*,*) A(i,i)
  enddo
  write(*,*) '**********************************************************'
  deallocate(mat,A,V)
  !*********************Part 1: Problem 3********************
  myFileName = 'Cmat.dat'
  open(10,file = myFileName)
  read(10,*) msize,nsize
  close(10)
  allocate(mat(msize,nsize),A(msize,nsize),xVec(msize,1))
  ! Initialize matrices as 0
  mat = 0.0
  A = 0.0
  xVec = 0.0
  call readMat(myFileName)
  A = mat
  write(*,*) '*********************Part 1: Problem 3********************'
  ! Print matrix A
  write(*,*) 'Matrix A ='
  ! Function call for printing matrix
  call printMat(A,msize,nsize)
  mu = (/-8.1,7.9,5.6,-1.6/)
  do i = 1,size(mu)
    j = mu(i)
    write(*,*) '******************After Iteration Algorithm***************'
    ! Function call for Iteration algorithm
    call inverseIter(A,msize,nsize,j,xVec)
    ! Print eigenvalues (diagonal entries of matrix A)
    write(*,*) 'mu =', j
    ! Print eigenvectors (matrix V)
    write(*,*) 'eigenvector ='
    call printMat(xVec,msize,1)
  enddo
  write(*,*) '**********************************************************'
  deallocate(mat,A,xVec)
End Program Driver_LinAl