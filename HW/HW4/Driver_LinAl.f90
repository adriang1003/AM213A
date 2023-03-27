! File: Driver_LinAL.f90
! Author: Adrian Garcia
! Date: 02/26/2023
! Purpose: Driver file for Homework 4 Part 1: Code.
! Run Command: Type 'make' to command window to create executable 'LinAl.exe'. Type './LinAl.exe' to run executable.
Program Driver_LinAl
  use LinAl
  implicit none
  character(len = 100) :: myFileNameA
  integer :: i,degree
  real(kind = kind(0.d0)) :: enormX, fnormE
  real(kind = kind(0.d0)),dimension(21,21) :: IMat
  myFileNameA = 'atkinson.dat'
  open(10,file = myFileNameA)
  ! This was personally added to the 'atkinson.dat' file
  read(10,*) msize,nsize
  close(10)
  allocate(mat(msize,nsize),xVec(msize,1),B(msize,1))
  ! Initialize matrices as 0
  mat = 0.0
  xVec = 0.0
  B = 0.0
  call readMat(myFileNameA)
  xVec(:,1) = mat(:,1)
  B(:,1) = mat(:,2)
  !*********************Part 1: Problem 1********************
  write(*,*) '*********************Part 1: Problem 1********************'
  ! Loop over degree values
  do degree = 3,5,2
    allocate(A(msize,degree+1),X(degree+1,1),E(degree+1,1),ANormal(degree+1,degree+1),BNormal(degree+1,1))
    ! Initialize matrices as 0
    A = 0.0
    X = 0.0
    E = 0.0
    ANormal = 0.0
    BNormal = 0.0
    ! Function call to create a Vandermonde matrix
    call vanderMat(A,msize,degree+1,xVec)
    ! Function call to form the normal equation
    call normalEq(A,B,ANormal,BNormal)
    write(*,*) '**********************Normal Equation*********************'
    write(*,*) 'Degree = ',degree
    ! Print matrix A and B (before Cholesky Decomposition)
    write(*,*) 'Matrix A ='
    call printMat(ANormal,degree+1,degree+1)
    write(*,*) 'Matrix B ='
    call printMat(BNormal,degree+1,1)
    ! Function call to perform Cholesky Decomposition
    call ChoDecomp(ANormal,singular,SPD)
    ! Print matrix A (after Cholesky Decomposition)
    write(*,*) '****************After Cholesky Decomposition**************'
    write(*,*) 'Matrix A ='
    call printMat(ANormal,degree+1,degree+1)
    write(*,*) '******Back substitution after Cholesky Decomposition******'
    ! Function call for back substitution after Cholesky Decomposition
    call ChoBacksub(ANormal,degree+1,BNormal,X,singular)
    ! Print solution matrix X
    write(*,*) 'Matrix X ='
    call printMat(X,degree+1,1)
    ! Save solution matrix X
    if (degree .eq. 3) then 
      open(10,file = 'ChoDegree3.dat')
      do i = 1,size(X,dim=1)
        write(10,*) X(i,:)
      enddo
      close(10)
    else 
      open(10,file = 'ChoDegree5.dat')
      do i = 1,size(X,dim=1)
        write(10,*) X(i,:)
      enddo
      close(10)
    endif
    ! Calculate error matrix E
    E = matmul(A,X)-B
    write(*,*) '***************************Error**************************'
    ! Function call for Euclidean norm
    call enormVector(E(:,1),degree+1,enormX)
    ! Print Euclidean norm
    write(*,*) 'Error: ||x||_2 =',enormX
    write(*,*) '**********************************************************'
    deallocate(A,X,E,ANormal,BNormal)
  enddo
  !********************High Degree Polynomial****************
  degree = 8
  allocate(A(msize,degree+1),X(degree+1,1),E(degree+1,1),ANormal(degree+1,degree+1),BNormal(degree+1,1))
  ! Initialize matrices as 0
  A = 0.0
  X = 0.0
  E = 0.0
  ANormal = 0.0
  BNormal = 0.0
  ! Function call to create a Vandermonde matrix
  call vanderMat(A,msize,degree+1,xVec)
  ! Function call to form the normal equation
  call normalEq(A,B,ANormal,BNormal)
  write(*,*) '**********************Normal Equation*********************'
  write(*,*) '"Large" Degree = ',degree
  ! Print matrix A and B (before Cholesky Decomposition)
  write(*,*) 'Matrix A ='
  call printMat(ANormal,degree+1,degree+1)
  write(*,*) 'Matrix B ='
  call printMat(BNormal,degree+1,1)
  ! Function call to perform Cholesky Decomposition
  call ChoDecomp(ANormal,singular,SPD)
  ! Print matrix A (after Cholesky Decomposition)
  write(*,*) '****************After Cholesky Decomposition**************'
  write(*,*) 'Matrix A ='
  call printMat(ANormal,degree+1,degree+1)
  write(*,*) '******Back substitution after Cholesky Decomposition******'
  ! Function call for back substitution after Cholesky Decomposition
  call ChoBacksub(ANormal,degree+1,BNormal,X,singular)
  ! Print solution matrix X
  write(*,*) 'Matrix X ='
  call printMat(X,degree+1,1)
  ! Calculate error matrix E
  E = matmul(A,X)-B
  write(*,*) '***************************Error**************************'
  ! Function call for Euclidean norm
  call enormVector(E(:,1),degree+1,enormX)
  ! Print Euclidean norm
  write(*,*) 'Error: ||x||_2 =',enormX
  write(*,*) '**********************************************************'
  deallocate(A,X,E,ANormal,BNormal)
  ! *********************Part 1: Problem 2********************
  write(*,*) '*********************Part 1: Problem 2********************'
  ! Loop over degree values
  do degree = 3,5,2
    allocate(A(msize,degree+1),X(degree+1,1),E(msize,degree+1),Q(msize,msize),R(msize,degree+1))
    allocate(Qhat(msize,degree+1),Rhat(degree+1,degree+1))
    ! Initialize matrices as 0
    A = 0.0
    X = 0.0
    E = 0.0
    Q = 0.0
    R = 0.0
    Qhat = 0.0
    Rhat = 0.0
    ! Function call to create a Vandermonde matrix
    call vanderMat(A,msize,degree+1,xVec)
    write(*,*) 'Degree = ',degree
    ! Function call to perform QR factorization
    call HouseQR(A,msize,degree+1,Q,R)
    write(*,*) '*******************After QR Decomposition*****************'
    ! Calculate error matrix E
    E = A-matmul(Q,R)
    ! Print E
    write(*,*) 'A - QR ='
    call printMat(E,msize,msize)
    write(*,*) '***************************Error**************************'
    ! Function call for Frobenius norm
    call fnormMat(E,msize,degree+1,fnormE)
    ! Print Frobenious norm
    write(*,*) 'Error: ||A - QR||_F =',fnormE
    write(*,*) '**********************************************************'
    ! Reset E
    E = 0.0
    ! Initialize Identity matrix
    IMat = 0.0
    do i = 1,msize
      IMat(i,i) = 1
    end do
    ! Calculate error matrix E
    E = matmul(transpose(Q),Q)-IMat
    ! Print E
    write(*,*) 'Q^TQ - I ='
    call printMat(E,msize,msize)
    write(*,*) '***************************Error**************************'
    ! Function call for Frobenius norm
    call fnormMat(E,msize,degree+1,fnormE)
    ! Print Frobenious norm
    write(*,*) 'Error: ||Q^TQ - I||_F =',fnormE
    write(*,*) '**********************************************************'
    ! Create Qhat and Rhat
    Qhat = Q(:,1:degree+1)
    Rhat = R(1:degree+1,:)
    ! Function call for back substitution (using gaussBacksub since R is an upper triangular matrix)
    call gaussBacksub(Rhat,matmul(transpose(Qhat),B),X)
    ! Print solution matrix X
    write(*,*) 'Matrix X ='
    call printMat(X,degree+1,1)
    ! Save solution matrix X
    if (degree .eq. 3) then 
      open(10,file = 'QRDegree3.dat')
      do i = 1,size(X,dim=1)
        write(10,*) X(i,:)
      enddo
      close(10)
    else 
      open(10,file = 'QRDegree5.dat')
      do i = 1,size(X,dim=1)
        write(10,*) X(i,:)
      enddo
      close(10)
    endif
    ! Calculate error matrix E
    E = matmul(A,X)-B
    write(*,*) '***************************Error**************************'
    ! Function call for Euclidean norm
    call enormVector(E(:,1),degree+1,enormX)
    ! Print Euclidean norm
    write(*,*) 'Error: ||x||_2 =',enormX
    write(*,*) '**********************************************************'
    deallocate(Qhat,Rhat)
    deallocate(A,X,E,Q,R)
  enddo
  deallocate(mat,xVec,B)
End Program Driver_LinAl