! File: LinAL.f90
! Author: Adrian Garcia
! Date: 03/18/2023
! Purpose: Module file containing all relevant subroutines for AM213A
module LinAl
  implicit none
  integer,save :: msize,nsize
  real(kind = kind(0.d0)),dimension(:,:),allocatable,save :: mat,A,b,x
contains
  !********************************************************
  ! Input(s): Data file
  ! Output(s): Matrix
  ! Purpose: Reads a file containing the matrix A
  subroutine readMat(filename)
    implicit none
    integer :: i,j
    character(len = *) :: filename
    open(10,file = filename)
    ! Read the matrix dimensions
    read(10,*) i,j
    ! Read matrix
    do i = 1,msize
       read(10,*) (mat(i,j), j = 1,nsize)
    enddo
    close(10)
  end subroutine readMat
  !**********************************************************
  ! Input(s): Matrix A, first dimension of A
  ! Output(s): Trace of A
  ! Purpose: Calculates the Trace of A
  subroutine traceMat(A,m,Tr)
    implicit none
    real(kind = kind(0.d0)),dimension(:,:),intent(in) :: A
    real(kind = kind(0.d0)),intent(out) :: Tr
    integer :: i,m
    ! Initialize the Trace as 0
    Tr = 0.0
    ! Calculate the Trace
    do i = 1,m
      Tr = Tr+A(i,i)
    end do
  end subroutine traceMat
  !**********************************************************
  ! Input(s): Vector x, dimension of x
  ! Output(s): Euclidean norm of x
  ! Purpose: Calculates the Euclidean norm of x
  subroutine enormVector(x,m,enorm)
    implicit none
    integer,intent(in) :: m
    real(kind = kind(0.d0)),dimension(m),intent(in) :: x
    real(kind = kind(0.d0)),intent(out) :: enorm
    integer :: i
    ! Initialize the Euclidean norm as 0
    enorm = 0.0
    ! Calculate the Euclidean norm
    do i = 1,m
      enorm = enorm+x(i)**2
    end do
    enorm  = sqrt(enorm)
  end subroutine enormVector
  !**********************************************************
  ! Input(s): Matrix X, dimensions of X
  ! Output(s): Frobenius norm of X
  ! Purpose: Calculates the Frobenius norm of x
  subroutine fnormMat(X,m,n,fnorm)
    implicit none
    integer,intent(in) :: m,n
    real(kind = kind(0.d0)),dimension(m,n),intent(in) :: X
    real(kind = kind(0.d0)),intent(out) :: fnorm
    integer :: i,j
    ! Initialize the Euclidean norm as 0
    fnorm = 0.0
    ! Calculate the Frobenious norm
    do i = 1,m
      do j = 1,n
        fnorm = fnorm+X(i,j)**2
      enddo
    enddo
    fnorm = sqrt(fnorm)
  end subroutine fnormMat
  !**********************************************************
  ! Input(s): Matrix A, dimensions of A (m rows, n columns)
  ! Output(s): None
  ! Purpose: Prints matrix A and its dimensions in human-readable format
  subroutine printMat(A,m,n)
    implicit none
    integer,intent(in) :: m,n
    real(kind = kind(0.d0)),dimension(m,n),intent(in) :: A
    integer :: i,j
    ! Print matrix A and its dimensions in human-readable format
    write(*,1) m,n
    1 format(2i4)
    do i = 1,m
      write(*,*) (A(i,j) , j = 1,n)
    end do
    write(*,*)
  end subroutine printMat
  !**********************************************************
  ! Input(s): Matrix A, dimensions of A, matrix B, dimensions of B
  ! Output(s): Matrix A (modified: upper triangular), matrix B (modified), logical value indicating wheter the problem is singular or not
  ! Purpose: Performs Gaussian Elimination with partial pivoting
  subroutine gaussElim(A,Am,An,B,Bm,Bn,singular)
    implicit none
    integer,intent(in) :: Am,An,Bm,Bn
    logical,intent(out) :: singular
    real(kind = kind(0.d0)),dimension(Am,An),intent(inout) :: A
    real(kind = kind(0.d0)),dimension(Bm,Bn),intent(inout) :: B
    integer :: i,j,K
    real(kind = kind(0.d0)),dimension(An) :: temprowA
    real(kind = kind(0.d0)),dimension(Bn) :: temprowB
    ! Initialize singular flag as false 
    singular = .FALSE.
    ! Loop over columns
    do j = 1,Am-1
      ! Find index K and pivot p such that p = |a_{Kj}| = max_{k = j,...,m}|a_{kj}|
      K = maxloc(abs(A(j:,j)),dim = 1)+j-1
      if (K .ne. j) then
        ! Interchange rows K and j of matrix A
        temprowA = A(j,:)
        A(j,:) = A(K,:)
        A(K,:) = temprowA
        ! Interchange rows K and j of matrix B
        temprowB = B(j,:)
        B(j,:) = B(K,:)
        B(K,:) = temprowB
      endif
      if (A(j,j) .eq. 0) then
        ! Matrix is singular
        singular = .TRUE.
        write(*,*) 'ERROR: System is singular'
        exit
      endif
      ! Loop over rows below row j
      do i = j+1,Am
        ! Transformation of RHS vector(s)
        B(i,:) = B(i,:)-A(i,j)*B(j,:)/A(j,j)
        ! Transformation of remaining submatrix
        A(i,:) = A(i,:)-A(i,j)*A(j,:)/A(j,j)
      enddo
    enddo
  end subroutine gaussElim
  !**********************************************************
  ! Input(s): Matrix U (upper triangular), Matrix B
  ! Output(s): Matrix X
  ! Purpose: Performs back-substitution after Gaussian Elimination
  subroutine gaussBacksub(U,B,X)
    implicit none
    real(kind = kind(0.d0)),dimension(:,:),intent(in) :: U,B
    real(kind = kind(0.d0)),dimension(:,:),intent(out) :: X
    integer :: n,i,j,k,cols
    ! Find the number of columns of matrix B
    cols = size(B,dim=2)
    ! Loop over each column
    do k = 1,cols
      ! Find size of vector B(:,k)
      n = size(B(:,k))
      ! Solve last row first
      X(n,k) = B(n,k)/U(n, n)
      ! Work backwards through rows
      do i = n-1,1,-1
        ! Initialize X(i,k) as B(i,k)
        X(i,k) = B(i,k)
        do j = i+1,n
          ! Subtract contributions from previously found unknowns
          X(i,k) = X(i,k)-U(i,j)*X(j,k)
        enddo
        ! Divide by the coefficient of X(i,k) to find X(i,k)
        X(i,k) = X(i,k)/U(i,i)
      enddo
    end do
  end subroutine gaussBacksub
  !**********************************************************
  ! Input(s): Matrix A, first dimension of A
  ! Output(s): Matrix A (modified: LU Decomposition), permutation vector s, logical value indicating wheter the problem is singular or not
  ! Purpose: Performs LU Decomposition with partial pivoting
  subroutine LUDecomp(A,m,s,singular)
    implicit none
    integer,intent(in) :: m
    logical,intent(out):: singular
    integer,dimension(:),intent(out) :: s
    real(kind = kind(0.d0)),dimension(:,:),intent(inout):: A
    integer :: i,j,l,K,tempval
    real(kind = kind(0.d0)),dimension(m) :: temp
    ! Initialize singular flag as false
    singular = .FALSE.
    ! Initialize permutation vector s
    do j = 1,m 
      s(j) = j 
    enddo 
    ! Loop over columns
    do j = 1, m
      ! Find index K and pivot p such that p = |a_{Kj}| = max_{k = j,...,m}|a_{kj}|
      K = maxloc(abs(A(j:,j)),dim = 1)+j-1
      if (K .ne. j) then
        ! Interchange rows K and j of matrix A
        temp = A(j,:)
        A(j,:) = A(K,:)
        A(K,:) = temp
        ! Interchange K and j entries of vector s
        tempval = s(j)
        s(j) = s(K)
        s(K) = tempval
      endif
      if (A(j,j) .eq. 0) then
        ! Matrix is singular
        singular = .TRUE.
        write(*,*) 'ERROR: System is singular'
        exit
      endif
      ! Loop over rows below row j
      do i = j+1,m
        ! Create the l_{ij} and stores them in a_{ij}
        A(i,j) = A(i,j)/A(j,j)
        do l = j+1,m
          ! Update matrix A
          A(i,l) = A(i,l)-A(i,j)*A(j,l)
        enddo
      enddo
    enddo
  end subroutine LUDecomp
  !**********************************************************
  ! Input(s): Matrix A (LU Decomposed), first dimension of A, matrix B, permutation vector s
  ! Output(s): Matrix X, logical value indicating wheter the problem is singular or not
  ! Purpose: Performs back-substitution after LU Decomposition
  subroutine LUBacksub(A,m,B,s,X,singular)
    implicit none
    integer,intent(in) :: m
    integer,dimension(:),intent(in) :: s
    logical,intent(out):: singular
    real(kind = kind(0.d0)),dimension(:,:),intent(in):: A, B
    real(kind = kind(0.d0)),dimension(:,:),intent(out) :: X
    integer :: i,j,k,l,cols
    real(kind = kind(0.d0)),dimension(m) :: y
    real(kind = kind(0.d0)) :: sum
    ! Find the number of columns of matrix B
    cols = size(B,dim=2)
    ! Loop over each column
    do l = 1,cols
      ! Initialize singular flag as false
      singular = .FALSE.
      ! Initialize y with Pb
      do i = 1,m
        y(i) = B(s(i),l)
      enddo
      ! Forward substitution, y = L^{-1}Pb
      do j = 1,m-1
        ! Do y = M_{j}y
        do i = j+1,m
          y(i) = y(i)-y(j)*A(i,j)
        enddo
      enddo
      ! Backward substitution, Ux = y
      do i = m,1,-1
        ! Loop over rows from bottom to top
        if (A(i,i) .eq. 0) then
          ! Matrix is singular
          singular = .TRUE.
          write(*,*) 'ERROR: System is singular'
          exit
        endif
        sum = 0.0
        do k = i+1,m
          sum = sum+A(i,k)*X(k,l)
        enddo
        X(i,l) = (y(i)-sum)/A(i,i)
      enddo
    enddo
  end subroutine LUBacksub
  !**********************************************************
  ! Input(s): Matrix A
  ! Output(s): Matrix A (modified: Cholesky Decomposed), logical value indicating wheter the problem is singular or not and SDP or not
  ! Purpose: Performs Cholesky Decomposition
  subroutine choDecomp(A,singular,SPD)
    implicit none
    logical,intent(out):: singular, SPD
    real(kind = kind(0.d0)),dimension(:,:),intent(inout):: A
    integer :: i,j,k,m
    ! Initialize singular flag as false
    singular = .FALSE.
    ! Initialize SPD flag as false
    SPD = .FALSE.
    ! Find the number of columns of matrix A
    m = size(A,dim=2)
    ! Loop over each column
    do j = 1,m
      ! Calculate new diagonal element
      do k = 1,j-1
        A(j,j) = A(j,j)-A(j,k)*A(j,k)
        if (A(j,j) .eq. 0) then
          ! Matrix is singular
          singular = .TRUE.
          write(*,*) 'ERROR: System is singular'
          exit
        else if (A(j,j) < 0) then
          ! Matrix is not SPD
          SPD = .TRUE.
          write(*,*) 'ERROR: System is not SPD'
          exit
        endif
      enddo
      A(j,j) = sqrt(A(j,j)) 
      ! Calculate elements below diagonal
      do i = j+1,m
        do k = 1,j-1
          A(i,j) = A(i,j) - A(i,k)*A(j,k)
        enddo
        A(i,j) = A(i,j)/A(j,j)
        if (A(i,j) .eq. 0 .or. A(i,j) < 0) then
          ! Matrix is not SPD
          SPD = .TRUE.
          write(*,*) 'ERROR: System is not SPD'
          exit
        endif
      enddo
    enddo
  end subroutine choDecomp
  !**********************************************************
  ! Input(s): Matrix A (Cholesky Decomposed), first dimension of A, matrix B
  ! Output(s): Matrix X, logical value indicating wheter the problem is singular or not
  ! Purpose: Performs back-substitution after Cholesky Decomposition
  subroutine choBacksub(A,m,B,X,singular)
    implicit none
    integer,intent(in) :: m
    logical,intent(out):: singular
    real(kind = kind(0.d0)),dimension(:,:),intent(in):: A, B
    real(kind = kind(0.d0)),dimension(:,:),intent(out) :: X
    integer :: i,j,k,l,cols
    real(kind = kind(0.d0)) :: sum
    real(kind = kind(0.d0)),dimension(m) :: y
    ! Find the number of columns in matrix B
    cols = size(B,dim=2)
    ! Initialize singular flag as false
    singular = .FALSE.
    ! Loop over each column
    do l = 1,cols
      ! Initialize y with B
      do i = 1,m
        y(i) = B(i,l)
      enddo
      ! Forward substitution, Ly = b
      do i = 1,m
        sum = B(i,l)
        do j = 1, i-1 
          sum = sum-y(j)*A(i,j)
        enddo
        y(i) = sum/A(i,i)  
      enddo
      ! Backward substitution, L^T x = y
      do i = m,1,-1
        if (A(i,i) .eq. 0) then
          ! Matrix is singular
          singular = .TRUE.
          write(*,*) 'ERROR: System is singular'
          exit
        endif
        do k = i+1,m 
          y(i) = y(i)-A(k,i)*y(k) 
        enddo
        y(i) = y(i)/A(i,i) 
      enddo
      do i = 1,m
        X(i,l) = y(i)
      enddo
    enddo
  end subroutine choBacksub
  !**********************************************************
  ! Input(s): Dimensions of desired Vandermonde matrix, vector X
  ! Output(s): Vandermonde matrix A
  ! Purpose: Creates a Vandermonde matrix given a vector X
  subroutine vanderMat(A,m,n,X)
    implicit none
    integer,intent(in) :: m,n
    real(kind = kind(0.d0)),dimension(m,1),intent(in):: X
    real(kind = kind(0.d0)),dimension(m,n),intent(out) :: A
    integer :: i
    ! Initialize the Vandermonde matrix as 1
    A = 1.0
    ! Create Vandermonde matrix
    do i = 2,n
      A(:,i) = X(:,1)**(i-1) 
    enddo
  end subroutine vanderMat
  !**********************************************************
  ! Input(s): Martix A, martix B
  ! Output(s): Martix A (modified: normal equation), martix B (modified: normal equation)
  ! Purpose: Forms the normal equation
  subroutine normalEq(A,B,ANormal,BNormal) 
    implicit none
    real(kind = kind(0.d0)),dimension(:,:),intent(in) :: A,B
    real(kind = kind(0.d0)),dimension(:,:),intent(out) :: ANormal,BNormal
    ! LHS of the normal equation
    ANormal = matmul(transpose(A),A)
    ! RHS of the normal equation
    BNormal = matmul(transpose(A),B)
  end subroutine normalEq
  !**********************************************************
  ! Input(s): Matrix A, the dimensions of A
  ! Output(s): Matrix Q, Matrix R
  ! Purpose: Performs a Householder based QR Decomposition
  subroutine houseQRDecomp(A,m,n,Q,R)
    implicit none
    integer,intent(in) :: m,n
    real(kind = kind(0.d0)),dimension(m,n),intent(in) :: A
    real(kind = kind(0.d0)),dimension(m,m),intent(out) :: Q
    real(kind = kind(0.d0)),dimension(m,n),intent(out) :: R
    integer :: i,j
    real(kind = kind(0.d0)) :: s
    real(kind = kind(0.d0)),dimension(m,1) :: v
    real(kind = kind(0.d0)),dimension(m,m) :: IMat
    ! Initialize Identity matrix
    IMat = 0.0
    do i = 1,m
      IMat(i,i) = 1.0
    end do
    ! Initialize R as A
    R = A
    ! Loop over columns
    do j = 1,n
      ! Compute signed norm
      s = sign(norm2(R(j:,j)),R(j,j))
      ! Compute Householder vector and normalize it
      v = 0.0
      v(j,1) = R(j,j)+s
      v(j+1:,1)= R(j+1:,j)
      v = v/norm2(v(:,1))
      ! Update R and Q
      R = R-2*matmul(matmul(v,transpose(v)),R)
      if (j .eq. 1) then
        Q = IMat-2*matmul(v,transpose(v))
      else 
        Q = matmul(Q,(IMat-2*matmul(v,transpose(v))))
      endif
    end do
  end subroutine houseQRDecomp
  !**********************************************************
  ! Input(s): Matrix A, the dimensions of A
  ! Output(s): Martix A (modified: Hessenberg form)
  ! Purpose: Performs Hessenberg Decomposition
  subroutine hessenDecomp(A,m,n)
    implicit none
    integer,intent(in) :: m,n
    real(kind = kind(0.d0)),dimension(m,n),intent(inout) :: A
    integer :: j
    real(kind = kind(0.d0)) :: s
    real(kind = kind(0.d0)),dimension(m,1) :: v
    ! Loop over columns
    do j = 1,m-2
      ! Compute signed norm
      s = sign(norm2(A(j+1:,j)),A(j+1,j))
      ! Compute Householder vector and normalize it
      v = 0.0
      v(j+1,1) = A(j+1,j)+s
      v(j+2:,1)= A(j+2:,j)
      v = v/norm2(v(:,1))
      ! Update A from the left
      A = A-2*matmul(matmul(v,transpose(v)),A)
      ! Update A from the right
      A = A-2*matmul(A,matmul(v,transpose(v)))
    end do
  end subroutine hessenDecomp
  !**********************************************************
  ! Input(s): Matrix A, the dimesions of A
  ! Output(s): Matrix A(modified: diagonal matrix), Matrix V
  ! Purpose: Performs QR Decomposition without shifts to calculate the eigenvectors/values of matrix A
  subroutine QRNoShift(A,m,n,V)
    implicit none
    integer,intent(in) :: m,n
    real(kind = kind(0.d0)),dimension(m,n),intent(inout) :: A
    real(kind = kind(0.d0)),dimension(m,n),intent(out) :: V
    integer :: i
    real(kind = kind(0.d0)) :: err
    real(kind = kind(0.d0)), dimension(m) :: errv
    real(kind = kind(0.d0)),dimension(m,n) :: R
    real(kind = kind(0.d0)),dimension(m,m) :: Q
    ! Initialize V
    V = 0.0
    do i = 1,m
      V(i,i) = 1.0
    end do
    ! Initialize error
    err = 1.0
    do while (err > 1e-8)
      ! Compute QR factorization of A
      call houseQRDecomp(A,m,n,Q,R)
      ! Replace A by product RQ
      A = matmul(R,Q)
      ! Update the eigenvector matrix
      V = matmul(V,Q)
      ! Update error
      errv = 0.0
      do i = 1,n-1
        errv(i) = abs(A(i+1,i))
      enddo
      err = norm2(errv)
    enddo
  end subroutine QRNoShift
  !**********************************************************
  ! Input(s): Matrix A, the dimesions of A
  ! Output(s): Matrix A(modified: diagonal matrix), Matrix V
  ! Purpose: Performs QR Decomposition with shifts to calculate the eigenvectors/values of matrix A
  subroutine QRShift(A,m,n,V)
    implicit none
    integer,intent(in) :: m,n
    real(kind = kind(0.d0)),dimension(:,:),intent(inout) :: A
    real(kind = kind(0.d0)),dimension(m,n),intent(out) :: V
    integer :: j
    real(kind = kind(0.d0)) :: err,mu
    real(kind = kind(0.d0)), dimension(m) :: errv
    real(kind = kind(0.d0)),dimension(m,n) :: R,I
    real(kind = kind(0.d0)),dimension(m,m) :: Q
    ! Initialize V and Identity matrix
    V = 0.0
    do j = 1,m
      V(j,j) = 1.0
    end do
    I = V
    ! Initialize error
    err = 1.0
    do while(err > 2e-2)
      ! Select shift
      mu = A(m,m)
      ! Compute QR factorization of A−mu*I
      call houseQRDecomp(A-mu*I,m,n,Q,R)
      ! Replace A by RQ+mu*I
      A = matmul(R,Q)+mu*I
      ! Update the eigenvector matrix
      V = matmul(V,Q)
      ! Update error
      errv = 0.0
      do j = 1,n-1
        errv(j) = A(j+1,j)
      enddo 
      err = norm2(errv)
    enddo 
  end subroutine QRShift
  !**********************************************************
  ! Input(s): Matrix A, the dimesions of A, scalar mu (estimated eigenvalue)
  ! Output(s): Vector x (eigenvector associated to estimated eigenvalue mu)
  ! Purpose: Performs an Inverse Iteration algorithm to calculate an eigenvector of matrix A
  subroutine inverseIter(A,m,n,mu,x)
    implicit none
    integer,intent(in) :: m,n
    real(kind = kind(0.d0)),intent(in) :: mu
    real(kind = kind(0.d0)),dimension(m,n),intent(inout) :: A
    real(kind = kind(0.d0)),dimension(m,1),intent(out) :: x
    integer :: j
    logical :: singular
    real(kind = kind(0.d0)) :: err
    real(kind = kind(0.d0)),dimension(m,1) :: r,y
    real(kind = kind(0.d0)),dimension(m,n) :: B,I,INV
    ! Initialize Identity matrix
    I = 0.0
    do j = 1,m
      I(j,j) = 1.0
    end do
    ! Initialize B: using Gaussian Elimination to find (A-mu*I)^(-1)
    INV = A-mu*I
    call gaussElim(INV,m,n,I,m,n,singular)
    call gaussBacksub(INV,I,B)
    ! Initialize x
    x = 0.0
    x(m,1) = 1.0
    ! Initialize error
    err = 1.0
    do while (err > 1e-8)
      ! Calculate new normalized vector
      y = matmul(B,x)/norm2(matmul(B,x))
      ! Calculate difference between old and new
      r = y-x
      ! Replace old by new
      x = y
      ! Update error
      err = norm2(r)
    enddo
  end subroutine inverseIter
  !*********************Part 2: Problem 1********************
  ! Input(s): Matrix A (square matrix), the first dimension of A, vector b, desired tolerance
  ! Output(s): Vector x
  ! Purpose: Solves the equation Ax = b using the Gauss-Jacobi algorithm
  subroutine gaussJacobi(A,m,b,tol,x)
    implicit none
    integer,intent(in) :: m
    real(kind = kind(0.d0)),intent(in) :: tol
    real(kind = kind(0.d0)),dimension(m,m),intent(in) :: A
    real(kind = kind(0.d0)),dimension(m,1),intent(in) :: b
    real(kind = kind(0.d0)),dimension(m,1),intent(out) :: x
    integer :: i,j,iter
    character(len = 100) :: name
    character(len = 10) :: str
    real(kind = kind(0.d0)) :: err
    real(kind = kind(0.d0)),dimension(m,1) :: rVec,y
    real(kind = kind(0.d0)),dimension(m,m) :: DInv,R
    ! Initialize DInv
    DInv = 0.0
    do i = 1,m
      DInv(i,i) = 1/A(i,i)
    end do
    ! Initialize R
    R = 0.0
    do i = 1,m
      do j = 1,m
        if (i.NE.j) then
          R(i,j) = A(i,j)
        end if
      end do
    end do
    ! Initialize x
    x = 0.0
    x(m,1) = 1.0
    ! Initialize error and iter
    err = 1.0
    iter = 0
    do while (err > tol)
      ! If the algorithm does not converge after 1000 iterations, then exit
      if (iter .eq. 1000) then
        exit
      end if
      ! Calculate new vector
      y = matmul(DInv,(b-matmul(R,x)))
      ! Calculate a new residual
      rVec = b-matmul(A,x)
      ! Replace old by new
      x = y
      ! Update error and iter
      err = norm2(rVec)
      iter = iter+1
      ! Save to error data file
      write(str,'(4I0)') int(A(1,1))
      name = trim('GJError')//trim(str)//trim('.dat')
      open(10,file = name,action = 'write',position = 'append')
      write(10,*) int(A(1,1)), iter, err
      close(10)
    enddo
  end subroutine gaussJacobi
  !**********************************************************
  ! Input(s): Matrix A (square matrix), the first dimension of A, vector b, desired tolerance
  ! Output(s): Vector x
  ! Purpose: Solves the equation Ax = b using the Gauss-Seidel algorithm
  subroutine gaussSeidel(A,m,b,tol,x)
    implicit none
    integer,intent(in) :: m
    real(kind = kind(0.d0)),intent(in) :: tol
    real(kind = kind(0.d0)),dimension(m,m),intent(in) :: A
    real(kind = kind(0.d0)),dimension(m,1),intent(in) :: b
    real(kind = kind(0.d0)),dimension(m,1),intent(out) :: x
    integer :: i,j,iter
    logical :: singular
    character(len = 100) :: name
    character(len = 10) :: str
    real(kind = kind(0.d0)) :: err
    real(kind = kind(0.d0)),dimension(m,1) :: r
    real(kind = kind(0.d0)),dimension(m,m) :: D,L,U,DLInv,INV,IMat,ICopy
    ! Initialize Identity matrix
    IMat = 0.0
    do j = 1,m
      IMat(j,j) = 1.0
    end do
    ! Initialize D
    D = 0.0
    do i = 1,m
      D(i,i) = A(i,i)
    end do
    ! Initialize L
    L = 0.0
    do j = 1, m
      do i = j+1,m
        L(i,j) = A(i,j)
      enddo
    enddo
    ! Initialize U
    U = A-L-D
    ! Initialize x
    x = 0.0
    x(m,1) = 1.0
    ! Initialize error and iter
    err = 1.0
    iter = 0
    do while (err > tol)
      ! If the algorithm does not converge after 1000 iterations, then exit
      if (iter .eq. 1000) then
        exit
      end if
      ! Initialize (D+L)^-1: using Gaussian Elimination
      ICopy = IMat
      INV = D+L
      call gaussElim(INV,m,m,ICopy,m,m,singular)
      call gaussBacksub(INV,ICopy,DLInv)
      ! Calculate new vector
      x = matmul(DLInv,(b-matmul(U,x)))
      ! Calculate a new residual
      r = b-matmul(A,x)
      ! Update error and iter
      err = norm2(r)
      iter = iter+1
      ! Save to error data file
      write(str,'(4I0)') int(A(1,1))
      name = trim('GSError')//trim(str)//trim('.dat')
      open(10,file = name,action = 'write',position = 'append')
      write(10,*) int(A(1,1)), iter, err
      close(10)
    enddo
  end subroutine gaussSeidel
  !*********************Part 2: Problem 2********************
  ! Input(s): Matrix A (square matrix), the first dimension of A, vector b, desired tolerance
  ! Output(s): Vector x
  ! Purpose: Performs the smart conjugate gradient algorithm without pre-conditioning to solve the problem Ax = b. 
  subroutine CGNoCond(A,m,b,tol,x)
    implicit none
    integer,intent(in) :: m
    real(kind = kind(0.d0)),intent(in) :: tol
    real(kind = kind(0.d0)),dimension(m,m),intent(in) :: A
    real(kind = kind(0.d0)),dimension(m,1),intent(in) :: b
    real(kind = kind(0.d0)),dimension(m,1),intent(out) :: x
    integer :: i,j,iter
    character(len = 100) :: name
    character(len = 10) :: str1,str2
    real(kind = kind(0.d0)) :: err,errNew
    real(kind = kind(0.d0)) :: beta,alpha
    real(kind = kind(0.d0)),dimension(m,1) :: r,p,y
    ! Check to see if matrix A is symmetric
    do i =1,m
      do j=1,m
        if (A(i,j) .ne. A(j,i)) then
          ! Matrix is not symmetric
          write(*,*) 'ERROR: Matrix is not symmetric'
          exit
        endif
      end do
    end do
    ! Initialize x
    x = 0.0
    ! Initialize r
    r = b-matmul(A,x)
    ! Initialize p
    p = r
    ! Initialize error and iter
    err = norm2(r)**2
    iter = 0
    do while (sqrt(err) > tol)
      ! If the algorithm does not converge after 1000 iterations, then exit
      if (iter .eq. 1000) then
        exit
      end if
      ! Calculate temporary vector
      y = matmul(A,p)
      ! Calculate alpha
      alpha = err/dot_product(p(:,1),y(:,1))
      ! Update x
      x = x+alpha*p
      ! We use this formula for r (which uses y) to avoid calculating Ax
      r = r-alpha*y
      ! Calculate E_new
      errNew = norm2(r)**2
      ! Calculate beta
      beta = errNew/err
      ! Update p
      p = r+beta*p
      ! Update error and iter
      err = errNew
      iter = iter+1
      ! Save to error data file
      write(str1,'(4I0)') int(m)
      write(str2,'(4I0)') int(A(1,1))
      name = trim('CGNoCond')//trim(str1)//trim(str2)//trim('.dat')
      open(10,file = name,action = 'write',position = 'append')
      write(10,*) int(A(1,1)), iter, err
      close(10)
    enddo
  end subroutine CGNoCond
  !**********************************************************
  ! Input(s): Matrix A (square matrix), the first dimension of A, vector b, desired tolerance
  ! Output(s): Vector x
  ! Purpose: Performs the smart conjugate gradient algorithm with pre-conditioning to solve the problem Ax = b.
  ! NOTE: All new lines of code have been labeled as "! NEW:..."
  subroutine CGPreCond(A,m,b,tol,x)
    implicit none
    integer,intent(in) :: m
    real(kind = kind(0.d0)),intent(in) :: tol
    real(kind = kind(0.d0)),dimension(m,m),intent(in) :: A
    real(kind = kind(0.d0)),dimension(m,1),intent(in) :: b
    real(kind = kind(0.d0)),dimension(m,1),intent(out) :: x
    integer :: i,j,iter
    character(len = 100) :: name
    character(len = 10) :: str1,str2
    real(kind = kind(0.d0)) :: err,errNew
    real(kind = kind(0.d0)) :: beta,alpha
    real(kind = kind(0.d0)),dimension(m,1) :: r,p,y,z
    real(kind = kind(0.d0)),dimension(m,m) :: MInv
    do i =1,m
      do j=1,m
        if (A(i,j) .ne. A(j,i)) then
          write(*,*) 'ERROR: Matrix is not symmetric'
          exit
        endif
      end do
    end do
    x = 0.0
    r = b-matmul(A,x)
    ! NEW: Initialize MInv
    MInv = 0.0
    do i = 1,m
      MInv(i,i) = 1/A(i,i)
    end do
    ! NEW: Initialize z
    z = matmul(MInv,r)
    ! NEW: Initialize p
    p = z
    ! NEW: Initialize error and iter
    err = dot_product(r(:,1),z(:,1))
    iter = 0
    do while (sqrt(err) > tol)
      if (iter .eq. 1000) then
        exit
      end if
      y = matmul(A,p)
      alpha = err/dot_product(p(:,1),y(:,1))
      x = x+alpha*p
      r = r-alpha*y
      ! NEW: Update z
      z = matmul(MInv,r)
      ! NEW: Calculate E_new
      errNew = dot_product(r(:,1),z(:,1))
      beta = errNew/err
      ! NEW: Update p
      p = z+beta*p
      err = errNew
      iter = iter+1
      write(str1,'(4I0)') int(m)
      write(str2,'(4I0)') int(A(1,1))
      name = trim('CGPreCond')//trim(str1)//trim(str2)//trim('.dat')
      open(10,file = name,action = 'write',position = 'append')
      write(10,*) int(A(1,1)), iter, err
      close(10)
    enddo
  end subroutine CGPreCond
end module LinAl