! File: LinAL.f90
! Author: Adrian Garcia
! Date: 02/10/2023
! DISCLAIMER: The original template code was given by AM213A Professor Dongwook Lee
module LinAl
  implicit none
  integer,save :: msize,nsize,Am,An,Bm,Bn,Xm,Xn,Em,En
  integer,dimension(:),allocatable,save :: s
  real(kind = kind(0.d0)),dimension(:,:),allocatable,save :: mat,A,ALU,B,X,E,L,U
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
  !********************Part 1: Problem 1*******************
  ! Input(s): Matrix A, first dimension of A
  ! Output(s): Trace of A
  ! Purpose: Calculates the Trace of A
  subroutine traceMat(A,m,Tr)
    implicit none
    real(kind = kind(0.d0)),dimension(:,:),intent(in) :: A
    real(kind = kind(0.d0)),intent(out) :: Tr
    integer :: i,m
    ! Initiallize the Trace as 0
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
    real(kind = kind(0.d0)),dimension(:),intent(in) :: x
    real(kind = kind(0.d0)),intent(out) :: enorm
    integer :: i
    ! Initiallize the Euclidean norm as 0
    enorm = 0.0
    ! Calculate the Euclidean norm
    do i = 1,m
      enorm = enorm+x(i)**2
    end do
    enorm  = sqrt(enorm)
  end subroutine enormVector
  !**********************************************************
  ! Input(s): Matrix A, dimensions of A (m rows, n columns)
  ! Output(s): None
  ! Purpose: Prints matrix A and its dimensions in human-readable format
  subroutine printMat(A,m,n)
    implicit none
    integer,intent(in) :: m,n
    real(kind = kind(0.d0)),dimension(:,:),intent(in) :: A
    integer :: i,j
    ! Print matrix A and its dimensions in human-readable format
    write(*,1) m,n
    1 format(2i4)
    do i = 1,m
      write(*,*) (A(i,j) , j = 1,n)
    end do
    write(*,*)
  end subroutine printMat
  !*******************Part 1: Problem 2*******************
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
    ! Initiallize singular flag as false 
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
        stop
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
  !*******************Part 1: Problem 3*******************
  ! Input(s): Matrix A, first dimension of A
  ! Output(s): Matrix A (modified: LU decomposition), permutation vector s, logical value indicating wheter the problem is singular or not
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
        stop
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
  ! Input(s): Matrix A (LU Decomposed), first dimension of A, matrix B, permutation vector s, logical value indicating wheter the problem is singular or not
  ! Output(s): Matrix X
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
          stop
        endif
        sum = 0.0
        do k = i+1,m
          sum = sum+A(i,k)*X(k,l)
        enddo
        X(i,l) = (y(i)-sum)/A(i,i)
      enddo
    enddo
  end subroutine LUBacksub
end module LinAl