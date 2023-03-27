! File: Driver_LinAL.f90
! Author: Adrian Garcia
! Date: 03/18/2023
Program Driver_LinAl
  use LinAl
  implicit none
  integer :: i,choice
  real(kind = kind(0.d0)) :: D,tol
  ! *********************Part 2: Problem 1********************
  msize = 10
  allocate(A(msize,msize),b(msize,1),x(msize,1))
  ! Initialize matrices
  A = 1.0
  b = 0.0
  x = 0.0
  tol = 1e-5
  ! Ask for D
  write(*,*) '*********************Part 1: Problem 1********************'
  write(*,*) 'D = '
  read(*,*) D
  ! Create matrix A
  do i = 1,msize
    A(i,i) = D
    b(i,1) = i
  end do
  ! Ask for choice
  write(*,*) 'Type 1, 2, or 3 to solve Ax = b using a Gauss-Jacobi (1), Gauss-Seidel (2), or Conjugate Gradient algorithm (3):'
  read(*,*) choice
  if (choice .eq. 1) then
    ! Function call for gaussJacobi
    call gaussJacobi(A,msize,b,tol,x)
    ! Print vector x
    write(*,*) 'Vector x = '
    call printMat(x,msize,1)
  else if (choice .eq. 2) then
    ! Function call for gaussSeidel
    call gaussSeidel(A,msize,b,tol,x)
    ! Print vector x
    write(*,*) 'Vector x = '
    call printMat(x,msize,1)
  else
    !*********************Part 2: Problem 2********************
    ! Function call for CGNoCond
    call CGNoCond(A,msize,b,tol,x)
    ! Print vector x
    write(*,*) 'Vector x = '
    call printMat(x,msize,1)
  end if
  write(*,*) '**********************************************************'
  deallocate(A,b,x)
  !****************Part 2: Problem 1 (MODIFIED)**************
  ! Ask for choice
  write(*,*) 'Do you want to run the a_ii = i case for problem 1? (yes (1) or no (0)):'
  read(*,*) choice
  if (choice .eq. 1) then
    allocate(A(msize,msize),b(msize,1),x(msize,1))
    ! Initialize matrices
    A = 1.0
    b = 0.0
    x = 0.0
    ! Create matrix A
    do i = 1,msize
      A(i,i) = i
      b(i,1) = i
    end do
    !***Function call for both algorithms for case a_ii = i****
    ! Function call for gaussJacobi
    call gaussJacobi(A,msize,b,tol,x)
    ! Print vector x
    write(*,*) 'Vector x = '
    call printMat(x,msize,1)
    write(*,*) '**********************************************************'
    ! Function call for gaussSeidel
    call gaussSeidel(A,msize,b,tol,x)
    ! Print vector x
    write(*,*) 'Vector x = '
    call printMat(x,msize,1)
    deallocate(A,b,x)
    write(*,*) '**********************************************************'
  end if
  !***************Part 2: Problem 2 (MODIFIED 1)*************
  ! Ask for choice
  write(*,*) 'Do you want to run the a_ii = i case for problem 2 (w/o conditioning)? (yes (1) or no (0)):'
  read(*,*) choice
  if (choice .eq. 1) then
    allocate(A(msize,msize),b(msize,1),x(msize,1))
    ! Initialize matrices
    A = 1.0
    b = 0.0
    x = 0.0
    ! Create matrix A
    do i = 1,msize
      A(i,i) = i
      b(i,1) = i
    end do
    !****Function call for a 10x10 matrix for case a_ii = i****
    ! Function call for CGNoCond
    call CGNoCond(A,msize,b,tol,x)
    ! Print vector x
    write(*,*) 'Vector x = '
    call printMat(x,msize,1)
    deallocate(A,b,x)
    write(*,*) '**********************************************************'
    msize = 100
    allocate(A(msize,msize),b(msize,1),x(msize,1))
    ! Initialize matrices
    A = 1.0
    b = 0.0
    x = 0.0
    ! Create matrix A
    do i = 1,msize
      A(i,i) = i
      b(i,1) = i
    end do
    !****Function call for a 100x100 matrix for case a_ii = i****
    ! Function call for CGNoCond
    call CGNoCond(A,msize,b,tol,x)
    ! Print vector x
    write(*,*) 'Vector x = '
    call printMat(x,msize,1)
    deallocate(A,b,x)
  end if
  !***************Part 2: Problem 2 (MODIFIED 2)*************
  ! Ask for choice
  write(*,*) 'Do you want to run the a_ii = i case for problem 2 (w/ conditioning)? (yes (1) or no (0)):'
  read(*,*) choice
  if (choice .eq. 1) then
    msize = 10
    allocate(A(msize,msize),b(msize,1),x(msize,1))
    ! Initialize matrices
    A = 1.0
    b = 0.0
    x = 0.0
    ! Create matrix A
    do i = 1,msize
      A(i,i) = i
      b(i,1) = i
    end do
    !****Function call for a 10x10 matrix for case a_ii = i****
    ! Function call for CGPreCond
    call CGPreCond(A,msize,b,tol,x)
    ! Print vector x
    write(*,*) 'Vector x = '
    call printMat(x,msize,1)
    deallocate(A,b,x)
    write(*,*) '**********************************************************'
    msize = 100
    allocate(A(msize,msize),b(msize,1),x(msize,1))
    ! Initialize matrices
    A = 1.0
    b = 0.0
    x = 0.0
    ! Create matrix A
    do i = 1,msize
      A(i,i) = i
      b(i,1) = i
    end do
    !****Function call for a 100x100 matrix for case a_ii = i****
    ! Function call for CGPreCond
    call CGPreCond(A,msize,b,tol,x)
    ! Print vector x
    write(*,*) 'Vector x = '
    call printMat(x,msize,1)
    deallocate(A,b,x)
  end if
End Program Driver_LinAl