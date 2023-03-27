! File: EX15.f90
! Author: Adrian Garcia
! Date: 1/16/2023
! Purpose: Exercise 15 from Homework 1 for AM213A
program EX15
    implicit none
    integer, parameter :: dimmat = 10000
    real, dimension(dimmat,dimmat) :: a,b,c
    integer :: i!,j
    ! This creates the matrices a and b
    a = 0.0
    b = 0.0
    a(1,2) = 1.0
    do i=2,dimmat-1
        a(i,i+1) = 1.0
        b(i,i-1) = 1.0
    enddo
    b(dimmat,dimmat-1) = 1.0
    a = 2*a
    ! This adds the matrices a and b
    c = a+b
    ! This prints c
    ! do i=1,dimmat
    !     write(*,*) (c(i,j),j=1,dimmat)
    ! enddo
end program EX15