! File: EX12.f90
! Author: Adrian Garcia
! Date: 1/16/2023
! Purpose: Exercise 12 from Homework 1 for AM213A
program addmats
    implicit none
    integer, parameter :: dimmat = 3
    real, dimension(dimmat,dimmat) :: a,b,c
    integer :: i,j
    ! This creates the matrices.
    a(1,2) = 2.0
    do i=2,dimmat-1
        a(i,i+1) = 2.0
        b(i,i-1) = 1.0
    enddo
    b(dimmat,dimmat-1) = 1.0
    ! This adds the matrices a and b
    do i=1,dimmat
        do j=1,dimmat
        c(i,j) = a(i,j)+b(i,j)
        enddo
        ! This prints c
        write(*,*) (c(i,j),j = 1,dimmat)
    enddo
end program addmats