! File: geometric.f90
! Author: Adrian Garcia
! Date: 1/16/2023
! Purpose: Beginner Fortran 90 tutorial
program geometric
    ! Exercise 6
    implicit none
    integer :: iter
    real :: a0,r
    write(*,*) 'What is the value of a0?'
    read(*,*) a0
    write(*,*) 'What is the value of r?'
    read(*,*) r
    write(*,*) 'How many iterations should I run?'
    read(*,*) iter
    open(10,file = 'geom_output.dat')
    do iter = 1,10
    write(10,*) iter,a0
    a0 = a0*r
    enddo
    close(10)
end program geometric