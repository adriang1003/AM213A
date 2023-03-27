! File: solbybisection.f90
! Author: Adrian Garcia
! Date: 1/16/2023
! Purpose: Driver file for Exercise 16 from Homework 1 for AM213A
program solbybisection
    implicit none
    real, parameter :: xmin = 0,xmax = 3.141,xmin1 = 2,xmax1 = 4
    integer, parameter :: iter = 100
    real :: sol,err,fcosx,ftanx
    external fcosx,ftanx
    ! Function call for cos(x) on interval [0,pi]
    call bisect(xmin,xmax,fcosx,sol,iter,err)
    write(*,*) sol,err
    ! Function call for tan(x) on interval [2,4]
    call bisect(xmin1,xmax1,ftanx,sol,iter,err)
    write(*,*) sol,err
end program solbybisection