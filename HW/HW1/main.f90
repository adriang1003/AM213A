! File: main.f90
! Author: Adrian Garcia
! Date: 1/16/2023
! Purpose: Driver file for Problem 2 from Homework 1 for AM213A
program main
    implicit none
    real :: thres
    real(kind = kind(0.d0)) :: diff,appx
    integer :: N
    thres = 1.e-16 ! = 1.e-4, 1.e-8, 1.e-12, 1.e-16
    ! Function call
    call pi_appx(thres,appx,diff,N)
    write(*,*) N,appx,diff
end program main