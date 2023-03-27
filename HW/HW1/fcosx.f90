! File: fcosx.f90
! Author: Adrian Garcia
! Date: 1/16/2023
! Purpose: Cosine function for Exercise 16 from Homework 1 for AM213A
function fcosx(x,a)
    implicit none
    real :: fcosx,x,a
    fcosx = cos(x+a)
end function fcosx