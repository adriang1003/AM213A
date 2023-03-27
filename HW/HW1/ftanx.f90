! File: fcosx.f90
! Author: Adrian Garcia
! Date: 1/16/2023
! Purpose: Tangent function for Exercise 16 from Homework 1 for AM213A
function ftanx(x)
    implicit none
    real :: ftanx,x
    ftanx = tan(x)
end function ftanx