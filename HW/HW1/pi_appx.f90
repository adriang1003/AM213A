! File: pi_appx.f90
! Author: Adrian Garcia
! Date: 1/16/2023
! Purpose: Subroutine for Problem 2 from Homework 1 for AM213A
subroutine pi_appx(thres,appx,diff,N)
    real :: thres
    real(kind = kind(0.d0)) :: diff,true,appx
    integer :: N
    true  = acos(-1.d0)
    appx = 0.d0
    N = 0
    diff = abs(appx-true) 
    do while (diff >= thres)
        ! Check if threshold is met
        diff = abs(appx-true)
        ! Increase accuracy of appx by adding another term of the summation
        appx = appx + (16.d0**(-N) * (4.d0/(8.d0*N+1.d0) - 2.d0/(8.d0*N+4.d0) - 1.d0/(8.d0*N+5.d0) - 1.d0/(8.d0*N+6.d0)))
        ! Note what term we are on
        N = N + 1
    end do
end subroutine pi_appx