! File: plotfunction.f90
! Author: Adrian Garcia
! Date: 1/16/2023
! Purpose: Beginner Fortran 90 tutorial
program plotfunction
    ! Exercise 7-10
    implicit none
    integer :: i
    real :: f,x
    real, parameter :: xmin = 0.,xmax=10.,a=-2.
    open(10,file = 'myplot.dat')
    do i = 1,100
    x = xmin + xmax*(i-1.0)/(100.0-1.0)
    write(10,*) x,f(x,a)
    enddo
    close(10)
    call system('gnuplot -p data_plot.plt')
    contains
end program plotfunction
! function f(x,a)
!     implicit none
!     real :: f,x,a
!     f = cos(x+a)
! end function f