! File: EX13.f90
! Author: Adrian Garcia
! Date: 1/30/2023
! Purpose: Exercise 13 from Homework 2 for AM213A
program EX13
    implicit none
    real :: i
    real(kind = kind(0.d0)) :: f,g,x
    real(kind = kind(0.d0)), parameter :: xmin = 1.920,xmax = 2.080
    open(10,file = 'myplot.dat')
    i = 0.001
    do while (x <= xmax)
    x = xmin + i
    i = i + 0.001
    write(10,*) x,f(x),g(x)
    enddo
    close(10)
    call system('gnuplot -p data_plot.plt')
    contains
end program EX13
function f(x)
    implicit none
    real(kind = kind(0.d0)) :: f,x
    f = (x - 2)**9
end function f
function g(x)
    implicit none
    real(kind = kind(0.d0)) :: g,x
    g = x**9 - 18*x**8 + 144*x**7 - 672*x**6 + 2016*x**5 - 4032*x**4 + 5376*x**3 - 4608*x**2 + 2304*x - 512
end function g