! File: myprogram.f90
! Author: Adrian Garcia
! Date: 1/16/2023
! Purpose: Beginner Fortran 90 tutorial
program myprogram
    ! Exercise 1
    implicit none
    integer :: i,j,k
    real :: x,y,z
    ! Exercise 4
    ! write(*,*) 'What is the value of x?'
    ! read(*,*) x
    ! write(*,*) 'What is the value of i?'
    ! read(*,*) i
    ! Exercise 5
    open(11,file = 'input.dat')
    read(11,*) x
    read(11,*) i
    close(11)
    y = cos(x)
    z = x + y
    i = 3
    j = i**2
    k = i-j
    ! Exercise 2
    write(*,*) 'i = ',i,',j = ',j,',k = ',k,',x = ',x,',y = ',y,', and z = ',z
    ! Exercise 3
    open(10,file = 'mydata.dat')
    write(10,*) 'i = ',i,',j = ',j,',k = ',k,',x = ',x,',y = ',y,', and z = ',z
    close(10)
end program myprogram