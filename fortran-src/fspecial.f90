!
! Wrapper for special functions
! WJH 25 Aug 2017
!
! For use with f2py or jupyter + "fortran magic"
!
! Example usage:
! 
! %%fortran --opt='-O3' --extra '-L/Users/will/Dropbox/Teaching/Estelar/Tarea-03-GreyAtmos -lspecial_functions.o' -vvv
!

module fspecial
  implicit none
  !# Exponential integrals (first and second)
  !# Based on from special_functions.f90 by Shanjie Zhang, Jianming Jin 
  !# http://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html
  external E1XA
  integer, parameter :: DP = kind(1.0d0), SP = kind(1.0)
contains
  function e1(x) result(rslt)
    ! First exponential integral
    ! Use E1XA routine from special_functions.f90 
    real, intent(in) :: x
    real :: rslt
    real(DP) :: drslt
    call E1XA(real(x, DP), drslt)
    rslt = real(drslt, SP)
  end function e1

  function e2(x) result(rslt)
    ! Second exponential integral
    ! Use recurrence relation from Hubeny & Mihalas, Ch 11, p 365
    real, intent(in) :: x
    real :: rslt
    rslt = exp(-x) - x*e1(x)
  end function e2
  
end module fspecial

