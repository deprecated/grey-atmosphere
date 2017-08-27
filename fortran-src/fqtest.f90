!
! Test of fquad module, using normal fortran only (no
! f2py or jupyter)
!

module qintegrand
  ! Function and context for use with integration routine
  implicit none
  integer :: PAR_n = 1
contains
  function integrand(x) result(rslt)
    ! QUADPACK requires function of one argument to be integrated, so
    ! any additional arguments must be passed via module variables
    ! (for example, PAR_n)
    real, intent(in) :: x
    real :: rslt
    rslt = x**PAR_n
  end function integrand

  function exact(a, b) result(rslt)
    ! Exact solution to integral of integrand(x) between x=a and x=b
    real, intent(in) :: a, b
    real :: rslt
    integer :: m
    m = PAR_n + 1
    if (m == 0) then
       rslt = log(b/a)
    else
       rslt = (b**m - a**m) / real(m)
    end if
  end function exact

  function exact_infinity(a) result(rslt)
    ! Exact solution to integral of integrand(x) between x=a and x=infty
    real, intent(in) :: a
    real :: rslt
    integer :: m
    m = PAR_n + 1
    if (m >= 0) then
       print *, 'ERROR: Integral diverges with n = ', PAR_n
       rslt = 0.0
    else
       rslt = a**m / real(-m)
    end if
  end function exact_infinity

  
end module qintegrand

program fqtest
  use fquad, only: integrate, integrate_infinity
  use qintegrand, only: integrand, exact, exact_infinity, PAR_n
  implicit none
  real :: a = 0.1, b = 5.0
  real :: integral

  print '(2(a,f5.2))', 'Testing integration of x**n from a = ', a, ' to b = ', b
  do PAR_n = -2, 5
     integral = integrate(integrand, a, b)
     print '(a,i2.1, 2(a,f8.2))', &
          & 'n = ', PAR_n, ', integral =', integral, ', exact = ', exact(a, b) 
  end do
  print *

  print '(a,f5.2,a)', 'Testing integration of x**n from a = ', a, ' to b = infty'
  do PAR_n = -6, -2
     integral = integrate_infinity(integrand, a)
     print '(a,i2.1, 2(a,f8.2))', &
          & 'n = ', PAR_n, ', integral =', integral, ', exact = ', exact_infinity(a) 
  end do
  print *
  
end program fqtest

