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
  
end module qintegrand

program fqtest
  use fquad, only: integrate
  use qintegrand, only: integrand, exact, PAR_n
  implicit none
  real :: a = 1.0, b = 5.0
  real :: integral

  do PAR_n = 0, 5
     integral = integrate(integrand, a, b)
     print '(a,i1, 2(a,f8.2))', &
          & 'n = ', PAR_n, ', integral =', integral, ', exact = ', exact(a, b) 
  end do
end program fqtest

