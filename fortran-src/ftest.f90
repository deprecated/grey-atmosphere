!
! Vanilla test of fspecial module, using normal fortran only (no
! f2py or jupyter)
!
program ftest
  use fspecial, only: e1, e2
  implicit none
  integer, parameter :: N = 101
  real, parameter :: x0 = 0.0001, x1 = 5.0001, dx = (x1 - x0)/real(N-1)
  integer :: i
  real, dimension(N) :: x = (/(real(i-1)*dx + x0, i = 1, N)/)
  real, dimension(N) :: y
  do i = 1, N
     y = (/(e2(x(i)), i = 1, N)/)
  end do

  print '(10(f6.3))', x
  print '(10(f6.3))', y
end program ftest
