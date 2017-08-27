!
! Test of fgrey set of modules, using normal fortran only (no
! f2py or jupyter)
!

program fgtest
  use planck_grey, only: p, planck
  use milne_flux_integral, only: flux
  implicit none
  real, allocatable, dimension(:) :: tau, alpha
  real, allocatable, dimension(:, :) :: ttau, aalpha

  ! Use f2003 syntax: square brackets and no need for explicit ALLOCATE
  tau = [0.0, 1.0, 2.0, 4.0, 8.0, 16.0]
  alpha = [1.0, 3.0, 5.5, 9.0, 15.0]

  print *, 'Test of inverse T function for vector argument'
  print '(a9, 6(f6.2))', 'tau = ', tau
  print '(a9, 6(f6.3))', 'p(tau) = ', p(tau)
  print *

  print *, 'Test of Planck function for vector tau and scalar alpha'
  print '(a18, 6(f6.2))', 'tau = ', tau
  print '(a18, 6(f6.3))', 'planck(3, tau) = ', planck(3.0, tau)
  print *
  
  print *, 'Test of Planck function for scalar tau and vector alpha'
  print '(a19, 5(f6.2))', 'alpha = ', alpha
  print '(a19, 5(f6.3))', 'planck(alpha, 1) = ', planck(alpha, 1.0)
  print *

  ! Make 2-d versions of alpha (x-axis, across columns) and tau
  ! (y-axis, down rows)
  aalpha = spread(alpha, dim=2, ncopies=size(tau))
  ttau = spread(tau, dim=1, ncopies=size(alpha))

  print *, 'Test of Planck function for 2-d tau and 2-d alpha'
  print '(a)', 'alpha = '
  print '(5(f6.2))', aalpha
  print '(a)', 'tau = '
  print '(5(f6.2))', ttau
  print '(a)', 'planck(alpha, tau) = '
  print '(5(f6.3))', planck(aalpha, ttau)
  print *

  print *, 'Test of flux integral for scalar args'
  print '(a23, 5(f6.3))', 'flux(alpha=3, tau=1) =', flux(3.0, 1.0)
  print *

  print *, 'Test of flux integral for vector tau and scalar alpha'
  print '(a18, 6(f6.2))', 'tau = ', tau
  print '(a18, 6(f6.3))', 'flux(5, tau) = ', flux(5.0, tau)
  print *
  
  print *, 'Test of flux integral for scalar tau and vector alpha'
  print '(a19, 5(f6.2))', 'alpha = ', alpha
  print '(a19, 5(f6.3))', 'flux(alpha, 1) = ', flux(alpha, 1.0)
  print *

  print *, 'Test of flux integral for vector tau and vector alpha'
  print '(a)', 'alpha = '
  print '(5(f6.2))', alpha
  print '(a)', 'tau = '
  print '(1(f6.2))', tau
  print '(a)', 'flux(alpha, tau) = '
  print '(5(f6.3))', flux(alpha, tau)
  print *
  
end program fgtest

