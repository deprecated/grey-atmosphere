module fquad
  ! Wrapper for quadrature routines from QUADPACK
  use qp_interface, only: QAG, QAGI
  implicit none

  ! Input variables that control the QUADPACK routines:
  !
  ! Required precision, absolute and relative
  real :: QP_epsabs = 0.0, QP_epsrel = 0.001
  ! Order of the integration rule
  integer :: QP_key = 6
  ! Type of semi-infinite integral for QAGI (1 is a->infty, -1 is
  ! -infty->a, 2 is -infty->infty)
  integer :: QP_inf = 1
  
  contains
    function integrate(func, a, b) result(rslt)
      ! Integrate func between limits a -> b
      real, external :: func
      real, intent(in) :: a, b
      real :: rslt
      ! Output arguments from QAG
      real :: QP_result, QP_abserr
      integer :: QP_neval, QP_ier
      call QAG(func, a, b, QP_epsabs, QP_epsrel, QP_key, &
           & QP_result, QP_abserr, QP_neval, QP_ier)
      if (QP_ier == 0) then
         rslt = QP_result
      else
         print *, 'Error in QAG: ', QP_ier, QP_neval, QP_abserr
         rslt = 0.0
      end if
    end function integrate

    function integrate_infinity(func, a) result(rslt)
      ! Integrate func between limits a -> infty
      real, external :: func
      real, intent(in) :: a
      real :: rslt
      ! Output arguments from QAGI
      real :: QP_result, QP_abserr
      integer :: QP_neval, QP_ier
      call QAGI(func, a, QP_inf, QP_epsabs, QP_epsrel, &
           & QP_result, QP_abserr, QP_neval, QP_ier)
      if (QP_ier == 0) then
         rslt = QP_result
      else
         print *, 'Error in QAGI: ', QP_ier, QP_neval, QP_abserr
         rslt = 0.0
      end if
    end function integrate_infinity

end module fquad
