module fquad
  implicit none
  ! Wrapper for quadrature routines from QUADPACK
  external QAG, QAGI

  ! Input variables that control the QUADPACK routines:
  !
  ! Required precision, absolute and relative
  real :: QP_epsabs = 0.0, QP_epsrel = 0.001
  ! Order of the integration rule
  integer :: QP_key = 6
  
  contains
    function integrate(func, a, b) result(rslt)
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

end module fquad
