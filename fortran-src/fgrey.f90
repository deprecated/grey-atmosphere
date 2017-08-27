module planck_grey
  ! Planck spectrum as function of dimensionless frequency
  ! (alpha = h nu / k T_eff) and optical depth (tau) for a
  ! LTE grey model atmosphere in the Eddington approximation. 
  implicit none
  real :: bigC = 0.15398972357101687
contains
  elemental function p(tau) result(rslt)
    real, intent(in) :: tau
    real :: rslt
    rslt = (4.0 / (3.0*tau + 2.0))**0.25
  end function p

  elemental function planck(alpha, tau) result(rslt)
    real, intent(in) :: alpha, tau
    real :: rslt
    rslt = bigC * alpha**3 / (exp(alpha*p(tau)) - 1.0)
  end function planck
    
!!$  function planck_11(alpha, tau) result(rslt)
!!$    real, intent(in) :: alpha(:), tau(:)
!!$    real :: rslt(size(alpha), size(tau))
!!$    real :: alpha_11(size(alpha), size(tau)), p_11(size(alpha), size(tau))
!!$    alpha_11 = spread(alpha, dim=2, ncopies=size(tau))
!!$    p_11 = spread(p(tau), dim=1, ncopies=size(alpha))
!!$    rslt = bigC * alpha_11**3 / (exp(alpha_11*p_11) - 1.0)
!!$  end function planck_11

end module planck_grey

module milne_context
  ! Integrand function and associated data for performing the Milne
  ! integral for the flux
  real :: PAR_tau, PAR_alpha
contains
  function milne_integrand(t) result(rslt)
    ! The integrand must be a function of 1 argument to meet the
    ! signature requirements of QUADPACK.  The extra arguments are
    ! sneaked in as module variables (see the PAR_* defined above)
    use fspecial, only: e2
    use planck_grey, only: planck
    real, intent(in) :: t
    real :: rslt
    rslt = 2.0 * e2(abs(t - PAR_tau)) * planck(PAR_alpha, t)
  end function milne_integrand
end module milne_context


module milne_flux_integral
  use fquad, only: integrate, integrate_infinity
  use milne_context, only: milne_integrand, PAR_tau, PAR_alpha
  implicit none
  interface flux
     ! Last 3 digits specify respective ranks of alpha, tau, result
     module procedure flux000, flux101, flux011, flux112
  end interface flux
contains
  function downward(alpha, tau) result(rslt)
    ! Integrate from 0->tau to find downward stream
    real, intent(in) :: alpha, tau
    real :: rslt
    PAR_tau = tau
    PAR_alpha = alpha
    rslt = integrate(milne_integrand, 0.0, tau)
  end function downward

  function upward(alpha, tau) result(rslt)
    ! Integrate from tau->infty to find upward stream
    real, intent(in) :: alpha, tau
    real :: rslt
    PAR_tau = tau
    PAR_alpha = alpha
    rslt = integrate_infinity(milne_integrand, tau)
  end function upward

  function flux000(alpha, tau) result(rslt)
    ! Flux is difference between upward and downward streams
    real, intent(in) :: alpha, tau
    real :: rslt
    rslt = upward(alpha, tau)
    if (tau > 0.0) then
       rslt = rslt - downward(alpha, tau)
    end if
  end function flux000

  function flux101(alpha, tau) result(rslt)
    ! (vector alpha, scalar tau) -> (vector flux)
    real, intent(in) :: alpha(:), tau
    real :: rslt(size(alpha))
    integer :: i
    rslt = [(flux(alpha(i), tau), i = 1, size(alpha))]
  end function flux101

  function flux011(alpha, tau) result(rslt)
    ! (scalar alpha, vector tau) -> (vector flux)
    real, intent(in) :: alpha, tau(:)
    real :: rslt(size(tau))
    integer :: i
    rslt = [(flux(alpha, tau(i)), i = 1, size(tau))]
  end function flux011

  function flux112(alpha, tau) result(rslt)
    ! (vector alpha, vector tau) -> (rank-2 flux)
    real, intent(in) :: alpha(:), tau(:)
    real :: rslt(size(alpha), size(tau))
    integer :: i, j
    ! Array literals can only be rank-1, so we construct the rank-2
    ! arrary by (1) making a flat rank-1 array by repeating the
    ! flux101 function ntau times, and (2) reshaping it to the desired
    ! [nalpha, ntau] shape
    rslt = reshape( source = [(flux(alpha, tau(j)), j = 1, size(tau))], &
         &           shape = [size(alpha), size(tau)] )
  end function flux112
  
end module milne_flux_integral
