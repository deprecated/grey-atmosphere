subroutine fm_e2(t, rslt)
  use fspecial, only: e2
  real, intent(in) :: t(:)
  real, intent(out) :: rslt(size(t))
  integer :: i
  rslt = (/(e2(t(i)), i = 1, size(t))/)
end subroutine

subroutine fm_p(tau, rslt)
  use planck_grey, only: p
  real, intent(in) :: tau(:)
  real, intent(out) :: rslt(size(tau))
  rslt = p(tau)
end subroutine

subroutine fm_planck(alpha, tau, rslt)
  use planck_grey, only: planck
  real, intent(in) :: alpha(:), tau(:)
  real, intent(out) :: rslt(size(alpha), size(tau))
  real :: aalpha(size(alpha), size(tau)), ttau(size(alpha), size(tau))
  aalpha = spread(alpha, dim=2, ncopies=size(tau))
  ttau = spread(tau, dim=1, ncopies=size(alpha))
  rslt = planck(aalpha, ttau)
end subroutine
