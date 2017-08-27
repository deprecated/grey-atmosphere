module qp_interface
  implicit none
  !
  ! Idealised interfaces for the QUADPACK routines that I actually use
  !
  ! WJH 26 Aug 2017
  !
  ! The full source for these routines is in quadpack.f90, from which
  ! I have extracted these signatures, but with 2 improvements:
  !
  ! 1. I have added intent(in/out) qualifiers for the arguments
  !
  ! 2. I have added a nested interface for the user-supplied function
  !    "f", following the description in the comment block, although
  !    this is actually declared as "external" in the source
  !
  interface
     subroutine QAG ( f, a, b, epsabs, epsrel, key, result, abserr, neval, ier )
       implicit none
       interface
          function f(x) result(rslt)
            real ( kind = 4 ), intent(in) :: x
            real ( kind = 4 ) :: rslt
          end function f
       end interface
       ! real ( kind = 4 ), external :: f
       real ( kind = 4 ), intent(in) :: a, b, epsabs, epsrel
       integer ( kind = 4 ), intent(in) :: key
       real ( kind = 4 ), intent(out) :: result, abserr
       integer ( kind = 4 ), intent(out) :: neval, ier
     end subroutine QAG

     subroutine QAGI ( f, bound, inf, epsabs, epsrel, result, abserr, neval, ier )
       implicit none
       interface 
          function f(x) result(rslt)
            real ( kind = 4 ), intent(in) :: x
            real ( kind = 4 ) :: rslt
          end function f
       end interface
       ! real ( kind = 4 ), external :: f
       real ( kind = 4 ), intent(in) :: bound, epsabs, epsrel
       integer ( kind = 4 ), intent(in) :: inf
       real ( kind = 4 ), intent(out) :: result, abserr
       integer ( kind = 4 ), intent(out) :: neval, ier
     end subroutine QAGI
  end interface
  
end module qp_interface
