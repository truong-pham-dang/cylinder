!======================================================================!
  LOGICAL FUNCTION Approx(A,B,tol)
!----------------------------------------------------------------------!
!   Returns true if A~B, false otherwise.                              !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  REAL          :: A,B 
  REAL,OPTIONAL :: tol
!-------------------------------[Locals]-------------------------------!
  REAL :: tolerance
!--------------------------------[CVS]---------------------------------!
!  $Id: Approx.f90,v 1.4 2002/10/30 16:29:31 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Library/Approx.f90,v $  
!======================================================================!

  if( .not. present(tol) ) then
    tolerance = 1.e-6
  else
    tolerance = tol
  end if

  if( (A  < (B + tolerance)) .and. (A  > (B - tolerance)) ) then
    approx = .TRUE.
  else
    approx = .FALSE.
  end if

  END FUNCTION Approx
