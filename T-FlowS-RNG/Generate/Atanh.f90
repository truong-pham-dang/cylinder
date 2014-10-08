!======================================================================!
  REAL FUNCTION Atanh(x)
!----------------------------------------------------------------------!
!   Calculates inverse of hyperbolic tangens.                          !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  REAL :: x
!--------------------------------[CVS]---------------------------------!
!  $Id: Atanh.f90,v 1.6 2002/10/30 16:29:20 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Generate/Atanh.f90,v $  
!======================================================================!
  if(x  > 1.0) then
    write(*,*) 'Error message from atanh: bad argument'
    stop
  end if 

  atanh=log( sqrt( (1.0+x)/(1.0-x) ) )

  END FUNCTION Atanh
