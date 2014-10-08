!======================================================================!
  SUBROUTINE IniDiv 
!----------------------------------------------------------------------!
!   Initialize processor numbers to 1                                  !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE par_mod 
  USE div_mod 
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c
!--------------------------------[CVS]---------------------------------!
!  $Id: IniDiv.f90,v 1.6 2002/10/30 16:29:18 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Divide/IniDiv.f90,v $   
!======================================================================!

!------------------------------!
!     Initialize parameters    !
!------------------------------!
  COORDINATE = 1
  INERTIAL   = 2

!------------------------------!
!     Initialize processors    !
!------------------------------!
  do c=1,NC
    proces(c)=1
  end do

  END SUBROUTINE IniDiv