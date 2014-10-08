!======================================================================!
  SUBROUTINE UserForce(var) 
!----------------------------------------------------------------------!
!   Discretizes and solves momentum conservation equations             !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: var ! 1 -> U,  2 -> V,  3 -> W 
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c
!--------------------------------[CVS]---------------------------------!
!  $Id: UserForce.f90,v 1.9 2002/10/30 16:30:02 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserForce.f90,v $    
!----------------------------------------------------------------------!
! Description:                                                         !
! ~~~~~~~~~~~~                                                         !
!   Adds bouyancy terms to the right hand side of velocity equation.   !
!                                                                      !
! Note:                                                                !
! ~~~~~                                                                !
!   Relies on two assumtions:                                          !
!   1. gravitational constant is assumed to be 1                       !
!   2. refference temperature is assumed to be 0                       !
!======================================================================!

  return

  if(var == 3) then  ! only for W velocity component
    do c=1,NC
      b(c)=b(c) + T % n(c)*volume(c)
    end do
  end if

  END SUBROUTINE UserForce
