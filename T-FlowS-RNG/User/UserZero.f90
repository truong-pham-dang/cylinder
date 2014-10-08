!======================================================================!
  SUBROUTINE UserZero(var) 
!----------------------------------------------------------------------!
! Description:                                                         !
! ~~~~~~~~~~~~                                                         !
!   Adds bouyancy terms to the right hand side of velocity equation    !
!   for zero Pr cavity L x H = 4 x 1.                                  !
!                                                                      ! 
! Note:                                                                !
! ~~~~~                                                                !
!   It should be used with zero.* problem in Test directory.           !  
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
!  $Id: UserZero.f90,v 1.4 2002/10/30 16:30:04 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserZero.f90,v $    
!======================================================================!

  if(var == 3) then  ! only for W velocity component
    do c=1,NC
      b(c)=b(c) + 0.25*(xc(c)-2.0)*volume(c)
    end do
  end if

  END SUBROUTINE UserZero
