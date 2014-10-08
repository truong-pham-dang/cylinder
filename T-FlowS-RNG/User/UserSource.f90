!======================================================================!
  SUBROUTINE UserSource() 
!----------------------------------------------------------------------!
!   This subroutine extracts the heat in case when periodic boundaries !
!   are used in order to keep energy balans                            !
!   This source is part of the convection                              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c
!--------------------------------[CVS]---------------------------------!
!  $Id: UserSource.f90,v 1.5 2002/11/25 10:32:50 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserSource.f90,v $    
!----------------------------------------------------------------------!
! Description:                                                         !
! ~~~~~~~~~~~~                                                         !
!   Adds source to the energy equation.                                !
!======================================================================!
! 
!
!  /
! |
! |dT/dx * Ux * volume 
! |
!/
! For channel flow: dT/dx = 2.0*Q*B/FLUXoX  
! where Q is heat flux through the wall, B is
! channel width and FLUXoX is mean mass flux in  
! x direction of the channel, 2 takes into accound
! upper and lower channel wall
!
! For pipe flow: dT/dz = D*pi*Q/FLUXoZ
!  
!  b(c)=b(c) - Qflux * W % n(c) / ( Wbulk(material(c)) + TINY ) * volume(c)
! 
!  In order to be consistent we did not use ideal d*pi but Area/L where
!  Area is total wall surface of pipe and L is lenght of pipe
!  AreaZ(1) is surface oposite to the flow stream
!  Qflux is calculated in CalBou

!return


  do c=1,NC
!    b(c)=b(c) -   2.0*3.1415926*(0.1/12) * W % n(c) / (FLUXoZ(material(c))) * volume(c)
!    b(c)=b(c) -   1.0*0.1 * U % n(c) / (FLUXoX(material(c))) * volume(c)
    b(c)=b(c) -  2.0*3.1415926*0.005*W % n(c) / ( FLUXoZ(1) + TINY ) * volume(c)
!    b(c)=b(c) + 0.02*volume(c)
  end do

  return

  END SUBROUTINE UserSource
