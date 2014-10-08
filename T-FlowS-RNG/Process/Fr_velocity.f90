!======================================================================!
  SUBROUTINE Fr_velocity()
!----------------------------------------------------------------------!
!  Purpose:                                                            !
!  Calculate Yplus in the near wall cells in order to perform swiching !
!  from  wall function approach to up to the wall approach             !
!                                                                      !
!  Authors: Muhamed Hadziabdic and Bojan Niceno                        !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER c1,c2, s 
  REAL :: UtotSq, Unor, UnorSq, Utan
!--------------------------------[CVS]---------------------------------!
!  $Id: Fr_velocity.f90,v 1.7 2008/12/10 15:32:46 IUS\mhadziabdic Exp $
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/Fr_velocity.f90,v $
!======================================================================!

  do s=1,NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)

    if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
      if(TypeBC(c2) == WALL .or. TypeBC(c2) == WALLFL) then
        UtotSq = U % n(c1) * U % n(c1) &
               + V % n(c1) * V % n(c1) &
               + W % n(c1) * W % n(c1)
        Unor = ( U % n(c1) * Sx(s)     &
               + V % n(c1) * Sy(s)     &
               + W % n(c1) * Sz(s) )   &
               /(Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s))**0.5
        UnorSq = Unor*Unor
  
        if( UtotSq   >  UnorSq) then
          Utan = sqrt(UtotSq - UnorSq)
        else
          Utan = TINY
        end if        

        Uf(c1)  = (CmuD*Kin%n(c1)*v_2%n(c1))**0.25  
        Ynd(c1) = Uf(c1)*WallDs(c1)/VISc
      end if
    end if
  end do 

  call Exchng(VISt)

  RETURN 
  END SUBROUTINE Fr_velocity
