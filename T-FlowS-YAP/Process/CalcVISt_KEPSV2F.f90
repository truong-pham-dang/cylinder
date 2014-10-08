!======================================================================!
  SUBROUTINE CalcVISt_KepsV2F() 
!----------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                  !
!                                                                      !
!   In the domain:                                                     !
!   ~~~~~~~~~~~~~~                                                     !
!   For k-eps model :                                                  !
!                       2                                              !
!   VISt = Cmu * rho * K  * Eps                                        ! 
!                                                                      !
!   On the boundary (wall viscosity):                                  !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                   !
!            +          kappa                                          !
!   VIStw = y  * VISt ----------                                       ! 
!                     E * ln(y+)                                       !
!                                                                      !
!    For k-eps-v2f model :                                             !
!                                                                      !
!    VISt = CmuD * rho * Tsc  * vv                                     !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------[Calling]-------------------------------!
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, c1, c2, s
  REAL    :: Utan, UnorSq, Unor, UtotSq, Cmu1
  REAL    :: lf, Gblend, Ustar, Ck, yPlus, Uplus, EBF

!--------------------------------[CVS]---------------------------------!
!  $Id: CalcVISt_KEPSV2F.f90,v 1.4 2008/12/10 14:14:53 IUS\mhadziabdic Exp $
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CalcVISt_KEPSV2F.f90,v $
!======================================================================!

!======================================================================!
!  K-EPS-V2F model
!======================================================================! 
  call Scale()

  if(SIMULA == K_EPS_VV) then
    do c = 1,NC
      VISt(c) = CmuD*v_2%n(c)*Tsc(c)
    end do
  else if(SIMULA == ZETA) then
    do c = 1,NC
      VISt(c) = CmuD*v_2%n(c)*Kin % n(c)*Tsc(c)
    end do
  else if(SIMULA == HYB_ZETA) then
    do c = 1,NC
      VISt(c) = CmuD*v_2%n(c)*Kin % n(c)*Tsc(c)
      VISt_eff(c) = max(VISt(c),VISt_sgs(c))
    end do
  call Exchng(VISt_eff)  
  end if

    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
    
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
!----- Compute tangential velocity component
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


          Uf(c1)  = Cmu**0.25*Kin%n(c1)**0.5
          Ynd(c1) = WallDs(c1)*Uf(c1)/VISc 
          Ck = sqrt(TauWall(c1))
          Ynd(c1) = WallDs(c1)*Ck/VISc
          Gblend  = 0.01*Ynd(c1)**4.0/(1.0+5.0*Ynd(c1))

          Yplus = max(Ynd(c1),1.1)
          Uplus = log(Yplus*Elog)/(kappa)
          
          if(Ck*WallDs(c1)/VISc< 5.0) then
            VISwall(c1) = VISt(c1) + VISc  
          else
            VISwall(c1) = Ynd(c1)*VISc/(Yplus*exp(-1.0*Gblend)+Uplus*exp(-1.0/Gblend) + TINY)
          end if

        end if  ! TypeBC(c2)==WALL or WALLFL
      end if    ! c2 < 0
    end do
!  end if  
  RETURN
 
  call Exchng(VISt)  
  RETURN

  END SUBROUTINE CalcVISt_KepsV2F  
