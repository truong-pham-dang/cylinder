!======================================================================!
  SUBROUTINE CalcVISt_EBM() 
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
  REAL    :: CK, Yplus, Fmu, Ret, yStar, Cmu_mod                                        
!--------------------------------[CVS]---------------------------------!
! '$Id: CalcVISt_KEps.f90,v 1.4 2009/06/30 11:56:12 IUS\mhadziabdic Exp $'
! '$Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CalcVISt_KEps.f90,v $'
!======================================================================!

!======================================================================!
!  K-EPS model
!======================================================================!
  Ret   = 0.0
  Fmu   = 0.0 
  yPlus = 0.0 
    call CalcShear(U % n, V % n, W % n, Shear)
    do c=1,NC

      Kin % n(c) = 0.5*max(uu%n(c)+vv%n(c)+ww%n(c),1.0e-12)

      Cmu_mod = max(-(uu%n(c)*Ux(c)+vv%n(c)*Vy(c)+ww%n(c)*Wz(c)+&
                  uv%n(c)*(Vx(c)+Uy(c))+uw%n(c)*(Uz(c)+&
                  Wx(c))+vw%n(c)*(Vz(c)+Wy(c)))/max(Kin%n(c)*Kin%n(c)/Eps%n(c)*Shear(c)*Shear(c),1.0e-12),0.0)

      Cmu_mod = min(0.09,Cmu_mod) 
      if(MODE==HYB) then
        VISt(c) = Cmu_mod*DENc(material(c)) * Kin%n(c) * Kin%n(c) / Eps % n(c)
      else
        VISt(c) = Cmu*DENc(material(c)) * Kin%n(c) * Kin%n(c) / Eps % n(c)
      end if         
    end do

  call Exchng(VISt)  
  RETURN

  END SUBROUTINE CalcVISt_EBM 
