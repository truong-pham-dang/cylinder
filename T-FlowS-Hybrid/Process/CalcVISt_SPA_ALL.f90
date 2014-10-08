!======================================================================!
  SUBROUTINE CalcVISt_SPA_ALL(n) 
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
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, n 
  REAL    :: Xrat, Fv1, lf, Cs 
!--------------------------------[CVS]---------------------------------!
!  $Id: CalcVISt_SPA_ALL.f90,v 1.4 2008/12/10 14:15:46 IUS\mhadziabdic Exp $
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CalcVISt_SPA_ALL.f90,v $
!======================================================================!

!======================================================================!
!  SPA_ALL model
!======================================================================!
 
  if(SIMULA == DES_SPA) then
    do c = 1,NC
      Xrat     = VIS % n(c)/VISc
      Fv1      = Xrat**3.0/(Xrat**3.0 + Cvis1**3.0)
      VISt(c)  = DENc(material(c)) * Fv1 * VIS % n(c)
    end do
  end if

  if(SIMULA == SPA_ALL) then
    do c = 1,NC
      Xrat     = VIS % n(c)/VISc
      Fv1      = Xrat**3.0/(Xrat**3.0 + Cvis1**3.0)
      VISt(c)  = DENc(material(c)) * Fv1 * VIS % n(c)
    end do
  end if

  call Exchng(VISt)  

  RETURN

  END SUBROUTINE CalcVISt_SPA_ALL
