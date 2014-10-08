!======================================================================!
  SUBROUTINE CalcVISt_Keps() 
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
  REAL    :: CK, Yplus, Fmu, Ret, yStar                                        
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
  if(MODE == HRe) then
    do c=1,NC
      VISt(c) = Cmu * DENc(material(c)) * Kin%n(c) * Kin%n(c) / (Eps % n(c)+1.0e-14)
    end do
    do s=1,NS
      c1 = SideC(1,s)
      c2 = SideC(2,s)
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then  
        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
          Ck = sqrt(TauWall(c1))
          yPlus = DENc(material(c1))*Ck*WallDs(c1)/VISc 
          VISwall(c1) = yPlus*VISc*kappa/LOG(Elog*yPlus)
        end if
      end if
    end do
  end if   ! end if mode == stan

  if(MODE==LRe) then
    do c = 1, NC 
      Ret = Kin % n(c)*Kin % n(c)/(VISc*Eps % n(c))
      Fmu = exp(-2.5/(1.0 + Ret/50.0))
      VISt(c) = Fmu * Cmu * DENc(material(c)) * Kin%n(c) * Kin%n(c) / Eps % n(c)
    end do
  end if

  if(SIMULA==HYB_PITM) then
    do c = 1, NC

      Ret = Kin % n(c)*Kin % n(c)/(VISc*Eps % n(c))

      yStar = (VISc * Eps % n(c))**0.25 * WallDs(c)/VISc

      Fmu = (1.0 - exp(-yStar/14.0))**2.0*(1.0                              &
            + 5.0*exp(-(Ret/200.0)*(Ret/200.0))/Ret**0.75)


      Fmu = Fmu / ( 1.0 + exp(-yStar/5.0)**1.5/0.06 )
      Fmu = min(1.0,Fmu)

      VISt(c) = Fmu * Cmu * DENc(material(c)) * Kin%n(c) * Kin%n(c) / Eps % n(c)
    end do
  end if


  call Exchng(VISt)  
  RETURN

  END SUBROUTINE CalcVISt_Keps 
