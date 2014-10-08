!======================================================================!
  SUBROUTINE CalcSGS_hybrid()
!----------------------------------------------------------------------!
!   Calculates SGS stresses and turbulent viscosity for LES.           !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE rans_mod
!----------------------------------------------------------------------!
! near(c) is the number of corresponding cell on the nearest wall.
! In case that, in parallel executions, the subdomain does not have 
! any nearwall cells, the near(c) is zero.
! near(c) is calculated in NearWallCells.f90, only ones in the beginig
! of a simulation.
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, s, c1, c2 
  REAL    :: Nx, Ny, Nz
  REAL    :: Cs, R
  REAL    :: Stot, lf, UtauL, Uff 
  REAL    :: Utot, Unor, Utan, Apow, Bpow, nu, dely, yPlus 
  REAL    :: fun
!--------------------------------[CVS]---------------------------------!
!  $Id: CalcSGS.f90,v 1.11 2008/12/10 14:09:43 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CalcSGS.f90,v $  
!======================================================================!

  
!===================!
!                   !
!     SGS terms     !
!                   !
!===================!

    do c=1,NC
      lf =  (volume(c)**0.333333)    
      VISt_sgs(c) = DENc(material(c))    &
                    * (lf*lf)              &          ! delta^2 
                    * Cdyn(c)              &          ! Cdynamic   
                    * Shear(c)      
    end do

    call Exchng(VISt_sgs)

  END SUBROUTINE CalcSGS_hybrid
