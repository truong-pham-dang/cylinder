!======================================================================!
  SUBROUTINE CalcWALE()
!----------------------------------------------------------------------!
!  
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  REAL    :: SijdSijd(-NbC:NC), She(-NbC:NC), Vor(-NbC:NC)
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, i
  REAL    :: S11, S22, S33, S12, S13, S23, S21, S31, S32
  REAL    :: S11d, S22d, S33d, S12d, S13d, S23d, S21d, S31d, S32d
  REAL    :: V11, V22, V33, V12, V13, V23, V21, V31, V32
 
!--------------------------------[CVS]---------------------------------!
!  $Id: CalcShear.f90,v 1.9 2008/11/19 14:51:08 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CalcShear.f90,v $  
!======================================================================!

!===================!
!                   !
!     SGS terms     !
!                   !
!===================!
  do c=1,NC
    S11 = Ux(c)
    S22 = Vy(c)
    S33 = Wz(c)
    S12 = 0.5*(Vx(c) + Uy(c))
    S13 = 0.5*(Uz(c) + Wx(c))
    S23 = 0.5*(Vy(c) + Wy(c))
    S21 = S12
    S31 = S13
    S32 = S23

    V11 = 0.0   
    V22 = 0.0  
    V33 = 0.0  
    V12 = 0.5*(Vx(c) - Uy(c))
    V13 = 0.5*(Uz(c) - Wx(c))
    V23 = 0.5*(Vy(c) - Wy(c))
    V21 = -V12
    V31 = -V13
    V32 = -V23

    She(c) = 0.5 * Shear(c) * Shear(c)
    Vor(c) = 0.5 * Vort(c) * Vort(c)

    S11d = S11*S11 + S12*S12 + S13*S13 - (V11*V11 + V12*V12 + V13*V13) - 1.0/3.0 * (She(c) - Vor(c))
    S22d = S12*S12 + S22*S22 + S23*S23 - (V12*V12 + V22*V22 + V23*V23) - 1.0/3.0 * (She(c) - Vor(c))
    S33d = S13*S13 + S23*S23 + S33*S33 - (V13*V13 + V23*V23 + V33*V33) - 1.0/3.0 * (She(c) - Vor(c))

    S12d = S11*S12 + S12*S22 + S13*S32 + (V11*V12 + V12*V22 + V13*V32) 
    S13d = S11*S13 + S12*S23 + S13*S33 + (V11*V13 + V12*V23 + V13*V33) 
    S23d = S21*S13 + S22*S23 + S23*S33 + (V21*V13 + V22*V23 + V23*V33) 
    S21d = S12d
    S31d = S13d
    S32d = S23d

    SijdSijd(c) = S11d*S11d + S22d*S22d + S33d*S33d + S12d*S12d + S13d*S13d + S23d*S23d
    
    WALEv(c) = sqrt(abs(SijdSijd(c)*SijdSijd(c)*SijdSijd(c)))/(sqrt(abs(She(c)*She(c)*She(c)*She(c)*She(c))) + &
               sqrt(sqrt(abs(SijdSijd(c)*SijdSijd(c)*SijdSijd(c)*SijdSijd(c)*SijdSijd(c)))) + tiny)
  end do 

  END SUBROUTINE CalcWALE
