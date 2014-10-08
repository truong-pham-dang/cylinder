!======================================================================!
  SUBROUTINE CalcShear(Ui, Vi, Wi, She)
!----------------------------------------------------------------------!
!  Computes the magnitude of the shear stress.                         !
!----------------------------------------------------------------------!
!  She = sqrt( 2 * Sij * Sij )                                       !
!  Sij = 1/2 ( dUi/dXj + dUj/dXi )                                     !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  REAL          :: Ui(-NbC:NC), Vi(-NbC:NC), Wi(-NbC:NC)
  REAL          :: She(-NbC:NC)
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, i
  REAL    :: Sii, Sjk 
!--------------------------------[CVS]---------------------------------!
!  $Id: CalcShear.f90,v 1.9 2008/11/19 14:51:08 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CalcShear.f90,v $  
!======================================================================!

  call Exchng(Ui)
  call Exchng(Vi)
  call Exchng(Wi)

!===================!
!                   !
!     SGS terms     !
!                   !
!===================!
  call GraPhi(Ui,1,Ux, .TRUE.)  ! dU/dx
  call GraPhi(Ui,2,Uy, .TRUE.)  ! dU/dy
  call GraPhi(Ui,3,Uz, .TRUE.)  ! dU/dz
  call GraPhi(Vi,1,Vx, .TRUE.)  ! dV/dx
  call GraPhi(Vi,2,Vy, .TRUE.)  ! dV/dy
  call GraPhi(Vi,3,Vz, .TRUE.)  ! dV/dz
  call GraPhi(Wi,1,Wx, .TRUE.)  ! dW/dx
  call GraPhi(Wi,2,Wy, .TRUE.)  ! dW/dy
  call GraPhi(Wi,3,Wz, .TRUE.)  ! dW/dz

  do c=1,NC
    She(c) = Ux(c)*Ux(c) + Vy(c)*Vy(c) + Wz(c)*Wz(c) + &
             0.5*(Vz(c) + Wy(c))*(Vz(c) + Wy(c)) + & 
             0.5*(Uz(c) + Wx(c))*(Uz(c) + Wx(c)) + & 
             0.5*(Vx(c) + Uy(c))*(Vx(c) + Uy(c)) 

     Vort(c) =  -(0.5*(Vz(c) - Wy(c))*(Vz(c) - Wy(c)) + &
                  0.5*(Uz(c) - Wx(c))*(Uz(c) - Wx(c)) + &
                  0.5*(Vx(c) - Uy(c))*(Vx(c) - Uy(c)))

  end do 

  She = sqrt(2.0 * She)
  Vort = sqrt(2.0*abs(Vort))

  END SUBROUTINE CalcShear
