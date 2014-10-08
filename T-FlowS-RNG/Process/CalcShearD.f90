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
  do i=1,3        

    if(i == 1) then
      call GraPhi(Ui,1,PHIx, .TRUE.)  ! dU/dx
      call GraPhi(Wi,2,PHIy, .TRUE.)  ! dW/dy
      call GraPhi(Vi,3,PHIz, .TRUE.)  ! dV/dz
    end if
    if(i == 2) then
      call GraPhi(Wi,1,PHIx, .TRUE.)  ! dW/dx
      call GraPhi(Vi,2,PHIy, .TRUE.)  ! dV/dy
      call GraPhi(Ui,3,PHIz, .TRUE.)  ! dU/dz
    end if
    if(i == 3) then
      call GraPhi(Vi,1,PHIx, .TRUE.)  ! dV/dx
      call GraPhi(Ui,2,PHIy, .TRUE.)  ! dU/dy
      call GraPhi(Wi,3,PHIz, .TRUE.)  ! dW/dz
    end if

    do c=1,NC
      if(i == 1) then
	Sii = PHIx(c)               ! Sxx:           dU/dx     
	Sjk = 0.5*(PHIy(c)+PHIz(c)) ! Syz=Szy:  .5(dV/dz+dW/dy)
	She(c) = ( Sii*Sii + 2.0*Sjk*Sjk )  
      end if
      if(i == 2) then
	Sii = PHIy(c)               ! Syy:           dV/dy     
	Sjk = 0.5*(PHIx(c)+PHIz(c)) ! Sxz=Szx:  .5(dU/dz+dW/dx)
	She(c) = She(c) + ( Sii*Sii + 2.0*Sjk*Sjk ) 
      end if
      if(i == 3) then
	Sii = PHIz(c)               ! Szz:           dW/dz     
	Sjk = 0.5*(PHIx(c)+PHIy(c)) ! Sxy=Syx:  .5(dV/dx+dU/dy)
	She(c) = She(c) + ( Sii*Sii + 2.0*Sjk*Sjk )
      end if
    end do 

  end do  ! i

  She = sqrt(2.0 * She)

  END SUBROUTINE CalcShear
