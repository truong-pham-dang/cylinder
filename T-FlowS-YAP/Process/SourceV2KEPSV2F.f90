!======================================================================!
  SUBROUTINE SourceV2KEPSV2F()
!----------------------------------------------------------------------!
! Purpose:                                                             !
! Calculate source terms in equation for v2                            !
! Term which is negativ is put on left hand side in diagonal of        !
! martix of coefficient                                                !
!                                                                      !  
! Authors: Muhamed Hadziabdic and Bojan Niceno                         !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c
!--------------------------------[CVS]---------------------------------!
!  $Id: SourceV2KEPSV2F.f90,v 1.4 2008/11/19 14:57:57 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/SourceV2KEPSV2F.f90,v $  
!======================================================================!


!----------------------------------------------------------------------!
!    In transport equation for v2 two                                  !
!    source terms exist which have form:                               !
!                                                                      !
!    /                 /                                               !
!   |                 |                                                !
!   | f22*Kin*dV  -   |(vi2*Eps/Kin)*dV                                !
!   |                 |                                                !
!  /                 /                                                 !
!                                                                      !
!  First term can appire as pozitiv and as negative as well so         !
!  depend of sign of term , it is placed on left or right hand side.   !
!  Second, negative source term is added to main diagonal left hand    !
!  side coefficient matrix in order to increase stability of solver    !
!----------------------------------------------------------------------!      

  if(SIMULA == ZETA.or.SIMULA==HYB_ZETA) then
!---- Positive source term 
    do c = 1, NC
      b(c) = b(c) + max(0.0,f22%n(c)*volume(c))
      Aval(Adia(c))= Aval(Adia(c))                                       &  
                   + max(0.0,-f22%n(c)*volume(c)/(v_2%n(c) + tiny))    
    end do

!---- Negative source term
  
    do c = 1, NC
      Aval(Adia(c)) =  Aval(Adia(c)) + volume(c)*Pk(c)/(Kin%n(c)+tiny)
    end do   
  else if(SIMULA == K_EPS_VV) then
    do c = 1, NC
      b(c) = b(c) + max(0.0,f22%n(c)*Kin%n(c)*volume(c))
      Aval(Adia(c))= Aval(Adia(c))                                       &
                     + max(0.0,-f22%n(c)*Kin%n(c)*volume(c)/(v_2%n(c) + tiny))
    end do
    do c = 1, NC  
      Aval(Adia(c)) =  Aval(Adia(c)) + volume(c)*Eps%n(c)/(Kin%n(c)+tiny) 
    end do
  end if
  RETURN
  END SUBROUTINE SourceV2KEPSV2F
