!======================================================================!
  SUBROUTINE CalcG(Boundary)
!----------------------------------------------------------------------!
!   Calculates gradient matrix.                                        !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  LOGICAL :: Boundary
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, c1, c2, s
  REAL    :: Dxc1, Dyc1, Dzc1, Dxc2, Dyc2, Dzc2
  REAL    :: Jac, Ginv(6)
!--------------------------------[CVS]---------------------------------!
!  $Id: CalcG.f90,v 1.11 2008/12/05 13:07:34 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CalcG.f90,v $  
!======================================================================!

!+++++++++++++++++++++++++++++++++++++++++!
!    Calculate matrix G for all cells     !
!+++++++++++++++++++++++++++++++++++++++++!
  do c=1,NC
    G(1,c) = 0.0
    G(2,c) = 0.0
    G(3,c) = 0.0
    G(4,c) = 0.0
    G(5,c) = 0.0
    G(6,c) = 0.0
  end do

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s) 

    Dxc1=Dx(s)
    Dyc1=Dy(s)
    Dzc1=Dz(s)
    Dxc2=Dx(s)
    Dyc2=Dy(s)
    Dzc2=Dz(s)

!---- Take care of material interfaces         ! 2mat
    if( StateMat(material(c1))==FLUID .and. &  ! 2mat
        StateMat(material(c2))==SOLID       &  ! 2mat 
        .or.                                &  ! 2mat
        StateMat(material(c1))==SOLID .and. &  ! 2mat
        StateMat(material(c2))==FLUID ) then   ! 2mat
      Dxc1 = xsp(s)-xc(c1)                     ! 2mat
      Dyc1 = ysp(s)-yc(c1)                     ! 2mat
      Dzc1 = zsp(s)-zc(c1)                     ! 2mat 
      Dxc2 = xsp(s)-xc(c2)                     ! 2mat
      Dyc2 = ysp(s)-yc(c2)                     ! 2mat 
      Dzc2 = zsp(s)-zc(c2)                     ! 2mat
    end if                                     ! 2mat

!---- With boundary cells, velocities, temperatures
    if(Boundary) then
      G(1,c1)=G(1,c1) + Dxc1*Dxc1  ! 1,1
      G(2,c1)=G(2,c1) + Dyc1*Dyc1  ! 2,2
      G(3,c1)=G(3,c1) + Dzc1*Dzc1  ! 3,3
      G(4,c1)=G(4,c1) + Dxc1*Dyc1  ! 1,2  &  2,1
      G(5,c1)=G(5,c1) + Dxc1*Dzc1  ! 1,3  &  3,1
      G(6,c1)=G(6,c1) + Dyc1*Dzc1  ! 2,3  &  3,2
      if(c2  > 0) then  ! this is enough even for parallel
        G(1,c2)=G(1,c2) + Dxc2*Dxc2  ! 1,1
        G(2,c2)=G(2,c2) + Dyc2*Dyc2  ! 2,2
        G(3,c2)=G(3,c2) + Dzc2*Dzc2  ! 3,3
        G(4,c2)=G(4,c2) + Dxc2*Dyc2  ! 1,2  &  2,1
        G(5,c2)=G(5,c2) + Dxc2*Dzc2  ! 1,3  &  3,1
        G(6,c2)=G(6,c2) + Dyc2*Dzc2  ! 2,3  &  3,2
      end if
!---- Without boundary cells => pressure
    else ! Don't use Boundary
      if(c2 > 0 .or. c2 < 0 .and. TypeBC(c2) == BUFFER) then  
        G(1,c1)=G(1,c1) + Dxc1*Dxc1  ! 1,1
        G(2,c1)=G(2,c1) + Dyc1*Dyc1  ! 2,2
        G(3,c1)=G(3,c1) + Dzc1*Dzc1  ! 3,3
        G(4,c1)=G(4,c1) + Dxc1*Dyc1  ! 1,2  &  2,1
        G(5,c1)=G(5,c1) + Dxc1*Dzc1  ! 1,3  &  3,1
        G(6,c1)=G(6,c1) + Dyc1*Dzc1  ! 2,3  &  3,2
      end if
      if(c2 > 0) then
        G(1,c2)=G(1,c2) + Dxc2*Dxc2  ! 1,1
        G(2,c2)=G(2,c2) + Dyc2*Dyc2  ! 2,2
        G(3,c2)=G(3,c2) + Dzc2*Dzc2  ! 3,3
        G(4,c2)=G(4,c2) + Dxc2*Dyc2  ! 1,2  &  2,1
        G(5,c2)=G(5,c2) + Dxc2*Dzc2  ! 1,3  &  3,1
        G(6,c2)=G(6,c2) + Dyc2*Dzc2  ! 2,3  &  3,2
      end if
    end if ! Boundary
  end do

!--------------------------------------!
!     Find the inverse of matrix G     !
!--------------------------------------!
  do c=1,NC
    Jac  =         G(1,c) * G(2,c) * G(3,c)                         &
	   -       G(1,c) * G(6,c) * G(6,c)                         &
	   -       G(4,c) * G(4,c) * G(3,c)                         &
	   + 2.0 * G(4,c) * G(5,c) * G(6,c)                         &
	   -       G(5,c) * G(5,c) * G(2,c)

    Ginv(1) = +( G(2,c)*G(3,c) - G(6,c)*G(6,c) ) / (Jac+TINY)
    Ginv(2) = +( G(1,c)*G(3,c) - G(5,c)*G(5,c) ) / (Jac+TINY)
    Ginv(3) = +( G(1,c)*G(2,c) - G(4,c)*G(4,c) ) / (Jac+TINY)
    Ginv(4) = -( G(4,c)*G(3,c) - G(5,c)*G(6,c) ) / (Jac+TINY)
    Ginv(5) = +( G(4,c)*G(6,c) - G(5,c)*G(2,c) ) / (Jac+TINY)
    Ginv(6) = -( G(1,c)*G(6,c) - G(4,c)*G(5,c) ) / (Jac+TINY)

    G(1,c) = Ginv(1) 
    G(2,c) = Ginv(2)
    G(3,c) = Ginv(3)
    G(4,c) = Ginv(4)
    G(5,c) = Ginv(5)
    G(6,c) = Ginv(6)
  end do 

  RETURN 

  END SUBROUTINE CalcG
