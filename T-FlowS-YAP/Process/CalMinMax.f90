!======================================================================!
  SUBROUTINE CalMinMax(PHI)
!----------------------------------------------------------------------!
!   Extrapoloate variables on the boundaries where needed              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  REAL    :: PHI(-NbC:NC) 
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c1, c2, s
!--------------------------------[CVS]---------------------------------!
!  $Id: CalMinMax.f90,v 1.3 2008/12/05 13:02:37 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CalMinMax.f90,v $  
!======================================================================!

  PHImax = PHI 
  PHImin = PHI 

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

    if( (c2>0) .or. (c2<0 .and. TypeBC(c2)==BUFFER) ) then
      PHImax(c1) = max(PHImax(c1),PHI(c2))
      PHImin(c1) = min(PHImin(c1),PHI(c2))
      PHImax(c2) = max(PHImax(c2),PHI(c1))
      PHImin(c2) = min(PHImin(c2),PHI(c1))
    end if

  end do

  call exchng(PHImax)
  call exchng(PHImin)

  RETURN 

  END SUBROUTINE CalMinMax
