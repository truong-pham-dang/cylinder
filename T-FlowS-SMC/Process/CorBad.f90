!======================================================================!
  SUBROUTINE CorBad(PHIi)
!----------------------------------------------------------------------!
!   Corrects the pressure gradients in the cells where they cannot     !
!   be computed, the so called "bad" cells.                            !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  REAL :: PHIi(-NbC:NC)
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, c1, c2, s
!--------------------------------[CVS]---------------------------------!
!  $Id: CorBad.f90,v 1.3 2008/12/10 14:19:22 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CorBad.f90,v $  
!======================================================================!

  do c=1,NC
    if(BadForG(c)) then
      PHIi(c) = 0.0
    end if
  end do 

  do s=1,NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)
     
    if(c2 > 0 .or. c2 < 0 .and. TypeBC(c2) == BUFFER) then
      if(BadForG(c1)) PHIi(c1) = PHIi(c1) + 0.5*PHIi(c2) 
      if(BadForG(c2)) PHIi(c2) = PHIi(c2) + 0.5*PHIi(c1) 
    end if
  end do

  call Exchng(PHIi)

  END SUBROUTINE CorBad 
