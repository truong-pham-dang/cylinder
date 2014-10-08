!======================================================================!
  SUBROUTINE FindBad()
!----------------------------------------------------------------------!
! Searches for cells which are "bad" for calculation of pressure       !
! gradients.                                                           !
!                                                                      !
! Practically, these are the tetrahedronal cells with two faces on the !
! boundary and two in the domain.                                      ! 
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: s, c, c1, c2, NumBad
!--------------------------------[CVS]---------------------------------!
!  $Id: FindBad.f90,v 1.3 2008/12/10 14:36:21 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/FindBad.f90,v $      
!======================================================================!

  BadForG = .FALSE. 
  NumGood = 0

  NumBad = 0

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2 < 0) then
      NumGood(c1) = NumGood(c1) + 1 
    end if  
  end do

  do c=1,NC
    if(NumGood(c)==2) then
      BadForG(c) = .TRUE.
      NumBad = NumBad + 1
    end if
  end do 

  NumGood = 0
  
  call IGlSum(NumBad)

  if(THIS < 2) write(*,*) '# There are ', NumBad, &
                          ' bad cells for gradients.'

  END SUBROUTINE FindBad
