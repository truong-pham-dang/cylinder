!======================================================================!
  SUBROUTINE UserSavRaw
!----------------------------------------------------------------------!
! Writes: NAME.rawdata                                                 !
! ~~~~~~~                                                              !
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: c, NCtot
  CHARACTER :: namOut*80
!--------------------------------[CVS]---------------------------------!
!  $Id: UserSavRaw.f90,v 1.5 2002/10/30 16:30:03 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserSavRaw.f90,v $  
!======================================================================!

!------------------------------!
!     Create raw data file     !
!------------------------------!
  namOut = name
  call NamFil(this, namOut, '.rawdata', len_trim('.rawdata') )
  open(9, FILE=namOut)
  write(6, *) 'Now creating the file:', namOut

!---- Total number of cells
  NCtot = NC
  call IGlSum(NCtot)
  if(this < 2) then
    write(9,'(I9)') NCtot
  end if 

!---- Write raw data file 
  do c=1,NC
    write(9,'(6E16.6)') xc(c),yc(c),zc(c),U % n(c),V % n(c),W % n(c)
  end do

  close(9)

  END SUBROUTINE UserSavRaw
