!======================================================================!
  SUBROUTINE BCelLoa
!----------------------------------------------------------------------!
! Reads:  NAME.faces.gmv  NAME.shadow.gmv                              !
! ~~~~~~                                                               !                                                                           
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod 
  USE div_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER      :: c1, c2, n, s, dum1
  CHARACTER*80 :: dummy, nameIn
!--------------------------------[CVS]---------------------------------!
!  $Id: BCelLoa.f90,v 1.10 2002/10/30 16:29:18 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Divide/BCelLoa.f90,v $  
!======================================================================!

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Read the GMV file with      !
!        cell face nodes          !
!                                 !
!     (it has been done by        !
!      copying and modifying      ! 
!          GenSav.f90)            !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  nameIn = name
  nameIn(len_trim(name)+1:len_trim(name)+10) = '.faces.gmv'
  write(*,*) '# Now reading the file:', nameIn
  open(9, FILE=nameIn)

!---------------!
!     start     !
!---------------!
  read(9,'(A80)') dummy 

!---------------!
!     nodes     !
!---------------!
  read(9,'(A80)') dummy 
  do n=1,3*NN
    read(9,'(A80)') dummy 
  end do  

!----------------------!
!     cell section     !
!----------------------!
  read(9,'(A80)') dummy 
  do s=1,NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)
    read(9,*) dummy, dum1
    if(dummy == 'tri') then 
      SideN(s,0) = 3
      read(9,*) &
        SideN(s,1), SideN(s,2), SideN(s,3)
    else if(dummy == 'quad') then
      SideN(s,0) = 4
      read(9,*) &
        SideN(s,1), SideN(s,2), SideN(s,3), SideN(s,4)  
    else
      write(*,*) 'Unsupported cell-face type:', dummy
      write(*,*) 'Exiting'
      stop
    end if
  end do  

  close(9)

!>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                          !
!     read shadow file     !
!                          !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  call NamFil(0, nameIn, '.shadow.gmv', len_trim('.shadow.gmv'))
  open(9, FILE=nameIn)
  write(6, *) 'Now reading the file:', nameIn

  do s=NS+1,NS+NSsh
    read(9,*) SideN(s,0)
    if(SideN(s,0)==3) then
      read(9,*) SideN(s,1), SideN(s,2), SideN(s,3),             &
                SideC(1,s), SideC(2,s)
    else if(SideN(s,0)==4) then
      read(9,*) SideN(s,1), SideN(s,2), SideN(s,3), SideN(s,4), &
                SideC(1,s), SideC(2,s)
    end if
  end do

  close(9)

  END SUBROUTINE BCelLoa
