!======================================================================!
  SUBROUTINE NearWallCell()
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER          :: k,  c, nearest  
  REAL             :: DISTnew, DISTold
  REAL             :: Dist

!--------------------------------[CVS]---------------------------------!
  character*80 rcs1,rcs2
  data rcs1/                                                        &
  '$Id: NearWallCell.f90,v 1.3 2008/12/10 14:47:15 IUS\mhadziabdic Exp $'/
  data rcs2/                                                        &
  '$Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/NearWallCell.f90,v $'/
!======================================================================!  

!--------------------------------------------------------------------------------!
! Purpose: This program is finding for every inside domain cell corresponding cell
! ~~~~~~~~ on the nearest wall. This is needed for calculation y+  
!--------------------------------------------------------------------------------!

  if(this  < 2)                                                     &
  write(*,*) '# Now searching for corresponding wall cells!'

  nearest = 0
  near = 0
  DISTold = HUGE
  do c = 1, NC
    DISTold = HUGE
    do k = 1, NC
      if(IsNearWall(k)) then
        DISTnew = Dist(xc(k),yc(k),zc(k),xc(c),yc(c),zc(c))
        if(DISTnew <= DISTold) then
          nearest =  k
          DISTold = DISTnew
        end if 
      end if
    end do
    near(c) = nearest 
  end do

  if(this < 2) write(*,*) '# Searching finished'

  END SUBROUTINE NearWallCell

