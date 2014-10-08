!======================================================================!
  SUBROUTINE Read9Ln(string, tn, ts, te) 
!----------------------------------------------------------------------!
! Reads a line from a file (unit 9) and parses it.                     ! 
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  CHARACTER :: string*300
  INTEGER   :: tn, ts(300), te(300)
!-------------------------------[Locals]-------------------------------!
  INTEGER :: i
!--------------------------------[CVS]---------------------------------!
  character*80 rcs1,rcs2
  data rcs1/                                                        &
  '$Id: Read9Ln.f90,v 1.4 2004/06/24 14:09:23 muhamed Exp $'/ 
  data rcs2/                                                        &
  '$Source: /home/muhamed/.CVSROOT/T-Rex/Neu2Trex/Read9Ln.f90,v $'/
!======================================================================!

  read(9,'(A300)') string

!---- Parse tokens. This is somehow cool !
  tn = 0
  if(string(1:1) >= '!') then
    tn = 1
    ts(1)=1
  end if
  do i=1,298
    if( string(i:i)  < '!' .and. string(i+1:i+1) >= '!') then
      tn = tn + 1
      ts(tn) = i+1
    end if
    if( string(i:i) >= '!' .and. string(i+1:i+1)  < '!') then
      te(tn) = i
    end if
  end do

  END SUBROUTINE Read9Ln
