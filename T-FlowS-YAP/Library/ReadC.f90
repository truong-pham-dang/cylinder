!======================================================================!
  SUBROUTINE ReadC(un, string, tn, ts, te) 
!----------------------------------------------------------------------!
!   Reads a line from a file (unit 9) and discards if it is comment.   !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  CHARACTER :: string*300
  INTEGER   :: un, tn, ts(300), te(300)
!-------------------------------[Locals]-------------------------------!
  INTEGER :: i
!--------------------------------[CVS]---------------------------------!
!  $Id: ReadC.f90,v 1.5 2002/10/30 16:29:33 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Library/ReadC.f90,v $            
!======================================================================!
!   A comment is each line which begins with one of the following      !
!   characters: ! # $ % / & [ ] { } : @ < > * (Quite a lot, huh ?)     !
!   Input line must not exceed the lenght of 300 characters.           !
!   Note: not very carefully checked, but so far it works.             !
!----------------------------------------------------------------------!

1 read(un,'(A300)') string

  if( string  ==  '' ) goto 1 ! see: man ascii

  if( string(1:1) == '!' .or.                                       &
      string(1:1) == '#' .or.                                       &
      string(1:1) == '$' .or.                                       &
      string(1:1) == '%' .or.                                       &
      string(1:1) == '/' .or.                                       &
      string(1:1) == '&' .or.                                       &
      string(1:1) == '[' .or.                                       &
      string(1:1) == ']' .or.                                       &
      string(1:1) == '{' .or.                                       &
      string(1:1) == '}' .or.                                       &
      string(1:1) == ':' .or.                                       &
      string(1:1) == ';' .or.                                       &
      string(1:1) == '@' .or.                                       &
      string(1:1) == '<' .or.                                       &
      string(1:1) == '>' .or.                                       &
      string(1:1) == '*' ) goto 1

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

  END SUBROUTINE ReadC
