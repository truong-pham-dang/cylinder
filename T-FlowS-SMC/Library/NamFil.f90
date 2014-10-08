!======================================================================!
  SUBROUTINE NamFil(sub, namOut, ext, lext)
!----------------------------------------------------------------------!
!   Creates the file name depending on the subdomain and file type.    !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER       :: sub, lext
  CHARACTER*(*) :: ext
  CHARACTER*(*) :: namOut
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: c
  CHARACTER :: numb*4
!--------------------------------[CVS]---------------------------------!
!  $Id: NamFil.f90,v 1.7 2002/10/30 16:29:33 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Library/NamFil.f90,v $  
!======================================================================!

  namOut = name

  if(sub == 0) then
    namOut(len_trim(name)+1:len_trim(name)+1+lext-1) = ext(1:lext) 
  else
    write(numb,'(I4)') sub
    write(namOut(len_trim(name)+1:len_trim(name)+5),'(A5)') '-0000' 
    do c=1,4
      if( numb(c:c) >= '0' .and. numb(c:c) <= '9' )                 &
	namOut(len_trim(name)+1+c:len_trim(name)+1+c) = numb(c:c)
    end do
    namOut(len_trim(name)+6:len_trim(name)+6+lext-1) = ext(1:lext) 
  end if 

  END SUBROUTINE NamFil
