!======================================================================!
  SUBROUTINE ToUppr(string)
!----------------------------------------------------------------------!
!   Transforms string to uppercase.                                    !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  CHARACTER*(*) :: string  
!-------------------------------[Locals]-------------------------------!
  INTEGER :: i, value
!--------------------------------[CVS]---------------------------------!
!  $Id: ToUppr.f90,v 1.7 2002/10/30 16:29:33 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Library/ToUppr.f90,v $  
!======================================================================!

  do i=1,len_trim(string)
    value = ichar(string(i:i))
    if (value  >=  97 .and. value  <=  122) then 
      string(i:i) = char(value-32) 
    end if
  end do

  END SUBROUTINE ToUppr
