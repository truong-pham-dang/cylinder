!======================================================================!
  SUBROUTINE Set_BC 
!----------------------------------------------------------------------!
!   Extrapoloate variables on the boundaries where needed              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c1, c2, s
  REAL    :: R
!--------------------------------[CVS]---------------------------------!
!  $Id: CalcConvect.f90,v 1.2 2002/10/30 16:29:48 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Process/CalcConvect.f90,v $  
!======================================================================!
 
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

    if(c2  < 0) then
      if(bcmark(c2) == 2) then
        R = sqrt(xc(c2)**2 + yc(c2)**2) 
        if(R <= 1.0) then
          bcmark(c2) = 3
        else if(R >= 2.0) then
          bcmark(c2) = 4
        end if
      end if
    end if
  end do

  RETURN 

  END SUBROUTINE Set_BC 
