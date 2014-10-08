!======================================================================!
  SUBROUTINE ForAPF 
!----------------------------------------------------------------------!
!   Forms the pressure system matrix for the fractional step method.   !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  REAL    :: A12
  INTEGER :: c, c1, c2, s 
!--------------------------------[CVS]---------------------------------!
!  $Id: ForAPF.f90,v 1.8 2008/12/10 14:36:38 IUS\mhadziabdic Exp $
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/ForAPF.f90,v $
!======================================================================!

  do c=1,Acol(NC+1) ! to je broj nonzero + 1
    Aval(c) = 0.0
  end do

!+++++++++++++++++++++++++++++++++!
!     Calculate system matrix     ! 
!+++++++++++++++++++++++++++++++++!
    do s=1,NS    

      c1=SideC(1,s)
      c2=SideC(2,s)

      if(c2  > 0) then
	A12 = dt * Scoef(s) 
	Aval(SidAij(1,s))  = -A12
	Aval(SidAij(2,s))  = -A12
	Aval(Adia(c1)) =                                            &
	Aval(Adia(c1)) +  A12
	Aval(Adia(c2)) =                                            &
	Aval(Adia(c2)) +  A12
      else
	if(TypeBC(c2) == BUFFER) then
	  A12 = dt * Scoef(s)
	  Aval(Adia(c1)) =                                          &
	  Aval(Adia(c1)) +  A12
	  Abou(c2) = -A12
	end if
      end if 

    end do ! through sides

  END SUBROUTINE ForAPF
