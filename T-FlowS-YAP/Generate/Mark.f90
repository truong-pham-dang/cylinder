!======================================================================!
  SUBROUTINE Mark
!----------------------------------------------------------------------!
!   Mark the region of the domain for local refinement.                !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, lev, regio, n1, n2, n3, n4, n5, n6, n7, n8
  REAL    :: x1, y1, z1, x8, y8, z8, x0, y0, z0
!--------------------------------[CVS]---------------------------------!
!  $Id: Mark.f90,v 1.7 2002/10/30 16:29:22 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Generate/Mark.f90,v $  
!======================================================================!

  do c=-MAXB,MAXN
    CelMar(c) = 0
  end do 

  do lev = 1,NRL
    do regio = 1,NR(lev) 
      x1=FRegio(lev,regio,1)
      y1=FRegio(lev,regio,2)
      z1=FRegio(lev,regio,3)
      x8=FRegio(lev,regio,4)
      y8=FRegio(lev,regio,5)
      z8=FRegio(lev,regio,6)

      do c=1,NC
	n1=CellN(c,1)
	n2=CellN(c,2)
	n3=CellN(c,3)
	n4=CellN(c,4)
	n5=CellN(c,5)
	n6=CellN(c,6)
	n7=CellN(c,7)
	n8=CellN(c,8)

	x0=1.25e-1*(x(n1)+x(n2)+x(n3)+x(n4)+x(n5)+x(n6)+x(n7)+x(n8))
	y0=1.25e-1*(y(n1)+y(n2)+y(n3)+y(n4)+y(n5)+y(n6)+y(n7)+y(n8))
	z0=1.25e-1*(z(n1)+z(n2)+z(n3)+z(n4)+z(n5)+z(n6)+z(n7)+z(n8))

!->>>   write(*,*) x0, y0, z0

	if(FRegio(lev,regio,0) == ELIPSO) then
	  if(  ( ((x1-x0)/x8)**2 +                                  &
		 ((y1-y0)/y8)**2 +                                  &
		 ((z1-z0)/z8)**2)  < 1.0 ) then
	    CelMar(c) = -1
	  end if
	else if(FRegio(lev,regio,0) == RECTAN) then 
	  if( (x1  < x0) .and. (x0  < x8) .and.                     &
	      (y1  < y0) .and. (y0  < y8) .and.                     &
	      (z1  < z0) .and. (z0  < z8) ) then
	    CelMar(c) = -1
	  endif
	else if(FRegio(lev,regio,0) == PLANE) then 
	  if( (x0-x1)*x8+(y0-y1)*y8+(z0-z1)*z8   >  0.0 ) then
	    CelMar(c) = -1
	  endif
	end if 
      end do   ! =-> cells

    end do   ! =-> regio
    call refine(lev)

    do c=-MAXB,MAXN
       CelMar(c) = 0
    enddo 

  end do  ! =-> lev

  END SUBROUTINE Mark
