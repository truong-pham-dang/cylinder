!======================================================================!
  SUBROUTINE Resid(N, NB, NONZ, A, Acol,Arow,Ab,x,r1) 
!----------------------------------------------------------------------!
!   Calculates residuals.                                              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER  :: N, NB, NONZ      

  REAL     :: A(NONZ),Ab(-NB:-1)
  INTEGER  :: Acol(N)
  INTEGER  :: Arow(NONZ)
  REAL     :: x(-NB:N), r1(N)             !  [A]{x}={r1}
!-------------------------------[Locals]-------------------------------!
  INTEGER  :: i,j,k,sub
!--------------------------------[CVS]---------------------------------!
!  $Id: Resid.f90,v 1.8 2008/12/10 14:52:43 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/Resid.f90,v $  
!======================================================================!

!+++++++++++++++++++++++++!
!     r = b - Ax          !
!     => Parallelized     ! 
!+++++++++++++++++++++++++!
  do i=1,N
    do j=Acol(i),Acol(i+1)-1     
      k = Arow(j)                 
      r1(i) = r1(i) - A(j) * x(k)  
    end do
  end do
!      call exchange(x) 
  do sub=1,Npro
    if(NBBe(sub)  <=  NBBs(sub)) then
      do k=NBBs(sub),NBBe(sub),-1
	i=BufInd(k)
	r1(i) = r1(i) - Ab(k)*x(k)
      end do
    end if
  end do

  END SUBROUTINE Resid
