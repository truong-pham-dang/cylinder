!======================================================================!
  SUBROUTINE INSort(X,indx,N)
!----------------------------------------------------------------------!
!   Sorts int. array X according to indx.                              !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: N,X(N),indx(N)
!-------------------------------[Locals]-------------------------------!
  INTEGER             :: i
  INTEGER,ALLOCATABLE :: work(:)
!--------------------------------[CVS]---------------------------------!
!  $Id: INSort.f90,v 1.6 2002/10/30 16:29:33 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Library/INSort.f90,v $  
!======================================================================!

  allocate(work(N)); work=0

  do i=1,N
    work(indx(i))=X(i)
  end do

  do i=1,N
    X(i)=work(i)
  end do

  deallocate(work)

  END SUBROUTINE INSort
