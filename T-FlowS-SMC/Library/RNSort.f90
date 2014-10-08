!======================================================================!
  SUBROUTINE RNSort(X,indx,N)
!----------------------------------------------------------------------!
!   Sorts real array X according to indx.                              !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: N,indx(N)
  REAL    :: X(N)
!-------------------------------[Locals]-------------------------------!
  INTEGER          :: i
  REAL,ALLOCATABLE :: work(:)
!--------------------------------[CVS]---------------------------------!
!  $Id: RNSort.f90,v 1.3 2002/10/30 16:29:33 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Library/RNSort.f90,v $  
!======================================================================!

  allocate(work(N)); work=0

  do i=1,N
    work(indx(i))=X(i)
  end do

  do i=1,N
    X(i)=work(i)
  end do

  deallocate(work)

  END SUBROUTINE RNSort
