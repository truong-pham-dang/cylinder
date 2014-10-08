!======================================================================!
  SUBROUTINE IGlSum(PHI) 
!----------------------------------------------------------------------!
!   Estimates global summ among all processors.                        !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------[Include]-------------------------------!
  INCLUDE 'mpif.h'
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: PHI
!-------------------------------[Locals]-------------------------------!
  INTEGER :: PHInew
  INTEGER :: error
!--------------------------------[CVS]---------------------------------!
!  $Id: IGlSum.f90,v 1.1 2002/11/01 15:12:16 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Parallel/Single/IGlSum.f90,v $  
!======================================================================!

!================================================
      call MPI_ALLREDUCE      &               
!-----------------------------------+------------
	     (PHI,            & ! send buffer
	      PHInew,         & ! recv buffer 
	      1,              & ! length     
	      MPI_INTEGER,    & ! datatype  
	      MPI_SUM,        & ! operation 
	      MPI_COMM_WORLD, &             
	      error) 
!================================================

  PHI = PHInew

  END SUBROUTINE IGlSum
