!======================================================================!
  SUBROUTINE Exchng(PHI) 
!----------------------------------------------------------------------!
!   Exchanges the values between the processors.                       !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE par_mod
  USE pro_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------[Include]-------------------------------!
  INCLUDE 'mpif.h'
!-----------------------------[Parameters]-----------------------------!
  REAL    :: PHI(-NbC:NC)
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c1, c2, sub, rtag, stag, length, error
  INTEGER :: status(MPI_STATUS_SIZE)
!--------------------------------[CVS]---------------------------------!
!  $Id: Exchng.f90,v 1.1 2002/11/01 15:12:16 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Parallel/Single/Exchng.f90,v $  
!======================================================================!

!///// fill the buffers with new values
  do sub=1,NPro
    if( NBBe(sub)  <=  NBBs(sub) ) then  
      do c2=NBBs(sub),NBBe(sub),-1
	c1 = BufInd(c2)
	PHI(c2) = PHI(c1)
      end do
    end if
  end do   

!///// exchange the values
  do sub=1,NPro
    if( NBBe(sub)  <=  NBBs(sub) ) then  

      length = NBBs(sub) - NBBe(sub) + 1
      stag = MAXPRO*this + sub    ! tag for sending
      rtag = MAXPRO*sub + this    ! tag for receivinging

!===============================================
      call MPI_SENDRECV_REPLACE & 
!-------------------------------------+---------
	     (PHI(NBBe(sub)),   & ! buffer  
	      length,           & ! length   
	      MPI_REAL,         & ! datatype  
!-------------------------------------+---------
	      (sub-1),          & ! dest,      
	      stag,             & ! sendtag,    
!-------------------------------------+---------
	      (sub-1),          & ! source,      
	      rtag,             & ! recvtag,      
!-------------------------------------+---------
	      MPI_COMM_WORLD,   &
	      status,           &
	      error) 
!===============================================

    end if  !  NBBe(sub)  /=  NBBs(sub)
  end do

  END SUBROUTINE Exchng
