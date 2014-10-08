!======================================================================!
  SUBROUTINE BufLoa 
!----------------------------------------------------------------------!
! Reads: NAME.buf                                                      !
! ~~~~~~                                                               !
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: c, dummy 
  INTEGER   :: sub, subo, NNsub,NCsub,NSsub,NBCsub,NBFsub
  CHARACTER :: nameIn*80
!--------------------------------[CVS]---------------------------------!
!  $Id: BufLoa.f90,v 1.1 2002/11/01 15:12:16 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Parallel/Single/BufLoa.f90,v $             
!======================================================================!
!  Each subdomain needs two buffers: a send buffer and a receive buffer.
!  A receive buffer will be stored as aditional boundary cells for each
!  subdomain. So each subdomain will have NBC physical boundary faces
!  and NBBC-NBC buffer bounndary cells. It is handy to do it that way,
!  because most of the algorythms can remain the same as they are now.
!  They won't even "know" that they use values from other processors.
!  On the other hand, a sending buffer has to be allocated in a new 
!  separate array called simply buffer(). An additional array is needed 
!  to keep track of all the indexes. That one is called BufInd().
!  BufInd() has stored cell numbers from it's own subdomain so that
!  later they can be copied with (well, something like that):
!  do i=1,BUFFSIZ
!    buffer(i) = U(BufInd(i))
!  end do
!----------------------------------------------------------------------!

  if(Npro == 0) return

  call NamFil(this, nameIn, '.buf', len_trim('.buf'))
  open(9, FILE=nameIn)
  if(this < 2) write(*,*) '# Now reading the file:', nameIn

!///// number of physical boundary cells
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) NBCsub

!///// initialize 
  do sub=0,NPRO
    NBBs(sub) = -(NBCsub) 
    NBBe(sub) = -(NBCsub)
  end do

!///// fill the indexes and the buffers
  do sub=1,NPRO
    if(sub  /=  this) then

!----- connections with subdomain          
      call ReadC(9,inp,tn,ts,te)
      read(inp,*) subo 

!----- number of local connections with subdomain sub 
      call ReadC(9,inp,tn,ts,te)
      read(inp,*) NBBe(sub)

      NBBs(sub) = NBBe(sub-1) - 1  
      NBBe(sub) = NBBs(sub) - NBBe(sub) + 1

      do c=NBBs(sub),NBBe(sub),-1
	call ReadC(9,inp,tn,ts,te)
	read(inp,*) dummy, BufInd(c) 
      end do 
    else
      NBBs(sub) = NBBe(sub-1)-1  ! just to become "sloppy" 
      NBBe(sub) = NBBe(sub-1)    ! this will be needed for next 
    end if
  end do   ! through subdomains

  close(9)

!///// correct the "sloppy" indexes
  do sub=1,NPRO
    if(NBBe(sub)  > NBBs(sub)) then  
      NBBs(sub) = -1 
      NBBe(sub) = 0 
    end if
  end do 

  call wait

!->>>  write(*,*) 'PE',this, '#===================#' 
!->>>  write(*,*) 'PE',this, '# Check connections #' 
!->>>  write(*,*) 'PE',this, '#-------------------#' 
!->>>  do sub=1,NPRO
!->>>    write(*,'(A2,I2,3I7)') 'PE',this, sub, NBBs(sub), NBBe(sub)
!->>>  end do   ! through subdomains

  END SUBROUTINE BufLoa
