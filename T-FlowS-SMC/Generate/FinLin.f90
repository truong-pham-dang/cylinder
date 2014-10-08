!======================================================================!
  SUBROUTINE FinLin(n1, n2, res) 
!----------------------------------------------------------------------!
!   Searches for a smallest block where the line defined by n1-n2 is.  !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: n1, n2, res
!-------------------------------[Locals]-------------------------------!
  INTEGER :: b, l1, l2
!--------------------------------[CVS]---------------------------------!
!  $Id: FinLin.f90,v 1.6 2002/10/30 16:29:21 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Generate/FinLin.f90,v $  
!======================================================================!

  do b=1,Nbloc
    do l1=1,8
      do l2=1,8
	if( (BlkPnt(b,l1) == n1) .and. &
	    (BlkPnt(b,l2) == n2) ) then
	  if( iabs(l2-l1) == 1 ) res=BlkRes(b,1) 
	  if( iabs(l2-l1) == 2 ) res=BlkRes(b,2) 
	  if( iabs(l2-l1) == 4 ) res=BlkRes(b,3) 
	  goto 1
	end if 
      end do
    end do     
  end do 

  write(6,*) 'ERROR MESSAGE FROM GENX'
  write(6,*) 'You tried to define the line', n1, n2, ' but it'
  write(6,*) 'doesn''t exists in the block specifications.'
  write(6,*) 'Exiting !'
  stop

1 return

  END SUBROUTINE FinLin
