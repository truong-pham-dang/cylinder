!======================================================================!
  INTEGER FUNCTION WchNod(c, n) 
!----------------------------------------------------------------------!
!   Returns the local number (1-8) of node n in cell c.                !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: c, n
!-------------------------------[Locals]-------------------------------!
  INTEGER :: i, j
!--------------------------------[CVS]---------------------------------!
!  $Id: WchNod.f90,v 1.6 2002/10/30 16:29:23 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Generate/WchNod.f90,v $  
!======================================================================!

  WchNod=0
  if (c  < 0) then 
    write(*,*) 'Which node: Cell non existent !'
    return
  endif

!----- try the node only 
  do i=1,8
    if( CellN(c,i)  ==  n) then
      goto 10
    end if 
  end do     

!----- if it failed try his twins also
  do j=1,TwinN(n,0)
    do i=1,8
      if( CellN(c,i)  ==  TwinN(n,j)) then
	goto 10
      end if 
    end do     
  end do

  WchNod=0
  write(*,*) 'Which node: Trouble, node not found !'
  write(*,*) 'x,y,z = ', x(n), y(n), z(n)
  write(*,*) 'cell  = ', c, level(c)
  return

10   WchNod=i
  return

  END FUNCTION WchNod
