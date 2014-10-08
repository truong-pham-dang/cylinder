!======================================================================!
  SUBROUTINE TopolM 
!----------------------------------------------------------------------!
!   Determines the topology of the system matrix.                      !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER             :: c, s, j, n
  INTEGER             :: c1, c2
  INTEGER,ALLOCATABLE :: stencw(:)
!--------------------------------[CVS]---------------------------------!
!  $Id: TopolM.f90,v 1.11 2008/11/19 14:54:06 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/TopolM.f90,v $  
!======================================================================!
!   Relies only on SideC structure. Try to keep it that way.
!----------------------------------------------------------------------!
		  
!---- memory allocation
  allocate(stencw(NC)); stencw=1

!====================================================================!
!                                                                    !
!     Determine the topology of unstructured system of equations     !
!                                                                    !
!====================================================================!

  if(this < 2) write(*,*) '# Determining the topology of the system matrix.'

!---- 
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2  > 0) then
      stencw(c1)=stencw(c1)+1
      stencw(c2)=stencw(c2)+1
    end if      
  end do

!---- Count the nonzero entries and allocate the memory for the array 
  n = 0  
  do c=1,NC
    n = n + stencw(c)
  end do   
  NONZERO = n+1
  allocate(Aval(n+1)); Aval=0 ! I think it reffers to Acol+1 somewhere
  allocate(Arow(n+1)); Arow=0 ! I think it reffers to Acol+1 somewhere
!---- 
  Acol(1)=1
  do c=1,NC
    Acol(c+1)=Acol(c)+stencw(c)
    Arow(Acol(c)) = c   ! sam sebi je prvi
    stencw(c)=1
  end do 
!---- 
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2  > 0) then
      Arow(Acol(c1)+stencw(c1)) = c2
      Arow(Acol(c2)+stencw(c2)) = c1
      stencw(c1)=stencw(c1)+1
      stencw(c2)=stencw(c2)+1
    end if      
  end do
!---- 
  do c=1,NC
    call isort(Arow(Acol(c)),                                       &
	       Arow(Acol(c)),stencw(c),1)
    do j=Acol(c),Acol(c+1)-1
      if(Arow(j) == c) Adia(c)=j
    end do
  end do 

!----
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2  > 0) then
!----- where is A(c1,c2) and ...
      do c=Acol(c1),Acol(c1+1)-1 
	if(Arow(c)  ==  c2) SidAij(1,s)=c
      end do
!----- where is A(c2,c1) and ...
      do c=Acol(c2),Acol(c2+1)-1 
	if(Arow(c)  ==  c1) SidAij(2,s)=c
      end do
    end if
  end do

  if(this < 2) write(*,*) '# Finished !'
 
  deallocate(stencw)

  RETURN

  END SUBROUTINE TopolM
