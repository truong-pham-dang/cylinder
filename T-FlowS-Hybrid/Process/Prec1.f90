!======================================================================!
  SUBROUTINE Prec1(N,NONZ,A,Acol,Arow,Adia,D,prec) 
!----------------------------------------------------------------------!
!   Solves the linear systems of equations by a precond. CG Method.    !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER  :: N, NONZ      

  REAL     :: A(NONZ)
  INTEGER  :: Acol(N),Adia(N)
  INTEGER  :: Arow(NONZ)
  REAL     :: D(N) 

  INTEGER  :: prec
!-------------------------------[Locals]-------------------------------!
  REAL     :: sum1
  INTEGER  :: i, j, k
!--------------------------------[CVS]---------------------------------!
!  $Id: Prec1.f90,v 1.9 2008/12/10 15:29:38 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/Prec1.f90,v $  
!======================================================================!
		 
!->>>
!      INTEGER c 
!      do c=1,N
!        write(*,*) 'Cell: ', c
!        write(*,*) 'Width: ', Acol(c+1)-Acol(c)
!        write(*,'(3I7)') Acol(c), Adia(c), Acol(c+1)-1
!        write(*,*) 'Diag: ', A(Adia(c))
!        write(*,'(25F15.9)') ( A(j),     j=Acol(c),Acol(c+1)-1 )
!        write(*,'(25I7)') ( Arow(j), j=Acol(c),Acol(c+1)-1 )
!        write(*,*) '- - - - - - - - - - - - - - - - - - - - - - -'
!      end do

!+++++++++++++++++++++++++!
!     preconditioning     !
!+++++++++++++++++++++++++!

!----- 1) diagonal preconditioning -----!
!     => Parallelized                   ! 
  if(prec == 1) then        
    do i=1,N                     
      D(i)=A(Adia(i))           
    end do                      

!----- 2) incomplete cholesky preconditioning -----!
  else if(prec == 2) then   
    do i=1,N
      sum1=A(Adia(i))          
      do j=Acol(i), Adia(i)-1         ! only lower traingular
	k=Arow(j)                    
	sum1= sum1- D(k) * A(j)*A(j)  
      end do
      D(i) = 1.0 / sum1                 ! BUG ?
    end do

!----- .) no preconditioning -----!
!     => Parallelized             ! 
  else                          
    do i=1,N
      D(i)=1.0
    end do
  end if 

  END SUBROUTINE Prec1
