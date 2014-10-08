!======================================================================!
  SUBROUTINE cg(N, NB, NONZ, A, Acol,Arow,Adia,Ab,x,r1,             &
		prec,niter,tol,                                     &
		IniRes,FinRes)
!----------------------------------------------------------------------!
!   Solves the linear systems of equations by a precond. CG Method.    !
!----------------------------------------------------------------------!
!   The structure of the system matrix is described by Acol, Arow,     !
!   and Adia vectors.                                                  !
!                                                                      !
!   Allows preconditioning of the system by:                           !
!     1. Diagonal preconditioning                                      !
!     2. Incomplete Cholesky preconditioning                           !
!                                                                      !
!   The type of precondtioning is chosen by setting the variable prec  !
!   to 0 (no preconditioning), 1 (diagonal preconditioning) or 2       !
!   (incomplete Cholesky preconditioning)                              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE sol_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER  :: N, NB, NONZ

  REAL     :: A(NONZ),Ab(-NB:-1)
  INTEGER  :: Acol(N+1),Adia(N)
  INTEGER  :: Arow(NONZ)
  REAL     :: x(-NB:N), r1(N)            !  [A]{x}={r1}

  INTEGER  :: prec,  niter               !  preconditioning
  REAL     :: tol                        !  tolerance
  REAL     :: IniRes, FinRes             !  residual
!--------------------------------[CVS]---------------------------------!
!  $Id: CG.f90,v 1.11 2008/11/19 14:48:29 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CG.f90,v $  
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
  call Prec1(N, NONZ, A, Acol,Arow,Adia,D,prec)

!????????????????????????????????????!
!      This is quite tricky point.   !
!     What if bnrm2 is very small ?  !
!     => Parallelized                ! 
!????????????????????????????????????!
  bnrm2=0.0
  do i=1,N
    bnrm2=bnrm2+r1(i)*r1(i)
  end do  
  call glosum(bnrm2) 
  bnrm2=sqrt(bnrm2)

  if(bnrm2 < tol) then 
    iter=0
    goto 1
  end if  

!+++++++++++++++++++++++++!
!     r = b - Ax          !
!     => Parallelized     ! 
!+++++++++++++++++++++++++!
  call Resid(N,NB,NONZ,A,Acol,Arow,Ab,x,r1) 

!+++++++++++++++++++++++++!
!     p = r               !
!     => Parallelized     ! 
!+++++++++++++++++++++++++!
  do i=1,N
    p1(i)=r1(i) 
  end do

!++++++++++++++++++++++++++++++++++++!
!     calculate initial residual     !
!     => Parallelized                ! 
!++++++++++++++++++++++++++++++++++++!
  error=0.0
  do i=1,N
    error=error + r1(i)*r1(i)
  end do
  call glosum(error) 
  error  = sqrt(error)  

!-------------------------------------------------------------------!
!     residual after the correction and before the new solution     !
!-------------------------------------------------------------------!
  IniRes=error

  if(error < tol) then
    iter=0
    goto 1
  end if  

!>>>>>>>>>>>>>>>>>>>!
!     MAIN LOOP     !
!<<<<<<<<<<<<<<<<<<<!
  do iter=1, niter

!++++++++++++++++++++++++!  
!       solve Mz = r     !
!     (q instead of z)   !
!++++++++++++++++++++++++!

!-------------------------------------! 
!     1) diagonal preconditioning     !
!     => Parallelized                 !
!-------------------------------------!
    if(prec == 1) then
      do i=1,N
	q1(i)=r1(i)/D(i)
      end do

!------------------------------------------------! 
!     2) incomplete cholesky preconditioning     !
!------------------------------------------------!
   else if(prec == 2) then

!----- forward substitutionn
      do i=1,N
	sum1=r1(i)
	do j=Acol(i),Adia(i)-1  ! only the lower triangular
	  k=Arow(j)             
	  sum1= sum1- A(j)*q1(k)  
	end do
	q1(i) = sum1* D(i)         ! BUG ?
      end do

      do i=1,N
	q1(i) = q1(i) / ( D(i) + TINY )
      end do

!----- backward substitution
      do i=N,1,-1
	sum1=q1(i)
	do j = Adia(i)+1, Acol(i+1)-1 ! upper triangular 
	  k=Arow(j)                  
	  sum1= sum1- A(j)*q1(k)      
	end do
	q1(i) = sum1* D(i)               ! BUG ?
      end do

!-------------------------------!
!     .) no preconditioning     !
!     => Parallelized           !
!-------------------------------!
    else
      do i=1,N
	q1(i)=r1(i)
      end do
    end if

!+++++++++++++++++++++++++!
!     rho = (r,z)         !
!     => Parallelized     !
!+++++++++++++++++++++++++!
    rho=0.0
    do i=1,N
      rho=rho+r1(i)*q1(i)
    end do
    call glosum(rho)

    if(iter == 1) then
      do i=1,N
	p1(i)=q1(i)
      end do        
    else
      beta=rho/rhoold
      do i=1,N
	p1(i) = q1(i) + beta*p1(i)
      end do
    end if

!++++++++++++++++++++++++++++!
!     q    = Ap              !     
!     => Parallelized        !
!++++++++++++++++++++++++++++!
    do i=1,N
      q1(i) = 0.0                    
      do j=Acol(i), Acol(i+1)-1  
	k=Arow(j)                
	q1(i) = q1(i) + A(j) * p1(k) 
      end do
    end do
    call Exchng(p1)
    do sub=1,Npro
      if(NBBe(sub)  <=  NBBs(sub)) then
	do k=NBBs(sub),NBBe(sub),-1
	  i=BufInd(k)
	  q1(i) = q1(i) + Ab(k)*p1(k)
	end do
      end if
    end do

!++++++++++++++++++++++++++++!
!     alfa = (r,z)/(p,q)     !
!     => Parallelized        !
!++++++++++++++++++++++++++++!
    alfa=0.0
    do i=1,N
      alfa=alfa+p1(i)*q1(i)
    end do
    call glosum(alfa)       
    alfa=rho/alfa

!+++++++++++++++++++++++++!
!     x = x + alfa p      !
!     r = r - alfa Ap     !
!     => Parallelized     !
!+++++++++++++++++++++++++!
    do i=1,N
      x(i)=x(i)   + alfa*p1(i)
      r1(i)=r1(i) - alfa*q1(i)
    end do

!???????????????????????????!
!     check convergence     !
!     => Parallelized       !
!???????????????????????????!
    error=0.0
    do i=1,N
      error=error+r1(i)*r1(i)
    end do  
    call glosum(error)       
    error=sqrt(error)  

    if(error < tol) goto 1

    rhoold=rho

  end do                ! iter 

1 FinRes = error
  niter = iter

  RETURN

  END SUBROUTINE cg
