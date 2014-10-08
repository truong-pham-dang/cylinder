!======================================================================!
  SUBROUTINE TestU()
!----------------------------------------------------------------------!
!   Discretizes and solves momentum conservation equations             !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: n,c,s,c1,c2,niter
  REAL    :: VISeff,FUex,FUim,PHIxS,PHIyS,PHIzS,A0,A12,A21,error
!--------------------------------[CVS]---------------------------------!
!  $Id: TestU.f90,v 1.17 2008/11/19 14:54:27 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/TestU.f90,v $      
!----------------------------------------------------------------------!
!     
!  The form of equations which I am solving:    
!     
!     /              /            
!    |     du       |             
!    | rho -- dV  = | mu DIV u dS 
!    |     dt       |             
!   /              /              
!
!======================================================================!

  do n=1,Acol(NC+1) ! to je broj nonzero + 1
    Aval(n) = 0.0
  end do

  do c=1,NC
    b(c)=0.0
  end do

!-----------------------------------------! 
!     Initialize variables and fluxes     !
!-----------------------------------------! 

!----- old values (o and oo)
  if(ini  < 2) then
    do c=1,NC
!----- U
      U % oo(c)  = U % o(c)
      U % o (c)  = U % n(c)
      U % Doo(c) = U % Do(c)
      U % Do (c) = 0.0 
      U % Xoo(c) = U % Xo(c)
      U % Xo (c) = U % X(c) 
    end do
  end if

!----- new values
  do c=1,NC
    U % X(c)    = 0.0
  end do

!----------------------------------------------------!
!     Browse through all the faces, where else ?     !
!----------------------------------------------------!
  do s=1,NS      

    c1=SideC(1,s)
    c2=SideC(2,s) 

!----- velocities on "orthogonal" cell centers 
    if(c2  > 0) then
      PHIxS = f(s)*PHIx(c1) + (1.0-f(s))*PHIx(c2)
      PHIyS = f(s)*PHIy(c1) + (1.0-f(s))*PHIy(c2)
      PHIzS = f(s)*PHIz(c1) + (1.0-f(s))*PHIz(c2)
    else
      PHIxS = PHIx(c1) 
      PHIyS = PHIy(c1) 
      PHIzS = PHIz(c1) 
    endif

    VISeff = VISc 
		      
!==================!
!                  !
!     Difusion     !
!                  !
!==================!

!--------------------------------!
!     Spatial Discretization     !
!--------------------------------!

!---- total (exact) diffusive flux
    FUex = VISeff*(PHIxS*Sx(s)+PHIyS*Sy(s)+PHIzS*Sz(s))

!---- implicit diffusive flux
    FUim = VISeff*Scoef(s)*                                         &
	    (  PHIxS*Dx(s)                                          &
	     + PHIyS*Dy(s)                                          &
	     + PHIzS*Dz(s) )

    A0 = VISeff * Scoef(s)

!---- straight diffusion part 
    if(ini  < 2) then
      if(c2  > 0) then
	U % Do(c1) = U % Do(c1) + VISeff*Scoef(s)*(U % n(c2)-U % n(c1)) 
	U % Do(c2) = U % Do(c2) - VISeff*Scoef(s)*(U % n(c2)-U % n(c1))   
      else
	if(TypeBC(c2) /= SYMMETRY) then 
	  U % Do(c1) = U % Do(c1) + VISeff*Scoef(s)*(U % n(c2)-U % n(c1))   
	end if
      end if 
    end if

!---- cross diffusion part
    U % X(c1) = U % X(c1) + FUex - FUim 
    if(c2  > 0) then
      U % X(c2) = U % X(c2) - FUex + FUim 
    end if 

!----- calculate the coefficients for the sysytem matrix
    if( (DIFFUS == CN) .or. (DIFFUS == FI) ) then     

      if(DIFFUS  ==  CN) then       ! Crank Nicholson
	A12 = 0.5 * A0 
	A21 = 0.5 * A0 
      end if

      if(DIFFUS  ==  FI) then       ! Fully implicit
	A12 = A0 
	A21 = A0
      end if

!----- fill the system matrix
      if(c2  > 0) then
	Aval(SidAij(1,s)) = Aval(SidAij(1,s)) - A12
	Aval(Adia(c1))    = Aval(Adia(c1))    + A12
	Aval(SidAij(2,s)) = Aval(SidAij(2,s)) - A21
	Aval(Adia(c2))    = Aval(Adia(c2))    + A21
      else if(c2  < 0) then
	if( (TypeBC(c2) == INFLOW).or.                                &
	    (TypeBC(c2) == WALL).or.                                  &
	    (TypeBC(c2) == WALLFL).or.                                &
	    (TypeBC(c2) == OUTFLOW).or.                               &
	    (TypeBC(c2) == BUFFER) ) then  
!
! Outflow is included because of the flux corrections which
! also affects velocities
!
	  Aval(Adia(c1)) = Aval(Adia(c1)) + A12
	  b(c1)  = b(c1)  + A12 * U % n(c2)
!->>>         write(*,'(3I7,6F12.6)') s, c1,c2,A12, U(c2)
	  Abou(c2) = -A12  ! cool new stuff
	endif
      end if       

    end if

  end do  ! through sides

!---------------------------------!
!     Temporal discretization     !
!---------------------------------!

!----- Adams-Bashfort scheeme for diffusion fluxes
  if(DIFFUS == AB) then 
    do c=1,NC
      b(c)  = b(c) + 1.5 * U % Do(c) - 0.5 * U % Doo(c)
    end do  
  end if

!----- Crank-Nicholson scheme for difusive terms
  if(DIFFUS == CN) then 
    do c=1,NC
      b(c)  = b(c) + 0.5 * U % Do(c)
    end do  
  end if
		  
!----- Fully implicit treatment for difusive terms
!      is handled via the linear system of equations 

!----- Adams-Bashfort scheeme for cross diffusion 
  if(CROSS == AB) then
    do c=1,NC
      b(c)  = b(c) + 1.5 * U % Xo(c) - 0.5 * U % Xoo(c)
    end do 
  end if

!----- Crank-Nicholson scheme for cross difusive terms
  if(CROSS == CN) then
    do c=1,NC
      b(c)  = b(c) + 0.5 * U % X(c) + 0.5 * U % Xo(c)
    end do 
  end if

!----- Fully implicit treatment for cross difusive terms
  if(CROSS == FI) then
    do c=1,NC
      if(U % X(c) >= 0.0) then
        b(c)  = b(c) + U % X(c)
      else
        Aval(Adia(c)) = Aval(Adia(c)) - U % X(c)/(U % n(c))
      end if
    end do 
  end if

!========================!
!                        !
!     Inertial terms     !
!                        !
!========================!

!----- Two time levels; Linear interpolation
  if(INERT == LIN) then
    do c=1,NC
      A0 = DENc(material(c))*volume(c)/dt
      Aval(Adia(c)) = Aval(Adia(c)) + A0
      b(c)  = b(c) + A0 * U % o(c)
    end do
  end if

!----- Three time levels; parabolic interpolation
  if(INERT == PAR) then
    do c=1,NC
      A0 = DENc(material(c))*volume(c)/dt
      Aval(Adia(c)) = Aval(Adia(c)) + 1.5 * A0
      b(c)  = b(c) + 2.0*A0 * U % o(c) - 0.5*A0 * U % oo(c)
    end do
  end if

!->>> take a look at the system of equations
!->>> do c=1,NC
!->>>   write(*,*) 'Width: ', Acol(c+1)-Acol(c)
!->>>   write(*,'(3I7)') Acol(c), Adia(c), Acol(c+1)-1
!->>>   write(*,*) 'Diag: ', Aval(Adia(c))
!->>>   write(*,'(F5.2)') ( Aval(j),  j=Acol(c),Acol(c+1)-1 )
!->>>   write(*,*) '- - - - - - - - - - - - - - - - - - - - - - -'
!->>> end do  

!===================================!
!                                   !
!     Solve the equations for U     !
!                                   !    
!===================================!
  if( (DIFFUS == AB) .and. (CONVEC == AB) ) then
    do c=1,NC
      U % n(c) = b(c)/Aval(Adia(c))
    end do
  else
    niter=25
    call cg(NC, NbC, NONZERO, Aval,Acol,Arow,Adia,Abou, &
            U % n, b, PREC,                             &
            niter, U % STol, res(1), error)
    write(LineRes(17:28), '(1PE12.3)') res(1)
    write(LineRes(77:80), '(I4)')      niter                                                              
  endif

  END SUBROUTINE TestU
