!======================================================================!
  SUBROUTINE UserCalcS(var, PHI, dPHIdx, dPHIdy, dPHIdz)
!----------------------------------------------------------------------!
! Purpose:                                                             !
! ~~~~~~~~                                                             !
!   Discretizes and solves energy equations, but under special         !
!   circumstances: It removes periodicity in the streamwise direction  !
!   of the flow. It was developed for the "Matrix of cubes" test case  !
!   but might be usefull in similar problems.                          !
!                                                                      !
! Assumptions:                                                         !
! ~~~~~~~~~~~~                                                         !
!   1. The fluid flows in the positive x direction.                    !
!   2. On the "inflow" zero temperature is prescribed.                 !
!   3. On the "outflow" zero gradient is prescribed.                   !
!                                                                      !
!  Authors: Adriaan Dronkers and Bojan Niceno                          !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Parameters]---------------------------!
  INTEGER       :: var
  TYPE(Unknown) :: PHI
  REAL          :: dPHIdx(-NbC:NC), dPHIdy(-NbC:NC), dPHIdz(-NbC:NC)
!-------------------------------[Locals]-------------------------------!
  INTEGER :: n,c,s,c1,c2,niter,miter
  REAL    :: A0,A12,A21,error
  REAL    :: CONeff1, FUex1, FUim1, PHIxS1, PHIyS1, PHIzS1
  REAL    :: CONeff2, FUex2, FUim2, PHIxS2, PHIyS2, PHIzS2
  REAL    :: Stot, PHIs1, PHIs2
!--------------------------------[CVS]---------------------------------!
!  $Id: UserCalcSc.f90,v 1.13 2002/10/29 13:11:23 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserCalcSc.f90,v $  
!----------------------------------------------------------------------!
!     
!  The form of equations which are solved:    
!     
!     /             /              /            
!    |     dT      |              | lambda      
!    | rho -- dV   | rho u T dS = | ------  DIV T dS 
!    |     dt      |              |   Cp        
!   /             /              /              
!
!
!  Dimension of the system under consideration
!
!     [A]{T} = {b}   [kg K/s]   
!
!  Dimensions of certain variables:
!
!     Cp     [J/kg K]
!     lambda [W/m K]
!
!     A      [kg/s]
!     T      [K]
!     b      [kg K/s] 
!     Flux   [kg/s]
!     CT*,   [kg K/s] 
!     DT*,   [kg K/s] 
!     XT*,   [kg K/s]
! 
!======================================================================!

  do n=1,Acol(NC+1) ! to je broj nonzero + 1
    Aval(n) = 0.0
  end do

  do c=1,NC
    b(c)=0.0
  end do

!----- This is important for "copy" boundary conditions. Find out why !
!----- (I just coppied this from NewUVW.f90. I don't have a clue if it
!----- is of any importance at all. Anyway, I presume it can do no harm.)
  do c=-NbC,-1
    Abou(c)=0.0
  end do

!-----------------------------------------! 
!     Initialize variables and fluxes     !
!-----------------------------------------! 

!----- old values (o and oo)
  if(ini.lt.2) then
    do c=1,NC
      PHI % oo(c)  = PHI % o(c)
      PHI % o (c)  = PHI % n(c)
      PHI % Coo(c) = PHI % Co(c)
      PHI % Co (c) = 0.0
      PHI % Doo(c) = PHI % Do(c)
      PHI % Do (c) = 0.0 
      PHI % Xoo(c) = PHI % Xo(c)
      PHI % Xo (c) = PHI % X(c) 
    end do
  end if

!----- new values
  do c=1,NC
    PHI % C(c)    = 0.0
    PHI % X(c)    = 0.0  ! use PHI % X for upwind convective fluxes
  end do

!----------------------------------------------------!
!     Browse through all the faces, where else ?     !
!----------------------------------------------------!
  do s=1,NS

    c1=SideC(1,s)
    c2=SideC(2,s)

!====================!
!                    !
!     Convection     !
!                    !
!====================!
  
!----- velocities on "orthogonal" cell centers
    if(c2>0 .or. c2<0.and.TypeBC(c2) == BUFFER) then
      PHIs1=f(s)*PHI % n(c1) + (1.0-f(s))*PHI % n(c2)
      PHIs2=PHIs1

!-> .......... Matrix ..........
      if( xc(c2) > xc(c1)+0.5 ) then  ! c2 is "outflow"; c1 "inflow"
        PHIs1 = 0.0
        PHIs2 = PHI % n(c2) 
      end if
      if( xc(c1) > xc(c2)+0.5 ) then  ! c1 is "outflow"; c2 "inflow"
        PHIs1 = PHI % n(c1) 
        PHIs2 = 0.0 
      end if
!-> .......... Matrix ..........   

!---- Central differencing for convection
      if(ini.eq.1) then
        if(c2.gt.0) then
          PHI % Co(c1)=PHI % Co(c1)-Flux(s)*PHIs1
          PHI % Co(c2)=PHI % Co(c2)+Flux(s)*PHIs2
        else
          PHI % Co(c1)=PHI % Co(c1)-Flux(s)*PHIs1
        endif
      end if
      if(c2.gt.0) then
        PHI % C(c1)=PHI % C(c1)-Flux(s)*PHIs1
        PHI % C(c2)=PHI % C(c2)+Flux(s)*PHIs2
      else
        PHI % C(c1)=PHI % C(c1)-Flux(s)*PHIs1
      endif
 
!---- Upwind
      if(BLEND(material(c1)) /= NO .or. BLEND(material(c2)) /= NO) then
        if(Flux(s) < 0.0) then  ! from c2 to c1
          PHIs1=PHI % n(c2)
          PHIs2=PHI % n(c2)
        else
          PHIs1=PHI % n(c1)
          PHIs2=PHI % n(c1)
        end if

!-> .......... Matrix ..........
        if( xc(c2) > xc(c1)+0.5 ) then  ! c2 is "outflow"; c1 "inflow"
          PHIs1 = 0.0
          PHIs2 = PHI % n(c2)
        end if
        if( xc(c1) > xc(c2)+0.5 ) then  ! c1 is "outflow"; c2 "inflow"
          PHIs1 = PHI % n(c1)
          PHIs2 = 0.0 
        end if
!-> .......... Matrix ..........   

        if(Flux(s).lt.0) then   ! from c2 to c1
          PHI % X(c1)=PHI % X(c1)-Flux(s)*PHIs1
          if(c2.gt.0) then
            PHI % X(c2)=PHI % X(c2)+Flux(s)*PHIs2
          endif
        else
          PHI % X(c1)=PHI % X(c1)-Flux(s)*PHIs1
          if(c2.gt.0) then
            PHI % X(c2)=PHI % X(c2)+Flux(s)*PHIs2
          endif
        end if
      end if   ! BLEND
    else       ! c2 < 0
      PHI % C(c1) = PHI % C(c1)-Flux(s)*PHI % n(c2)
!---- Full upwind
      if(TypeBC(c2).eq.OUTFLOW) then
        if(BLEND(material(c1)).eq.YES) then
          if(Flux(s).lt.0) then   ! from c2 to c1
            PHI % X(c1)=PHI % X(c1)-Flux(s)*PHI % n(c2)
          else
            PHI % X(c1)=PHI % X(c1)-Flux(s)*PHI % n(c1)
          end if
        end if ! BLEND = YES
      end if   ! OUTFLOW
    end if     ! c2 > 0

  end do  ! through sides

!---------------------------------!
!     Temporal discretization     !
!---------------------------------!
   
!----- Adams-Bashforth scheeme for convective fluxes
  if(CONVEC.eq.AB) then
    do c=1,NC
      b(c) = b(c) + URFC(material(c)) * &
                    (1.5*PHI % Co(c) - 0.5*PHI % Coo(c) - PHI % X(c))
    end do
  endif

!----- Crank-Nicholson scheeme for convective fluxes
  if(CONVEC.eq.CN) then
    do c=1,NC
      b(c) = b(c) + URFC(material(c)) * & 
                    (0.5 * ( PHI % C(c) + PHI % Co(c) ) - PHI % X(c))
    end do
  endif

!----- Fully implicit treatment of convective fluxes
  if(CONVEC.eq.FI) then
    do c=1,NC
      b(c) = b(c) + URFC(material(c)) * & 
                    ( PHI % C(c) - PHI % X(c) )
    end do
  end if

!==================!
!                  !
!     Difusion     !
!                  !
!==================!

!----- Set PHI % X back to zero 
  do c=1,NC
    PHI % X(c) = 0.0  
  end do

!--------------------------------!
!     Spatial Discretization     !
!--------------------------------!

  do s=1,NS

    c1=SideC(1,s)
    c2=SideC(2,s)

!----- gradients on the cell face 
    if(c2  > 0 .or. c2  < 0.and.TypeBC(c2) == BUFFER) then
      if(material(c1) == material(c2)) then
        PHIxS1 = f(s)*dPHIdx(c1) + (1.0-f(s))*dPHIdx(c2)
        PHIyS1 = f(s)*dPHIdy(c1) + (1.0-f(s))*dPHIdy(c2)
        PHIzS1 = f(s)*dPHIdz(c1) + (1.0-f(s))*dPHIdz(c2)
        PHIxS2 = PHIxS1 
        PHIyS2 = PHIyS1 
        PHIzS2 = PHIzS1 
        CONeff1 = f(s) * CONc(material(c1)) + (1.-f(s)) * CONc(material(c2)) 
        CONeff2 = CONeff1 
      else 
        PHIxS1 = dPHIdx(c1) 
        PHIyS1 = dPHIdy(c1) 
        PHIzS1 = dPHIdz(c1) 
        PHIxS2 = dPHIdx(c2) 
        PHIyS2 = dPHIdy(c2) 
        PHIzS2 = dPHIdz(c2) 
        CONeff1 = CONc(material(c1)) 
        CONeff2 = CONc(material(c2)) 
      end if
    else
      PHIxS1 = dPHIdx(c1) 
      PHIyS1 = dPHIdy(c1) 
      PHIzS1 = dPHIdz(c1) 
      PHIxS2 = PHIxS1 
      PHIyS2 = PHIyS1 
      PHIzS2 = PHIzS1 
      CONeff1 = CONc(material(c1)) 
      CONeff2 = CONeff1 
    endif

!---- total (exact) diffusive flux
    FUex1 = CONeff1*(PHIxS1*Sx(s)+PHIyS1*Sy(s)+PHIzS1*Sz(s))
    FUex2 = CONeff2*(PHIxS2*Sx(s)+PHIyS2*Sy(s)+PHIzS2*Sz(s))

!---- implicit diffusive flux
    FUim1 = CONeff1*Scoef(s)*    &
            (  PHIxS1*Dx(s)      &
             + PHIyS1*Dy(s)      &
             + PHIzS1*Dz(s) )
    FUim2 = CONeff2*Scoef(s)*    &
            (  PHIxS2*Dx(s)      &
             + PHIyS2*Dy(s)      &
             + PHIzS2*Dz(s) )

!---- straight diffusion part 
    if(ini.lt.2) then
      if(c2.gt.0) then
        if(material(c1) == material(c2)) then
          PHI % Do(c1) = PHI % Do(c1) + CONeff1*Scoef(s)*(PHI % n(c2)-PHI % n(c1)) 
          PHI % Do(c2) = PHI % Do(c2) - CONeff2*Scoef(s)*(PHI % n(c2)-PHI % n(c1))   
        else
          PHI % Do(c1) = PHI % Do(c1) + 2.*CONc(material(c1))*Scoef(s)*(PHIside(s)-PHI % n(c1)) 
          PHI % Do(c2) = PHI % Do(c2) - 2.*CONc(material(c2))*Scoef(s)*(PHI % n(c2)-PHIside(s))   
        end if
      else
        if(TypeBC(c2).ne.SYMMETRY) then 
          if(material(c1) == material(c2)) then
            PHI % Do(c1) = PHI % Do(c1) + CONeff1*Scoef(s)*(PHI % n(c2)-PHI % n(c1))   
          else
            PHI % Do(c1) = PHI % Do(c1) + 2.*CONc(material(c1))*Scoef(s)*(PHIside(s)-PHI % n(c1)) 
          end if
        end if
      end if 
    end if

!---- cross diffusion part
    PHI % X(c1) = PHI % X(c1) + FUex1 - FUim1 
    if(c2.gt.0) then
      PHI % X(c2) = PHI % X(c2) - FUex2 + FUim2 
    end if 

!----- calculate the coefficients for the sysytem matrix
    if( (DIFFUS.eq.CN) .or. (DIFFUS.eq.FI) ) then

      if(DIFFUS .eq. CN) then       ! Crank Nicholson
        if(material(c1) == material(c2)) then
          A12 = .5*CONeff1*Scoef(s)  
          A21 = .5*CONeff2*Scoef(s)  
        else
          A12 = CONc(material(c1))*Scoef(s)  
          A21 = CONc(material(c2))*Scoef(s)  
        end if
      end if

      if(DIFFUS .eq. FI) then       ! Fully implicit
        if(material(c1) == material(c2)) then
          A12 = CONeff1*Scoef(s)  
          A21 = CONeff2*Scoef(s)  
        else
          A12 = 2.*CONc(material(c1))*Scoef(s)  
          A21 = 2.*CONc(material(c2))*Scoef(s)  
        end if
      end if

      if(BLEND(material(c1)) /= NO .or. BLEND(material(c2)) /= NO) then
        A12 = A12  - min(Flux(s), 0.0)
        A21 = A21  + max(Flux(s), 0.0)
      endif
                
!----- fill the system matrix
      if(c2.gt.0) then
	Aval(Adia(c1)) = Aval(Adia(c1)) + A12
	Aval(Adia(c2)) = Aval(Adia(c2)) + A21
        if(material(c1) == material(c2)) then
          Aval(SidAij(1,s)) = Aval(SidAij(1,s)) - A12
          Aval(SidAij(2,s)) = Aval(SidAij(2,s)) - A21
        else
          b(c1) = b(c1) + A12*PHIside(s)
          b(c2) = b(c2) + A21*PHIside(s)
        end if
!-> .......... Matrix ..........
        if( xc(c2) > xc(c1)+0.5 ) then  ! c2 is "outflow"; c1 "inflow"
          Aval(SidAij(2,s)) = 0.0
          Aval(SidAij(1,s)) = 0.0
          Aval(Adia(c2))    = Aval(Adia(c2))    - A21
        end if
        if( xc(c1) > xc(c2)+0.5 ) then  ! c1 is "outflow"; c2 "inflow"
          Aval(SidAij(1,s)) = 0.0
          Aval(SidAij(2,s)) = 0.0
          Aval(Adia(c1))    = Aval(Adia(c1))    - A12
        end if
!-> .......... Matrix ..........   
      else if(c2.lt.0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Outflow is included because of the flux corrections which
! also affects velocities
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if( (TypeBC(c2).eq.INFLOW).or.    &
            (TypeBC(c2).eq.WALL).or.      &
            (TypeBC(c2).eq.OUTFLOW) ) then    
          Aval(Adia(c1)) = Aval(Adia(c1)) + A12
          b(c1)  = b(c1)  + A12 * PHI % n(c2)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Buffer: System matrix and parts belonging to other
! subdomains are filled here.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        else if(TypeBC(c2).eq.BUFFER) then
          Aval(Adia(c1)) = Aval(Adia(c1)) + A12
          if(material(c1) == material(c2)) then
            Abou(c2) = -A12  ! cool parallel stuff
          else
            b(c1) = b(c1) + A12*PHIside(s)
          end if
!-> .......... Matrix ..........
          if( xc(c2) > xc(c1) + 0.5 ) then
            Abou(c2) = 0
          end if
          if( xc(c1) > xc(c2) + 0.5 ) then
            Aval(Adia(c1)) = Aval(Adia(c1)) - A12
            Abou(c2) = 0
          end if
!-> .......... Matrix ..........     
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! In case of wallflux 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        else if(TypeBC(c2).eq.WALLFL) then
          Stot  = sqrt(Sx(s)*Sx(s)+Sy(s)*Sy(s)+Sz(s)*Sz(s))
          b(c1) = b(c1) + Stot * PHI % q(c2)
        endif
      end if

    end if

  end do  ! through sides

!---------------------------------*
!     Temporal discretization     *
!---------------------------------*

!----- Adams-Bashfort scheeme for diffusion fluxes
  if(DIFFUS.eq.AB) then 
    do c=1,NC
      b(c)  = b(c) + 1.5 * PHI % Do(c) - 0.5 * PHI % Doo(c)
    end do  
  end if

!----- Crank-Nicholson scheme for difusive terms
  if(DIFFUS.eq.CN) then 
    do c=1,NC
      b(c)  = b(c) + 0.5 * PHI % Do(c)
    end do  
  end if

!----- Fully implicit treatment for difusive terms
!      is handled via the linear system of equations 

!----- Adams-Bashfort scheeme for cross diffusion 
  if(CROSS.eq.AB) then
    do c=1,NC
      b(c)  = b(c) + 1.5*PHI % Xo(c) - 0.5*PHI % Xoo(c)
    end do 
  end if

!----- Crank-Nicholson scheme for cross difusive terms
  if(CROSS.eq.CN) then
    do c=1,NC
      if( (PHI % X(c)+PHI % Xo(c))  >= 0) then
        b(c)  = b(c) + 0.5*(PHI % X(c) + PHI % Xo(c))
      else
        Aval(Adia(c)) = Aval(Adia(c)) &
             - 0.5 * (PHI % X(c) + PHI % Xo(c)) / (PHI % n(c)+1.e-6)
      end if
    end do 
  end if

!----- Fully implicit treatment for cross difusive terms
  if(CROSS.eq.FI) then
    do c=1,NC
      if(PHI % X(c) >= 0) then
        b(c)  = b(c) + PHI % X(c)
      else
        Aval(Adia(c)) = Aval(Adia(c)) - PHI % X(c)/(PHI % n(c)+1.e-6)
      end if
    end do 
  end if

!========================*
!                        *
!     Inertial terms     *
!                        *
!========================*

!----- Two time levels; Linear interpolation
  if(INERT.eq.LIN) then
    do c=1,NC
      A0 = volume(c)/dt
      Aval(Adia(c)) = Aval(Adia(c)) + A0
      b(c)  = b(c) + A0 * PHI % o(c)
    end do
  end if

!----- Three time levels; parabolic interpolation
  if(INERT.eq.PAR) then
    do c=1,NC
      A0 = volume(c)/dt
      Aval(Adia(c)) = Aval(Adia(c)) + 1.5 * A0
      b(c)  = b(c) + 2.0 * A0 * PHI % o(c) - 0.5 * A0 * PHI % oo(c)
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
!     Solve the equations for PHI   !
!                                   !    
!===================================!
  do c=1,NC
    b(c) = b(c) + Aval(Adia(c))*(1.0 - T % URF) * PHI % n(c) / T % URF
    Aval(Adia(c)) = Aval(Adia(c)) / T % URF
    Asave(c) = Aval(Adia(c))
  end do  

  if(ALGOR == SIMPLE)   miter=10
  if(ALGOR == FRACT)    miter=5 

  niter=miter
  call bicg(NC, Nbc, NONZERO, Aval,Acol,Arow,Adia,Abou,  &
            PHI % n, b, PREC,                            &
            niter,PHI % STol, res(var), error)
  write(*,*) 'niter =', niter
 
  call Exchng(PHI)

  END SUBROUTINE UserCalcS
