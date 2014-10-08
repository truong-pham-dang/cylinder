!======================================================================!
  SUBROUTINE NewUVW(var, Ui,                                &
		    dUidi, dUidj, dUidk,                            &
		    Si, Sj, Sk,                                     &
		    Di, Dj, Dk,                                     &
		    Pi, dUjdi, dUkdi)
!----------------------------------------------------------------------!
!   Discretizes and solves momentum conservation equations             !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE rans_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER       :: var
  TYPE(Unknown) :: Ui
  REAL          :: dUidi(-NbC:NC), dUidj(-NbC:NC), dUidk(-NbC:NC)
  REAL          :: Si(NS), Sj(NS), Sk(NS) 
  REAL          :: Di(NS), Dj(NS), Dk(NS) 
  REAL          :: Pi(-NbC:NC), dUjdi(-NbC:NC), dUkdi(-NbC:NC) 
  REAL          :: uuS, vvS, wwS, uvS, uwS, vwS
!-------------------------------[Locals]-------------------------------!
  INTEGER :: s, c, c1, c2, niter, miter, mat
  REAL    :: Fex, Fim 
  REAL    :: Uis
  REAL    :: A0, A12, A21
  REAL    :: error
  REAL    :: VISeff, VIStS, Fstress 
  REAL    :: dUidiS,dUidjS,dUidkS,dUjdiS,dUkdiS
!--------------------------------[CVS]---------------------------------!
!  $Id: NewUVW.f90,v 1.44 2009/06/30 12:10:38 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/NewUVW.f90,v $    
!----------------------------------------------------------------------!
!
!  Stress tensor on the face s:
!
!    T = mi * [    2*dU/dx     dU/dy+dV/dx   dU/dz+dW/dx  ]
!             [  dU/dy+dV/dx     2*dV/dy     dV/dz+dW/dy  ]
!             [  dU/dz+dW/dx   dV/dz+dW/dy     2*dW/dz    ]
!
!  The forces, acting on the cell face are:
!  
!    Fx = T11*Sx + T12*Sy + T13*Sz
!    Fy = T21*Sx + T22*Sy + T23*Sz
!    Fz = T31*Sx + T32*Sy + T33*Sz
!
!  which could also be written in the compact form:
!
!    {F} = [T]{S}
!
!  U component:
!
!    Fx = Txx*Sx + Txy*Sy + Txz*Sz 
!
!  V component:    
!
!    Fy = Tyx*Sx + Tyy*Sy + Tyz*Sz   
!
!  W component:
!
!    Fz = Tzx*Sx + Tzy*Sy + Tzz*Sz
!
!----------------------------------------------------------------------!
!     
!  The form of equations which I am solving:    
!     
!     /             /              /               /
!    |     du      |              |               |
!    | rho -- dV + | rho u u dS = | mu DIV u dS - | p dS
!    |     dt      |              |               |
!   /             /              /               /
!
!  Dimension of the system under consideration
!
!     [A]{u} = {b}   [kgm/s^2]   [N]
!
!  Dimensions of certain variables
!
!     A              [kg/s]
!     U, V, W        [m/s]
!     bU, bV, bW     [kgm/s^2]   [N]
!     P, PP          [kg/ms^2]   [N/m^2]
!     Flux           [kg/s]
!     CU*, CV*, CW*  [kgm/s^2]   [N]
!     DU*, DV*, DW*  [kgm/s^2]   [N]
!     XU*, XV*, XW*  [kgm/s^2]   [N]
!
!======================================================================!

  b = 0.0
  Aval = 0.0

!----- This is important for "copy" boundary conditions. Find out why !
  do c=-NbC,-1
    Abou(c)=0.0
  end do

!-----------------------------------------! 
!     Initialize variables and fluxes     !
!-----------------------------------------! 

!----- old values (o and oo)
  if(ini == 1) then
    do c=1,NC
      Ui % oo(c)  = Ui % o(c)
      Ui % o (c)  = Ui % n(c)
      Ui % Coo(c) = Ui % Co(c)
      Ui % Co (c) = 0.0 
      Ui % Doo(c) = Ui % Do(c)
      Ui % Do (c) = 0.0 
      Ui % Xoo(c) = Ui % Xo(c)
      Ui % Xo (c) = Ui % X(c) 
    end do
  end if

!====================!
!                    !
!     Convection     !
!                    !
!====================!

!----- Compute PHImax and PHImin
  do mat=1,Nmat
    if(BLEND(mat) /= NO) then
      call CalMinMax(Ui % n)  ! or Ui % o ???
      goto 1
    end if
  end do

!----- new values
1 do c=1,NC
    Ui % C(c)    = 0.0
    Ui % X(c)    = 0.0
  end do

!--------------------------------!
!     Spatial Discretization     !
!--------------------------------!

  do s=1,NS

    c1=SideC(1,s)
    c2=SideC(2,s) 

!---- Central differencing
    Uis=f(s)*Ui % n(c1) + (1.0-f(s))*Ui % n(c2)

    if(BLEND(material(c1)) /= NO .or. BLEND(material(c2)) /= NO) then
      call ConvScheme(Uis, s, Ui % n, dUidi, dUidj, dUidk, Di, Dj, Dk, &
                           max(BLEND(material(c1)),BLEND(material(c2))) ) 
    end if 

!---- Central differencing for convection 
    if(ini == 1) then 
      if(c2  > 0) then
        Ui % Co(c1)=Ui % Co(c1)-Flux(s)*Uis
        Ui % Co(c2)=Ui % Co(c2)+Flux(s)*Uis
      else
        Ui % Co(c1)=Ui % Co(c1)-Flux(s)*Uis
      endif 
    end if
    if(c2  > 0) then
      Ui % C(c1)=Ui % C(c1)-Flux(s)*Uis
      Ui % C(c2)=Ui % C(c2)+Flux(s)*Uis
    else
      Ui % C(c1)=Ui % C(c1)-Flux(s)*Uis
    endif 

!---- Upwind 
    if(BLEND(material(c1)) /= NO .or. BLEND(material(c2)) /= NO) then
      if(Flux(s)  < 0) then   ! from c2 to c1
        Ui % X(c1)=Ui % X(c1)-Flux(s)*Ui % n(c2)
        if(c2  > 0) then
          Ui % X(c2)=Ui % X(c2)+Flux(s)*Ui % n(c2)
        endif
      else 
        Ui % X(c1)=Ui % X(c1)-Flux(s)*Ui % n(c1)
        if(c2  > 0) then
          Ui % X(c2)=Ui % X(c2)+Flux(s)*Ui % n(c1)
        endif
      end if
    end if   ! BLEND 
  end do    ! through sides

!---------------------------------!
!     Temporal discretization     !
!---------------------------------!

!----- Adams-Bashforth scheeme for convective fluxes
  if(CONVEC == AB) then
    do c=1,NC
      b(c) = b(c) + URFC(material(c)) * & 
                    (1.5*Ui % Co(c) - 0.5*Ui % Coo(c) - Ui % X(c))
    end do  
  endif

!----- Crank-Nicholson scheeme for convective fluxes
  if(CONVEC == CN) then
    do c=1,NC
      b(c) = b(c) + URFC(material(c)) * & 
                    (0.5 * ( Ui % C(c) + Ui % Co(c) ) - Ui % X(c))
    end do  
  endif

!----- Fully implicit treatment of convective fluxes 
  if(CONVEC == FI) then
    do c=1,NC
      b(c) = b(c) + URFC(material(c)) * & 
                    (Ui % C(c) - Ui % X(c))
    end do  
  end if     
	    
!----------------------------------------------------!
!     Browse through all the faces, where else ?     !
!----------------------------------------------------!

!----- new values
  do c=1,NC
    Ui % X(c) = 0.0
  end do

!==================!
!                  !
!     Difusion     !
!                  !
!==================!

!--------------------------------!
!     Spatial Discretization     !
!--------------------------------!
  do s=1,NS       

    c1=SideC(1,s)
    c2=SideC(2,s)   

    VISeff = fF(s)*VISt(c1)+(1.0-fF(s))*VISt(c2) + VISc

    if(SIMULA==HYB_ZETA) VISeff = fF(s)*VISt_eff(c1)+(1.0-fF(s))*VISt_eff(c2) + VISc

    if(c2 < 0 .and. SIMULA == LES) then
      if(TypeBC(c2) == WALL .or. TypeBC(c2) == WALLFL) then
        VISeff = VISwall(c1)
      end if
    end if 

    if(SIMULA == ZETA.or.SIMULA==K_EPS_VV.or.(SIMULA == K_EPS.and.MODE==HRe).or.SIMULA==HYB_ZETA) then
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
        if(TypeBC(c2) == WALL .or. TypeBC(c2) == WALLFL) then
          VISeff = VISwall(c1)
        end if
      end if
    end if

!---- Add influence of Re stresses for EBM
    if(SIMULA == EBM) then
      if(var == 1) then
        uuS = fF(s)*uu % n(c1)+(1.0-fF(s))*uu % n(c2)
        uvS = fF(s)*uv % n(c1)+(1.0-fF(s))*uv % n(c2)
        uwS = fF(s)*uw % n(c1)+(1.0-fF(s))*uw % n(c2)
       Fstress = - (uuS*Sx(s)+uvS*Sy(s)+uwS*Sz(s))  
      else if(var == 2) then 
        uvS = fF(s)*uv % n(c1)+(1.0-fF(s))*uv % n(c2)
        vvS = fF(s)*vv % n(c1)+(1.0-fF(s))*vv % n(c2)
        vwS = fF(s)*vw % n(c1)+(1.0-fF(s))*vw % n(c2)
       Fstress =  - (uvS*Sx(s)+vvS*Sy(s)+vwS*Sz(s))  
      else if(var == 3) then 
        uwS = fF(s)*uw % n(c1)+(1.0-fF(s))*uw % n(c2)
        vwS = fF(s)*vw % n(c1)+(1.0-fF(s))*vw % n(c2)
        wwS = fF(s)*ww % n(c1)+(1.0-fF(s))*ww % n(c2)
       Fstress =  - (uwS*Sx(s)+vwS*Sy(s)+wwS*Sz(s))  
      end if 
    end if



    dUidiS = fF(s)*dUidi(c1) + (1.0-fF(s))*dUidi(c2)
    dUidjS = fF(s)*dUidj(c1) + (1.0-fF(s))*dUidj(c2)
    dUidkS = fF(s)*dUidk(c1) + (1.0-fF(s))*dUidk(c2)
    dUjdiS = fF(s)*dUjdi(c1) + (1.0-fF(s))*dUjdi(c2)
    dUkdiS = fF(s)*dUkdi(c1) + (1.0-fF(s))*dUkdi(c2)

!---- total (exact) diffusive flux
    Fex=VISeff*( 2.0*dUidiS*Si(s)                                 &
	          + (dUidjS+dUjdiS)*Sj(s)                         &
        	  + (dUidkS+dUkdiS)*Sk(s) )

    A0 = VISeff * Scoef(s)

!---- implicit diffusive flux
!.... this is a very crude approximation: Scoef is not
!.... corrected at interface between materials
    Fim=( dUidiS*Di(s)                                      &
	 +dUidjS*Dj(s)                                      &
	 +dUidkS*Dk(s))*A0

!---- this is yet another crude approximation:
!.... A0 is calculated approximatelly
!    if( StateMat(material(c1))==FLUID .and.  &  ! 2mat
!        StateMat(material(c2))==SOLID        &  ! 2mat
!        .or.                                 &  ! 2mat 
!        StateMat(material(c1))==SOLID .and.  &  ! 2mat
!        StateMat(material(c2))==FLUID ) then    ! 2mat
!      A0 = A0 + A0                              ! 2mat
!    end if                                      ! 2mat

!---- straight diffusion part 
    if(ini == 1) then
      if(c2  > 0) then
	Ui % Do(c1) = Ui % Do(c1) + (Ui % n(c2)-Ui % n(c1))*A0   
	Ui % Do(c2) = Ui % Do(c2) - (Ui % n(c2)-Ui % n(c1))*A0    
      else
	if(TypeBC(c2) /= SYMMETRY) then
	  Ui % Do(c1) = Ui % Do(c1) + (Ui % n(c2)-Ui % n(c1))*A0   
	end if 
      end if 
    end if

!---- cross diffusion part
    Ui % X(c1) = Ui % X(c1) + Fex - Fim + Fstress
    if(c2  > 0) then
      Ui % X(c2) = Ui % X(c2) - Fex + Fim - Fstress
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
    
      if(BLEND(material(c1)) /= NO .or. BLEND(material(c2)) /= NO) then
	A12 = A12  - min(Flux(s), real(0.0)) 
	A21 = A21  + max(Flux(s), real(0.0))
      endif

!----- fill the system matrix
      if(c2  > 0) then
	Aval(SidAij(1,s)) = Aval(SidAij(1,s)) - A12
	Aval(Adia(c1))    = Aval(Adia(c1))    + A12
	Aval(SidAij(2,s)) = Aval(SidAij(2,s)) - A21
	Aval(Adia(c2))    = Aval(Adia(c2))    + A21
      else if(c2  < 0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Outflow is not included because it was causing problems     !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -! 
	if( (TypeBC(c2) == INFLOW).or.                                &
	    (TypeBC(c2) == WALL).or.                                  &
	    (TypeBC(c2) == CONVECT).or.                               &
	    (TypeBC(c2) == WALLFL) ) then                                
!---------  (TypeBC(c2) == OUTFLOW) ) then   
	  Aval(Adia(c1)) = Aval(Adia(c1)) + A12
	  b(c1) = b(c1) + A12 * Ui % n(c2)
	else if( TypeBC(c2) == BUFFER ) then  
	  Aval(Adia(c1)) = Aval(Adia(c1)) + A12
	  Abou(c2) = - A12  ! cool parallel stuff
	endif
      end if     

    end if

  end do  ! through sides

!---------------------------------!
!     Temporal discretization     !
!---------------------------------!
!---------------------------------!
!
! Add Re stress influence on momentum
!
  if(SIMULA == EBM.and.MODE/=HYB) then
    if(var == 31) then
      call GraPhi(uu%n,1,VAR2x,.TRUE.)
      call GraPhi(uv%n,2,VAR2y,.TRUE.)
      call GraPhi(uw%n,3,VAR2z,.TRUE.)
      do c = 1, NC
        b(c) = b(c) - (VAR2x(c)+VAR2y(c)+VAR2z(c))*volume(c)
      end do
    else if(var == 32) then
      call GraPhi(uv%n,1,VAR2x,.TRUE.)
      call GraPhi(vv%n,2,VAR2y,.TRUE.)
      call GraPhi(vw%n,3,VAR2z,.TRUE.)
      do c = 1, NC
        b(c) = b(c) - (VAR2x(c)+VAR2y(c)+VAR2z(c))*volume(c)
      end do
    else if(var == 33) then
      call GraPhi(uw%n,1,VAR2x,.TRUE.)
      call GraPhi(vw%n,2,VAR2y,.TRUE.)
      call GraPhi(ww%n,3,VAR2z,.TRUE.)
      do c = 1, NC
        b(c) = b(c) - (VAR2x(c)+VAR2y(c)+VAR2z(c))*volume(c)
      end do
    end if
!
! Here we clean up momentum from the false diffusion
!
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)

      VIStS = (fF(s)*VISt(c1)+(1.0-fF(s))*VISt(c2))
      A0 = Scoef(s)*VIStS 
      VISeff = VIStS

      dUidiS = fF(s)*dUidi(c1) + (1.0-fF(s))*dUidi(c2)
      dUidjS = fF(s)*dUidj(c1) + (1.0-fF(s))*dUidj(c2)
      dUidkS = fF(s)*dUidk(c1) + (1.0-fF(s))*dUidk(c2)
      dUjdiS = fF(s)*dUjdi(c1) + (1.0-fF(s))*dUjdi(c2)
      dUkdiS = fF(s)*dUkdi(c1) + (1.0-fF(s))*dUkdi(c2)

      Fex=VISeff*( 2.0*dUidiS*Si(s)                                 &
             + (dUidjS+dUjdiS)*Sj(s)                         &
             + (dUidkS+dUkdiS)*Sk(s) )

      Fim=( dUidiS*Di(s)                                      &
          +dUidjS*Dj(s)                                      &
          +dUidkS*Dk(s))*VISeff*Scoef(s)

      b(c1) = b(c1) - VISeff*(Ui%n(c2)-Ui%n(c1))*Scoef(s)- Fex + Fim
      if(c2  > 0) then
        b(c2) = b(c2) + VISeff*(Ui%n(c2)-Ui%n(c1))*Scoef(s)+ Fex - Fim
      end if
    end do
  end if

!----- Adams-Bashfort scheeme for diffusion fluxes
  if(DIFFUS == AB) then 
    do c=1,NC
      b(c) = b(c) + 1.5 * Ui % Do(c) - 0.5 * Ui % Doo(c)
    end do  
  end if

!----- Crank-Nicholson scheme for difusive terms
  if(DIFFUS == CN) then 
    do c=1,NC
      b(c) = b(c) + 0.5 * Ui % Do(c)
    end do  
  end if
		 
!----- Fully implicit treatment for difusive terms
!      is handled via the linear system of equations 

!----- Adams-Bashfort scheeme for cross diffusion 
  if(CROSS == AB) then
    do c=1,NC
      b(c) = b(c) + 1.5 * Ui % Xo(c) - 0.5 * Ui % Xoo(c)
    end do 
  end if

!----- Crank-Nicholson scheme for cross difusive terms
  if(CROSS == CN) then
    do c=1,NC
      b(c) = b(c) + 0.5 * Ui % X(c) + 0.5 * Ui % Xo(c)
    end do 
  end if

!----- Fully implicit treatment for cross difusive terms
  if(CROSS == FI) then
    do c=1,NC
      b(c) = b(c) + Ui % X(c)
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
      b(c) = b(c) + A0 * Ui % o(c)
    end do
  end if

!----- Three time levels; parabolic interpolation
  if(INERT == PAR) then
    do c=1,NC
      A0 = DENc(material(c))*volume(c)/dt
      Aval(Adia(c)) = Aval(Adia(c)) + 1.5 * A0
      b(c) = b(c) + 2.0*A0 * Ui % o(c) - 0.5*A0 * Ui % oo(c)
    end do
  end if

!=====================================!
!                                     !
!     Pressure term contributions     !
!                                     !
!=====================================!

!------------------------------!
!     Global pressure drop     !
!------------------------------!
  if(var == 1) then
    do c=1,NC
      b(c) = b(c)  + PdropX(material(c)) * volume(c)
    end do
  else if(var == 2) then
    do c=1,NC
      b(c) = b(c)  + PdropY(material(c)) * volume(c)
    end do
  else if(var == 3) then
    do c=1,NC
      b(c) = b(c)  + PdropZ(material(c)) * volume(c)
    end do
  end if

!-------------------------------------!
!     Local pressure distribution     !
!-------------------------------------!
  do c=1,NC
    b(c) = b(c) - Pi(c)*volume(c)
  end do

!--------------------------------------------!
!     All other terms defined by the user    !
!--------------------------------------------!
  if(HOT == YES) call UserForce(var)

!  do c=1,NC                                         ! 2mat
!    if(StateMat(material(c))==SOLID) then           ! 2mat
!      b(c)          = 0.0                           ! 2mat
!    end if                                          ! 2mat
!  end do                                            ! 2mat

!---- Disconnect the SOLID cells from FLUID system  ! 2mat
!  do s=1,NS                                         ! 2mat
!    c1=SideC(1,s)                                   ! 2mat
!    c2=SideC(2,s)                                   ! 2mat
!    if(c2>0 .or. c2<0.and.TypeBC(c2)==BUFFER) then  ! 2mat
!      if(c2 > 0) then ! => not buffer               ! 2mat
!        if(StateMat(material(c1)) == SOLID) then    ! 2mat
!          Aval(SidAij(1,s)) = 0.0                   ! 2mat 
!          Aval(SidAij(2,s)) = 0.0                   ! 2mat
!        end if                                      ! 2mat 
!        if(StateMat(material(c2)) == SOLID) then    ! 2mat
!          Aval(SidAij(2,s)) = 0.0                   ! 2mat
!          Aval(SidAij(1,s)) = 0.0                   ! 2mat 
!        end if                                      ! 2mat 
!      else            ! => buffer region            ! 2mat 
!        if(StateMat(material(c1)) == SOLID .or.  &  ! 2mat
!           StateMat(material(c2)) == SOLID) then    ! 2mat
!          Abou(c2) = 0.0                            ! 2mat
!        end if                                      ! 2mat 
!      end if                                        ! 2mat
!    end if                                          ! 2mat
!  end do                                            ! 2mat

!=======================================!
!                                       !
!     Solve the equations for U,V,W     !
!                                       !    
!=======================================!
  do c=1,NC
    Asave(c) = Aval(Adia(c))
    b(c) = b(c) + Aval(Adia(c)) * (1.0-U % URF)*Ui % n(c) / U % URF
    Aval(Adia(c)) = Aval(Adia(c)) / U % URF
  end do

  if(ALGOR == SIMPLE)   miter=10
  if(ALGOR == FRACT)    miter=5

  niter=miter

  call cg(NC, Nbc, NONZERO, Aval,Acol,Arow,Adia,Abou,   &
	    Ui % n, b, PREC,                              &
	    niter,U % STol, res(var), error)

  if(var == 1) then
    write(LineRes(17:28), '(1PE12.3)') res(var) 
    write(LineRes(77:80), '(I4)')      niter 
  end if
  if(var == 2) then
    write(LineRes(29:40), '(1PE12.3)') res(var) 
    write(LineRes(81:84), '(I4)')      niter 
  end if
  if(var == 3) then
    write(LineRes(41:52), '(1PE12.3)') res(var) 
    write(LineRes(85:88), '(I4)')      niter 
  end if

  call Exchng(Ui % n)

  END SUBROUTINE NewUVW
