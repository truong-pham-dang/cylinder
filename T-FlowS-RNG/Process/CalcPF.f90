!======================================================================!
  SUBROUTINE CalcPF
!----------------------------------------------------------------------!
!   Forms and solves pressure equation for the fractional step method. !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: s, c, c1, c2, niter
  REAL    :: Pmax, Pmin
  REAL    :: error
  REAL    :: Us, Vs, Ws, DENs
  REAL    :: A12
!--------------------------------[CVS]---------------------------------!
!  $Id: CalcPF.f90,v 1.19 2008/12/05 13:22:27 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CalcPF.f90,v $  
!----------------------------------------------------------------------!
!     
!  The form of equations which I am solving:    
!     
!     /               /            
!    |               |             
!    | rho u dS = dt | GRAD pp dS
!    |               |             
!   /               /              
!
!  Dimension of the system under consideration
!   
!     [App] {pp} = {bpp}               [kg/s]
!   
!  Dimensions of certain variables
!
!     APP            [ms]
!     PP,            [kg/ms^2]
!     b              [kg/s]
!     Flux           [kg/s]
!   
!======================================================================!

  Aval = 0.0

!------------------------------------------!
!      Initialize the source term for      !
!     the pressure correction equation     !
!------------------------------------------!
  b = 0.0 

!-----------------------------------------------------!
!     Calculate the mass fluxes on the cell faces     !
!-----------------------------------------------------!
  do s=1, NS
    c1=SideC(1,s)
    c2=SideC(2,s)

!---- handle two materials
    if( StateMat(material(c1))==FLUID .and.      &
        StateMat(material(c2))==FLUID) then
      DENs =      f(s)  * DENc(material(c1))     &
           + (1.0-f(s)) * DENc(material(c2))
    else if( StateMat(material(c1))==FLUID .and. &
             StateMat(material(c2))==SOLID) then
      DENs = DENc(material(c1)) 
    else if( StateMat(material(c1))==SOLID .and. &
             StateMat(material(c2))==FLUID) then
      DENs = DENc(material(c2)) 
    else
      DENs =      f(s)  * DENc(material(c1))     &
           + (1.0-f(s)) * DENc(material(c2))
    end if  

!---- side is inside the domain
    if( c2  > 0 .or. c2  < 0 .and. TypeBC(c2) == BUFFER) then 

!---- extract the "centred" pressure terms from cell velocities
      Us = f(s)*      (U % n(c1)+Px(c1)*volume(c1)/Asave(c1))       &
         + (1.0-f(s))*(U % n(c2)+Px(c2)*volume(c2)/Asave(c2))

      Vs = f(s)*      (V % n(c1)+Py(c1)*volume(c1)/Asave(c1))       &
	 + (1.0-f(s))*(V % n(c2)+Py(c2)*volume(c2)/Asave(c2))

      Ws = f(s)*      (W % n(c1)+Pz(c1)*volume(c1)/Asave(c1))       &
	 + (1.0-f(s))*(W % n(c2)+Pz(c2)*volume(c2)/Asave(c2))

!---- add the "staggered" pressure terms to face velocities
      Us= Us + (P % n(c1)-P % n(c2))*Sx(s)*                         &
	 ( f(s)/Asave(c1) + (1.0-f(s))/Asave(c2) )
      Vs=Vs+(P % n(c1)-P % n(c2))*Sy(s)*                            &
	 ( f(s)/Asave(c1) + (1.0-f(s))/Asave(c2) )
      Ws=Ws+(P % n(c1)-P % n(c2))*Sz(s)*                            &
	 ( f(s)/Asave(c1) + (1.0-f(s))/Asave(c2) )

!---- now calculate the flux through cell face
      Flux(s) = DENs * ( Us*Sx(s) + Vs*Sy(s) + Ws*Sz(s) )

      A12=DENs*(Sx(s)*Sx(s)+Sy(s)*Sy(s)+Sz(s)*Sz(s))
      A12=A12*(f(s)/Asave(c1)+(1.-f(s))/Asave(c2))

      if(c2  > 0) then 
	Aval(SidAij(1,s))  = -A12
	Aval(SidAij(2,s))  = -A12
	Aval(Adia(c1)) = Aval(Adia(c1)) +  A12
	Aval(Adia(c2)) = Aval(Adia(c2)) +  A12
      else
	Abou(c2) = -A12
	Aval(Adia(c1)) = Aval(Adia(c1)) +  A12
      endif

      b(c1)=b(c1)-Flux(s)
      if(c2  > 0) b(c2)=b(c2)+Flux(s)

!----- side is on the boundary
    else
      Us = U % n(c2)
      Vs = V % n(c2)
      Ws = W % n(c2)

      Flux(s) = DENs * ( Us*Sx(s) + Vs*Sy(s) + Ws*Sz(s) )

      b(c1) = b(c1)-Flux(s)
    end if

  end do

!--------------------------------------------!
!     Initialize the pressure correction     !
!--------------------------------------------!
  PP % n = 0.0 

  errmax=0.0
  do c=1,NC
!        errmax=max(errmax, abs(b(c)))
    errmax=errmax + abs(b(c))
  end do
  call glosum(errmax)                       
!>>>>>write(*,*) 'debalans before: ', errmax
!>>>>>write(*,*) 'mass before', errmax !, u(cm), w(cm), p(cm)

!------------------------------------------------!
!     Solve the pressure correction equation     !
!------------------------------------------------!
!>>>>> Aval(1) = Aval(1) * 2.0

!---- Give the "false" flux back and set it to zero ! 2mat
  do s=1,NS                                         ! 2mat
    c1=SideC(1,s)                                   ! 2mat
    c2=SideC(2,s)                                   ! 2mat
    if(c2>0 .or. c2<0.and.TypeBC(c2)==BUFFER) then  ! 2mat
      if(StateMat(material(c1))==SOLID .or. &       ! 2mat
         StateMat(material(c2))==SOLID) then        ! 2mat
        b(c1) = b(c1) + Flux(s)                     ! 2mat
        if(c2 > 0) b(c2) = b(c2) - Flux(s)          ! 2mat
        Flux(s) = 0.0                               ! 2mat
      end if                                        ! 2mat
    end if                                          ! 2mat
  end do                                            ! 2mat

!---- Disconnect the SOLID cells from FLUID system  ! 2mat
  do s=1,NS                                         ! 2mat
    c1=SideC(1,s)                                   ! 2mat
    c2=SideC(2,s)                                   ! 2mat
    if(c2>0 .or. c2<0.and.TypeBC(c2)==BUFFER) then  ! 2mat 
      if(c2 > 0) then ! => not BUFFER               ! 2mat
        if(StateMat(material(c1)) == SOLID) then    ! 2mat
          A12 = -Aval(SidAij(2,s))                  ! 2mat
          Aval(SidAij(1,s)) = 0.0                   ! 2mat
          Aval(SidAij(2,s)) = 0.0                   ! 2mat
          if(StateMat(material(c2)) == FLUID) then  ! 2mat
	    Aval(Adia(c2)) = Aval(Adia(c2)) -  A12  ! 2mat
          endif                                     ! 2mat
        end if                                      ! 2mat
        if(StateMat(material(c2)) == SOLID) then    ! 2mat
          A12 = -Aval(SidAij(1,s))                  ! 2mat
          Aval(SidAij(2,s)) = 0.0                   ! 2mat
          Aval(SidAij(1,s)) = 0.0                   ! 2mat
          if(StateMat(material(c1)) == FLUID) then  ! 2mat
	    Aval(Adia(c1)) = Aval(Adia(c1)) -  A12  ! 2mat
          endif                                     ! 2mat
        end if                                      ! 2mat
      else            ! => BUFFER                   ! 2mat
        if(StateMat(material(c1)) == SOLID  .or. &  ! 2mat
           StateMat(material(c2)) == SOLID) then    ! 2mat
          A12 = -Abou(c2)                           ! 2mat
          Abou(c2) = 0.0                            ! 2mat
          if(StateMat(material(c1)) == FLUID) then  ! 2mat
	    Aval(Adia(c1)) = Aval(Adia(c1)) -  A12  ! 2mat
          endif                                     ! 2mat
        end if                                      ! 2mat
      end if                                        ! 2mat
    end if                                          ! 2mat
  end do                                            ! 2mat

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!     Don't solve the pressure corection too accurate.      !
!     Value 1.e-18 blows the solution.                      !
!     Value 1.e-12 keeps the solution stable                !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!  
  if(ALGOR == FRACT)  niter = 200
  if(ALGOR == SIMPLE) niter =  15
  call cg(NC, NbC, NONZERO, Aval, Acol,Arow,Adia,Abou,   &
	  PP % n, b, PREC, niter, PP % STol,             &
	  res(4), error) 
  write(LineRes(53:64),  '(1PE12.3)') res(4)
  write(LineRes(89:92),  '(I4)')      niter

!-----------------------------------!
!     Update the pressure field     !
!-----------------------------------!
  P % n  =  P % n  +  P % URF  *  PP % n

!--------------------------------------!
!     Normalize the pressure field     !
!--------------------------------------!
  Pmax  = maxval(P % n(1:NC))
  Pmin  = minval(P % n(1:NC))

  call glomax(Pmax) 
  call glomin(Pmin) 

  P % n  =  P % n  -  0.5 * (Pmax+Pmin)

  call Exchng(PP % n) 

  RETURN 

  END SUBROUTINE CalcPF
