!======================================================================!
  SUBROUTINE IniCoDoXo(PHI)
!----------------------------------------------------------------------!
!  Initializes convective, diffusive and cross-diffusive terms from    !
!  previous calculations.                                              !
!----------------------------------------------------------------------!
!                                                                      ! 
!  The form of equations which are solved:                             !   
!                                                                      !
!     /               /                /                     /         !
!    |     dPHI      |                | mu_eff              |          !
!    | rho ---- dV + | rho u PHI dS = | ------ DIV PHI dS + | G dV     !
!    |      dt       |                |  sigma              |          !
!   /               /                /                     /           !
!                                                                      !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE par_mod
  USE les_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  TYPE(Unknown) :: PHI
!-------------------------------[Locals]-------------------------------!
  INTEGER :: s, c, c1, c2
  REAL    :: dPHIdx(-NbC:NC), dPHIdy(-NbC:NC), dPHIdz(-NbC:NC)
  REAL    :: PHIs, Fex, Fim, A0, VISeff, Sfac
  REAL    :: dPHIdxS, dPHIdyS, dPHIdzS
  REAL    :: xc1, xc2, yc1, yc2, zc1, zc2
  REAL    :: DeltaX, DeltaY, DeltaZ, Weigth
!--------------------------------[CVS]---------------------------------!
!  $Id: CalcTurb.f90,v 1.16 2003/01/21 11:29:37 mirza Exp $  
!  $Source: /home/mirza/.CVSROOT_mirza/T-Rex/Process/CalcTurb.f90,v $    
!======================================================================!


!-----Gradients
  call GraPhi(PHI % n, 1, dPHIdx, .TRUE.)
  call GraPhi(PHI % n, 2, dPHIdy, .TRUE.)
  call GraPhi(PHI % n, 3, dPHIdz, .TRUE.)

!-----Old values (o and oo)
!  if(ini == 1) then
!    do c=1,NC
!      PHI % oo(c)  = PHI % o(c)
!      PHI % o (c)  = PHI % n(c)
!!      PHI % Coo(c) = PHI % Co(c)
!      PHI % Co (c) = 0.0 
!!      PHI % Doo(c) = PHI % Do(c)
!      PHI % Do (c) = 0.0 
!!      PHI % Xoo(c) = PHI % Xo(c)
!      PHI % Xo (c) = PHI % X(c) 
!    end do
!  end if


!-----New values
  do c=1,NC
    PHI % C(c) = 0.0
    PHI % X(c) = 0.0
    PHI % Co(c) = 0.0
    PHI % Xo(c) = 0.0
    PHI % Do(c) = 0.0
  end do
 
!-----Browse through faces
  do s=1,NS
  c1=SideC(1,s)
  c2=SideC(2,s)

!--delta
    xc1=xc(c1)
    yc1=yc(c1)
    zc1=zc(c1) 

    xc2=xc(c2) + Dx(s)
    yc2=yc(c2) + Dy(s)
    zc2=zc(c2) + Dz(s)

    DeltaX = xc2-xc1  !Dx(s)
    DeltaY = yc2-yc1  !Dy(s)
    DeltaZ = zc2-zc1  !Dz(s)

!--weigth
    Weigth = f(s)                                 ! 2mat
                                                  ! 2mat 
    if( StateMat(material(c1))==FLUID .and.  &    ! 2mat
        StateMat(material(c2))==SOLID ) then      ! 2mat
      Weigth = 1.0                                ! 2mat
    end if                                        ! 2mat
                                                  ! 2mat
    if( StateMat(material(c1))==SOLID .and.  &    ! 2mat 
        StateMat(material(c2))==FLUID ) then      ! 2mat 
      Weigth = 0.0                                ! 2mat
    end if                                        ! 2mat
                                                  ! 2mat
    if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then    ! 2mat
      Weigth = 1.0                                ! 2mat
    end if                                        ! 2mat


!====================!
!     Convection     !
!====================!

!-----Velocities on "orthogonal" cell centers 
    if(c2 > 0 .or. c2 < 0.and.TypeBC(c2) == BUFFER) then
      PHIs = Weigth*PHI % n(c1) + (1.0-Weigth)*PHI % n(c2)

!-----Central differencing for convection 
      if(c2 > 0) then
        PHI % C(c1) = PHI % C(c1) - Flux(s)*PHIs
        PHI % C(c2) = PHI % C(c2) + Flux(s)*PHIs
        PHI % Co(c1) = PHI % Co(c1) - Flux(s)*PHIs
        PHI % Co(c2) = PHI % Co(c2) + Flux(s)*PHIs
      else
        PHI % C(c1) = PHI % C(c1) - Flux(s)*PHIs
        PHI % Co(c1) = PHI % Co(c1) - Flux(s)*PHIs
      end if 
    else       ! c2 < 0
!=================================================
    end if     ! c2 > 0 


!==================!
!     Difusion     !
!==================!

    VISeff = VISc + (Weigth*VISt(c1) + (1.0-Weigth)*VISt(c2)) !/ PHI%Sigma  !izbaceno

    if(SIMULA==SPA_ALL.or.SIMULA==DES_SPA)                                &
    VISeff = VISc + (Weigth*VIS%n(c1) + (1.0-Weigth)*VIS%n(c2)) !/ PHI%Sigma

    dPHIdxS = Weigth*dPHIdx(c1) + (1.0-Weigth)*dPHIdx(c2)
    dPHIdyS = Weigth*dPHIdy(c1) + (1.0-Weigth)*dPHIdy(c2)
    dPHIdzS = Weigth*dPHIdz(c1) + (1.0-Weigth)*dPHIdz(c2)

!-----This implements zero gradient for k
    if(c2 < 0 .and. (TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL)) then
      if(SIMULA==K_EPS .or. SIMULA==HYB_PITM .or. &
        (SIMULA==K_EPS_VV).or.(SIMULA==ZETA)) then
        dPHIdxS = 0.0 
        dPHIdyS = 0.0
        dPHIdzS = 0.0
        VISeff  = 0.0
      end if
    end if

!-----Total (exact) diffusive flux
    Fex = ( dPHIdxS * Sx(s)          &
          + dPHIdyS * Sy(s)          &
          + dPHIdzS * Sz(s) ) * VISeff

    Sfac = (Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s))  &
         / (DeltaX*Sx(s) + DeltaY*Sy(s) + DeltaZ*Sz(s))

    A0 = VISeff * Sfac

!---- implicit diffusive flux
    Fim = ( dPHIdxS * DeltaX      &
          + dPHIdyS * DeltaY      &
          + dPHIdzS * DeltaZ ) * A0

!---- two materials
    if( StateMat(material(c1))==FLUID .and.  &  ! 2mat
        StateMat(material(c2))==SOLID        &  ! 2mat
        .or.                                 &  ! 2mat 
        StateMat(material(c1))==SOLID .and.  &  ! 2mat
        StateMat(material(c2))==FLUID ) then    ! 2mat
      A0 = A0 + A0                              ! 2mat
    end if                                      ! 2mat

!-----Straight diffusion part 
    if(c2  > 0) then
      PHI % Do(c1) = PHI % Do(c1) + (PHI % n(c2)-PHI % n(c1))*A0   
      PHI % Do(c2) = PHI % Do(c2) - (PHI % n(c2)-PHI % n(c1))*A0    
    else
      if(TypeBC(c2) /= SYMMETRY) then
        PHI % Do(c1) = PHI % Do(c1) + (PHI % n(c2)-PHI % n(c1))*A0   
      end if 
    end if 

!---- cross diffusion part
    PHI % X(c1) = PHI % X(c1) + Fex - Fim 
    PHI % Xo(c1) = PHI % Xo(c1) + Fex - Fim 
    if(c2  > 0) then
      PHI % X(c2) = PHI % X(c2) - Fex + Fim 
      PHI % Xo(c2) = PHI % Xo(c2) - Fex + Fim 
    end if 
  end do

  
!  call Exchng(PHI)       !%n)  !IZMJENA this caused trouble

  RETURN
  END SUBROUTINE IniCoDoXo
