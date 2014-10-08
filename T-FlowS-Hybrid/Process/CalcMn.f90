!======================================================================!
  SUBROUTINE CalcMn(n0, n1,restart)   
!----------------------------------------------------------------------!
!   Calculates time averaged velocity and velocity fluctuations.       !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE les_mod
  USE pro_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: n0, n1
  LOGICAL :: restart
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, n,c1,c2,s
  REAL    :: Fpx, Fpy, Fvx, Fvy, Fvz, FtotalX, FtotalY
!--------------------------------[CVS]---------------------------------!
!  $Id: CalcMn.f90,v 1.17 2008/12/05 13:22:06 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CalcMn.f90,v $   
!======================================================================!

  Fpx = 0.0
  Fpy = 0.0
  Fvx = 0.0
  Fvy = 0.0
  Fvz = 0.0
  FtotalX = 0.0
  FtotalY = 0.0
  

  n=n1-n0

  if(n  > -1) then

  if(SIMULA == ZETA) then

   if (this < 2 .and. ( .not. restart)) write(765,'(A9,7A15)') "n1","Fpx","Fpy","Fvx","Fvy","Fvz","FtotalX","FtotalY"

    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      if(c2  < 0) then
        if(bcmark(c2) == 1) then !cylinders wall
            Fpx = Fpx + P%n(c1)*Sx(s)
            Fpy = Fpy + P%n(c1)*Sy(s)

            Fvx = Fvx + (VISc+VISt(c1))*(Ux(c1)*Sx(s)+Uy(c1)*Sy(s))
            Fvy = Fvy + (VISc+VISt(c1))*(Vx(c1)*Sx(s)+Vy(c1)*Sy(s))
            Fvz = Fvz + (VISc+VISt(c1))*(Wx(c1)*Sx(s)+Wy(c1)*Sy(s))

 
        end if
      end if
    end do



    FtotalX = 0.0 + Fpx + Fvx
    FtotalY = 0.0 + Fpy + Fvy


    call Wait
    call GloSum(Fpx)
    call GloSum(Fpy)
    call GloSum(Fvx)
    call GloSum(Fvy)
    call GloSum(Fvz)
    call GloSum(FtotalX)
    call GloSum(FtotalY)

    !if(this<2) write(*,*) 'in CalcForce; Fpy = ', Fpy

    if(this < 2 .and. mod(n1,10)==0) then !write every 10 steps
      write(765,'(I9,7E15.7)') n1, Fpx, Fpy, Fvx, Fvy, Fvz, FtotalX, FtotalY
    endif
   
    call Wait

    !!call CalcShear( U % mean, V % mean, W % mean, ShearMean) !sqrt(2<Sij><Sij>)

  end if

  do c=1,NC
!-----------------------!
!      mean values      !
!-----------------------!
    VISt_mean(c) = ( VISt_mean(c) * (1.*n) + VISt(c)   ) / ( 1. *(n+1) ) !nu_t
    VISt_sgs_mean(c) = ( VISt_sgs_mean(c) * (1.*n) + VISt_sgs(c)   ) / ( 1. *(n+1) ) !nu_t
    VISt_eff_mean(c) = ( VISt_eff_mean(c) * (1.*n) + VISt_eff(c)   ) / ( 1. *(n+1) ) !nu_t

    U   % mean(c) = ( U  % mean(c) * (1.*n) + U % n(c) ) / ( 1. *(n+1) )
    V   % mean(c) = ( V  % mean(c) * (1.*n) + V % n(c) ) / ( 1. *(n+1) )
    W   % mean(c) = ( W  % mean(c) * (1.*n) + W % n(c) ) / ( 1. *(n+1) )
    P   % mean(c) = ( P  % mean(c) * (1.*n) + P  % n(c) ) / ( 1. *(n+1) )
    Ls  % mean(c) = ( Ls % mean(c) * (1.*n) + ( ( Eps % n(c) / ( Shear(c) / ( 2**0.5 ) ) ** 3 )**0.5 ) ) / ( 1. *(n+1) )
    Kin % mean(c) = ( Kin% mean(c) * (1.*n) + Kin % n(c) ) / ( 1. *(n+1) )
    Eps % mean(c) = ( Eps% mean(c) * (1.*n) + Eps % n(c) ) / ( 1. *(n+1) )

!----------------------------------------------------------------------------------------------------------------------!
!!    mean 2nd order сovariance components  ( arithmetical mean )
!----------------------------------------------------------------------------------------------------------------------!

    !normal
    UU % mean(c) = ( UU % mean(c) * (1.*n) + U % n(c) * U % n(c) ) / ( 1. *(n+1) )!<UU>
    VV % mean(c) = ( VV % mean(c) * (1.*n) + V % n(c) * V % n(c) ) / ( 1. *(n+1) )!<VV>
    WW % mean(c) = ( WW % mean(c) * (1.*n) + W % n(c) * W % n(c) ) / ( 1. *(n+1) )!<WW>
    !shear
    UV % mean(c) = ( UV % mean(c) * (1.*n) + U % n(c) * V % n(c) ) / ( 1. *(n+1) )!<UV>
    UW % mean(c) = ( UW % mean(c) * (1.*n) + U % n(c) * W % n(c) ) / ( 1. *(n+1) )!<UW>
    VW % mean(c) = ( VW % mean(c) * (1.*n) + V % n(c) * W % n(c) ) / ( 1. *(n+1) )!<VW>
  end do
!----------------------------------------------------------------------------------------------------------------------!
!!    Tracking local minimum of Fpy to determine the phase of Uy
!----------------------------------------------------------------------------------------------------------------------!
  if (SIMULA == ZETA) then
   !call GraPhi(U % n(c),1,grad_U_x, .TRUE.) ! grad(U)_x
   !call GraPhi(U % n(c),2,grad_U_y, .TRUE.) ! grad(U)_y
   !call GraPhi(U % n(c),3,grad_U_z, .TRUE.) ! grad(U)_z
   !call GraPhi(V % n(c),1,grad_V_x, .TRUE.) ! grad(V)_x
   !call GraPhi(V % n(c),2,grad_V_y, .TRUE.) ! grad(V)_y
   !call GraPhi(V % n(c),3,grad_V_z, .TRUE.) ! grad(V)_z
   !call GraPhi(W % n(c),1,grad_W_x, .TRUE.) ! grad(W)_x
   !call GraPhi(W % n(c),2,grad_W_y, .TRUE.) ! grad(W)_y
   !call GraPhi(W % n(c),3,grad_W_z, .TRUE.) ! grad(W)_z


    if (this<2) write(*,*) "Fpy=",Fpy
    if (this<2) write(*,*) "MinFound=",MinFoundOnThisCycle
    if (this<2) write(*,*) "Fpy_prev=",Fpy_prev
    if ( (mod(n1,2)==0) .and. (Fpy < 0) ) then ! every second one. working somewhere below zero
      !searching minimum 
      if (     (Fpy < Fpy_prev) .and. (MinFoundOnThisCycle==0) ) then !descending on curve, no minimum yet
        Fpy_min  = Fpy
      elseif ( (Fpy > Fpy_prev) .and. (MinFoundOnThisCycle==0) ) then !global minimum CANDIDATE
        MinFoundOnThisCycle=-1 !candidate state
        CounterFromGlobalMinimumCandidate=0
        do c=-NbC,NC
          UcohMINCanditate(c) = U % n(c)
          VcohMINCanditate(c) = V % n(c)
          WcohMINCanditate(c) = W % n(c)
          
          VIStcohMINCanditate(c) =  VISt(c)
          VortcohMINCanditate(c) =  Vx(c) - Uy(c)


          !UUcohMINCanditate(c) = U % n(c) * U % n(c) !- U % mean(c)
          !VVcohMINCanditate(c) = V % n(c) * V % n(c)
          !WWcohMINCanditate(c) = W % n(c) * W % n(c)
          !UVcohMINCanditate(c) = U % n(c) * V % n(c)
          !UWcohMINCanditate(c) = U % n(c) * W % n(c)
          !VWcohMINCanditate(c) = V % n(c) * W % n(c)
          !2SijVISt - turbulent viscosity
          !VISTxxcohMINCanditate(c) = 2*VISt(c)*( Ux(c) )
          !VISTyycohMINCanditate(c) = 2*VISt(c)*( Vy(c) )
          !VISTzzcohMINCanditate(c) = 2*VISt(c)*( Wz(c) )
          !VISTxycohMINCanditate(c) =   VISt(c)*( Uy(c) + Vx(c) )
          !VISTxzcohMINCanditate(c) =   VISt(c)*( Uz(c) + Wx(c) )
          !VISTyzcohMINCanditate(c) =   VISt(c)*( Vz(c) + Wy(c) )
          
        end do
      elseif ( (Fpy > Fpy_prev) .and. (MinFoundOnThisCycle==-1) .and. (CounterFromGlobalMinimumCandidate<=6) ) then !counting step from 
        CounterFromGlobalMinimumCandidate=CounterFromGlobalMinimumCandidate+1
      if (this<2) write(*,*) "Points from possible minimum=",CounterFromGlobalMinimumCandidate
      elseif ( (Fpy < Fpy_prev) .and. (MinFoundOnThisCycle==-1) .and. (CounterFromGlobalMinimumCandidate<=6) ) then !local minimum
        MinFoundOnThisCycle=0
      elseif ( (Fpy > Fpy_prev) .and. (MinFoundOnThisCycle==-1) .and. (CounterFromGlobalMinimumCandidate>6) ) then !global minimum CANDIDATE
        MinFoundOnThisCycle=1
        do c=-NbC,NC
          UcohMIN(c) = ( UcohMIN(c)*(1.*kMIN) + UcohMINCanditate(c) ) / (1.*(kMIN+1))
          VcohMIN(c) = ( VcohMIN(c)*(1.*kMIN) + VcohMINCanditate(c) ) / (1.*(kMIN+1))
          WcohMIN(c) = ( WcohMIN(c)*(1.*kMIN) + WcohMINCanditate(c) ) / (1.*(kMIN+1))
          VIStcohMIN(c) = ( VIStcohMIN(c)*(1.*kMIN) + VIStcohMINCanditate(c) ) / (1.*(kMIN+1))
          VortcohMIN(c) = ( VortcohMIN(c)*(1.*kMIN) + VortcohMINCanditate(c) ) / (1.*(kMIN+1))

  !       if (this<2) write(*,*) "Cohmin submitted"
          !!so, kMINRead needs to be read in processor.f90 as = the row counf of cohMIN.dat file
          !UUcohMIN(c) = ( UUcohMIN(c)*(1.*kMIN) + UUcohMINCanditate(c) ) / (1.*(kMIN+1))
          !VVcohMIN(c) = ( VVcohMIN(c)*(1.*kMIN) + VVcohMINCanditate(c) ) / (1.*(kMIN+1))
          !WWcohMIN(c) = ( WWcohMIN(c)*(1.*kMIN) + WWcohMINCanditate(c) ) / (1.*(kMIN+1))
          !UVcohMIN(c) = ( UVcohMIN(c)*(1.*kMIN) + UVcohMINCanditate(c) ) / (1.*(kMIN+1))
          !UWcohMIN(c) = ( UWcohMIN(c)*(1.*kMIN) + UWcohMINCanditate(c) ) / (1.*(kMIN+1))
          !VWcohMIN(c) = ( VWcohMIN(c)*(1.*kMIN) + VWcohMINCanditate(c) ) / (1.*(kMIN+1))
          !! same scheme is here ^
          !VISTxxcohMIN(c) = ( VISTxxcohMIN(c)*(1.*kMIN) + VISTxxcohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTyycohMIN(c) = ( VISTyycohMIN(c)*(1.*kMIN) + VISTyycohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTzzcohMIN(c) = ( VISTzzcohMIN(c)*(1.*kMIN) + VISTzzcohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTxycohMIN(c) = ( VISTxycohMIN(c)*(1.*kMIN) + VISTxycohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTxzcohMIN(c) = ( VISTxzcohMIN(c)*(1.*kMIN) + VISTxzcohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTyzcohMIN(c) = ( VISTyzcohMIN(c)*(1.*kMIN) + VISTyzcohMINCanditate(c) ) / (1.*(kMIN+1))

        end do
        kMIN = kMIN+1 !start value - zero

        if(this<2) then
          write(*,*) 'Minimum detected at n=',n !formally, n-1, Fpy_prev
          open(766,FILE='cohMIN.dat',POSITION='APPEND')
          write(766,'(I9,E15.7,I9)') n1, Fpy, kMIN
          close(766)
        endif    

      end if !minimum recorded


      !now, searching 2/3 of minimum, ascending on curve
      if (    (Fpy > Fpy_prev) .and. (MinFoundOnThisCycle==1) .and. &
           (Fpy<(0.6666666-0.03)*Fpy_min) .and. (Fpy>(0.6666666+0.03)*Fpy_min) ) then !2/3 phase +- 3%
        MinFoundOnThisCycle=-23 !candidate state
        do c=-NbC,NC
          UcohMINCanditate(c) = U % n(c)
          VcohMINCanditate(c) = V % n(c)
          WcohMINCanditate(c) = W % n(c)
          VIStcohMINCanditate(c) =  VISt(c)
          VortcohMINCanditate(c) =  Vx(c) - Uy(c)


          !UUcohMINCanditate(c) = U % n(c) * U % n(c)
          !VVcohMINCanditate(c) = V % n(c) * V % n(c)
          !WWcohMINCanditate(c) = W % n(c) * W % n(c)          
          !UVcohMINCanditate(c) = U % n(c) * V % n(c)          
          !UWcohMINCanditate(c) = U % n(c) * W % n(c)          
          !VWcohMINCanditate(c) = V % n(c) * W % n(c) 
          !!2SijVISt - turbulent viscosity
          !VISTxxcohMINCanditate(c) = 2*VISt(c)*( Ux(c) )
          !VISTyycohMINCanditate(c) = 2*VISt(c)*( Vy(c) )
          !VISTzzcohMINCanditate(c) = 2*VISt(c)*( Wz(c) )
          !VISTxycohMINCanditate(c) =   VISt(c)*( Uy(c) + Vx(c) )
          !VISTxzcohMINCanditate(c) =   VISt(c)*( Uz(c) + Wx(c) )
          !VISTyzcohMINCanditate(c) =   VISt(c)*( Vz(c) + Wy(c) )
        end do
      elseif ( (Fpy > Fpy_prev) .and. (MinFoundOnThisCycle==-23) .and. &
           (Fpy<(0.6666666-0.03)*Fpy_min) .and. (Fpy>(0.6666666+0.03)*Fpy_min) ) then
        if ( abs(Fpy_prev-2*Fpy_min/3)>abs(Fpy-2*Fpy_min/3) ) then !new Fpy is closer to 2/3 min
          do c=-NbC,NC !better candidate
            UcohMINCanditate(c) = U % n(c)
            VcohMINCanditate(c) = V % n(c)
            WcohMINCanditate(c) = W % n(c)
            VIStcohMINCanditate(c) =  VISt(c)
            VortcohMINCanditate(c) =  Vx(c) - Uy(c)

            !UUcohMINCanditate(c) = U % n(c) * U % n(c)
            !VVcohMINCanditate(c) = V % n(c) * V % n(c)
            !WWcohMINCanditate(c) = W % n(c) * W % n(c)          
            !UVcohMINCanditate(c) = U % n(c) * V % n(c)          
            !UWcohMINCanditate(c) = U % n(c) * W % n(c)          
            !VWcohMINCanditate(c) = V % n(c) * W % n(c)
            !!2SijVISt - turbulent viscosity
            !VISTxxcohMINCanditate(c) = 2*VISt(c)*( Ux(c) )
            !VISTyycohMINCanditate(c) = 2*VISt(c)*( Vy(c) )
            !VISTzzcohMINCanditate(c) = 2*VISt(c)*( Wz(c) )
            !VISTxycohMINCanditate(c) =   VISt(c)*( Uy(c) + Vx(c) )
            !VISTxzcohMINCanditate(c) =   VISt(c)*( Uz(c) + Wx(c) )
            !VISTyzcohMINCanditate(c) =   VISt(c)*( Vz(c) + Wy(c) )
          end do
        end if
      end if

      !now, claiming 2/3 of minimum, ascending on curve
      if (    (Fpy > Fpy_prev) .and. (MinFoundOnThisCycle==-23) .and. &
           (Fpy_prev<(0.6666666-0.03)*Fpy_min) .and. (Fpy_prev>(0.6666666+0.03)*Fpy_min) &
           .and. (Fpy>(0.6666666-0.03)*Fpy_min)) then ! перевал

        MinFoundOnThisCycle=1
        do c=-NbC,NC
          Ucoh23MIN(c) = ( Ucoh23MIN(c)*(1.*k23MIN) + UcohMINCanditate(c) ) / (1.*(k23MIN+1))
          Vcoh23MIN(c) = ( Vcoh23MIN(c)*(1.*k23MIN) + VcohMINCanditate(c) ) / (1.*(k23MIN+1))
          Wcoh23MIN(c) = ( Wcoh23MIN(c)*(1.*k23MIN) + WcohMINCanditate(c) ) / (1.*(k23MIN+1))
          !so, kMINRead needs to be read in processor.f90 as = the row counf of cohMIN.dat file
          VIStcoh23MIN(c) = ( VIStcoh23MIN(c)*(1.*k23MIN) + VIStcohMINCanditate(c) ) / (1.*(k23MIN+1))
          Vortcoh23MIN(c) = ( Vortcoh23MIN(c)*(1.*k23MIN) + VortcohMINCanditate(c) ) / (1.*(k23MIN+1))
          !UUcoh23MIN(c) = ( UUcoh23MIN(c)*(1.*k23MIN) + UUcohMINCanditate(c) ) / (1.*(k23MIN+1))
          !VVcoh23MIN(c) = ( VVcoh23MIN(c)*(1.*k23MIN) + VVcohMINCanditate(c) ) / (1.*(k23MIN+1))
          !WWcoh23MIN(c) = ( WWcoh23MIN(c)*(1.*k23MIN) + WWcohMINCanditate(c) ) / (1.*(k23MIN+1))
          !UVcoh23MIN(c) = ( UVcoh23MIN(c)*(1.*k23MIN) + UVcohMINCanditate(c) ) / (1.*(k23MIN+1))
          !UWcoh23MIN(c) = ( UWcoh23MIN(c)*(1.*k23MIN) + UWcohMINCanditate(c) ) / (1.*(k23MIN+1))
          !VWcoh23MIN(c) = ( VWcoh23MIN(c)*(1.*k23MIN) + VWcohMINCanditate(c) ) / (1.*(k23MIN+1))

          !VISTxxcoh23MIN(c) = ( VISTxxcoh23MIN(c)*(1.*kMIN) + VISTxxcohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTyycoh23MIN(c) = ( VISTyycoh23MIN(c)*(1.*kMIN) + VISTyycohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTzzcoh23MIN(c) = ( VISTzzcoh23MIN(c)*(1.*kMIN) + VISTzzcohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTxycoh23MIN(c) = ( VISTxycoh23MIN(c)*(1.*kMIN) + VISTxycohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTxzcoh23MIN(c) = ( VISTxzcoh23MIN(c)*(1.*kMIN) + VISTxzcohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTyzcoh23MIN(c) = ( VISTyzcoh23MIN(c)*(1.*kMIN) + VISTyzcohMINCanditate(c) ) / (1.*(kMIN+1))

        end do

        k23MIN = k23MIN+1
        if(this<2) then
          write(*,*) '2/3 of minimum detected at n=',n
          open(767,FILE='coh23MIN.dat',POSITION='APPEND')
          write(767,'(I9,E15.7,I9)') n1, Fpy, k23MIN
          close(767)
        endif      
      end if   

      !now, searching 1/3 of minimum, ascending on curve
      if (    (Fpy > Fpy_prev) .and. (MinFoundOnThisCycle==1) .and. &
           (Fpy<(0.3333333-0.03)*Fpy_min) .and. (Fpy>(0.3333333+0.03)*Fpy_min) ) then !2/3 phase +- 3%
        MinFoundOnThisCycle=-13 !candidate state
        do c=-NbC,NC
          UcohMINCanditate(c) = U % n(c)
          VcohMINCanditate(c) = V % n(c)
          WcohMINCanditate(c) = W % n(c)
          VIStcohMINCanditate(c) =  VISt(c)
          VortcohMINCanditate(c) =  Vx(c) - Uy(c)


          !UUcohMINCanditate(c) = U % n(c) * U % n(c)
          !VVcohMINCanditate(c) = V % n(c) * V % n(c)
          !WWcohMINCanditate(c) = W % n(c) * W % n(c)
          !UVcohMINCanditate(c) = U % n(c) * V % n(c)
          !UWcohMINCanditate(c) = U % n(c) * W % n(c)
          !VWcohMINCanditate(c) = V % n(c) * W % n(c)
          !!2SijVISt - turbulent viscosity
          !VISTxxcohMINCanditate(c) = 2*VISt(c)*( Ux(c) )
          !VISTyycohMINCanditate(c) = 2*VISt(c)*( Vy(c) )
          !VISTzzcohMINCanditate(c) = 2*VISt(c)*( Wz(c) )
          !VISTxycohMINCanditate(c) =   VISt(c)*( Uy(c) + Vx(c) )
          !VISTxzcohMINCanditate(c) =   VISt(c)*( Uz(c) + Wx(c) )
          !VISTyzcohMINCanditate(c) =   VISt(c)*( Vz(c) + Wy(c) )
        end do
      elseif ( (Fpy > Fpy_prev) .and. (MinFoundOnThisCycle==-13) .and. &
           (Fpy<(0.3333333-0.03)*Fpy_min) .and. (Fpy>(0.3333333+0.03)*Fpy_min) ) then
        if ( abs(Fpy_prev-1*Fpy_min/3)>abs(Fpy-1*Fpy_min/3) ) then !new Fpy is closer to 2/3 min
          do c=-NbC,NC !better candidate
            UcohMINCanditate(c) = U % n(c)
            VcohMINCanditate(c) = V % n(c)
            WcohMINCanditate(c) = W % n(c)
            VIStcohMINCanditate(c) =  VISt(c)
            VortcohMINCanditate(c) =  Vx(c) - Uy(c)



            !UUcohMINCanditate(c) = U % n(c) * U % n(c)
            !VVcohMINCanditate(c) = V % n(c) * V % n(c)
            !WWcohMINCanditate(c) = W % n(c) * W % n(c)
            !UVcohMINCanditate(c) = U % n(c) * V % n(c)
            !UWcohMINCanditate(c) = U % n(c) * W % n(c)
            !VWcohMINCanditate(c) = V % n(c) * W % n(c) 
            !!2SijVISt - turbulent viscosity
            !VISTxxcohMINCanditate(c) = 2*VISt(c)*( Ux(c) )
            !VISTyycohMINCanditate(c) = 2*VISt(c)*( Vy(c) )
            !VISTzzcohMINCanditate(c) = 2*VISt(c)*( Wz(c) )
            !VISTxycohMINCanditate(c) =   VISt(c)*( Uy(c) + Vx(c) )
            !VISTxzcohMINCanditate(c) =   VISt(c)*( Uz(c) + Wx(c) )
            !VISTyzcohMINCanditate(c) =   VISt(c)*( Vz(c) + Wy(c) )
          end do
        end if
      end if

      !now, claiming 1/3 of minimum, ascending on curve
      if (    (Fpy > Fpy_prev) .and. (MinFoundOnThisCycle==-13) .and. &
           (Fpy_prev<(0.3333333-0.03)*Fpy_min) .and. (Fpy_prev>(0.3333333+0.03)*Fpy_min) &
           .and. (Fpy>(0.3333333-0.03)*Fpy_min)) then ! перевал
      
        MinFoundOnThisCycle=1
        do c=-NbC,NC
          Ucoh13MIN(c) = ( Ucoh13MIN(c)*(1.*k13MIN) + UcohMINCanditate(c) ) / (1.*(k13MIN+1))
          Vcoh13MIN(c) = ( Vcoh13MIN(c)*(1.*k13MIN) + VcohMINCanditate(c) ) / (1.*(k13MIN+1))
          Wcoh13MIN(c) = ( Wcoh13MIN(c)*(1.*k13MIN) + WcohMINCanditate(c) ) / (1.*(k13MIN+1))
          !!so, kMINRead needs to be read in processor.f90 as = the row counf of cohMIN.dat file
          VIStcoh13MIN(c) = ( VIStcoh13MIN(c)*(1.*k13MIN) + VIStcohMINCanditate(c) ) / (1.*(k13MIN+1))
          Vortcoh13MIN(c) = ( Vortcoh13MIN(c)*(1.*k13MIN) + VortcohMINCanditate(c) ) / (1.*(k13MIN+1))
          !UUcoh13MIN(c) = ( UUcoh13MIN(c)*(1.*k13MIN) + UUcohMINCanditate(c) ) / (1.*(k13MIN+1))
          !VVcoh13MIN(c) = ( VVcoh13MIN(c)*(1.*k13MIN) + VVcohMINCanditate(c) ) / (1.*(k13MIN+1))
          !WWcoh13MIN(c) = ( WWcoh13MIN(c)*(1.*k13MIN) + WWcohMINCanditate(c) ) / (1.*(k13MIN+1))
          !UVcoh13MIN(c) = ( UVcoh13MIN(c)*(1.*k13MIN) + UVcohMINCanditate(c) ) / (1.*(k13MIN+1))
          !UWcoh13MIN(c) = ( UWcoh13MIN(c)*(1.*k13MIN) + UWcohMINCanditate(c) ) / (1.*(k13MIN+1))
          !VWcoh13MIN(c) = ( VWcoh13MIN(c)*(1.*k13MIN) + VWcohMINCanditate(c) ) / (1.*(k13MIN+1))

          !VISTxxcoh13MIN(c) = ( VISTxxcoh13MIN(c)*(1.*kMIN) + VISTxxcohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTyycoh13MIN(c) = ( VISTyycoh13MIN(c)*(1.*kMIN) + VISTyycohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTzzcoh13MIN(c) = ( VISTzzcoh13MIN(c)*(1.*kMIN) + VISTzzcohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTxycoh13MIN(c) = ( VISTxycoh13MIN(c)*(1.*kMIN) + VISTxycohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTxzcoh13MIN(c) = ( VISTxzcoh13MIN(c)*(1.*kMIN) + VISTxzcohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTyzcoh13MIN(c) = ( VISTyzcoh13MIN(c)*(1.*kMIN) + VISTyzcohMINCanditate(c) ) / (1.*(kMIN+1))
        end do

        k13MIN = k13MIN+1
        if(this<2) then
          write(*,*) '1/3 of minimum detected at n=',n
          open(768,FILE='coh13MIN.dat',POSITION='APPEND')
          write(768,'(I9,E15.7,I9)') n1, Fpy, k13MIN
          close(768)
        endif 
      end if

      Fpy_prev = Fpy

    elseif ( (mod(n1,2)==0) .and. (Fpy > 0.0) .and. (Fpy_prev < 0.0) .and. (MinFoundOnThisCycle==1) )  then ! zero point

      !now, searching 0, ascending on curve
        MinFoundOnThisCycle=0
        do c=-NbC,NC
          Ucoh0(c) = ( Ucoh0(c)*(1.*k0) + U % n(c) ) / (1.*(k0+1))
          Vcoh0(c) = ( Vcoh0(c)*(1.*k0) + V % n(c) ) / (1.*(k0+1))
          Wcoh0(c) = ( Wcoh0(c)*(1.*k0) + W % n(c) ) / (1.*(k0+1))
          !so, kMINRead needs to be read in processor.f90 as = the row counf of cohMIN.dat file
          VIStcoh0(c) = ( VIStcoh0(c)*(1.*k0) + VISt(c) ) / (1.*(k0+1))
          Vortcoh0(c) = ( Vortcoh0(c)*(1.*k0) + Vx(c) - Uy(c) ) / (1.*(k0+1))
          !UUcoh0MIN(c) = ( UUcoh0MIN(c)*(1.*k0MIN) + UUcohMINCanditate(c) ) / (1.*(k0MIN+1))
          !VVcoh0MIN(c) = ( VVcoh0MIN(c)*(1.*k0MIN) + VVcohMINCanditate(c) ) / (1.*(k0MIN+1))
          !WWcoh0MIN(c) = ( WWcoh0MIN(c)*(1.*k0MIN) + WWcohMINCanditate(c) ) / (1.*(k0MIN+1))
          !UVcoh0MIN(c) = ( UVcoh0MIN(c)*(1.*k0MIN) + UVcohMINCanditate(c) ) / (1.*(k0MIN+1))
          !UWcoh0MIN(c) = ( UWcoh0MIN(c)*(1.*k0MIN) + UWcohMINCanditate(c) ) / (1.*(k0MIN+1))
          !VWcoh0MIN(c) = ( VWcoh0MIN(c)*(1.*k0MIN) + VWcohMINCanditate(c) ) / (1.*(k0MIN+1))

          !VISTxxcoh0MIN(c) = ( VISTxxcoh0MIN(c)*(1.*kMIN) + VISTxxcohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTyycoh0MIN(c) = ( VISTyycoh0MIN(c)*(1.*kMIN) + VISTyycohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTzzcoh0MIN(c) = ( VISTzzcoh0MIN(c)*(1.*kMIN) + VISTzzcohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTxycoh0MIN(c) = ( VISTxycoh0MIN(c)*(1.*kMIN) + VISTxycohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTxzcoh0MIN(c) = ( VISTxzcoh0MIN(c)*(1.*kMIN) + VISTxzcohMINCanditate(c) ) / (1.*(kMIN+1))
          !VISTyzcoh0MIN(c) = ( VISTyzcoh0MIN(c)*(1.*kMIN) + VISTyzcohMINCanditate(c) ) / (1.*(kMIN+1))
        end do
        k0 = k0+1

        if(this<2) then
          write(*,*) 'zero detected at n=',n
          open(769,FILE='coh0.dat',POSITION='APPEND')
          write(769,'(I9,E15.7,I9)') n1, Fpy, k0
          close(769)
        endif         

    end if  !mod 2
  end if !tracking forces
  end if ! n > -1

call Wait
RETURN 

END SUBROUTINE CalcMn
