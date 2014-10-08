!======================================================================!
  SUBROUTINE ReaCom(restar) 
!----------------------------------------------------------------------!
!   Reads the T-Rex.cmn file, with commands to T-Rex.                  ! 
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  LOGICAL   :: restar
!------------------------------[Calling]-------------------------------!
  REAL      :: Dist
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: i, j, l, m
  REAL      :: Mres(MAXP), MresT, dummy
  REAL      :: xm(MAXP), ym(MAXP), zm(MAXP)
  CHARACTER :: answer*80, nammon*80
!--------------------------------[CVS]---------------------------------!
!  $Id: ReaCom.f90,v 1.37 2009/06/30 12:11:49 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/ReaCom.f90,v $  
!======================================================================!

  call Wait   

!----- Angular velocity vector
  if(ROT == YES) then
    if(this  < 2)  write(*,*) '# Angular velocity vector: '
    call ReadC(7,inp,tn,ts,te)
    read(inp,*)  omegaX, omegaY, omegaZ
    if(this  < 2)  write(*,*) '  ', omegaX, omegaY, omegaZ
  end if
!----- The number of time steps
  if(this  < 2)  write(*,*) '# Enter the number of time steps:'
    !write(*,*) '# (type 0 if you just want to analyse results)'
  call ReadC(7,inp,tn,ts,te)
  read(inp,*)  Ndt
  if(this  < 2)  write(*,*) '  ', Ndt

!----- Starting time step for statistics 
  if(this  < 2)  write(*,*) '# Starting time step for statistics:'
  call ReadC(7,inp,tn,ts,te)
  read(inp,*)  Nstat
  if(this  < 2)  write(*,*) '  ', Nstat
  if(BUDG == YES) then
    read(inp(ts(2):te(2)),*) Nbudg
    if(this  < 2)  write(*,*) '  ', Nbudg
  end if



  if(this  < 2) write(*,*) '# Number of monitoring points:'
  call ReadC(7,inp,tn,ts,te)
  read(inp,*) Nmon 
  if(this  < 2)  write(*,*) '  ', Nmon
  if(this  < 2)  write(*,*) '# Enter the coordinates of monitoring point(s)'
  do i=1,Nmon
    call ReadC(7,inp,tn,ts,te)
    read(inp,*)  xm(i), ym(i), zm(i)
    if(this  < 2)  write(*,*) '  ',  xm(i), ym(i), zm(i)
  end do

!----- Find the monitoring cells
  nammon=name 
  nammon(len_trim(name)+1:len_trim(name)+10)="-monit.000"
  l=len_trim(nammon) 
  do j=1,Nmon
    Mres(j)=HUGE
    do i=1,NC
      if(dist(xm(j),ym(j),zm(j),  &
	      xc(i),yc(i),zc(i))  < Mres(j)) then
	Cm(j)=i
	Mres(j)=dist(xm(j),ym(j),zm(j),xc(i),yc(i),zc(i))
      end if
    end do
    MresT=Mres(j)
    call GloMin(MresT)
    if(MresT /= Mres(j)) then ! there is a cell which is nearer
      Cm(j) = 0               ! so erase this monitoring point     
    end if 
  end do

  do j=1,Nmon
    if(Cm(j)  > 0) then
      if(j  <  10) then
	write(nammon(l  :l),'(I1)') j
      else if(j  < 100) then
	write(nammon(l-1:l),'(I2)') j
      else
	write(nammon(l-2:l),'(I3)') j
      end if
      if(Ndtt == 0) then
        open(10+j,FILE=nammon)
      else
        open(10+j,FILE=nammon,POSITION='APPEND')
      endif

      write(10+j,'(A24,3F16.6)') '# Monitoring point:',xc(Cm(j)),yc(Cm(j)),zc(Cm(j))
    end if
  end do

!----- Plane for calcution of overall mass fluxes
  do m=1,Nmat
    if(this  < 2)  write(*,*) '# Enter the coordinates of monitoring plane:'
    call ReadC(7,inp,tn,ts,te)
    read(inp,*) xp(m), yp(m), zp(m)
    if(this  < 2)  write(*,*) '  ', xp(m), yp(m), zp(m)
  end do

!----- Kind of simulation
  if(this  < 2)  write(*,*) '# Type of simulation: '
  call ReadC(7,inp,tn,ts,te)
  read(inp(ts(1):te(1)),'(A)')  answer
  call ToUppr(answer)
  if(answer == 'DNS') then
    SIMULA = DNS
    if(this  < 2) write(*,*) '  DNS -> Direct Numerical Simulation'
  else if(answer == 'LES') then
    SIMULA = LES
    if(this  < 2) write(*,*) '  LES -> Large Eddy Simulation'
  else if(answer == 'K_EPS') then
    SIMULA = K_EPS
    if(this  < 2) write(*,*) '  K_EPS -> High Reynolds k-eps model.' 
  else if(answer == 'K_EPS_VV') then
    SIMULA = K_EPS_VV
    if(this  < 2) write(*,*) '  K_EPS_VV -> Durbin`s model.' 
  else if(answer == 'SPA_ALL') then
    SIMULA = SPA_ALL
    if(this  < 2) write(*,*) '  SPA_ALL -> Spalart-Allmaras model.' 
  else if(answer == 'DES_SPA') then
    SIMULA = DES_SPA
    if(this  < 2) write(*,*) '  DES_SPA -> ???' 
  else if(answer == 'ZETA') then
    SIMULA = ZETA
    if(this  < 2) write(*,*) '  ZETA -> k-eps-zeta-f model.' 
  else if(answer == 'HYB_PITM') then
    SIMULA = HYB_PITM
    if(this  < 2) write(*,*) '  HYB_PITM -> ???' 
  else if(answer == 'HYB_ZETA') then
    SIMULA = HYB_ZETA
    if(this  < 2) write(*,*) '  HYB_ZETA -> Hybrid ZETA and LES' 
  else if(answer == 'EBM') then
    SIMULA = EBM
    if(this  < 2) write(*,*) '  EBM -> Second Momentum Closure' 
  else if(answer == 'HJ') then
    SIMULA = HJ
    if(this  < 2) write(*,*) '  HJ -> ???' 
  else
    if(this  < 2) write(*,*) '  Error in input ! Exiting'
    stop
  endif



  if(SIMULA==EBM.or.SIMULA==HJ) then
    read(inp(ts(2):te(2)),'(A8)') answer
    call ToUppr(answer)
    if(answer == 'HYB') then
      MODE = HYB 
    end if
  end if 
  if(SIMULA==K_EPS) then
    read(inp(ts(2):te(2)),'(A8)') answer
    call ToUppr(answer)
    if(answer == 'HRE') then
      MODE = HRE
    else if(answer == 'LRE') then
      MODE = LRE 
    else
      if(this  < 2) write(*,*) '  Error in input ! Exiting'
      stop
    end if
  end if

  if(SIMULA==LES.or.SIMULA==HYB_ZETA) then
    read(inp(ts(2):te(2)),'(A8)') answer
    call ToUppr(answer)
    if(answer == 'SMAG') then
      MODE = SMAG 
    else if(answer == 'DYN') then
      MODE = DYN 
    else if(answer == 'WALE') then
      MODE = WALE
    else if(answer == 'MIX') then
      MODE = MIX
    else
      if(this  < 2) write(*,*) '  Error in input ! Exiting'
      stop
    end if
  end if

  if(SIMULA  ==  LES.and.MODE == SMAG) then
    if(this  < 2)  write(*,*) '# C Smagorinsky = '
    read(inp(ts(3):te(3)),*) Cs0
    if(this  < 2)  write(*,*) '  ', Cs0 
  endif 

    
  do m=1,Nmat
    if(SIMULA  ==  LES .or. SIMULA == DNS .or. SIMULA == DES_SPA &
      .or.SIMULA  ==  HYB_PITM .or. SIMULA == HYB_ZETA) then
      if(this  < 2) write(*,*) '# Do you want to shake the velocity field ?'
      call ReadC(7,inp,tn,ts,te)
      read(inp(ts(1):te(1)),'(A)')  answer
      call ToUppr(answer)
      if(answer == 'YES') then
        if(this  < 2)  write(*,*) '  YES -> shake'
        SHAKE(m) = YES
      else if(answer == 'NO') then
        SHAKE(m) = NO
        if(this  < 2) write(*,*) '  NO -> don''t shake'
      else
        if(this  < 2) write(*,*) '  Error in input ! Exiting'
        stop
      endif
      if(SHAKE(m) == YES) then
        if(this < 2) write(*,*) '# For how many time steps you want to shake ?'
        call ReadC(7,inp,tn,ts,te)
        read(inp,*) SHAKE_PER(m) 
        if(this  < 2)  write(*,*) '  ', SHAKE_PER(m)
        if(this < 2) write(*,*) '# Interval for shaking:'
        call ReadC(7,inp,tn,ts,te)
        read(inp,*) SHAKE_INT(m)
        if(this  < 2)  write(*,*) '  ', SHAKE_INT(m)
      end if
    endif 
  end do

  if(.not. restar) call UnkAloc

  
  if(SIMULA  ==  K_EPS.or.SIMULA  ==  HYB_PITM) then
    Ce1 = 1.44
    Ce2 = 1.92
    Cmu = 0.09
    Cmu25 = sqrt(sqrt(Cmu))
    Cmu75 = Cmu25**3 
    kappa = 0.41
    Elog  = 8.342
    Kin % Sigma = 1.0
    Eps % Sigma = 1.3
    if(MODE == LRe) then
      Ce1 = 1.55
      Ce2 = 2.0
    end if
  endif
 
  if(SIMULA  ==  K_EPS_VV) then
    Ce1   = 1.4
    Ce2   = 1.9
    Cmu   = 0.09
    CmuD  = 0.22
    Cmu25 = sqrt(sqrt(Cmu))
    Cmu75 = Cmu25**3
    kappa = 0.41
    Elog  = 8.342
    Cl    = 0.22
    Ct    = 6.0
    Cni   = 85.0
    alpha = 0.045
    cf1   = 0.697
    cf2   = 4.4  
    cf3   = 1.6 
    Cf_1  = 1.4
    Cf_2  = 0.3
    Lim   = 11.0
    Kin % Sigma = 1.0
    Eps % Sigma = 1.3
    v_2 % Sigma = 1.0
  endif

  if(SIMULA  ==  EBM) then
    Ce1   = 1.44
    Ce2   = 1.83
    Ce3   = 0.55
    CmuD  = 0.21
    Cmu25 = sqrt(sqrt(Cmu))
    Cmu75 = Cmu25**3
    Cl    = 0.161
    Ct    = 6.0
    Cni   = 80.0
    Kin % Sigma = 1.0
    Eps % Sigma = 1.15
    uu % Sigma = 1.0
    vv % Sigma = 1.0
    ww % Sigma = 1.0
    uv % Sigma = 1.0
    uw % Sigma = 1.0
    vw % Sigma = 1.0
    g1    = 3.4
    g1_star = 1.8
    g2    = 4.2
    g3    = 0.8
    g3_star = 1.3
    g4    = 1.25
    g5    = 0.4
  endif

  if(SIMULA  ==  HJ) then
    Ce1   = 1.44
    Ce2   = 1.83
    Ce3   = 0.55
    CmuD  = 0.21
    Cmu25 = sqrt(sqrt(Cmu))
    Cmu75 = Cmu25**3
    Cl    = 0.161
    Ct    = 6.0
    Cni   = 80.0
    Kin % Sigma = 1.0
    Eps % Sigma = 1.15
    uu % Sigma = 1.0
    vv % Sigma = 1.0
    ww % Sigma = 1.0
    uv % Sigma = 1.0
    uw % Sigma = 1.0
    vw % Sigma = 1.0
    g1    = 3.4
    g1_star = 1.8
    g2    = 4.2
    g3    = 0.8
    g3_star = 1.3
    g4    = 1.25
    g5    = 0.4
  endif

  if(SIMULA  ==  ZETA.or.SIMULA  ==  HYB_ZETA) then
    Ce1   = 1.4
    Ce2   = 1.9
    Cmu   = 0.09 
    CmuD  = 0.22
    Cmu25 = sqrt(sqrt(Cmu))
    Cmu75 = Cmu25**3
    kappa = 0.41
    Elog  = 8.342
    Cl    = 0.36
    Ct    = 6.0
    Cni   = 85.0
    alpha = 0.012
    cf1   = 0.697
    cf2   = 4.4
    cf3   = 1.6
    Cf_1  = 1.4
    Cf_2  = 0.3
    Lim   = 11.0
    Kin % Sigma = 1.0
    Eps % Sigma = 1.3
    v_2 % Sigma = 1.2
  endif

  if(SIMULA == SPA_ALL .or. SIMULA==DES_SPA) then
    kappa  = 0.41
    Cb1    = 0.1355
    Cb2    = 0.622
    Cvis1  = 7.1
    Cw2    = 0.3
    Cw3    = 2.0
    VIS % Sigma = 2.0/3.0
    Cw1    = Cb1/kappa**2.0 + (1+Cb2)/VIS % Sigma
    SIGMAv = 2.0/3.0
  end if

!----- Time stepping scheme
  if(this  < 2) write(*,*) '# Algorythm for time-integration: '
  call ReadC(7,inp,tn,ts,te)
  read(inp(ts(1):te(1)),'(A)')  answer
  call ToUppr(answer)
  if(answer == 'SIMPLE') then
    if(this  < 2) write(*,*) '  SIMPLE [Nini] -> S. I. M. P. L. E.'
    ALGOR = SIMPLE
    Nini  = 10
    if(tn==2) read(inp(ts(2):te(2)),*) Nini
  else if(answer == 'FRACTION') then
    if(this  < 2) write(*,*) '  FRACTION -> Fractional step method'
    ALGOR = FRACT
    Nini  = 1
  else
    if(this  < 2) write(*,*) 'Error in input ! Exiting'
    stop
  endif

  if(ALGOR == SIMPLE) then
    if(this < 2) write(*,*) '# Under Relaxation Factor for velocity:'
    call ReadC(7,inp,tn,ts,te)
    read(inp,*)  U % URF
    if(this < 2) write(*,*) '  ', U % URF
    if(this < 2) write(*,*) '# Under Relaxation Factor for pressure:'
    call ReadC(7,inp,tn,ts,te)
    read(inp,*)  P % URF
    if(this < 2) write(*,*) '  ', P % URF
    if(HOT == YES) then
      if(this < 2) write(*,*) '# Under Relaxation Factor for temperature:'
      call ReadC(7,inp,tn,ts,te)
      read(inp,*)  T % URF
      if(this < 2) write(*,*) '  ', T % URF
    end if
    if(SIMULA /= LES .and. SIMULA /= DNS) then
      if(this < 2) write(*,*) '# Under Relaxation Factor for turbulent variables:'
      call ReadC(7,inp,tn,ts,te)
      read(inp,*)  URFT
      if(this < 2) write(*,*) '  ', T % URF
    end if
  endif
  Kin % URF   = URFT
  Eps % URF   = URFT
  v_2 % URF   = URFT
  f22 % URF   = URFT 
  VIS % URF   = URFT
  uu  % URF   = URFT
  vv  % URF   = URFT
  ww  % URF   = URFT
  uv  % URF   = URFT
  uw  % URF   = URFT
  vw  % URF   = URFT


  if(this  < 2) write(*,*) '# Integration of inertial terms: '
  call ReadC(7,inp,tn,ts,te)
  read(inp(ts(1):te(1)),'(A)')  answer
  call ToUppr(answer)
  if(answer == 'LIN') then
    if(this < 2) write(*,*) '  LIN -> Linear'
    INERT = LIN
  else if(answer == 'PAR') then
    if(this < 2) write(*,*) '  PAR -> Parabolic'
    INERT = PAR
  else
    if(this  < 2) write(*,*) '  Error in input ! Exiting1'
    stop
  endif

  if(this  < 2) then
    write(*,*) '# Integration of convective terms: '
  endif 
  call ReadC(7,inp,tn,ts,te)
  read(inp(ts(1):te(1)),'(A)')  answer
  call ToUppr(answer)
  if(answer == 'AB') then
    CONVEC = AB
    if(this < 2) write(*,*) '  AB -> Adams-Bashforth'
  else if(answer == 'CN') then
    CONVEC = CN
    if(this < 2) write(*,*) '  CN -> Crank-Nicholson'
  else if(answer == 'FI') then
    CONVEC = FI
    if(this < 2) write(*,*) '  FI -> Fully Implicit'
  else
    if(this  < 2) write(*,*) '  Error in input ! Exiting2'
    stop
  endif

  if(this  < 2) write(*,*) '# Integration of diffusive terms: '
  call ReadC(7,inp,tn,ts,te)
  read(inp(ts(1):te(1)),'(A)')  answer
  call ToUppr(answer)
  if(answer == 'AB') then
    DIFFUS = AB
    if(this < 2) write(*,*) '  AB -> Adams-Bashforth'
  else if(answer == 'CN') then
    DIFFUS = CN
    if(this < 2) write(*,*) '  CN -> Crank-Nicholson'
  else if(answer == 'FI') then
    DIFFUS = FI
    if(this < 2) write(*,*) '  FI -> Fully Implicit'
  else
    if(this  < 2) write(*,*) '  Error in input ! Exiting3'
    stop
  endif

  if(this  < 2) write(*,*) '# Integration of cross-diffusive terms: '
  call ReadC(7,inp,tn,ts,te)
  read(inp(ts(1):te(1)),'(A)')  answer
  call ToUppr(answer)
  if(answer == 'AB') then
    CROSS = AB
    if(this < 2) write(*,*) '  AB -> Adams-Bashforth'
  else if(answer == 'CN') then
    CROSS = CN
    if(this < 2) write(*,*) '  CN -> Crank-Nicholson'
  else if(answer == 'FI') then
    if(this < 2) write(*,*) '  FI -> Fully Implicit'
    CROSS = FI
  else
    if(this  < 2) write(*,*) '  Error in input ! Exiting'
    stop
  endif

!----- Upwind blending
  do m=1,Nmat
    URFC(m) = 1.0
    if(this  < 2) write(*,*) '# Convetive schemes for momentum equation:'
    !if(this  < 2) write(*,*) '  # Do you want to use upwind blending: '
    call ReadC(7,inp,tn,ts,te)
    read(inp(ts(1):te(1)),'(A)')  answer
    call ToUppr(answer)
    if(answer == 'BLEND_CDS_UDS') then
      BLEND(m) = YES
      if(this  < 2) write(*,*) '  YES -> use blening'
      if(tn==2) read(inp(ts(2):te(2)),*) URFC(m)
    else if(answer == 'NO') then
      BLEND(m) = NO
      if(this  < 2) write(*,*) '  NO -> don''t use blending'
    else if(answer == 'UDS') then
      BLEND(m) = YES
      if(this  < 2) write(*,*) '  UDS -> cross-differencing'
      URFC(m)  = 0.0
    else if(answer == 'CDS') then
      BLEND(m) = CDS 
      if(this  < 2) write(*,*) '  CDS -> central differencing'
    else if(answer == 'LUDS') then
      BLEND(m) = LUDS
      if(this  < 2) write(*,*) '  LUDS -> linear upwind'
    else if(answer == 'QUICK') then
      BLEND(m) = QUICK
      if(this  < 2) write(*,*) '  QUICK -> self descriptive'
    else if(answer == 'MINMOD') then
      BLEND(m) = MINMOD
      if(this  < 2) write(*,*) '  MINMOD -> self descriptive'
    else if(answer == 'SMART') then
      BLEND(m) = SMART 
      if(this  < 2) write(*,*) '  SMART -> self descriptive'
    else if(answer == 'AVL_SMART') then
      BLEND(m) = AVL_SMART 
      if(this  < 2) write(*,*) '  AVL_SMART -> self descriptive'
    else if(answer == 'SUPERBEE') then
      BLEND(m) = SUPERBEE 
      if(this  < 2) write(*,*) '  SUPERBEE -> ???'
    else if(answer == 'GAMMA') then
      BLEND(m) = GAMMA 
      if(this  < 2) write(*,*) '  GAMMA -> ???'
    else
      if(this  < 2) write(*,*) '  Error in input in line:' 
      if(this  < 2) write(*,*) answer
      if(this  < 2) write(*,*) '  Exiting !'
      stop
    endif
  end do

  if(HOT==YES) then
    do m=1,Nmat
      URFC_Tem(m) = 1.0
      if(this  < 2) write(*,*) '# Convetive schemes for energy equation:'
      call ReadC(7,inp,tn,ts,te)
      read(inp(ts(1):te(1)),'(A)')  answer
      call ToUppr(answer)
      if(answer == 'BLEND_TEM_CDS_UDS') then
        BLEND_TEM(m) = YES
        if(this  < 2) write(*,*) '  YES -> use blening'
        if(tn==2) read(inp(ts(2):te(2)),*) URFC_Tem(m)
      else if(answer == 'NO') then
        BLEND_TEM(m) = NO
        if(this  < 2) write(*,*) '  NO -> don''t use blending'
      else if(answer == 'UDS') then
        BLEND_TEM(m) = YES
        if(this  < 2) write(*,*) '  UDS -> cross'
        URFC_Tem(m)  = 0.0
      else if(answer == 'CDS') then
        BLEND_TEM(m) = CDS 
        if(this  < 2) write(*,*) '  CDS -> central differencing'
      else if(answer == 'LUDS') then
        BLEND_TEM(m) = LUDS
        if(this  < 2) write(*,*) '  LUDS -> linear upwind'
      else if(answer == 'QUICK') then
        BLEND_TEM(m) = QUICK
        if(this  < 2) write(*,*) '  QUICK -> self descriptive'
      else if(answer == 'MINMOD') then
        BLEND_TEM(m) = MINMOD
        if(this  < 2) write(*,*) '  MINMOD -> self descriptive'
      else if(answer == 'SMART') then
        BLEND_TEM(m) = SMART 
        if(this  < 2) write(*,*) '  SMART -> self descriptive'
      else if(answer == 'AVL_SMART') then
        BLEND_TEM(m) = AVL_SMART 
        if(this  < 2) write(*,*) '  AVL_SMART -> self descriptive'
      else if(answer == 'SUPERBEE') then
        BLEND_TEM(m) = SUPERBEE 
        if(this  < 2) write(*,*) '  SUPERBEE -> ???'
      else if(answer == 'GAMMA') then
        BLEND_TEM(m) = GAMMA 
        if(this  < 2) write(*,*) '  GAMMA -> ???'
      else
       if(this  < 2) write(*,*) '  Error in input in line:' 
       if(this  < 2) write(*,*) answer
       if(this  < 2) write(*,*) '  Exiting !'
        stop
      endif
    end do
  end if

  if(SIMULA/=LES.and.SIMULA/=DNS) then
    do m=1,Nmat
      URFC_Tur(m) = 1.0
      if(this  < 2) write(*,*) '# Convetive schemes for transport equation:'
      call ReadC(7,inp,tn,ts,te)
      read(inp(ts(1):te(1)),'(A)')  answer
      call ToUppr(answer)
      if(answer == 'BLEND_TUR_CDS_UDS') then
        BLEND_TUR(m) = YES
        if(this  < 2) write(*,*) '  YES    -> use blening'
        if(tn==2) read(inp(ts(2):te(2)),*) URFC_Tur(m)
      else if(answer == 'NO') then
        BLEND_TUR(m) = NO
        if(this  < 2) write(*,*) '  NO -> don''t use blending'
      else if(answer == 'UDS') then
        BLEND_TUR(m) = YES
        URFC_Tur(m)  = 0.0
        if(this  < 2) write(*,*) '  UDS -> cross'
      else if(answer == 'CDS') then
        BLEND_TUR(m) = CDS 
        if(this  < 2) write(*,*) '  CDS -> central differencing'
      else if(answer == 'LUDS') then
        BLEND_TUR(m) = LUDS
        if(this  < 2) write(*,*) '  LUDS -> linear upwind'
      else if(answer == 'QUICK') then
        BLEND_TUR(m) = QUICK
        if(this  < 2) write(*,*) '  QUICK -> self descriptive'
      else if(answer == 'MINMOD') then
        BLEND_TUR(m) = MINMOD
        if(this  < 2) write(*,*) '  MINMOD -> self descriptive'
      else if(answer == 'SMART') then
        BLEND_TUR(m) = SMART 
        if(this  < 2) write(*,*) '  SMART -> self descriptive'
      else if(answer == 'AVL_SMART') then
        BLEND_TUR(m) = AVL_SMART 
        if(this  < 2) write(*,*) '  AVL_SMART -> self descriptive'
      else if(answer == 'SUPERBEE') then
        BLEND_TUR(m) = SUPERBEE 
        if(this  < 2) write(*,*) '  SUPERBEE -> ???'
      else if(answer == 'GAMMA') then
        BLEND(m) = GAMMA 
        if(this  < 2) write(*,*) '  GAMMA -> ???'
      else
        if(this  < 2) write(*,*) '  Error in input in line:' 
        if(this  < 2) write(*,*) answer
        if(this  < 2) write(*,*) '  Exiting !'
        stop
      endif
    end do
  end if

!----- Solver parameters
  if(this  < 2) write(*,*) '# Preconditioning of the system matrix: '
  call ReadC(7,inp,tn,ts,te)
  read(inp(ts(1):te(1)),'(A)')  answer
  call ToUppr(answer)
  if(answer == 'NO') then
    PREC = 0 
    if(this  < 2) write(*,*) '  NO -> No preconditioning'
  else if(answer == 'DI') then
    PREC = 1 
    if(this  < 2) write(*,*) '  DI -> Diagonal preconditioning'
  else if(answer == 'IC') then
    PREC = 2 
    if(this  < 2) write(*,*) '  IC -> Incomplete Cholesky'
  else
    if(this  < 2) write(*,*) '  Error in input ! Exiting'
    stop
  endif

  if(this  < 2) write(*,*) '# Tolerance for velocity solver:'
    call ReadC(7,inp,tn,ts,te)
    read(inp,*)    U % STol
    if(this  < 2) write(*,*) '  ' , U % STol
    V % Stol     = U % Stol
    W % Stol     = U % Stol
  if(this  < 2) write(*,*) '# Tolerance for pressure solver:'
    call ReadC(7,inp,tn,ts,te)
    read(inp,*)   PP % STol
    if(this  < 2) write(*,*) '  ' , PP % STol
    P % Stol = PP % Stol
  if(SIMULA/=LES.and.SIMULA/=DNS) then
    call ReadC(7,inp,tn,ts,te)
    read(inp,*)    Kin % STol
    Eps % Stol   = Kin % Stol
    v_2 % Stol   = Kin % Stol
    f22 % Stol   = Kin % Stol
    VIS % Stol   = Kin % Stol
    uu  % Stol   = Kin % Stol
    vv  % Stol   = Kin % Stol
    ww  % Stol   = Kin % Stol
    uv  % Stol   = Kin % Stol
    uw  % Stol   = Kin % Stol
    vw  % Stol   = Kin % Stol
  end if
  if(HOT == YES) then
    if(this  < 2) write(*,*) '# Tolerance for temperature solver:'
    call ReadC(7,inp,tn,ts,te)
    read(inp,*)    T % STol
    if(this  < 2) write(*,*) '  ' , T % STol
  end if
 
  if(ALGOR == SIMPLE) then
    if(this  < 2) write(*,*) '# Tolerance for SIMPLE:'
    call ReadC(7,inp,tn,ts,te)
    read(inp,*)   SIMTol
    if(this  < 2) write(*,*) '  ' , SIMTol
  endif     

!----- Time step
  if(this  < 2) write(*,*) '# Time step:'
    call ReadC(7,inp,tn,ts,te)
  read(inp,*)   dt
  if(this  < 2) write(*,*) '  ' , dt

!----- Wall velocity 
  do m=1,Nmat
    if(this  < 2) write(*,*) '# Enter Pdrop (x, y, z) for domain ', m
    call ReadC(7,inp,tn,ts,te) 
    if(.not. restar) then 
      read(inp,*)  PdropX(m), PdropY(m), PdropZ(m)
      UTau(m) = sqrt(abs(PdropX(m))) ! delta=1, nu=1 
      VTau(m) = sqrt(abs(PdropY(m))) ! delta=1, nu=1 
      WTau(m) = sqrt(abs(PdropZ(m))) ! delta=1, nu=1 
    else
      read(inp,*)  dummy, dummy, dummy 
      if(this  < 2) write(*,*) '  ' ,dummy,dummy,dummy 
    end if     
  end do

!----- Mass fluxes
  do m=1,Nmat
    if(this  < 2) then
      write(*,*) '# Enter the wanted mass flux through domain  (type 0.0 to keep the pressure drop constant)', m
    endif
    call ReadC(7,inp,tn,ts,te) 
    if(.not. restar) read(inp,*)  FLUXoX(m), FLUXoY(m), FLUXoZ(m)
    if(restar)       read(inp,*)  dummy, dummy, dummy 
  end do

  call Wait 

  RETURN

  END SUBROUTINE ReaCom

