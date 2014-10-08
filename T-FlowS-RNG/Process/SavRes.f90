!======================================================================!
  SUBROUTINE SavRes(namAut)
!----------------------------------------------------------------------!
! Writes: NAME.restart                                                 !
! ~~~~~~~                                                              !
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER             :: c, s, m
  CHARACTER           :: namOut*80, answer*80
  CHARACTER, OPTIONAL :: namAut*(*)
!--------------------------------[CVS]---------------------------------!
!  $Id: SavRes.f90,v 1.33 2008/12/10 14:58:16 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/SavRes.f90,v $  
!======================================================================!

  if(PRESENT(namAut)) then
!---- save the name
    answer = name
    name = namAut
  else
    if(this  < 2)                                                     &
      write(*,*) '# Output restart file name [skip cancels]:'
    call ReadC(7,inp,tn,ts,te)
!->>> write(*,*) inp(1:300)
    read(inp(ts(1):te(1)), '(A80)')  namOut
    answer=namOut
    call ToUppr(answer) 

    if(answer == 'SKIP') return 

!---- save the name
    answer = name
    name = namOut
  end if

!-----------------------------!
!     Create restart file     !
!-----------------------------!
  call NamFil(this, namOut, '.restart', len_trim('.restart') )
  open(9, FILE=namOut, FORM='UNFORMATTED')
  if(this  < 2) write(6, *) '# Now creating the file:', namOut

!---- version
  write(9) 0.0  ! version

!---- 60 INTEGER parameters -----------------------------------------
  write(9)        0,      NbC,       NC,       NS,     Ndtt,    Nstat
  write(9)       Cm,        0,        0,        0,        0,        0
  write(9)    ALGOR,    INERT,   CONVEC,    CROSS,   DIFFUS,   SIMULA
  write(9)   POSPRO,  CHANNEL,     TEST,    OTHER,      HOT,        0
  write(9)    BLEND,        0,        0,     MODE,     PIPE,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
!--------------------------------------------------------------------

!---- 60 REAL parameters --------------------------------------
  write(9)     0.0,    0.0,    0.0,     xp,     yp,     zp  
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0     
  write(9)   ReTau,    0.0,    Cs0,    0.0,    0.0,    0.0
  write(9)      dt,   Time,  Kflow,    0.0,    0.0,    0.0
  write(9)   U%URF,  P%URF,   URFC, SIMTol, U%Stol,    0.0 
  write(9) PP%Stol,  T%URF, T%STol,    0.0,    0.0,    0.0
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0
!--------------------------------------------------------------

  write(9) (U % n(c),   c=-NbC,NC)
  write(9) (V % n(c),   c=-NbC,NC)
  write(9) (W % n(c),   c=-NbC,NC)
  write(9) (U % o(c),   c=1,NC)
  write(9) (V % o(c),   c=1,NC)
  write(9) (W % o(c),   c=1,NC)

  write(9) (U % C(c),   c=1,NC)
  write(9) (V % C(c),   c=1,NC)
  write(9) (W % C(c),   c=1,NC)
  write(9) (U % Co(c),  c=1,NC)
  write(9) (V % Co(c),  c=1,NC)
  write(9) (W % Co(c),  c=1,NC)

  write(9) (U % Do(c),  c=1,NC)
  write(9) (V % Do(c),  c=1,NC)
  write(9) (W % Do(c),  c=1,NC)

  write(9) (U % X(c),   c=1,NC)
  write(9) (V % X(c),   c=1,NC)
  write(9) (W % X(c),   c=1,NC)
  write(9) (U % Xo(c),  c=1,NC)
  write(9) (V % Xo(c),  c=1,NC)
  write(9) (W % Xo(c),  c=1,NC)

  write(9) (P % n(c),   c=-NbC,NC)
  write(9) (PP % n(c),  c=-NbC,NC)

  write(9) (Px(c),   c=-NbC,NC)
  write(9) (Py(c),   c=-NbC,NC)
  write(9) (Pz(c),   c=-NbC,NC)

  write(9) (VISt(c), c=-NbC,NC)

!---- Fluxes 
  write(9) (Flux(s), s=1,NS)

!---- Temperature
  if(HOT == YES) then
    write(9) (T % n(c),    c=-NbC,NC)
    write(9) (T % q(c),    c=-NbC,-1)
    write(9) (T % o(c),    c=1,NC)
    write(9) (T % C(c),    c=1,NC)
    write(9) (T % Co(c),   c=1,NC)
    write(9) (T % Do(c),   c=1,NC)
    write(9) (T % X(c),    c=1,NC)
    write(9) (T % Xo(c),   c=1,NC)
    if(SIMULA == LES .or. SIMULA == DNS.or.SIMULA == HYB_PITM &
       .or.SIMULA == HYB_ZETA) then
      write(9) (T % mean(c), c=-NbC,NC)
      write(9) (TT % mean(c), c=-NbC,NC)
      write(9) (uT % mean(c), c=-NbC,NC)
      write(9) (vT % mean(c), c=-NbC,NC)
      write(9) (wT % mean(c), c=-NbC,NC)
    end if
  end if

!---- LES and DNS
  if(SIMULA == DNS .or. SIMULA == LES) then
    write(9) (U % mean(c),  c=-NbC,NC)
    write(9) (V % mean(c),  c=-NbC,NC)
    write(9) (W % mean(c),  c=-NbC,NC)

    write(9) (uu % mean(c), c=-NbC,NC)
    write(9) (vv % mean(c), c=-NbC,NC)
    write(9) (ww % mean(c), c=-NbC,NC)
    write(9) (uv % mean(c), c=-NbC,NC)
    write(9) (uw % mean(c), c=-NbC,NC)
    write(9) (vw % mean(c), c=-NbC,NC)

    write(9) (P % mean(c),  c=1,NC)
    if (SIMULA == LES) then
      write(9) (VISt_mean(c), c=1,NC)
      write(9) (near(c),   c=-NbC,NC)
      if(MODE == DYN) write(9) (Cdyn_mean(c), c=1,NC)
    end if 
  end if

  if(SIMULA == HYB_PITM) then
    write(9) (U % mean(c),  c=-NbC,NC)
    write(9) (V % mean(c),  c=-NbC,NC)
    write(9) (W % mean(c),  c=-NbC,NC)
    write(9) (uu % mean(c), c=-NbC,NC)
    write(9) (vv % mean(c), c=-NbC,NC)
    write(9) (ww % mean(c), c=-NbC,NC)
    write(9) (uv % mean(c), c=-NbC,NC)
    write(9) (uw % mean(c), c=-NbC,NC)
    write(9) (vw % mean(c), c=-NbC,NC)

    write(9) (P % mean(c),  c=1,NC)
    !write(9) (VISt_mean(c), c=1,NC)
  end if

  if(SIMULA == HYB_ZETA .or. SIMULA == ZETA) then
    write(9) (U % mean(c),  c=1,NC)
    write(9) (V % mean(c),  c=1,NC)
    write(9) (W % mean(c),  c=1,NC)

    write(9) (uu % mean(c), c=1,NC)
    write(9) (vv % mean(c), c=1,NC)
    write(9) (ww % mean(c), c=1,NC)
    write(9) (uv % mean(c), c=1,NC)
    write(9) (uw % mean(c), c=1,NC)
    write(9) (vw % mean(c), c=1,NC)

    write(9) (P % mean(c),  c=1,NC)
    write(9) (VISt_mean(c), c=1,NC)

    write(9) (VISt_sgs(c),      c=1,NC)
    write(9) (VISt_sgs_mean(c), c=1,NC)
    write(9) (VISt_eff(c),      c=1,NC)
    write(9) (VISt_eff_mean(c), c=1,NC)
    write(9) (Kin % mean(c),    c=1,NC)
    write(9) (Eps % mean(c),    c=1,NC)

  end if


!---- Rans models
  if(SIMULA==EBM.or.SIMULA==HJ) then
    write(9) (uu  % n(c),    c=-NbC,NC)
    write(9) (uu  % o(c),    c=1,NC)
    write(9) (uu  % C(c),    c=1,NC)
    write(9) (uu  % Co(c),   c=1,NC)
    write(9) (uu  % Do(c),   c=1,NC)
    write(9) (uu  % X(c),    c=1,NC)
    write(9) (uu  % Xo(c),   c=1,NC)

    write(9) (vv  % n(c),    c=-NbC,NC)
    write(9) (vv  % o(c),    c=1,NC)
    write(9) (vv  % C(c),    c=1,NC)
    write(9) (vv  % Co(c),   c=1,NC)
    write(9) (vv  % Do(c),   c=1,NC)
    write(9) (vv  % X(c),    c=1,NC)
    write(9) (vv  % Xo(c),   c=1,NC)
    
    write(9) (ww  % n(c),    c=-NbC,NC)
    write(9) (ww  % o(c),    c=1,NC)
    write(9) (ww  % C(c),    c=1,NC)
    write(9) (ww  % Co(c),   c=1,NC)
    write(9) (ww  % Do(c),   c=1,NC)
    write(9) (ww  % X(c),    c=1,NC)
    write(9) (ww  % Xo(c),   c=1,NC)
    
    write(9) (uv  % n(c),    c=-NbC,NC)
    write(9) (uv  % o(c),    c=1,NC)
    write(9) (uv  % C(c),    c=1,NC)
    write(9) (uv  % Co(c),   c=1,NC)
    write(9) (uv  % Do(c),   c=1,NC)
    write(9) (uv  % X(c),    c=1,NC)
    write(9) (uv  % Xo(c),   c=1,NC)
    
    write(9) (uw  % n(c),    c=-NbC,NC)
    write(9) (uw  % o(c),    c=1,NC)
    write(9) (uw  % C(c),    c=1,NC)
    write(9) (uw  % Co(c),   c=1,NC)
    write(9) (uw  % Do(c),   c=1,NC)
    write(9) (uw  % X(c),    c=1,NC)
    write(9) (uw  % Xo(c),   c=1,NC)
    
    write(9) (vw  % n(c),    c=-NbC,NC)
    write(9) (vw  % o(c),    c=1,NC)
    write(9) (vw  % C(c),    c=1,NC)
    write(9) (vw  % Co(c),   c=1,NC)
    write(9) (vw  % Do(c),   c=1,NC)
    write(9) (vw  % X(c),    c=1,NC)
    write(9) (vw  % Xo(c),   c=1,NC)
  
    if(SIMULA==EBM) then  
      write(9) (f22 % n(c),    c=-NbC,NC)
      write(9) (f22 % o(c),    c=1,NC)
      write(9) (f22 % Do(c),   c=1,NC)
      write(9) (f22 % X(c),    c=1,NC)
      write(9) (f22 % Xo(c),   c=1,NC)
    end if

    write(9) (Eps % n(c),    c=-NbC,NC)
    write(9) (Eps % o(c),    c=1,NC)
    write(9) (Eps % C(c),    c=1,NC)
    write(9) (Eps % Co(c),   c=1,NC)
    write(9) (Eps % Do(c),   c=1,NC)
    write(9) (Eps % X(c),    c=1,NC)
    write(9) (Eps % Xo(c),   c=1,NC)
    
    write(9) (Pk(c),         c=-NbC,NC)
    write(9) (Kin % n(c),    c=-NbC,NC)
    write(9) (VISt(c),       c=-NbC,NC)
  end if

  if(SIMULA == K_EPS .or. SIMULA == HYB_PITM) then
    write(9) (KIN % n(c),    c=-NbC,NC)
    write(9) (KIN % o(c),    c=1,NC)
    write(9) (KIN % C(c),    c=1,NC)
    write(9) (KIN % Co(c),   c=1,NC)
    write(9) (KIN % Do(c),   c=1,NC)
    write(9) (KIN % X(c),    c=1,NC)
    write(9) (KIN % Xo(c),   c=1,NC)
    write(9) (EPS % n(c),    c=-NbC,NC)
    write(9) (EPS % o(c),    c=1,NC)
    write(9) (EPS % C(c),    c=1,NC)
    write(9) (EPS % Co(c),   c=1,NC)
    write(9) (EPS % Do(c),   c=1,NC)
    write(9) (EPS % X(c),    c=1,NC)
    write(9) (EPS % Xo(c),   c=1,NC)
    
    write(9) (Pk(c),         c=-NbC,NC)
    write(9) (Uf(c),         c=-NbC,NC)
    write(9) (Ynd(c),        c=-NbC,NC)
    write(9) (VISwall(c),    c=-NbC,NC)
  end if

  if(SIMULA==K_EPS_VV.or.SIMULA==ZETA.or.SIMULA==HYB_ZETA) then

    write(9) (Kin%n(c),    c=-NbC,NC)
    write(9) (Eps%n(c),    c=-NbC,NC)
    write(9) (v_2%n(c),    c=-NbC,NC)
    write(9) (f22%n(c),    c=-NbC,NC)
 
    write(9) (Kin%o(c),   c=1,NC)
    write(9) (Eps%o(c),   c=1,NC)
    write(9) (v_2%o(c),   c=1,NC)
    write(9) (f22%o(c),   c=1,NC)
 
    write(9) (Kin%C(c),   c=1,NC)
    write(9) (Eps%C(c),   c=1,NC)
    write(9) (v_2%C(c),   c=1,NC)
 
    write(9) (Kin%Co(c),  c=1,NC)
    write(9) (Eps%Co(c),  c=1,NC)
    write(9) (v_2%Co(c),  c=1,NC)
 
    write(9) (Kin%Do(c),  c=1,NC)
    write(9) (Eps%Do(c),  c=1,NC)
    write(9) (v_2%Do(c),  c=1,NC)
    write(9) (f22%Do(c),  c=1,NC)
 
    write(9) (Kin%X(c),   c=1,NC)
    write(9) (Eps%X(c),   c=1,NC)
    write(9) (v_2%X(c),   c=1,NC)
    write(9) (f22%X(c),   c=1,NC)
 
    write(9) (Kin%Xo(c),  c=1,NC)
    write(9) (Eps%Xo(c),  c=1,NC)
    write(9) (v_2%Xo(c),  c=1,NC)
    write(9) (f22%Xo(c),  c=1,NC)
 
    write(9) (Pk(c),     c=-NbC,NC)
    write(9) (Tsc(c),    c=-NbC,NC)
    write(9) (Lsc(c),    c=-NbC,NC)

    write(9) (Uf(c),     c=-NbC,NC)
    write(9) (Ynd(c),    c=-NbC,NC)
    write(9) (VISwall(c),c=-NbC,NC)
  end if

if(SIMULA==ZETA .or. SIMULA == HYB_ZETA) then
!!!! new variables
!coherent variables
!for cylinder task
    write(9) (Ucoh0(c), c=1,NC)
    write(9) (Vcoh0(c), c=1,NC)
    write(9) (Wcoh0(c), c=1,NC)
    write(9) (VIStcoh0(c), c=1,NC)
    write(9) (Vortcoh0(c), c=1,NC)

    write(9) (Ucoh13MIN(c), c=1,NC)
    write(9) (Vcoh13MIN(c), c=1,NC)
    write(9) (Wcoh13MIN(c), c=1,NC)
    write(9) (VIStcoh13MIN(c), c=1,NC)
    write(9) (Vortcoh13MIN(c), c=1,NC)

    write(9) (Ucoh23MIN(c), c=1,NC)
    write(9) (Vcoh23MIN(c), c=1,NC)
    write(9) (Wcoh23MIN(c), c=1,NC)
    write(9) (VIStcoh23MIN(c), c=1,NC)
    write(9) (Vortcoh23MIN(c), c=1,NC)

    write(9) (UcohMIN(c), c=1,NC)
    write(9) (VcohMIN(c), c=1,NC)
    write(9) (WcohMIN(c), c=1,NC)
    write(9) (VIStcohMIN(c), c=1,NC)
    write(9) (VortcohMIN(c), c=1,NC)

    write(9) (UcohMINCanditate(c), c=1,NC)
    write(9) (VcohMINCanditate(c), c=1,NC)
    write(9) (WcohMINCanditate(c), c=1,NC)
    write(9) (VIStcohMINCanditate(c), c=1,NC)
    write(9) (VortcohMINCanditate(c), c=1,NC)

    write(9) (Fpy_min)
    write(9) (Fpy_prev)
    write(9) (kMIN)
    write(9) (k23MIN)
    write(9) (k13MIN)
    write(9) (k0)
    !currently code running without saving variables below(it will deny one cycle)
    write(9) (CounterFromGlobalMinimumCandidate)
    write(9) (MinFoundOnThisCycle)

end if

 if(SIMULA == SPA_ALL) then
    write(9) (VIS % n(c),    c=-NbC,NC)
    write(9) (VIS % o(c),    c=1,NC)
    write(9) (VIS % C(c),    c=1,NC)
    write(9) (VIS % Co(c),   c=1,NC)
    write(9) (VIS % Do(c),   c=1,NC)
    write(9) (VIS % X(c),    c=1,NC)
    write(9) (VIS % Xo(c),   c=1,NC)

    write(9) (Vort(c),       c=-NbC,NC)
  end if

 if(SIMULA == DES_SPA) then
    write(9) (VIS % n(c),    c=-NbC,NC)
    write(9) (VIS % o(c),    c=1,NC)
    write(9) (VIS % C(c),    c=1,NC)
    write(9) (VIS % Co(c),   c=1,NC)
    write(9) (VIS % Do(c),   c=1,NC)
    write(9) (VIS % X(c),    c=1,NC)
    write(9) (VIS % Xo(c),   c=1,NC)

    write(9) (Vort(c),       c=-NbC,NC)

    write(9) (U % mean(c),  c=-NbC,NC)
    write(9) (V % mean(c),  c=-NbC,NC)
    write(9) (W % mean(c),  c=-NbC,NC)
    write(9) (uu % mean(c), c=-NbC,NC)
    write(9) (vv % mean(c), c=-NbC,NC)
    write(9) (ww % mean(c), c=-NbC,NC)
    write(9) (uv % mean(c), c=-NbC,NC)
    write(9) (uw % mean(c), c=-NbC,NC)
    write(9) (vw % mean(c), c=-NbC,NC)

    write(9) (P % mean(c),  c=1,NC)
    write(9) (VISt_mean(c), c=1,NC)
  end if



!---- Pressure drops in each material (domain)
  do m=1,Nmat
    write(9) PdropX(m), PdropY(m), PdropZ(m)
    write(9) FLUXoX(m), FLUXoY(m), FLUXoZ(m)
    write(9) FLUXx(m),  FLUXy(m),  FLUXz(m)
    write(9) AreaX(m),  AreaY(m),  AreaZ(m)
    write(9) Ubulk(m),  Vbulk(m),  Wbulk(m)
  end do

  close(9)

!---- restore the name
  name = answer 

  END SUBROUTINE SavRes
