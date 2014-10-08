!======================================================================!
  SUBROUTINE UnkAloc
!----------------------------------------------------------------------!
! Allocates the memory for unknowns. It is called either from LoaRes   !
! or from Processor.                                                   !
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!--------------------------------[CVS]---------------------------------!
!  $Id: UnkAloc.f90,v 1.28 2008/11/19 14:55:56 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/UnkAloc.f90,v $  
!======================================================================!


  allocate (U % n(-NbC:NC)); U % n=0.
  allocate (V % n(-NbC:NC)); V % n=0.
  allocate (W % n(-NbC:NC)); W % n=0.

  allocate (U % mean(-NbC:NC));  U % mean =0.
  allocate (V % mean(-NbC:NC));  V % mean =0.
  allocate (W % mean(-NbC:NC));  W % mean =0.

  allocate (P % n(-NbC:NC));  P % n=0.
  allocate (P % mean(-NbC:NC));  P % mean=0.
  allocate (PP % n(-NbC:NC)); PP % n=0.

  allocate (Px(-NbC:NC)); Px=0.
  allocate (Py(-NbC:NC)); Py=0.
  allocate (Pz(-NbC:NC)); Pz=0.  

  allocate (Ux(-NbC:NC)); Ux=0.
  allocate (Uy(-NbC:NC)); Uy=0.
  allocate (Uz(-NbC:NC)); Uz=0.  

  allocate (Vx(-NbC:NC)); Vx=0.
  allocate (Vy(-NbC:NC)); Vy=0.
  allocate (Vz(-NbC:NC)); Vz=0.  

  allocate (Wx(-NbC:NC)); Wx=0.
  allocate (Wy(-NbC:NC)); Wy=0.
  allocate (Wz(-NbC:NC)); Wz=0.  

  allocate (Ps(NS)); Ps=0.;

  allocate (PHI1x(-NbC:NC)); PHI1x=0.
  allocate (PHI1y(-NbC:NC)); PHI1y=0.
  allocate (PHI1z(-NbC:NC)); PHI1z=0.

  allocate (PHI2x(-NbC:NC)); PHI2x=0.
  allocate (PHI2y(-NbC:NC)); PHI2y=0.
  allocate (PHI2z(-NbC:NC)); PHI2z=0.

  allocate (PHI3x(-NbC:NC)); PHI3x=0.
  allocate (PHI3y(-NbC:NC)); PHI3y=0.
  allocate (PHI3z(-NbC:NC)); PHI3z=0.

  allocate (PHI4x(-NbC:NC)); PHI4x=0.
  allocate (PHI4y(-NbC:NC)); PHI4y=0.
  allocate (PHI4z(-NbC:NC)); PHI4z=0.

  allocate (PHIx(-NbC:NC)); PHIx=0.
  allocate (PHIy(-NbC:NC)); PHIy=0.
  allocate (PHIz(-NbC:NC)); PHIz=0.
  allocate (PHIside(NS)); PHIside=0.

  allocate (PHImax(-NbC:NC)); PHImax=0.
  allocate (PHImin(-NbC:NC)); PHImin=0.

  allocate (G(6,NC)); G=0

  allocate (U % o(NC));  U % o=0.
  allocate (V % o(NC));  V % o=0.
  allocate (W % o(NC));  V % o=0.
  allocate (U % oo(NC)); U % oo=0.
  allocate (V % oo(NC)); V % oo=0.
  allocate (W % oo(NC)); W % oo=0.

  allocate (U % C(NC));   U % C=0.
  allocate (V % C(NC));   V % C=0.
  allocate (W % C(NC));   W % C=0.
  allocate (U % Co(NC));  U % Co=0.
  allocate (V % Co(NC));  V % Co=0.
  allocate (W % Co(NC));  W % Co=0.
  allocate (U % Coo(NC)); U % Coo=0.
  allocate (V % Coo(NC)); V % Coo=0.
  allocate (W % Coo(NC)); W % Coo=0.

  allocate (U % Do(NC));  U % Do=0.
  allocate (V % Do(NC));  V % Do=0.
  allocate (W % Do(NC));  W % Do=0.
  allocate (U % Doo(NC)); U % Doo=0.
  allocate (V % Doo(NC)); V % Doo=0.
  allocate (W % Doo(NC)); W % Doo=0.

  allocate (U % X(NC));   U % X=0.
  allocate (V % X(NC));   V % X=0.
  allocate (W % X(NC));   W % X=0.
  allocate (U % Xo(NC));  U % Xo=0.
  allocate (V % Xo(NC));  V % Xo=0.
  allocate (W % Xo(NC));  W % Xo=0.
  allocate (U % Xoo(NC)); U % Xoo=0.
  allocate (V % Xoo(NC)); V % Xoo=0.
  allocate (W % Xoo(NC)); W % Xoo=0.

  allocate (Flux(NS));     Flux=0.

  allocate (PdropX(Nmat)); PdropX=0.0
  allocate (PdropY(Nmat)); PdropY=0.0 
  allocate (PdropZ(Nmat)); PdropZ=0.0 

  allocate (Utau(Nmat));   Utau=0.0
  allocate (Vtau(Nmat));   Vtau=0.0
  allocate (Wtau(Nmat));   Wtau=0.0

  allocate (FLUXx(Nmat));  FLUXx=0.0
  allocate (FLUXy(Nmat));  FLUXy=0.0
  allocate (FLUXz(Nmat));  FLUXz=0.0

  allocate (FLUXoX(Nmat)); FLUXoX=0.0
  allocate (FLUXoY(Nmat)); FLUXoY=0.0
  allocate (FLUXoZ(Nmat)); FLUXoZ=0.0

  allocate (Ubulk(Nmat));  Ubulk=0.0
  allocate (Vbulk(Nmat));  Vbulk=0.0
  allocate (Wbulk(Nmat));  Wbulk=0.0

  allocate (MassIn(Nmat)); MassIn=0.0
  allocate (MasOut(Nmat)); MasOut=0.0

  allocate (BadForG(NC));  BadForG = .FALSE.
  allocate (NumGood(NC));  NumGood = 0          
  allocate (NumNeig(NC));  NumNeig = 0         

  allocate (near(-NbC:NC));  near  = 0.
  allocate (VISwall(-NbC:NC)); VISwall =0.0

!---- variables for temperature
  if(HOT==YES) then
    allocate (T % n(-NbC:NC)); T % n=0.
    allocate (T % o(NC));      T % o=0.
    allocate (T % oo(NC));     T % oo=0.
    allocate (T % C(NC));      T % C=0.
    allocate (T % Co(NC));     T % Co=0.
    allocate (T % Coo(NC));    T % Coo=0.
    allocate (T % Do(NC));     T % Do=0.
    allocate (T % Doo(NC));    T % Doo=0.
    allocate (T % X(NC));      T % X=0.
    allocate (T % Xo(NC));     T % Xo=0.
    allocate (T % Xoo(NC));    T % Xoo=0.
    allocate (T % q(-NbC:-1)); T % q=0.
  end if

  if(SIMULA==EBM.or.SIMULA==HJ) then
    allocate (VAR1x(-NbC:NC)); VAR1x=0.
    allocate (VAR1y(-NbC:NC)); VAR1y=0.
    allocate (VAR1z(-NbC:NC)); VAR1z=0.

    allocate (VAR2x(-NbC:NC)); VAR2x=0.
    allocate (VAR2y(-NbC:NC)); VAR2y=0.
    allocate (VAR2z(-NbC:NC)); VAR2z=0.

    allocate (uu % n(-NbC:NC)); uu % n=0.
    allocate (uu % o(NC));      uu % o=0.
    allocate (uu % oo(NC));     uu % oo=0.
    allocate (uu % C(NC));      uu % C=0.
    allocate (uu % Co(NC));     uu % Co=0.
    allocate (uu % Coo(NC));    uu % Coo=0.
    allocate (uu % Do(NC));     uu % Do=0.
    allocate (uu % Doo(NC));    uu % Doo=0.
    allocate (uu % X(NC));      uu % X=0.
    allocate (uu % Xo(NC));     uu % Xo=0.
    allocate (uu % Xoo(NC));    uu % Xoo=0.
    allocate (uu % mean(NC));   uu % mean=0.

    allocate (vv % n(-NbC:NC)); vv % n=0.
    allocate (vv % o(NC));      vv % o=0.
    allocate (vv % oo(NC));     vv % oo=0.
    allocate (vv % C(NC));      vv % C=0.
    allocate (vv % Co(NC));     vv % Co=0.
    allocate (vv % Coo(NC));    vv % Coo=0.
    allocate (vv % Do(NC));     vv % Do=0.
    allocate (vv % Doo(NC));    vv % Doo=0.
    allocate (vv % X(NC));      vv % X=0.
    allocate (vv % Xo(NC));     vv % Xo=0.
    allocate (vv % Xoo(NC));    vv % Xoo=0.
    allocate (vv % mean(NC));   vv % mean=0.

    allocate (ww % n(-NbC:NC)); ww % n=0.
    allocate (ww % o(NC));      ww % o=0.
    allocate (ww % oo(NC));     ww % oo=0.
    allocate (ww % C(NC));      ww % C=0.
    allocate (ww % Co(NC));     ww % Co=0.
    allocate (ww % Coo(NC));    ww % Coo=0.
    allocate (ww % Do(NC));     ww % Do=0.
    allocate (ww % Doo(NC));    ww % Doo=0.
    allocate (ww % X(NC));      ww % X=0.
    allocate (ww % Xo(NC));     ww % Xo=0.
    allocate (ww % Xoo(NC));    ww % Xoo=0.
    allocate (ww % mean(NC));   ww % mean=0.

    allocate (uv % n(-NbC:NC)); uv % n=0.
    allocate (uv % o(NC));      uv % o=0.
    allocate (uv % oo(NC));     uv % oo=0.
    allocate (uv % C(NC));      uv % C=0.
    allocate (uv % Co(NC));     uv % Co=0.
    allocate (uv % Coo(NC));    uv % Coo=0.
    allocate (uv % Do(NC));     uv % Do=0.
    allocate (uv % Doo(NC));    uv % Doo=0.
    allocate (uv % X(NC));      uv % X=0.
    allocate (uv % Xo(NC));     uv % Xo=0.
    allocate (uv % Xoo(NC));    uv % Xoo=0.
    allocate (uv % mean(NC));   uv % mean=0.

    allocate (uw % n(-NbC:NC)); uw % n=0.
    allocate (uw % o(NC));      uw % o=0.
    allocate (uw % oo(NC));     uw % oo=0.
    allocate (uw % C(NC));      uw % C=0.
    allocate (uw % Co(NC));     uw % Co=0.
    allocate (uw % Coo(NC));    uw % Coo=0.
    allocate (uw % Do(NC));     uw % Do=0.
    allocate (uw % Doo(NC));    uw % Doo=0.
    allocate (uw % X(NC));      uw % X=0.
    allocate (uw % Xo(NC));     uw % Xo=0.
    allocate (uw % Xoo(NC));    uw % Xoo=0.
    allocate (uw % mean(NC));   uw % mean=0.

    allocate (vw % n(-NbC:NC)); vw % n=0.
    allocate (vw % o(NC));      vw % o=0.
    allocate (vw % oo(NC));     vw % oo=0.
    allocate (vw % C(NC));      vw % C=0.
    allocate (vw % Co(NC));     vw % Co=0.
    allocate (vw % Coo(NC));    vw % Coo=0.
    allocate (vw % Do(NC));     vw % Do=0.
    allocate (vw % Doo(NC));    vw % Doo=0.
    allocate (vw % X(NC));      vw % X=0.
    allocate (vw % Xo(NC));     vw % Xo=0.
    allocate (vw % Xoo(NC));    vw % Xoo=0.
    allocate (vw % mean(NC));   vw % mean=0.

    allocate (Tsc(-NbC:NC));     Tsc   =0.0
    allocate (Lsc(-NbC:NC));     Lsc   =0.0
    allocate (Kin % n(-NbC:NC)); Kin % n =0.0
    allocate (Kin % mean(NC));   Kin % mean =0.0
    allocate (Pk(-NbC:NC));      Pk    =0.0

    allocate (f22 % n(-NbC:NC)); f22 % n=0.
    if(SIMULA==EBM) then
      allocate (f22 % o(NC));      f22 % o=0.
      allocate (f22 % oo(NC));     f22 % oo=0.
      allocate (f22 % Do(NC));     f22 % Do=0.
      allocate (f22 % Doo(NC));    f22 % Doo=0.
      allocate (f22 % X(NC));      f22 % X=0.
      allocate (f22 % Xo(NC));     f22 % Xo=0.
      allocate (f22 % Xoo(NC));    f22 % Xoo=0.
    end if

    allocate (Eps % n(-NbC:NC)); Eps % n=0.
    allocate (Eps % o(NC));      Eps % o=0.
    allocate (Eps % oo(NC));     Eps % oo=0.
    allocate (Eps % C(NC));      Eps % C=0.
    allocate (Eps % Co(NC));     Eps % Co=0.
    allocate (Eps % Coo(NC));    Eps % Coo=0.
    allocate (Eps % Do(NC));     Eps % Do=0.
    allocate (Eps % Doo(NC));    Eps % Doo=0.
    allocate (Eps % X(NC));      Eps % X=0.
    allocate (Eps % Xo(NC));     Eps % Xo=0.
    allocate (Eps % Xoo(NC));    Eps % Xoo=0.
    allocate (Eps % mean(NC));   Eps % mean=0.
    allocate (Ls  % mean(NC));   Ls  % mean=0.
    allocate (VISt_mean(NC)); VISt_mean = 0.0

  end if
!---- variables for Rans models
  if(SIMULA==K_EPS.or.SIMULA == HYB_PITM) then
    allocate (Kin % n(-NbC:NC)); Kin % n=0.
    allocate (Kin % o(NC));      Kin % o=0.
    allocate (Kin % oo(NC));     Kin % oo=0.
    allocate (Kin % C(NC));      Kin % C=0.
    allocate (Kin % Co(NC));     Kin % Co=0.
    allocate (Kin % Coo(NC));    Kin % Coo=0.
    allocate (Kin % Do(NC));     Kin % Do=0.
    allocate (Kin % Doo(NC));    Kin % Doo=0.
    allocate (Kin % X(NC));      Kin % X=0.
    allocate (Kin % Xo(NC));     Kin % Xo=0.
    allocate (Kin % Xoo(NC));    Kin % Xoo=0.
    allocate (Kin % mean(NC));   Kin % mean=0.


    allocate (Eps % n(-NbC:NC)); Eps % n=0.
    allocate (Eps % o(NC));      Eps % o=0.
    allocate (Eps % oo(NC));     Eps % oo=0.
    allocate (Eps % C(NC));      Eps % C=0.
    allocate (Eps % Co(NC));     Eps % Co=0.
    allocate (Eps % Coo(NC));    Eps % Coo=0.
    allocate (Eps % Do(NC));     Eps % Do=0.
    allocate (Eps % Doo(NC));    Eps % Doo=0.
    allocate (Eps % X(NC));      Eps % X=0.
    allocate (Eps % Xo(NC));     Eps % Xo=0.
    allocate (Eps % Xoo(NC));    Eps % Xoo=0.
    allocate (Eps % mean(NC));   Eps % mean=0.

    allocate (Uf(-NbC:NC));      Uf    =0.0
    allocate (Ufmean(-NbC:NC));  Ufmean=0.0
    allocate (Pk(-NbC:NC));      Pk    =0.0
    allocate (Ynd(-NbC:NC));     Ynd   =0.0
  end if

  if(SIMULA==K_EPS_VV.or.SIMULA==ZETA.or.SIMULA==HYB_ZETA) then
    allocate (Kin % n(-NbC:NC)); Kin % n=0.
    allocate (Kin % o(NC));      Kin % o=0.
    allocate (Kin % oo(NC));     Kin % oo=0.
    allocate (Kin % C(NC));      Kin % C=0.
    allocate (Kin % Co(NC));     Kin % Co=0.
    allocate (Kin % Coo(NC));    Kin % Coo=0.
    allocate (Kin % Do(NC));     Kin % Do=0.
    allocate (Kin % Doo(NC));    Kin % Doo=0.
    allocate (Kin % X(NC));      Kin % X=0.
    allocate (Kin % Xo(NC));     Kin % Xo=0.
    allocate (Kin % Xoo(NC));    Kin % Xoo=0.
    allocate (Kin % mean(NC));   Kin % mean=0.
 
 
    allocate (Eps % n(-NbC:NC)); Eps % n=0.
    allocate (Eps % o(NC));      Eps % o=0.
    allocate (Eps % oo(NC));     Eps % oo=0.
    allocate (Eps % C(NC));      Eps % C=0.
    allocate (Eps % Co(NC));     Eps % Co=0.
    allocate (Eps % Coo(NC));    Eps % Coo=0.
    allocate (Eps % Do(NC));     Eps % Do=0.
    allocate (Eps % Doo(NC));    Eps % Doo=0.
    allocate (Eps % X(NC));      Eps % X=0.
    allocate (Eps % Xo(NC));     Eps % Xo=0.
    allocate (Eps % Xoo(NC));    Eps % Xoo=0.
    allocate (Eps % mean(NC));   Eps % mean=0.

    allocate (v_2 % n(-NbC:NC)); v_2 % n=0.
    allocate (v_2 % o(NC));      v_2 % o=0.
    allocate (v_2 % oo(NC));     v_2 % oo=0.
    allocate (v_2 % C(NC));      v_2 % C=0.
    allocate (v_2 % Co(NC));     v_2 % Co=0.
    allocate (v_2 % Coo(NC));    v_2 % Coo=0.
    allocate (v_2 % Do(NC));     v_2 % Do=0.
    allocate (v_2 % Doo(NC));    v_2 % Doo=0.
    allocate (v_2 % X(NC));      v_2 % X=0.
    allocate (v_2 % Xo(NC));     v_2 % Xo=0.
    allocate (v_2 % Xoo(NC));    v_2 % Xoo=0.
    allocate (v_2 % mean(NC));   v_2 % mean=0.

    allocate (f22 % n(-NbC:NC)); f22 % n=0.
    allocate (f22 % o(NC));      f22 % o=0.
    allocate (f22 % oo(NC));     f22 % oo=0.
    allocate (f22 % C(NC));      f22 % C=0.
    allocate (f22 % Co(NC));     f22 % Co=0.
    allocate (f22 % Coo(NC));    f22 % Coo=0.
    allocate (f22 % Do(NC));     f22 % Do=0.
    allocate (f22 % Doo(NC));    f22 % Doo=0.
    allocate (f22 % X(NC));      f22 % X=0.
    allocate (f22 % Xo(NC));     f22 % Xo=0.
    allocate (f22 % Xoo(NC));    f22 % Xoo=0.
    allocate (f22 % mean(NC));   f22 % mean=0.
 
    allocate (Tsc(-NbC:NC));     Tsc   =0.0
    allocate (Lsc(-NbC:NC));     Lsc   =0.0
    allocate (Uf(-NbC:NC));      Uf    =0.0
    allocate (Ufmean(-NbC:NC));  Ufmean=0.0  
    allocate (Pk(-NbC:NC));      Pk    =0.0  
    allocate (Ynd(-NbC:NC));     Ynd   =0.0
  end if                    

  if(SIMULA == HYB_ZETA) then
    allocate (VISt_sgs(-NbC:NC));  VISt_sgs=0.
    allocate (VISt_eff(-NbC:NC));  VISt_eff=0.
  end if

  if(SIMULA == DES_SPA) then
    allocate (Ksgs(-NbC:NC));  Ksgs=0.
  end if

  if(SIMULA == SPA_ALL.or.SIMULA == DES_SPA) then
    allocate (VIS % n(-NbC:NC)); VIS % n=0.
    allocate (VIS % o(NC));      VIS % o=0.
    allocate (VIS % oo(NC));     VIS % oo=0.
    allocate (VIS % C(NC));      VIS % C=0.
    allocate (VIS % Co(NC));     VIS % Co=0.
    allocate (VIS % Coo(NC));    VIS % Coo=0.
    allocate (VIS % Do(NC));     VIS % Do=0.
    allocate (VIS % Doo(NC));    VIS % Doo=0.
    allocate (VIS % X(NC));      VIS % X=0.
    allocate (VIS % Xo(NC));     VIS % Xo=0.
    allocate (VIS % Xoo(NC));    VIS % Xoo=0.
  end if

  if(SIMULA == DES_SPA) then
    allocate (VIS % mean(NC));   VIS % mean=0.
  end if

!---- variables defined in les_mod.h90:
  if(SIMULA == LES.or.SIMULA==HYB_ZETA) then
    if(MODE == WALE) then 
      allocate (WALEv(-NbC:NC));  WALEv =0.
    end if
    if(MODE == DYN.or.MODE==MIX) then 
      allocate (U % filt(-NbC:NC));  U % filt =0.
      allocate (V % filt(-NbC:NC));  V % filt =0.
      allocate (W % filt(-NbC:NC));  W % filt =0.
   
      allocate (Cdyn(-NbC:NC)); Cdyn = 0
      allocate(UUf(NC));   UUf = 0.0
      allocate(VVf(NC));   VVf = 0.0
      allocate(WWf(NC));   WWf = 0.0
      allocate(UVf(NC));   UVf = 0.0
      allocate(UWf(NC));   UWf = 0.0
      allocate(VWf(NC));   VWf = 0.0

      allocate(M11f(NC));   M11f = 0.0
      allocate(M22f(NC));   M22f = 0.0
      allocate(M33f(NC));   M33f = 0.0
      allocate(M12f(NC));   M12f = 0.0
      allocate(M13f(NC));   M13f = 0.0
      allocate(M23f(NC));   M23f = 0.0
    end if 
    allocate(ShearTest(-NbC:NC));   ShearTest = 0.0
    allocate (Ksgs(-NbC:NC));  Ksgs=0.
    allocate (Cdyn_mean(-NbC:NC)); Cdyn_mean = 0
  end if

  if(SIMULA == LES.or.SIMULA==DNS.or.SIMULA==DES_SPA) then
    allocate (uu % mean(-NbC:NC)); uu % mean=0.
    allocate (vv % mean(-NbC:NC)); vv % mean=0.
    allocate (ww % mean(-NbC:NC)); ww % mean=0.
    allocate (uv % mean(-NbC:NC)); uv % mean=0.
    allocate (uw % mean(-NbC:NC)); uw % mean=0.
    allocate (vw % mean(-NbC:NC)); vw % mean=0.

    allocate (uuu % mean(-NbC:NC)); uuu % mean=0.
    allocate (uuv % mean(-NbC:NC)); uuv % mean=0.
    allocate (uuw % mean(-NbC:NC)); uuw % mean=0.

    allocate (vvu % mean(-NbC:NC)); vvu % mean=0.
    allocate (vvv % mean(-NbC:NC)); vvv % mean=0.
    allocate (vvw % mean(-NbC:NC)); vvw % mean=0.

    allocate (wwu % mean(-NbC:NC)); wwu % mean=0.
    allocate (wwv % mean(-NbC:NC)); wwv % mean=0.
    allocate (www % mean(-NbC:NC)); www % mean=0.

    allocate (uwu % mean(-NbC:NC)); uwu % mean=0.
    allocate (uwv % mean(-NbC:NC)); uwv % mean=0.
    allocate (uww % mean(-NbC:NC)); uww % mean=0.

    if(BUDG==YES) then
    allocate (uu % n(-NbC:NC)); uu % n=0.
    allocate (vv % n(-NbC:NC)); vv % n=0.
    allocate (ww % n(-NbC:NC)); ww % n=0.
    allocate (uv % n(-NbC:NC)); uv % n=0.
    allocate (uw % n(-NbC:NC)); uw % n=0.
    allocate (vw % n(-NbC:NC)); vw % n=0.

    allocate (U % fluc(-NbC:NC));  U % fluc =0.
    allocate (V % fluc(-NbC:NC));  V % fluc =0.
    allocate (W % fluc(-NbC:NC));  W % fluc =0.
    allocate (P % fluc(-NbC:NC));  P % fluc =0.

    allocate (Puu_mean(1:NC));  Puu_mean =0.
    allocate (Pvv_mean(1:NC));  Pvv_mean =0.
    allocate (Pww_mean(1:NC));  Pww_mean =0.
    allocate (Puv_mean(1:NC));  Puv_mean =0.
    allocate (Puw_mean(1:NC));  Puw_mean =0.
    allocate (Pvw_mean(1:NC));  Pvw_mean =0.


    allocate (Diss_uu_mean(1:NC)); Diss_uu_mean  =0.
    allocate (Diss_vv_mean(1:NC)); Diss_vv_mean  =0.
    allocate (Diss_ww_mean(1:NC)); Diss_ww_mean  =0.
    allocate (Diss_uv_mean(1:NC)); Diss_uv_mean  =0.
    allocate (Diss_uw_mean(1:NC)); Diss_uw_mean  =0.
    allocate (Diss_vw_mean(1:NC)); Diss_vw_mean  =0.

    allocate (Diss_sgs_mean(1:NC)); Diss_sgs_mean  =0.


  allocate (PHI5x(-NbC:NC)); PHI5x=0.
  allocate (PHI5y(-NbC:NC)); PHI5y=0.
  allocate (PHI5z(-NbC:NC)); PHI5z=0.

  allocate (PHI6x(-NbC:NC)); PHI6x=0.
  allocate (PHI6y(-NbC:NC)); PHI6y=0.
  allocate (PHI6z(-NbC:NC)); PHI6z=0.

  allocate (PHI7x(-NbC:NC)); PHI7x=0.
  allocate (PHI7y(-NbC:NC)); PHI7y=0.
  allocate (PHI7z(-NbC:NC)); PHI7z=0.

  allocate (PHI8x(-NbC:NC)); PHI8x=0.
  allocate (PHI8y(-NbC:NC)); PHI8y=0.
  allocate (PHI8z(-NbC:NC)); PHI8z=0.

  allocate (PHI9x(-NbC:NC)); PHI9x=0.
  allocate (PHI9y(-NbC:NC)); PHI9y=0.
  allocate (PHI9z(-NbC:NC)); PHI9z=0.

  allocate (PHI10x(-NbC:NC)); PHI10x=0.
  allocate (PHI10y(-NbC:NC)); PHI10y=0.
  allocate (PHI10z(-NbC:NC)); PHI10z=0.

  allocate (Kx(-NbC:NC)); Kx=0.

    if(HOT==YES) then
    allocate (Put_mean(1:NC));  Put_mean =0.
    allocate (Pvt_mean(1:NC));  Pvt_mean =0.
    allocate (Pwt_mean(1:NC));  Pwt_mean =0.
    allocate (Ptt_mean(1:NC));  Ptt_mean =0.
    allocate (Difv_ut_tot(1:NC));   Difv_ut_tot=0.
    allocate (Difv_vt_tot(1:NC));   Difv_vt_tot=0.
    allocate (Difv_wt_tot(1:NC));   Difv_wt_tot=0.
    allocate (Diss_ut_mean(1:NC)); Diss_ut_mean  =0.
    allocate (Diss_vt_mean(1:NC)); Diss_vt_mean  =0.
    allocate (Diss_wt_mean(1:NC)); Diss_wt_mean  =0.
    allocate (Diss_tt_mean(1:NC)); Diss_tt_mean  =0.
    allocate (Dift_ut_mean(1:NC));   Dift_ut_mean=0.
    allocate (Dift_vt_mean(1:NC));   Dift_vt_mean=0.
    allocate (Dift_wt_mean(1:NC));   Dift_wt_mean=0.
    allocate (Dift_tt_mean(1:NC));   Dift_tt_mean=0.
    allocate (Difv_ut_mean(1:NC));   Difv_ut_mean=0.
    allocate (Difv_vt_mean(1:NC));   Difv_vt_mean=0.
    allocate (Difv_wt_mean(1:NC));   Difv_wt_mean=0.
    allocate (Difv_tt_mean(1:NC));   Difv_tt_mean=0.
    allocate (C_ut_mean(1:NC));   C_ut_mean=0.
    allocate (C_vt_mean(1:NC));   C_vt_mean=0.
    allocate (C_wt_mean(1:NC));   C_wt_mean=0.
    allocate (C_tt_mean(1:NC));   C_tt_mean=0.
    allocate (PD_ut_mean(-NbC:NC));   PD_ut_mean=0.
    allocate (PD_vt_mean(-NbC:NC));   PD_vt_mean=0.
    allocate (PD_wt_mean(-NbC:NC));   PD_wt_mean=0.
    allocate (PR_ut_mean(-NbC:NC));   PR_ut_mean=0.
    allocate (PR_vt_mean(-NbC:NC));   PR_vt_mean=0.
    allocate (PR_wt_mean(-NbC:NC));   PR_wt_mean=0.
    allocate (T % fluc(-NbC:NC));  T % fluc=0.
    end if

    allocate (Dift_uu_mean(1:NC));   Dift_uu_mean=0.
    allocate (Dift_vv_mean(1:NC));   Dift_vv_mean=0.
    allocate (Dift_ww_mean(1:NC));   Dift_ww_mean=0.
    allocate (Dift_uv_mean(1:NC));   Dift_uv_mean=0.
    allocate (Dift_uw_mean(1:NC));   Dift_uw_mean=0.
    allocate (Dift_vw_mean(1:NC));   Dift_vw_mean=0.

    allocate (Difv_uu_mean(1:NC));   Difv_uu_mean=0.
    allocate (Difv_vv_mean(1:NC));   Difv_vv_mean=0.
    allocate (Difv_ww_mean(1:NC));   Difv_ww_mean=0.
    allocate (Difv_uv_mean(1:NC));   Difv_uv_mean=0.
    allocate (Difv_uw_mean(1:NC));   Difv_uw_mean=0.
    allocate (Difv_vw_mean(1:NC));   Difv_vw_mean=0.

    allocate (C_uu_mean(1:NC));   C_uu_mean=0.
    allocate (C_vv_mean(1:NC));   C_vv_mean=0.
    allocate (C_ww_mean(1:NC));   C_ww_mean=0.
    allocate (C_uv_mean(1:NC));   C_uv_mean=0.
    allocate (C_uw_mean(1:NC));   C_uw_mean=0.
    allocate (C_vw_mean(1:NC));   C_vw_mean=0.

    allocate (PD_uu_mean(-NbC:NC));   PD_uu_mean=0.
    allocate (PD_vv_mean(-NbC:NC));   PD_vv_mean=0.
    allocate (PD_ww_mean(-NbC:NC));   PD_ww_mean=0.
    allocate (PD_uv_mean(-NbC:NC));   PD_uv_mean=0.
    allocate (PD_uw_mean(-NbC:NC));   PD_uw_mean=0.
    allocate (PD_vw_mean(-NbC:NC));   PD_vw_mean=0.

    allocate (PR_uu_mean(-NbC:NC));   PR_uu_mean=0.
    allocate (PR_vv_mean(-NbC:NC));   PR_vv_mean=0.
    allocate (PR_ww_mean(-NbC:NC));   PR_ww_mean=0.
    allocate (PR_uv_mean(-NbC:NC));   PR_uv_mean=0.
    allocate (PR_uw_mean(-NbC:NC));   PR_uw_mean=0.
    allocate (PR_vw_mean(-NbC:NC));   PR_vw_mean=0.
    end if

    allocate(VISt_mean(NC)); VISt_mean = 0.0
    allocate (ShearMean(NC));  ShearMean=0.
    if(HOT==YES) then
      allocate (T % mean(-NbC:NC));  T % mean=0.
      allocate (TT % mean(-NbC:NC)); TT % mean=0.
      allocate (uT % mean(-NbC:NC)); uT % mean=0.
      allocate (vT % mean(-NbC:NC)); vT % mean=0.
      allocate (wT % mean(-NbC:NC)); wT % mean=0.
    end if
  end if

  if(SIMULA == HYB_ZETA.or.SIMULA==HYB_PITM) then
    allocate (uu % mean(-NbC:NC)); uu % mean=0.
    allocate (vv % mean(-NbC:NC)); vv % mean=0.
    allocate (ww % mean(-NbC:NC)); ww % mean=0.
    allocate (uv % mean(-NbC:NC)); uv % mean=0.
    allocate (uw % mean(-NbC:NC)); uw % mean=0.
    allocate (vw % mean(-NbC:NC)); vw % mean=0.


    allocate(VISt_mean(NC)); VISt_mean = 0.0
    allocate (ShearMean(NC));  ShearMean=0.
    if(HOT==YES) then
      allocate (T % mean(-NbC:NC));  T % mean=0.
      allocate (TT % mean(-NbC:NC)); TT % mean=0.
      allocate (uT % mean(-NbC:NC)); uT % mean=0.
      allocate (vT % mean(-NbC:NC)); vT % mean=0.
      allocate (wT % mean(-NbC:NC)); wT % mean=0.
    end if
  end if
 
  allocate (VISt(-NbC:NC)); VISt=0
  allocate (IsNearWall(NC)); IsNearWall = .FALSE.

  allocate (Vort(-NbC:NC));  Vort=0.
  allocate (Shear(-NbC:NC)); Shear=0.
  allocate (TauWall(NC));    TauWall=0.

!??????????????????????????????????????????!
!     Is there enough allocated memory     !
!??????????????????????????????????????????!
! Do something !  

  END SUBROUTINE UnkAloc
