!======================================================================!
  SUBROUTINE IniPar 
!----------------------------------------------------------------------!
!   Initialize various solver parameters.                              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: i
!--------------------------------[CVS]---------------------------------!
!  $Id: IniPar.f90,v 1.23 2009/06/30 12:08:38 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/IniPar.f90,v $  
!======================================================================!

!-------------------------------!
!     Initialize parameters     !
!-------------------------------!
  AB         =  1
  CN         =  2
  FI         =  3
  LIN        =  4
  PAR        =  5
  SIMPLE     =  6
  FRACT      =  7
  AVS        =  8
  GMV        =  9
  YES        = 10
  NO         = 11
  LES        = 12
  DNS        = 13
  K_EPS      = 14
  K_EPS_VV   = 15
  SPA_ALL    = 16
  DES_SPA    = 17
  ZETA       = 18
  EBM        = 19
!---- convective schemes
  CDS        = 20
  QUICK      = 21 
  LUDS       = 22
  MINMOD     = 23
  SMART      = 24
  AVL_SMART  = 25
  SUPERBEE   = 26
  GAMMA      = 27
!---- for wall functions
  MODE       = 28
  LRe        = 29
  HRe        = 30
!---- LowRe k-eps models
  J_L        = 31  
  S_L_Y      = 32 
  NAG        = 33
  RNG        = 34
  SMAG       = 35
  WALE       = 36
  DYN        = 37
  MIX        = 38
  HYB_ZETA   = 39
  HYB_PITM   = 40
  BUDG       = 40
  HJ         = 41




  ALGOR  = FRACT   ! FRACT, SIMPLE
  INERT  = LIN     ! LIN, PAR 
  CONVEC = AB      ! AB, CN, FI
  CROSS  = AB      ! AB, CN, FI
  DIFFUS = CN      ! AB, CN, FI
  SIMULA = DNS     ! LES, DNS  
  POSPRO = GMV     ! AVS, GMV
  CHANNEL = NO     ! YES, NO
  PIPE    = NO     ! YES, NO
  JET     = NO     ! YES, NO
  TGV    = NO      ! YES, NO
  TEST   = NO      ! YES, NO
  OTHER  = NO      ! YES, NO
  HOT    = NO      ! YES, NO
  ROT    = NO      ! YES, NO
  HOTini = NO      ! YES, NO
  SHAKE  = NO      ! YES, NO
  BLEND  = NO      ! YES, NO, QUICK,, LUDS, MINMOD, SMART, AVL_SMART
  SHAKE_PER = -10
  SHAKE_INT = -10
  XHOM   = NO
  YHOM   = NO
  ZHOM   = NO

!----- Under relaxation factors
  U % URF      = 1.0
  P % URF      = 1.0
  T % URF      = 1.0
  Kin % URF    = 1.0
  Eps % URF    = 1.0
  v_2 % URF    = 1.0
  VIS % URF    = 1.0
  URFC         = 0.9
  
  Ndyn         = 0
 
!----- Algorythm and solver tolerances
  SIMTol       = 1.e-4
  U % STol     = 1.e-6 
  T % STol     = 1.e-6 
  PP % STol    = 1.e-6
  Kin % STol   = 1.e-6
  Eps % STol   = 1.e-6
  v_2 % STol   = 1.e-6
  f22 % STol   = 1.e-6
  VIS % STol   = 1.e-6

!----- Set the default monitoring point
  Nmon=1
  do i=1,Nmon
    Cm(i) = 1
  end do

  END SUBROUTINE IniPar
