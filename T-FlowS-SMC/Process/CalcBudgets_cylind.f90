!======================================================================!
  SUBROUTINE CalcBudgets_cylind(n0,n1)  
!----------------------------------------------------------------------!
!   Calculates time averaged velocity and velocity fluctuations.       !
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
  INTEGER :: n0, n1, n2
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c1, c2, s, c, n, new, i 
  REAL    :: VolDom, Lf, Lmix, Esgs, vol1
  REAL    :: Utan, Urad, Utan_mean, Urad_mean, R, Flux_mean, Us, Vs, Ws   
  REAL    :: Fim1, Fex1, Fim2, Fex2, Fim3, Fex3, Fim4, Fex4, Fim5, Fex5, Fim6, Fex6
  REAL    :: Fim7, Fex7, Fim8, Fex8, Fim9, Fex9, Fim10, Fex10 
  REAL    :: uu_s, vv_s, ww_s, uv_s, uw_s, vw_s, A0, R2, R1, PHIs1, PHIs2,PHIs3 
  REAL    :: B1, B2, CONc_s

  REAL    :: G1r, G1t, G2r, G2t, G3r, G3t, G4r, G4t, G5r, G5t 
  REAL    :: G6r, G6t, G7r, G7t, G8r, G8t, G9r, G9t, Gr, Gt 
  REAL    :: Gz, G1z, G2z, G3z, G4z, G5z, G6z, G7z, G8z, G9z
  REAL    :: G10r, G10t, G10z

  REAL,ALLOCATABLE :: Puu(:), Pvv(:), Pww(:), Puv(:), Puw(:), Pvw(:)
  REAL,ALLOCATABLE :: dUdr(:), dUdt(:), dUdz(:)
  REAL,ALLOCATABLE :: dVdr(:), dVdt(:), dVdz(:)
  REAL,ALLOCATABLE :: dWdr(:), dWdt(:), dWdz(:) 
  REAL,ALLOCATABLE :: Diss_uu(:), Diss_vv(:), Diss_ww(:), Diss_sgs(:) 
  REAL,ALLOCATABLE :: Diss_uv(:), Diss_uw(:), Diss_vw(:) 
  REAL,ALLOCATABLE :: C_uu(:), C_vv(:), C_ww(:) 
  REAL,ALLOCATABLE :: C_uv(:), C_uw(:), C_vw(:) 
  REAL,ALLOCATABLE :: Difv_uu(:), Difv_vv(:), Difv_ww(:) 
  REAL,ALLOCATABLE :: Difv_uv(:), Difv_uw(:), Difv_vw(:) 
  REAL,ALLOCATABLE :: Dift_uu(:), Dift_vv(:), Dift_ww(:) 
  REAL,ALLOCATABLE :: Dift_uv(:), Dift_uw(:), Dift_vw(:) 
  REAL,ALLOCATABLE :: R_uu(:), R_vv(:), R_ww(:) 
  REAL,ALLOCATABLE :: R_uv(:), R_uw(:), R_vw(:) 
  REAL,ALLOCATABLE :: PD_uu(:), PD_vv(:), PD_ww(:) 
  REAL,ALLOCATABLE :: PD_uv(:), PD_uw(:), PD_vw(:) 
  REAL,ALLOCATABLE :: uu_sgs(:), vv_sgs(:), ww_sgs(:), uv_sgs(:), uw_sgs(:), vw_sgs(:)

!--------------------------------[CVS]---------------------------------!
!  $Id: CalcMn.f90,v 1.15 2002/10/30 16:29:48 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Process/CalcMn.f90,v $   
!======================================================================!

  new = n1 - n0

  if(new < 0) return 

  allocate (Puu(-NbC:NC)); Puu = 0.0
  allocate (Pvv(-NbC:NC)); Pvv = 0.0
  allocate (Pww(-NbC:NC)); Pww = 0.0
  allocate (Puv(-NbC:NC)); Puv = 0.0
  allocate (Puw(-NbC:NC)); Puw = 0.0
  allocate (Pvw(-NbC:NC)); Pvw = 0.0
  allocate (C_uu(-NbC:NC)); C_uu = 0.0
  allocate (C_vv(-NbC:NC)); C_vv = 0.0
  allocate (C_ww(-NbC:NC)); C_ww = 0.0
  allocate (C_uv(-NbC:NC)); C_uv = 0.0
  allocate (C_uw(-NbC:NC)); C_uw = 0.0
  allocate (C_vw(-NbC:NC)); C_vw = 0.0
  allocate (R_uu(-NbC:NC)); R_uu = 0.0
  allocate (R_vv(-NbC:NC)); R_vv = 0.0
  allocate (R_ww(-NbC:NC)); R_ww = 0.0
  allocate (R_uv(-NbC:NC)); R_uv = 0.0
  allocate (R_uw(-NbC:NC)); R_uw = 0.0
  allocate (R_vw(-NbC:NC)); R_vw = 0.0
  allocate (dUdr(-NbC:NC)); dUdr = 0.0
  allocate (dUdt(-NbC:NC)); dUdt = 0.0
  allocate (dUdz(-NbC:NC)); dUdz = 0.0
  allocate (dVdr(-NbC:NC)); dVdr = 0.0
  allocate (dVdt(-NbC:NC)); dVdt = 0.0
  allocate (dVdz(-NbC:NC)); dVdz = 0.0
  allocate (dWdr(-NbC:NC)); dWdr = 0.0
  allocate (dWdt(-NbC:NC)); dWdt = 0.0
  allocate (dWdz(-NbC:NC)); dWdz = 0.0
  allocate (Diss_uu(-NbC:NC)); Diss_uu = 0.0
  allocate (Diss_vv(-NbC:NC)); Diss_vv = 0.0
  allocate (Diss_ww(-NbC:NC)); Diss_ww = 0.0
  allocate (Diss_uv(-NbC:NC)); Diss_uv = 0.0
  allocate (Diss_uw(-NbC:NC)); Diss_uw = 0.0
  allocate (Diss_vw(-NbC:NC)); Diss_vw = 0.0
  allocate (Diss_sgs(-NbC:NC)); Diss_sgs = 0.0
  allocate (Difv_uu(-NbC:NC)); Difv_uu = 0.0
  allocate (Difv_vv(-NbC:NC)); Difv_vv = 0.0
  allocate (Difv_ww(-NbC:NC)); Difv_ww = 0.0
  allocate (Difv_uv(-NbC:NC)); Difv_uv = 0.0
  allocate (Difv_uw(-NbC:NC)); Difv_uw = 0.0
  allocate (Difv_vw(-NbC:NC)); Difv_vw = 0.0
  allocate (Dift_uu(-NbC:NC)); Dift_uu = 0.0
  allocate (Dift_vv(-NbC:NC)); Dift_vv = 0.0
  allocate (Dift_ww(-NbC:NC)); Dift_ww = 0.0
  allocate (Dift_uv(-NbC:NC)); Dift_uv = 0.0
  allocate (Dift_uw(-NbC:NC)); Dift_uw = 0.0
  allocate (Dift_vw(-NbC:NC)); Dift_vw = 0.0
  allocate (PD_uu(-NbC:NC)); PD_uu = 0.0
  allocate (PD_vv(-NbC:NC)); PD_vv = 0.0
  allocate (PD_ww(-NbC:NC)); PD_ww = 0.0
  allocate (PD_uv(-NbC:NC)); PD_uv = 0.0
  allocate (PD_uw(-NbC:NC)); PD_uw = 0.0
  allocate (PD_vw(-NbC:NC)); PD_vw = 0.0
  allocate (uu_sgs(-NbC:NC)); uu_sgs = 0.0
  allocate (vv_sgs(-NbC:NC)); vv_sgs = 0.0
  allocate (ww_sgs(-NbC:NC)); ww_sgs = 0.0
  allocate (uv_sgs(-NbC:NC)); uv_sgs = 0.0
  allocate (uw_sgs(-NbC:NC)); uw_sgs = 0.0
  allocate (vw_sgs(-NbC:NC)); vw_sgs = 0.0

  do c = 1, NC
    R     = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 
    R1      = 1.0/R 
    R2      = 1.0/(R*R) 

    Ux(c) = U % mean(c) * xc(c) * R1  + V % mean(c) * yc(c) * R1  
    Vx(c) = -U % mean(c) * yc(c) * R1  + V % mean(c) * xc(c) * R1 
    uu % n(c) = uu % mean(c) - Ux(c) * Ux(c)
    vv % n(c) = vv % mean(c) - Vx(c) * Vx(c)
    ww % n(c) = ww % mean(c) - W % mean(c) * W % mean(c)

    uv % n(c) = uv % mean(c) - Ux(c) * Vx(c)    
    uw % n(c) = uw % mean(c) - Ux(c) * W % mean(c)    
    vw % n(c) = vw % mean(c) - Vx(c) * W % mean(c)    

    U % fluc(c) = U % n(c) * xc(c) * R1  + V % n(c) * yc(c) * R1 - Ux(c)
    V % fluc(c) = -U % n(c) * yc(c) * R1  + V % n(c) * xc(c) * R1 - Vx(c)
    W % fluc(c) = W % n(c) - W % mean(c)
    P % fluc(c) = P % n(c) - P % mean(c)
  end do

if(mod(n1,1) == 0) then   !IZMJENA
!----------------------------------------!
!  PRODUCTION OF STRESSES COMPONENTS
!----------------------------------------!
    call GraPhi(Ux,1,PHIx, .TRUE.)  
    call GraPhi(Ux,2,PHIy, .TRUE.)
    call GraPhi(Ux,3,PHIz, .TRUE.)
    call GraPhi(Vx,1,PHI1x, .TRUE.)  
    call GraPhi(Vx,2,PHI1y, .TRUE.)  
    call GraPhi(Vx,3,PHI1z, .TRUE.)  
    call GraPhi(W % mean,1,PHI2x, .TRUE.) 
    call GraPhi(W % mean,2,PHI2y, .TRUE.) 
    call GraPhi(W % mean,3,PHI2z, .TRUE.) 

    do c = 1, NC
      R       = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 + tiny 
      R1      = 1.0/R 
      R2      = 1.0/(R*R) 

      dUdr(c) = PHIx(c) * xc(c) * R1  + PHIy(c) * yc(c) * R1 
      dUdt(c) = -PHIx(c) * yc(c) * R1  + PHIy(c) * xc(c) * R1
      dUdz(c) = PHIz(c)

      dVdr(c) = PHI1x(c) * xc(c) * R1  + PHI1y(c) * yc(c) * R1
      dVdt(c) = -PHI1x(c) * yc(c) * R1  + PHI1y(c) * xc(c) * R1
      dVdz(c) = PHI1z(c)

      dWdr(c) = PHI2x(c) * xc(c) * R1  + PHI2y(c) * yc(c) * R1
      dWdt(c) = -PHI2x(c) * yc(c) * R1  + PHI2y(c) * xc(c) * R1
      dWdz(c) = PHI2z(c)

      G1r = PHI3x(c) * xc(c) * R1  + PHI3y(c) * yc(c) * R1
      G1t = -PHI3x(c) * yc(c) * R1  + PHI3y(c) * xc(c) * R1
      G1z = PHI3z(c)


      Puu(c) = (2.0*uu % n(c)*dUdr(c) + 2.0*uv % n(c)*R1*dUdt(c) + 2.0*uw % n(c)*dUdz(c) & 
             - 2.0*Vx(c)*R1*uv % n(c))     
      Pvv(c) = (2.0*uv % n(c) * dVdr(c) + 2.0*vv % n(c) * R1 * dVdt(c) + 2.0*vw % n(c) * dVdz(c)  &
             + 2.0*Ux(c) * R1 * vv % n(c))        

      Pww(c) = (2.0*uw % n(c) * dWdr(c) + 2.0*vw % n(c) * R1 * dWdt(c) + 2.0*ww % n(c) * dWdz(c)) 

      Puv(c) = (uu % n(c) * dVdr(c) + uv%n(c)*R1 * dVdt(c) + uw%n(c)*dVdz(c) &
             + uv % n(c) * dUdr(c) + vv%n(c)*R1 * dUdt(c) + vw%n(c)*dUdz(c) - vv % n(c)*Vx(c)*R1 + uv%n(c)*Ux(c)*R1)

      Puw(c) = (uw % n(c) * dUdr(c) + vw % n(c) * R1 * dUdt(c) + ww % n(c) * dUdz(c) &
             + uu % n(c) * dWdr(c) + uv % n(c) * R1 * dWdt(c) + uw % n(c) * dWdz(c)   &
             - Vx(c) * R1 * vw % n(c))  
      Pvw(c) = (uw % n(c) * dVdr(c) + vw % n(c) * R1 * dVdt(c) + ww % n(c) * dVdz(c)  &
             + uv % n(c) * dWdr(c) + vv % n(c) * R1 * dWdt(c) + vw % n(c) * dWdz(c)   &
             + Ux(c) * R1 * vw % n(c))  

!      C_uu(c)=C_uu(c) + 2.0*R1*Vx(c)*uv % n(c)*volume(c) 
!      C_vv(c)=C_vv(c) - 2.0*R1*Vx(c)*uv % n(c)*volume(c)
!      C_uv(c)=C_uv(c) + R1*Vx(c)*vv % n(c)*volume(c) + R1*Vx(c)*uu % n(c)*volume(c)
!      C_uw(c)=C_uw(c) + R1*Vx(c)*vw % n(c)*volume(c)
!      C_vw(c)=C_vw(c) - R1*Vx(c)*uw % n(c)*volume(c)


      Puu_mean(c) = Puu(c) 
      Pvv_mean(c) = Pvv(c)
      Pww_mean(c) = Pww(c) 
      Puv_mean(c) = Puv(c) 
      Puw_mean(c) = Puw(c) 
      Pvw_mean(c) = Pvw(c) 

      Ksgs(c)  = (uv%n(c)*(dVdr(c)-V % mean(c)*R1) + uw%n(c)*dWdr(c))/(Shear(c)*Shear(c))
      uu_sgs(c) = -2.0 * VISt(c) * dUdr(c) + 2.0/3.0 * ksgs(c) 
      vv_sgs(c) = -2.0 * VISt(c) * dVdt(c) * R1 + 2.0/3.0 * ksgs(c)
      ww_sgs(c) = -2.0 * VISt(c) * dWdz(c) + 2.0/3.0 * ksgs(c)
      uv_sgs(c) = -VISt(c) * (dUdt(c) * R1 + dVdr(c))
      uw_sgs(c) = -VISt(c) * (dUdz(c) + dWdr(c)) 
      vw_sgs(c) = -VISt(c) * (dWdt(c) * R1 + dVdz(c))

    end do    ! c 
!----------------------------------------!
!  CONVECTION & VISCOUS DIFFUSION
!----------------------------------------!

    call GraPhi(uu % n,1,PHIx, .TRUE.) 
    call GraPhi(uu % n,2,PHIy, .TRUE.) 
    call GraPhi(uu % n,3,PHIz, .TRUE.) 

    call GraPhi(vv % n,1,PHI1x, .TRUE.) 
    call GraPhi(vv % n,2,PHI1y, .TRUE.) 
    call GraPhi(vv % n,3,PHI1z, .TRUE.) 

    call GraPhi(ww % n,1,PHI2x, .TRUE.) 
    call GraPhi(ww % n,2,PHI2y, .TRUE.) 
    call GraPhi(ww % n,3,PHI2z, .TRUE.) 

    call GraPhi(uv % n,1,PHI3x, .TRUE.) 
    call GraPhi(uv % n,2,PHI3y, .TRUE.) 
    call GraPhi(uv % n,3,PHI3z, .TRUE.) 

    call GraPhi(uw % n,1,PHI4x, .TRUE.) 
    call GraPhi(uw % n,2,PHI4y, .TRUE.) 
    call GraPhi(uw % n,3,PHI4z, .TRUE.) 

    call GraPhi(vw % n,1,PHI5x, .TRUE.) 
    call GraPhi(vw % n,2,PHI5y, .TRUE.) 
    call GraPhi(vw % n,3,PHI5z, .TRUE.) 

    do s=1,NS

      c1=SideC(1,s)
      c2=SideC(2,s)

      uu_s =f(s)*uu % n(c1) + (1.0-f(s))*uu % n(c2)
      vv_s =f(s)*vv % n(c1) + (1.0-f(s))*vv % n(c2)
      ww_s =f(s)*ww % n(c1) + (1.0-f(s))*ww % n(c2)
      uv_s =f(s)*uv % n(c1) + (1.0-f(s))*uv % n(c2)
      uw_s =f(s)*uw % n(c1) + (1.0-f(s))*uw % n(c2)
      vw_s =f(s)*vw % n(c1) + (1.0-f(s))*vw % n(c2)
      Us   =f(s)*U % mean(c1) + (1.0-f(s))*U % mean(c2)
      Vs   =f(s)*V % mean(c1) + (1.0-f(s))*V % mean(c2)
      Ws   =f(s)*W % mean(c1) + (1.0-f(s))*W % mean(c2)
      Flux_mean = Us*Sx(s) + Vs*Sy(s) + Ws*Sz(s)

      A0 = VISc * Scoef(s)

      G1r = fF(s)*PHIx(c1) + (1.0-fF(s))*PHIx(c2) 
      G1t = fF(s)*PHIy(c1) + (1.0-fF(s))*PHIy(c2) 
      G1z = fF(s)*PHIz(c1) + (1.0-fF(s))*PHIz(c2) 

      G2r = fF(s)*PHI1x(c1) + (1.0-fF(s))*PHI1x(c2) 
      G2t = fF(s)*PHI1y(c1) + (1.0-fF(s))*PHI1y(c2) 
      G2z = fF(s)*PHI1z(c1) + (1.0-fF(s))*PHI1z(c2) 

      G3r = fF(s)*PHI2x(c1) + (1.0-fF(s))*PHI2x(c2) 
      G3t = fF(s)*PHI2y(c1) + (1.0-fF(s))*PHI2y(c2) 
      G3z = fF(s)*PHI2z(c1) + (1.0-fF(s))*PHI2z(c2) 

      G4r = fF(s)*PHI3x(c1) + (1.0-fF(s))*PHI3x(c2) 
      G4t = fF(s)*PHI3y(c1) + (1.0-fF(s))*PHI3y(c2) 
      G4z = fF(s)*PHI3z(c1) + (1.0-fF(s))*PHI3z(c2) 

      G5r = fF(s)*PHI4x(c1) + (1.0-fF(s))*PHI4x(c2) 
      G5t = fF(s)*PHI4y(c1) + (1.0-fF(s))*PHI4y(c2) 
      G5z = fF(s)*PHI4z(c1) + (1.0-fF(s))*PHI4z(c2) 

      G6r = fF(s)*PHI5x(c1) + (1.0-fF(s))*PHI5x(c2) 
      G6t = fF(s)*PHI5y(c1) + (1.0-fF(s))*PHI5y(c2) 
      G6z = fF(s)*PHI5z(c1) + (1.0-fF(s))*PHI5z(c2) 

      G7r = fF(s)*PHI6x(c1) + (1.0-fF(s))*PHI6x(c2) 
      G7t = fF(s)*PHI6y(c1) + (1.0-fF(s))*PHI6y(c2) 
      G7z = fF(s)*PHI6z(c1) + (1.0-fF(s))*PHI6z(c2) 

      G8r = fF(s)*PHI7x(c1) + (1.0-fF(s))*PHI7x(c2) 
      G8t = fF(s)*PHI7y(c1) + (1.0-fF(s))*PHI7y(c2) 
      G8z = fF(s)*PHI7z(c1) + (1.0-fF(s))*PHI7z(c2) 

      G9r = fF(s)*PHI8x(c1) + (1.0-fF(s))*PHI8x(c2) 
      G9t = fF(s)*PHI8y(c1) + (1.0-fF(s))*PHI8y(c2) 
      G9z = fF(s)*PHI8z(c1) + (1.0-fF(s))*PHI8z(c2) 

      G10r = fF(s)*PHI9x(c1) + (1.0-fF(s))*PHI9x(c2) 
      G10t = fF(s)*PHI9y(c1) + (1.0-fF(s))*PHI9y(c2) 
      G10z = fF(s)*PHI9z(c1) + (1.0-fF(s))*PHI9z(c2) 

      Fim1 =  VISc*( G1r*Sx(s) + G1t*Sy(s) + G1z*Sz(s) )
      Fex1 = ( G1r*Dx(s) + G1t*Dy(s) + G1z*Dz(s))*A0 
 
      Fim2 =  VISc*( G2r*Sx(s) + G2t*Sy(s) + G2z*Sz(s) )
      Fex2 = ( G2r*Dx(s) + G2t*Dy(s) + G2z*Dz(s))*A0 

      Fim3 =  VISc*( G3r*Sx(s) + G3t*Sy(s) + G3z*Sz(s) )
      Fex3 = ( G3r*Dx(s) + G3t*Dy(s) + G3z*Dz(s))*A0 
 
      Fim4 =  VISc*( G4r*Sx(s) + G4t*Sy(s) + G4z*Sz(s) )
      Fex4 = ( G4r*Dx(s) + G4t*Dy(s) + G4z*Dz(s))*A0 
 
      Fim5 =  VISc*( G5r*Sx(s) + G5t*Sy(s) + G5z*Sz(s) )
      Fex5 = ( G5r*Dx(s) + G5t*Dy(s) + G5z*Dz(s))*A0 
 
      Fim6 =  VISc*( G6r*Sx(s) + G6t*Sy(s) + G6z*Sz(s) )
      Fex6 = ( G6r*Dx(s) + G6t*Dy(s) + G6z*Dz(s))*A0 
 
      if(c2  > 0) then
        C_uu(c1)=C_uu(c1)-Flux_mean*uu_s
        C_uu(c2)=C_uu(c2)+Flux_mean*uu_s

        C_vv(c1)=C_vv(c1)-Flux_mean*vv_s
        C_vv(c2)=C_vv(c2)+Flux_mean*vv_s

        C_ww(c1)=C_ww(c1)-Flux_mean*ww_s
        C_ww(c2)=C_ww(c2)+Flux_mean*ww_s

        C_uv(c1)=C_uv(c1)-Flux_mean*uv_s
        C_uv(c2)=C_uv(c2)+Flux_mean*uv_s

        C_uw(c1)=C_uw(c1)-Flux_mean*uw_s
        C_uw(c2)=C_uw(c2)+Flux_mean*uw_s

        C_vw(c1)=C_vw(c1)-Flux_mean*vw_s
        C_vw(c2)=C_vw(c2)+Flux_mean*vw_s

        Difv_uu(c1)=Difv_uu(c1) + (uu % n(c2)-uu % n(c1))*A0  + Fex1 - Fim1
        Difv_uu(c2)=Difv_uu(c2) - (uu % n(c2)-uu % n(c1))*A0  + Fex1 - Fim1

        Difv_vv(c1)=Difv_vv(c1) + (vv % n(c2)-vv % n(c1))*A0  + Fex2 - Fim2
        Difv_vv(c2)=Difv_vv(c2) - (vv % n(c2)-vv % n(c1))*A0  + Fex2 - Fim2

        Difv_ww(c1)=Difv_ww(c1) + (ww % n(c2)-ww % n(c1))*A0  + Fex3 - Fim3
        Difv_ww(c2)=Difv_ww(c2) - (ww % n(c2)-ww % n(c1))*A0  + Fex3 - Fim3

        Difv_uv(c1)=Difv_uv(c1) + (uv % n(c2)-uv % n(c1))*A0  + Fex4 - Fim4
        Difv_uv(c2)=Difv_uv(c2) - (uv % n(c2)-uv % n(c1))*A0  + Fex4 - Fim4
  
        Difv_uw(c1)=Difv_uw(c1) + (uw % n(c2)-uw % n(c1))*A0  + Fex5 - Fim5
        Difv_uw(c2)=Difv_uw(c2) - (uw % n(c2)-uw % n(c1))*A0  + Fex5 - Fim5

        Difv_vw(c1)=Difv_vw(c1) + (vw % n(c2)-vw % n(c1))*A0  + Fex6 - Fim6
        Difv_vw(c2)=Difv_vw(c2) - (vw % n(c2)-vw % n(c1))*A0  + Fex6 - Fim6
  
      else
        C_uu(c1)=C_uu(c1)-Flux_mean*uu_s
        C_vv(c1)=C_vv(c1)-Flux_mean*vv_s
        C_ww(c1)=C_ww(c1)-Flux_mean*ww_s
        C_uv(c1)=C_uv(c1)-Flux_mean*uv_s
        C_uw(c1)=C_uw(c1)-Flux_mean*uw_s
        C_vw(c1)=C_vw(c1)-Flux_mean*vw_s

      endif
    end do 

!-------------------------------------------------!
!  EXTRA TERMS IN VISCOUS DIFFUSION 
!-------------------------------------------------!

    do c = 1, NC
      R       = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 + tiny
      R1      = 1.0/R 
      R2      = 1.0/(R*R) 
      dUdr(c) = -PHI3x(c)  * yc(c) * R1  + PHI3y(c) * xc(c) * R1
      dUdt(c) = -PHI1x(c) * yc(c) * R1  + PHI1y(c) * xc(c) * R1 
      dVdr(c) = -PHIx(c) * yc(c) * R1  + PHIy(c) * xc(c) * R1
      dVdt(c) = -PHI5x(c) * yc(c) * R1  + PHI5y(c) * xc(c) * R1
      dWdt(c) = -PHI4x(c) * yc(c) * R1  + PHI4y(c) * xc(c) * R1

      Difv_uu(c)=Difv_uu(c) + VISc*volume(c)*(-2.0*R2*uu % n(c) + 2.0*R2*vv % n(c) - 4.0*R2*dUdr(c)) 
      Difv_vv(c)=Difv_vv(c) + VISc*volume(c)*(2.0*R2* uu % n(c) - 2.0*R2*vv % n(c) + 4.0*R2*dUdr(c)) 
      Difv_uv(c)=Difv_uv(c) + VISc*volume(c)*(-4.0*R2*uv % n(c) -2.0*R2*dUdt(c)+ 2.0*R2*dVdr(c))
      Difv_uw(c)=Difv_uw(c) + VISc*volume(c)*(-2.0*R2*dVdt(c) - R2*uw % n(c))
      Difv_vw(c)=Difv_vw(c) + VISc*volume(c)*(2.0*R2*dWdt(c) - R2*vw % n(c))


      C_uu_mean(c) = C_uu(c) - 2.0*R1*Vx(c)*uv % n(c)*volume(c) 
      C_vv_mean(c) = C_vv(c) + 2.0*R1*Vx(c)*uv % n(c)*volume(c)
      C_ww_mean(c) = C_ww(c) 
      C_uv_mean(c) = C_uv(c) + (R1*Vx(c)*uu % n(c) - R1*Vx(c)*vv %n(c))*volume(c)
      C_uw_mean(c) = C_uw(c) - R1*Vx(c)*vw % n(c)*volume(c)
      C_vw_mean(c) = C_vw(c) + R1*Vx(c)*uw % n(c)*volume(c)
      
      Difv_uu_mean(c) = Difv_uu(c)
      Difv_vv_mean(c) = Difv_vv(c) 
      Difv_ww_mean(c) = Difv_ww(c) 
      Difv_uv_mean(c) = Difv_uv(c) 
      Difv_uw_mean(c) = Difv_uw(c) 
      Difv_vw_mean(c) = Difv_vw(c) 

    end do
  end if   !end mod(1,0)
!----------------------------------------!
!  DISSIPATION & PRESSURE REDISTRIBUTION
!----------------------------------------!
  call GraPhi(U % fluc,1,PHIx, .TRUE.) 
  call GraPhi(U % fluc,2,PHIy, .TRUE.) 
  call GraPhi(U % fluc,3,PHIz, .TRUE.) 

  call GraPhi(V % fluc,1,PHI1x, .TRUE.)
  call GraPhi(V % fluc,2,PHI1y, .TRUE.) 
  call GraPhi(V % fluc,3,PHI1z, .TRUE.) 

  call GraPhi(W % fluc,1,PHI2x, .TRUE.) 
  call GraPhi(W % fluc,2,PHI2y, .TRUE.) 
  call GraPhi(W % fluc,3,PHI2z, .TRUE.) 

  do c = 1, NC
    R       = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 + tiny
    R1      = 1.0/R 
    R2      = 1.0/(R*R) 

    dUdr(c) =  PHIx(c) * xc(c) * R1   + PHIy(c) * yc(c) * R1
    dUdt(c) = -PHIx(c) * yc(c) * R1   + PHIy(c) * xc(c) * R1
    dUdz(c) =  PHIz(c)

    dVdr(c) =  PHI1x(c) * xc(c) * R1  + PHI1y(c) * yc(c) * R1
    dVdt(c) = -PHI1x(c) * yc(c) * R1  + PHI1y(c) * xc(c) * R1
    dVdz(c) =  PHI1z(c)       

    dWdr(c) =  PHI2x(c) * xc(c) * R1  + PHI2y(c) * yc(c) * R1
    dWdt(c) = -PHI2x(c) * yc(c) * R1  + PHI2y(c) * xc(c) * R1
    dWdz(c) =  PHI2z(c)       

    G1r = PHI3x(c) * xc(c) * R1  + PHI3y(c) * yc(c) * R1
    G1t = -PHI3x(c) * yc(c) * R1  + PHI3y(c) * xc(c) * R1
    G1z = PHI3z(c)

    PHI4x(c) = dUdr(c)   
    PHI4y(c) = dUdt(c)   
    PHI4z(c) = dUdz(c)   
    PHI5x(c) = dVdr(c)   
    PHI5y(c) = dVdt(c)   
    PHI5z(c) = dVdz(c)   
    PHI6x(c) = dWdr(c)   
    PHI6y(c) = dWdt(c)   
    PHI6z(c) = dWdz(c)   
    PHI7x(c) = G1r 
    PHI7y(c) = G1t
    PHI7z(c) = G1z

    Diss_uu(c) = Diss_uu(c) - 2*VISc*(dUdr(c)*dUdr(c) + dUdt(c)*dUdt(c) * R2  &
               + dUdz(c)*dUdz(c) + vv % n(c) * R2 - 2.0 * R2 * V % fluc(c)*dUdt(c) )

    Diss_vv(c) = Diss_vv(c) - 2*VISc*(dVdr(c)*dVdr(c) + dVdt(c)*dVdt(c) * R2  &
               + dVdz(c)*dVdz(c) + uu % n(c) * R2 + 2.0 * R2 * U % fluc(c)*dVdt(c) ) 
    
    Diss_ww(c) = Diss_ww(c) - 2*VISc*(dWdr(c)*dWdr(c) + dWdt(c)*dWdt(c)*R2 + dWdz(c)*dWdz(c))

    Diss_uv(c) = Diss_uv(c) - 2.0*VISc*( dUdr(c)*dVdr(c) + dUdt(c)*dVdt(c)*R2  &
               + dUdz(c)*dVdz(c) - uv % n(c)*R2 - V % fluc(c)*dVdt(c)*R2 + U % fluc(c)*dUdt(c)*R2)

    Diss_uw(c) = DiSs_uw(c) - 2.0*VISc*( dUdr(c)*dWdr(c) + dUdt(c)*dWdt(c)*R2  &
               + dUdz(c)*dWdz(c) - V % fluc(c)*dWdt(c)*R2 )

    Diss_vw(c) = Diss_vw(c) - 2.0*VISc*( dVdr(c)*dWdr(c) + dVdt(c)*dWdt(c)*R2  &
               + dVdz(c)*dWdz(c) + U % fluc(c)*dWdt(c)*R2 )

    R_uu(c)=R_uu(c) + 2.0*P % fluc(c)*dUdr(c) 
    R_vv(c)=R_vv(c) + 2.0*P % fluc(c)*dVdt(c)*R1 + 2.0*R1*P % fluc(c)*U % fluc(c) 
    R_ww(c)=R_ww(c) + 2.0*P % fluc(c)*dWdz(c) 
    R_uv(c)=R_uv(c) + P % fluc(c)*( R1*dUdt(c) + dVdr(c) - R1*V % fluc(c))
    R_uw(c)=R_uw(c) + P % fluc(c)*( dUdz(c) + dWdr(c) )
    R_vw(c)=R_vw(c) + P % fluc(c)*( dVdz(c) + R1*dWdt(c) ) 
    Diss_sgs(c) = Diss_sgs(c) +  (uu_sgs(c) * dUdr(c) + uv_sgs(c) * dUdt(c) * R1 + uw_sgs(c) &
                  * dUdz(c) - uv_sgs(c) * V % fluc(c) * R1 + &
                   uv_sgs(c) * dVdr(c) + vv_sgs(c) * dVdt(c) * R1 + vw_sgs(c) * dVdz(c) &
                    + vv_sgs(c) * U % fluc(c) * R1 + &
                   uw_sgs(c) * dWdr(c) + vw_sgs(c) * dWdt(c) * R1 + ww_sgs(c) * dWdz(c))

!!!
!!! EXTRA TERMS FOR VISCOUS DIFFUSION OF TURBULENT HEAT FLUX
!!! IZMJENA
  end do 

  call GraPhi(PHI4x,1,PHIx, .TRUE.) 
  call GraPhi(PHI4x,2,PHIy, .TRUE.) 

  call GraPhi(PHI4y,1,PHI1x, .TRUE.) 
  call GraPhi(PHI4y,2,PHI1y, .TRUE.) 
  call GraPhi(PHI4z,3,PHI1z, .TRUE.) 

  call GraPhi(PHI5x,1,PHI2x, .TRUE.) 
  call GraPhi(PHI5x,2,PHI2y, .TRUE.) 

  call GraPhi(PHI5y,1,PHI3x, .TRUE.) 
  call GraPhi(PHI5y,2,PHI3y, .TRUE.) 
  call GraPhi(PHI5z,3,PHI3z, .TRUE.) 

  call GraPhi(PHI6x,1,PHI8x, .TRUE.) 
  call GraPhi(PHI6x,2,PHI8y, .TRUE.) 

  call GraPhi(PHI6y,1,PHI9x, .TRUE.) 
  call GraPhi(PHI6y,2,PHI9y, .TRUE.) 
  call GraPhi(PHI6z,3,PHI9z, .TRUE.) 

  call GraPhi(PHI7x,1,PHI5x, .TRUE.) 
  call GraPhi(PHI7x,2,PHI5y, .TRUE.) 

  call GraPhi(PHI7y,1,PHI4x, .TRUE.) 
  call GraPhi(PHI7y,2,PHI4y, .TRUE.) 
  call GraPhi(PHI7z,3,PHI4z, .TRUE.) 


  do c = 1, NC
    R       = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 + tiny
    R1      = 1.0/R
    R2      = 1.0/(R*R)
    Gr =  PHIx(c) * xc(c) * R1   + PHIy(c) * yc(c) * R1
    Gt = -PHI1x(c) * yc(c) * R1   + PHI1y(c) * xc(c) * R1
    Gz =  PHI1z(c)

    G1r =  PHI2x(c) * xc(c) * R1   + PHI2y(c) * yc(c) * R1
    G1t = -PHI3x(c) * yc(c) * R1   + PHI3y(c) * xc(c) * R1
    G1z =  PHI3z(c)

    G2r =  PHI8x(c) * xc(c) * R1   + PHI8y(c) * yc(c) * R1
    G2t = -PHI9x(c) * yc(c) * R1   + PHI9y(c) * xc(c) * R1
    G2z =  PHI9z(c)

    G3r =  PHI5x(c) * xc(c) * R1   + PHI5y(c) * yc(c) * R1
    G3t = -PHI4x(c) * yc(c) * R1   + PHI4y(c) * xc(c) * R1
    G3z =  PHI4z(c)

  end do

!----------------------------------------!
!  TURBULENT DIFFUSION
!----------------------------------------!

!  do c = 1, NC
!    Ux(c)     = U % fluc(c) * U % fluc(c) * U % fluc(c)
!    Uy(c)     = U % fluc(c) * U % fluc(c) * V % fluc(c)
!    Uz(c)     = U % fluc(c) * U % fluc(c) * W % fluc(c)
!
!    Vx(c)     = U % fluc(c) * V % fluc(c) * V % fluc(c)
!    Vy(c)     = V % fluc(c) * V % fluc(c) * V % fluc(c)
!    Vz(c)     = V % fluc(c) * V % fluc(c) * W % fluc(c)
!
!    Wx(c)     = U % fluc(c) * W % fluc(c) * W % fluc(c)
!    Wy(c)     = V % fluc(c) * W % fluc(c) * W % fluc(c)
!    Wz(c)     = W % fluc(c) * W % fluc(c) * W % fluc(c)
!   
!    Kx(c)     = U % fluc(c) * V % fluc(c) * W % fluc(c)
!  end do
!  The commented part was original implementation. The both previous one and
!  the new one are correct and tested. The new implementation use uuu which is
!   averaged in CalcMn together with U,V,W. 

  call GraPhi(uuu % mean,1,PHIx, .TRUE.)  
  call GraPhi(uuu % mean,2,PHIy, .TRUE.)  
  call GraPhi(uuu % mean,3,PHIz, .TRUE.)  

  call GraPhi(uuv % mean,1,PHI1x, .TRUE.)  
  call GraPhi(uuv % mean,2,PHI1y, .TRUE.)  
  call GraPhi(uuv % mean,3,PHI1z, .TRUE.)  

  call GraPhi(uuw % mean,1,PHI2x, .TRUE.)  
  call GraPhi(uuw % mean,2,PHI2y, .TRUE.)  
  call GraPhi(uuw % mean,3,PHI2z, .TRUE.)  
 
  call GraPhi(vvu % mean,1,PHI3x, .TRUE.)  
  call GraPhi(vvu % mean,2,PHI3y, .TRUE.)  
  call GraPhi(vvu % mean,3,PHI3z, .TRUE.)  

  call GraPhi(vvv % mean,1,PHI4x, .TRUE.)  
  call GraPhi(vvv % mean,2,PHI4y, .TRUE.)  
  call GraPhi(vvv % mean,3,PHI4z, .TRUE.)  

  call GraPhi(vvw % mean,1,PHI5x, .TRUE.)  
  call GraPhi(vvw % mean,2,PHI5y, .TRUE.)  
  call GraPhi(vvw % mean,3,PHI5z, .TRUE.)  

  call GraPhi(wwu % mean,1,PHI6x, .TRUE.)  
  call GraPhi(wwu % mean,2,PHI6y, .TRUE.)  
  call GraPhi(wwu % mean,3,PHI6z, .TRUE.)  

  call GraPhi(wwv % mean,1,PHI7x, .TRUE.)  
  call GraPhi(wwv % mean,2,PHI7y, .TRUE.)  
  call GraPhi(wwv % mean,3,PHI7z, .TRUE.)  

  call GraPhi(www % mean,1,PHI8x, .TRUE.)  
  call GraPhi(www % mean,2,PHI8y, .TRUE.)  
  call GraPhi(www % mean,3,PHI8z, .TRUE.)  

  call GraPhi(uwv % mean,1,PHI9x, .TRUE.)  
  call GraPhi(uwv % mean,2,PHI9y, .TRUE.)  
  call GraPhi(uwv % mean,3,PHI9z, .TRUE.)  

!  call GraPhi(Ux,1,PHIx, .TRUE.)  
!  call GraPhi(Ux,2,PHIy, .TRUE.)  
!  call GraPhi(Ux,3,PHIz, .TRUE.)  

!  call GraPhi(Uy,1,PHI1x, .TRUE.)  
!  call GraPhi(Uy,2,PHI1y, .TRUE.)  
!  call GraPhi(Uy,3,PHI1z, .TRUE.)  

!  call GraPhi(Uz,1,PHI2x, .TRUE.)  
!  call GraPhi(Uz,2,PHI2y, .TRUE.)  
!  call GraPhi(Uz,3,PHI2z, .TRUE.)  
 
!  call GraPhi(Vx,1,PHI3x, .TRUE.)  
!  call GraPhi(Vx,2,PHI3y, .TRUE.)  
!  call GraPhi(Vx,3,PHI3z, .TRUE.)  

!  call GraPhi(Vy,1,PHI4x, .TRUE.)  
!  call GraPhi(Vy,2,PHI4y, .TRUE.)  
!  call GraPhi(Vy,3,PHI4z, .TRUE.)  

!  call GraPhi(Vz,1,PHI5x, .TRUE.)  
!  call GraPhi(Vz,2,PHI5y, .TRUE.)  
!  call GraPhi(Vz,3,PHI5z, .TRUE.)  

!  call GraPhi(Wx,1,PHI6x, .TRUE.)  
!  call GraPhi(Wx,2,PHI6y, .TRUE.)  
!  call GraPhi(Wx,3,PHI6z, .TRUE.)  

!  call GraPhi(Wy,1,PHI7x, .TRUE.)  
!  call GraPhi(Wy,2,PHI7y, .TRUE.)  
!  call GraPhi(Wy,3,PHI7z, .TRUE.)  

!  call GraPhi(Wz,1,PHI8x, .TRUE.)  
!  call GraPhi(Wz,2,PHI8y, .TRUE.)  
!  call GraPhi(Wz,3,PHI8z, .TRUE.)  

!  call GraPhi(Kx,1,PHI9x, .TRUE.)  
!  call GraPhi(Kx,2,PHI9y, .TRUE.)  
!  call GraPhi(Kx,3,PHI9z, .TRUE.)  


  do c = 1, NC
    R       = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 + tiny
    R1      = 1.0/R
    R2      = 1.0/(R*R)

    Gr =  PHIx(c) * xc(c) * R1   + PHIy(c) * yc(c) * R1
    Gt = -PHIx(c) * yc(c) * R1   + PHIy(c) * xc(c) * R1
    Gz =  PHIz(c)

    G1r =  PHI1x(c) * xc(c) * R1   + PHI1y(c) * yc(c) * R1
    G1t = -PHI1x(c) * yc(c) * R1   + PHI1y(c) * xc(c) * R1
    G1z =  PHI1z(c)

    G2r =  PHI2x(c) * xc(c) * R1   + PHI2y(c) * yc(c) * R1
    G2t = -PHI2x(c) * yc(c) * R1   + PHI2y(c) * xc(c) * R1
    G2z =  PHI2z(c)

    G3r =  PHI3x(c) * xc(c) * R1   + PHI3y(c) * yc(c) * R1
    G3t = -PHI3x(c) * yc(c) * R1   + PHI3y(c) * xc(c) * R1
    G3z =  PHI3z(c)

    G4r =  PHI4x(c) * xc(c) * R1   + PHI4y(c) * yc(c) * R1
    G4t = -PHI4x(c) * yc(c) * R1   + PHI4y(c) * xc(c) * R1
    G4z =  PHI4z(c)

    G5r =  PHI5x(c) * xc(c) * R1   + PHI5y(c) * yc(c) * R1
    G5t = -PHI5x(c) * yc(c) * R1   + PHI5y(c) * xc(c) * R1
    G5z =  PHI5z(c)

    G6r =  PHI6x(c) * xc(c) * R1   + PHI6y(c) * yc(c) * R1
    G6t = -PHI6x(c) * yc(c) * R1   + PHI6y(c) * xc(c) * R1
    G6z =  PHI6z(c)

    G7r =  PHI7x(c) * xc(c) * R1   + PHI7y(c) * yc(c) * R1
    G7t = -PHI7x(c) * yc(c) * R1   + PHI7y(c) * xc(c) * R1
    G7z =  PHI7z(c)

    G8r =  PHI8x(c) * xc(c) * R1   + PHI8y(c) * yc(c) * R1
    G8t = -PHI8x(c) * yc(c) * R1   + PHI8y(c) * xc(c) * R1
    G8z =  PHI8z(c)

    G9r =  PHI9x(c) * xc(c) * R1   + PHI9y(c) * yc(c) * R1
    G9t = -PHI9x(c) * yc(c) * R1   + PHI9y(c) * xc(c) * R1
    G9z =  PHI9z(c)

    Dift_uu(c)=Dift_uu(c) + Gr  + R1*G1t + G2z - 2.0*R1*vvu % mean(c) + R1*uuu % mean(c)  
    Dift_vv(c)=Dift_vv(c) + G3r + R1*G4t + G5z + 2.0*R1*vvu % mean(c) + R1*vvu % mean(c)
    Dift_ww(c)=Dift_ww(c) + G6r + R1*G7t + G8z +                        R1*wwu % mean(c)
    Dift_uv(c)=Dift_uv(c) + G1r + R1*G3t + G9z + 2.0*R1*uuv % mean(c) + R1*vvv % mean(c)
    Dift_uw(c)=Dift_uw(c) + G2r + R1*G9t + G6z -     R1*vvw % mean(c) + R1*uuw % mean(c)  
    Dift_vw(c)=Dift_vw(c) + G9r + R1*G5t + G7z + 2.0*R1*Kx(c) 

!    Dift_uu(c)=Dift_uu(c) + Gr  + R1*G1t + G2z - 2.0*R1*Vx(c) + R1*Ux(c)  
!    Dift_vv(c)=Dift_vv(c) + G3r + R1*G4t + G5z + 2.0*R1*Vx(c) + R1*Vx(c)
!    Dift_ww(c)=Dift_ww(c) + G6r + R1*G7t + G8z +                R1*Wx(c)
!    Dift_uv(c)=Dift_uv(c) + G1r + R1*G3t + G9z + 2.0*R1*Uy(c) + R1*Vy(c)
!    Dift_uw(c)=Dift_uw(c) + G2r + R1*G9t + G6z -     R1*Vz(c) + R1*Uz(c)  
!    Dift_vw(c)=Dift_vw(c) + G9r + R1*G5t + G7z + 2.0*R1*Kx(c) 

    Ux(c) = P % fluc(c) * U % fluc(c) 
    Uy(c) = P % fluc(c) * V % fluc(c) 
    Uz(c) = P % fluc(c) * W % fluc(c) 
  end do 

!----------------------------------------!
!  PRESSURE DIFFUSION
!----------------------------------------!
  call GraPhi(Ux,1,PHIx, .TRUE.)  
  call GraPhi(Ux,2,PHIy, .TRUE.)  
  call GraPhi(Ux,3,PHIz, .TRUE.)  

  call GraPhi(Uy,1,PHI1x, .TRUE.)  
  call GraPhi(Uy,2,PHI1y, .TRUE.)  
  call GraPhi(Uy,3,PHI1z, .TRUE.)  

  call GraPhi(Uz,1,PHI2x, .TRUE.)  
  call GraPhi(Uz,2,PHI2y, .TRUE.)  
  call GraPhi(Uz,3,PHI2z, .TRUE.)  
  
  do c = 1, NC
    R       = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 + tiny
    R1      = 1.0/R
    R2      = 1.0/(R*R)
    Gr =  PHIx(c) * xc(c) * R1   + PHIy(c) * yc(c) * R1
    Gt = -PHIx(c) * yc(c) * R1   + PHIy(c) * xc(c) * R1
    Gz =  PHIz(c)

    G1r =  PHI1x(c) * xc(c) * R1   + PHI1y(c) * yc(c) * R1
    G1t = -PHI1x(c) * yc(c) * R1   + PHI1y(c) * xc(c) * R1
    G1z =  PHI1z(c)

    G2r =  PHI2x(c) * xc(c) * R1   + PHI2y(c) * yc(c) * R1
    G2t = -PHI2x(c) * yc(c) * R1   + PHI2y(c) * xc(c) * R1
    G2z =  PHI2z(c)

    G3r =  PHI4x(c) * xc(c) * R1   + PHI4y(c) * yc(c) * R1
    G3t = -PHI4x(c) * yc(c) * R1   + PHI4y(c) * xc(c) * R1
    G3z =  PHI4z(c)

    PD_uu(c)=PD_uu(c) + 2.0*Gr
    PD_vv(c)=PD_vv(c) + 2.0*R1*G1t + 2.0*R1*P % fluc(c)*U % fluc(c)
    PD_ww(c)=PD_ww(c) + 2.0*G2z
    PD_uv(c)=PD_uv(c) + G1r + R1*Gt - R1*P % fluc(c)*V % fluc(c)
    PD_uw(c)=PD_uw(c) + G2r + Gz
    PD_vw(c)=PD_vw(c) + R1*G2t + G1z
  end do

!----------------------------------------!
!  TURBULENT DIFFUSION FOR HEAT FLUX
!----------------------------------------!
  if(new  > -1) then
    do c=1 ,NC
      R       = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 + tiny
      Diss_uu_mean(c) = ( Diss_uu_mean(c) * (1.*new) + Diss_uu(c) ) / (1.*(new+1))
      Diss_vv_mean(c) = ( Diss_vv_mean(c) * (1.*new) + Diss_vv(c) ) / (1.*(new+1))
      Diss_ww_mean(c) = ( Diss_ww_mean(c) * (1.*new) + Diss_ww(c) ) / (1.*(new+1))
      Diss_uv_mean(c) = ( Diss_uv_mean(c) * (1.*new) + Diss_uv(c) ) / (1.*(new+1))
      Diss_uw_mean(c) = ( Diss_uw_mean(c) * (1.*new) + Diss_uw(c) ) / (1.*(new+1))
      Diss_vw_mean(c) = ( Diss_vw_mean(c) * (1.*new) + Diss_vw(c) ) / (1.*(new+1))

      Diss_sgs_mean(c) = (Diss_sgs_mean(c)*(1.*new)+Diss_sgs(c) )/(1.*(new+1))

      Dift_uu_mean(c) = ( Dift_uu_mean(c) * (1.*new) + Dift_uu(c) ) / (1.*(new+1))
      Dift_vv_mean(c) = ( Dift_vv_mean(c) * (1.*new) + Dift_vv(c) ) / (1.*(new+1))
      Dift_ww_mean(c) = ( Dift_ww_mean(c) * (1.*new) + Dift_ww(c) ) / (1.*(new+1))
      Dift_uv_mean(c) = ( Dift_uv_mean(c) * (1.*new) + Dift_uv(c) ) / (1.*(new+1))
      Dift_uw_mean(c) = ( Dift_uw_mean(c) * (1.*new) + Dift_uw(c) ) / (1.*(new+1))
      Dift_vw_mean(c) = ( Dift_vw_mean(c) * (1.*new) + Dift_vw(c) ) / (1.*(new+1))

      PR_uu_mean(c) = ( PR_uu_mean(c) * (1.*new) + R_uu(c) ) / (1.*(new+1))
      PR_vv_mean(c) = ( PR_vv_mean(c) * (1.*new) + R_vv(c) ) / (1.*(new+1))
      PR_ww_mean(c) = ( PR_ww_mean(c) * (1.*new) + R_ww(c) ) / (1.*(new+1))
      PR_uv_mean(c) = ( PR_uv_mean(c) * (1.*new) + R_uv(c) ) / (1.*(new+1))
      PR_uw_mean(c) = ( PR_uw_mean(c) * (1.*new) + R_uw(c) ) / (1.*(new+1))
      PR_vw_mean(c) = ( PR_vw_mean(c) * (1.*new) + R_vw(c) ) / (1.*(new+1))

      PD_uu_mean(c) = ( PD_uu_mean(c) * (1.*new) + PD_uu(c) ) / (1.*(new+1))
      PD_vv_mean(c) = ( PD_vv_mean(c) * (1.*new) + PD_vv(c) ) / (1.*(new+1))
      PD_ww_mean(c) = ( PD_ww_mean(c) * (1.*new) + PD_ww(c) ) / (1.*(new+1))
      PD_uv_mean(c) = ( PD_uv_mean(c) * (1.*new) + PD_uv(c) ) / (1.*(new+1))
      PD_uw_mean(c) = ( PD_uw_mean(c) * (1.*new) + PD_uw(c) ) / (1.*(new+1))
      PD_vw_mean(c) = ( PD_vw_mean(c) * (1.*new) + PD_vw(c) ) / (1.*(new+1))

    end do
  end if

  deallocate (Puu)
  deallocate (Pvv)
  deallocate (Pww)
  deallocate (Puv)
  deallocate (Puw)
  deallocate (Pvw)
  deallocate (C_uu)
  deallocate (C_vv)
  deallocate (C_ww)
  deallocate (C_uv)
  deallocate (C_uw)
  deallocate (C_vw)
  deallocate (R_uu)
  deallocate (R_vv)
  deallocate (R_ww)
  deallocate (R_uv)
  deallocate (R_uw)
  deallocate (R_vw)
  deallocate (dUdr)
  deallocate (dUdt)
  deallocate (dUdz)
  deallocate (dVdr)
  deallocate (dVdt)
  deallocate (dVdz)
  deallocate (dWdr)
  deallocate (dWdt)
  deallocate (dWdz)
  deallocate (Diss_uu)
  deallocate (Diss_vv)
  deallocate (Diss_ww)
  deallocate (Diss_uv)
  deallocate (Diss_uw)
  deallocate (Diss_vw)
  deallocate (Diss_sgs)
  deallocate (Difv_uu)
  deallocate (Difv_vv)
  deallocate (Difv_ww)
  deallocate (Difv_uv)
  deallocate (Difv_uw)
  deallocate (Difv_vw)
  deallocate (Dift_uu)
  deallocate (Dift_vv)
  deallocate (Dift_ww)
  deallocate (Dift_uv)
  deallocate (Dift_uw)
  deallocate (Dift_vw)
  deallocate (PD_uu)
  deallocate (PD_vv)
  deallocate (PD_ww)
  deallocate (PD_uv)
  deallocate (PD_uw)
  deallocate (PD_vw)
  deallocate (uu_sgs)
  deallocate (vv_sgs)
  deallocate (ww_sgs)
  deallocate (uv_sgs)
  deallocate (uw_sgs)
  deallocate (vw_sgs)
  RETURN 

  END SUBROUTINE CalcBudgets_cylind
