!======================================================================!
  SUBROUTINE SourcesHJ(var)
!----------------------------------------------------------------------!
! Purpose:                                                             !
! Calculate source terms in equation for v2                            !
! Term which is negativ is put on left hand side in diagonal of        !
! martix of coefficient                                                !
!                                                                      !  
! Authors: Muhamed Hadziabdic and Bojan Niceno                         !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE rans_mod
  USE les_mod

!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, s, c1, c2, i, var,icont
  REAL    :: Prod, Diss, VAR_hom, VAR_wall, mag, VAR_tot, Esor
  REAL    :: a11, a22, a33, a12, a13, a21, a31, a23, a32
  REAL    :: n1,n2,n3,AA2,AA3,AA,Ret,feps,ff2,fd,FF1,CC,C1W,C2W,fw,uu_nn
  REAL    :: e11,e12,e13,e21,e22,e23,e31,e32,e33
  REAL    :: VAR_211,VAR_212,VAR_213,VAR_222,VAR_223,VAR_233,phi2_nn
  REAL    :: VAR_hm1,VAR_wal1,VAR_wal2
  REAL    :: fss,e11w,e12w,e13w,e22w,e23w,e33w,E2,E3,EE,CC1,CC2
!
!  REAL    :: b11, b22, b33, b12, b13, b21, b31, b23, b32
!  REAL    :: S11, S22, S33, S12, S13, S21, S31, S23, S32
!  REAL    :: V11, V22, V33, V12, V13, V21, V31, V23, V32
!  REAL    :: n1, n2, n3, AA, b_mn_b_mn, b_lk_s_lk, uiujn, Ce11, uu_nn 
!
  REAL    :: Uxx, Uyy, Uzz, Uxy, Uxz, Uyz, Uzy, Uzx, Uyx               
  REAL    :: Diss_wall, Diss_hom, r13, r23
  REAL,ALLOCATABLE :: Diss1(:)
!  REAL,ALLOCATABLE :: VAR1x(:),VAR1y(:),VAR1z(:)
!  REAL,ALLOCATABLE :: VAR2x(:),VAR2y(:),VAR2z(:)
!  REAL,ALLOCATABLE :: a11(:), a12(:), a13(:), a21(:), a22(:), a23(:)
!  REAL,ALLOCATABLE :: a31(:), a32(:), a33(:), AA(:), n1(:), n2(:), n3(:)
!--------------------------------[CVS]---------------------------------!
!  $Id: SourceV2KEPSV2F.f90,v 1.4 2008/11/19 14:57:57 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/SourceV2KEPSV2F.f90,v $  
!======================================================================!

  allocate (Diss1(1:NC)); Diss1 = 0.0

! calculate alpha as l=(k^3/2)/eps
! and normal as grad(alpha)
!
  do c=1,NC
    Kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)),1.0e-12)
    f22 % n(c)=  (Kin % n(c))**1.5/eps % n(c)
  end do
!
  do i=1,3
    if(i == 1) then
      call GraPhi(Ux,1,VAR2x, .TRUE.)  ! d2U/dxdx
      call GraPhi(Ux,2,VAR2y, .TRUE.)  ! d2U/dxdy
      call GraPhi(Ux,3,VAR2z, .TRUE.)  ! d2U/dxdz
      call GraPhi(Uy,2,VAR1x, .TRUE.) ! d2U/dydy
      call GraPhi(Uy,3,VAR1y, .TRUE.) ! d2U/dydz
      call GraPhi(Uz,3,VAR1z, .TRUE.) ! d2U/dzdz
    end if
    if(i == 2) then
      call GraPhi(Vx,1,VAR2x, .TRUE.)  ! d2V/dxdx
      call GraPhi(Vx,2,VAR2y, .TRUE.)  ! d2V/dxdy
      call GraPhi(Vx,3,VAR2z, .TRUE.)  ! d2V/dxdz
      call GraPhi(Vy,2,VAR1x, .TRUE.) ! d2V/dydy
      call GraPhi(Vy,3,VAR1y, .TRUE.) ! d2V/dydz
      call GraPhi(Vz,3,VAR1z, .TRUE.) ! d2V/dzdz
    end if
    if(i == 3) then
      call GraPhi(Wx,1,VAR2x, .TRUE.)  ! d2W/dxdx
      call GraPhi(Wx,2,VAR2y, .TRUE.)  ! d2W/dxdy
      call GraPhi(Wx,3,VAR2z, .TRUE.)  ! d2W/dxdz
      call GraPhi(Wy,2,VAR1x, .TRUE.) ! d2W/dydy
      call GraPhi(Wy,3,VAR1y, .TRUE.) ! d2W/dydz
      call GraPhi(Wz,3,VAR1z, .TRUE.) ! d2W/dzdz
    end if

    do c=1,NC
      if(i == 1) then
        Uxx = VAR2x(c)
        Uxy = VAR2y(c)
        Uyx = Uxy
        Uxz = VAR2z(c)
        Uzx = Uxz
        Uyy = VAR1x(c)
        Uyz = VAR1y(c)
        Uzy = Uyz
        Uzz = VAR1z(c)
        Diss1(c) = &
                uu % n(c)*(Uxx*Uxx+Uxy*Uxy+Uxz*Uxz)+&
                uv % n(c)*(Uxx*Uyx+Uxy*Uyy+Uxz*Uyz)+&
                uw % n(c)*(Uxx*Uzx+Uxy*Uzy+Uxz*Uzz)+&
                uv % n(c)*(Uyx*Uxx+Uyy*Uxy+Uyz*Uxz)+&
                vv % n(c)*(Uyx*Uyx+Uyy*Uyy+Uyz*Uyz)+&
                vw % n(c)*(Uyx*Uzx+Uyy*Uzy+Uyz*Uzz)+&
                uw % n(c)*(Uzx*Uxx+Uzy*Uxy+Uzz*Uxz)+&
                vw % n(c)*(Uzx*Uyx+Uzy*Uyy+Uzz*Uyz)+&
                ww % n(c)*(Uzx*Uzx+Uzy*Uzy+Uzz*Uzz)
      end if
      if(i == 2) then
        Uxx = VAR2x(c)
        Uxy = VAR2y(c)
        Uyx = Uxy
        Uxz = VAR2z(c)
        Uzx = Uxz
        Uyy = VAR1x(c)
        Uyz = VAR1y(c)
        Uzy = Uyz
        Uzz = VAR1z(c)
        Diss1(c) = Diss1(c) +&
                uu % n(c)*(Uxx*Uxx+Uxy*Uxy+Uxz*Uxz)+&
                uv % n(c)*(Uxx*Uyx+Uxy*Uyy+Uxz*Uyz)+&
                uw % n(c)*(Uxx*Uzx+Uxy*Uzy+Uxz*Uzz)+&
                uv % n(c)*(Uyx*Uxx+Uyy*Uxy+Uyz*Uxz)+&
                vv % n(c)*(Uyx*Uyx+Uyy*Uyy+Uyz*Uyz)+&
                vw % n(c)*(Uyx*Uzx+Uyy*Uzy+Uyz*Uzz)+&
                uw % n(c)*(Uzx*Uxx+Uzy*Uxy+Uzz*Uxz)+&
                vw % n(c)*(Uzx*Uyx+Uzy*Uyy+Uzz*Uyz)+&
                ww % n(c)*(Uzx*Uzx+Uzy*Uzy+Uzz*Uzz)
      end if
      if(i == 3) then
        Uxx = VAR2x(c)
        Uxy = VAR2y(c)
        Uyx = Uxy
        Uxz = VAR2z(c)
        Uzx = Uxz
        Uyy = VAR1x(c)
        Uyz = VAR1y(c)
        Uzy = Uyz
        Uzz = VAR1z(c)
        Diss1(c) = Diss1(c) +&
                uu % n(c)*(Uxx*Uxx+Uxy*Uxy+Uxz*Uxz)+&
                uv % n(c)*(Uxx*Uyx+Uxy*Uyy+Uxz*Uyz)+&
                uw % n(c)*(Uxx*Uzx+Uxy*Uzy+Uxz*Uzz)+&
                uv % n(c)*(Uyx*Uxx+Uyy*Uxy+Uyz*Uxz)+&
                vv % n(c)*(Uyx*Uyx+Uyy*Uyy+Uyz*Uyz)+&
                vw % n(c)*(Uyx*Uzx+Uyy*Uzy+Uyz*Uzz)+&
                uw % n(c)*(Uzx*Uxx+Uzy*Uxy+Uzz*Uxz)+&
                vw % n(c)*(Uzx*Uyx+Uzy*Uyy+Uzz*Uyz)+&
                ww % n(c)*(Uzx*Uzx+Uzy*Uzy+Uzz*Uzz)
      end if
    end do
  end do  ! i
!
  call GraPhi(f22 % n,1,VAR2x,.TRUE.)             ! df22/dx
  call GraPhi(f22 % n,2,VAR2y,.TRUE.)             ! df22/dy
  call GraPhi(f22 % n,3,VAR2z,.TRUE.)             ! df22/dz
!
!!!!  call Scale

  r23 = 2.0/3.0 
  do  c = 1, NC
!    Kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)),1.0e-12)
    Pk(c)      = max(-(uu % n(c)*Ux(c) + uv % n(c)*Uy(c) + uw % n(c)*Uz(c) +&
                   uv % n(c)*Vx(c) + vv % n(c)*Vy(c) + vw % n(c)*Vz(c) +&
                   uw % n(c)*Wx(c) + vw % n(c)*Wy(c) + ww % n(c)*Wz(c)),1.0e-12)                
  
    mag = max(1.0e-12,sqrt(VAR2x(c)*VAR2x(c)+VAR2y(c)*VAR2y(c)+VAR2z(c)*VAR2z(c)))       
    n1 = VAR2x(c)/mag 
    n2 = VAR2y(c)/mag 
    n3 = VAR2z(c)/mag 


    a11 = uu % n(c)/(Kin % n(c)) - r23 
    a22 = vv % n(c)/(Kin % n(c)) - r23
    a33 = ww % n(c)/(Kin % n(c)) - r23
    a12 = uv % n(c)/(Kin % n(c))   
    a21 = a12
    a13 = uw % n(c)/(Kin % n(c))    
    a31 = a13
    a23 = vw % n(c)/(Kin % n(c))    
    a32 = a23
    
    AA2=(a11**2)+(a22**2)+(a33**2)+2*((a12**2)+(a13**2)+(a23**2))
    AA3=a11*(a11*a11+a12*a21+a13*a31)+ &
     2*a12*(a21*a11+a22*a21+a23*a31)+ &
     2*a13*(a31*a11+a32*a21+a33*a31)+ &
       a22*(a21*a21+a22*a22+a23*a32)+ &
     2*a13*(a31*a11+a32*a21+a33*a31)+ &
       a33*(a31*a13+a32*a23+a33*a33)
    AA=1.0 - (9.0/8.0)*(AA2-AA3)
    AA=max(AA,0.0)
    AA=min(AA,1.0)
 
    Tsc(c)=Kin % n(c) / (Eps % n(c) + tiny)

    Ret= (Kin % n(c)**2)/(visc*eps % n(c)+tiny)
    Ret=max(Ret,tiny)
   
    feps=-((1.92-1.44)/1.92)*exp(-(Ret/6.0)**2)
    ff2=min((Ret/150)**1.5, 1.0)
    fd=1.0/(1.0+0.1*Ret)
    FF1=min(0.6, AA2)
    CC=2.5*AA*FF1**0.25*ff2
    C1W=max((1.0-0.7*C), 0.3)
    C2W=min (AA2,0.3)
    fw=min((Kin%n(c)**1.5)/(2.5*Eps%n(c)*WallDs(c)+tiny),1.4)
!     Write(*,*) 'c=',c
!     Write(*,*) 'aa-dep coeff',Pk(c),n1,n2,n3,aa,Tsc(c),ret,feps,ff2,fd,ff1,cc,c1w,c2w,fw

!    S11 = Ux(c) 
!    S22 = Vy(c) 
!    S33 = Wz(c) 
!    S12 = 0.5*(Uy(c)+Vx(c)) 
!    S21 = S12
!    S13 = 0.5*(Uz(c)+Wx(c)) 
!    S31 = S13 
!    S23 = 0.5*(Vz(c)+Wy(c)) 
!    S32 = S23 

!    V11 = 0.0
!    V22 = 0.0
!    V33 = 0.0
!    V12 = 0.5*(Uy(c)-Vx(c)) - omegaZ
!    V21 = -V12 + omegaZ
!    V13 = 0.5*(Uz(c)-Wx(c)) + omegaY
!    V31 = -V13 - omegaY
!    V23 = 0.5*(Vz(c)-Wy(c)) - omegaX
!    V32 = -V23 + omegaX

!    b_mn_b_mn = b11*b11 + b22*b22 + b33*b33 + 2.0*(b12*b12+b13*b13+b23*b23)
!    b_lk_s_lk = b11*S11 + b22*S22 + b33*S33 + 2.0*(b12*S12+b13*S13+b23*S23)
    uu_nn     = (uu % n(c)*n1*n1+uv % n(c)*n1*n2+uw % n(c)*n1*n3 &
               + uv % n(c)*n2*n1+vv % n(c)*n2*n2+vw % n(c)*n2*n3 &
               + uw % n(c)*n3*n1+vw % n(c)*n3*n2+ww % n(c)*n3*n3)

      e11w = (1.0/Tsc(c))* &
                  (uu % n(c) + (2*(uu % n(c)*n1*n1 + uv % n(c)*n1*n2 + uw % n(c) * n1*n3) + &
                   uu_nn)*fd)/(1.0 + 1.5*(uu_nn/Kin % n(c)*fd))
      e12w = (1.0/Tsc(c))* &
                  (uv % n(c) + (uu % n(c)*n2*n1 + uv % n(c)*n2*n2 + uw % n(c) * n2*n3 + &
                                uv % n(c)*n1*n1 + vv % n(c)*n1*n2 + vw % n(c) * n1*n3 + &
                   uu_nn)*fd)/(1.0 + 1.5*(uu_nn/Kin % n(c)*fd))
      e13w = (1.0/Tsc(c))* &
                  (uw % n(c) + (uu % n(c)*n3*n1 + uv % n(c)*n3*n2 + uw % n(c) * n3*n3 + &
                                uw % n(c)*n1*n1 + vw % n(c)*n1*n2 + ww % n(c) * n1*n3 + &
                   uu_nn)*fd)/(1.0 + 1.5*(uu_nn/Kin % n(c)*fd))
      e22w = (1.0/Tsc(c))* &
                  (vv % n(c) + (2*(uv % n(c)*n1*n2 + vv % n(c)*n2*n2 + vw % n(c) * n2*n3) + &
                   uu_nn)*fd)/(1.0 + 1.5*(uu_nn/Kin % n(c)*fd))
      e23w = (1.0/Tsc(c))* &
                  (vw % n(c) + (uv % n(c)*n3*n1 + vv % n(c)*n3*n2 + vw % n(c) * n3*n3 + &
                                uw % n(c)*n1*n2 + vw % n(c)*n2*n2 + ww % n(c) * n2*n3 + &
                   uu_nn)*fd)/(1.0 + 1.5*(uu_nn/Kin % n(c)*fd))
      e33w = (1.0/Tsc(c))* &
                  (ww % n(c) + (2*(uw % n(c)*n1*n3 + vw % n(c)*n2*n3 + ww % n(c) * n3*n3) + &
                   uu_nn)*fd)/(1.0 + 1.5*(uu_nn/Kin % n(c)*fd))
!!     Write(*,*) 'EE construction',e11w,e12w,e13w,e22w,e23w,e33w,uu_nn
      EE=AA
      do icont=1,10
        fss=1.0-(AA**0.5*EE**2.0)
        e11=(fss*e11w+(1.0-fss)*r23*Eps % n(c))/(eps%n(c))-r23
        e22=(fss*e22w+(1.0-fss)*r23*Eps % n(c))/(eps%n(c))-r23
        e33=(fss*e33w+(1.0-fss)*r23*Eps % n(c))/(eps%n(c))-r23
        e12=fss*e12w/(eps%n(c))
        e13=fss*e13w/(eps%n(c))
        e23=fss*e23w/(eps%n(c))
        e21=e12
        e31=e13
        e32=e23
        E2=(e11**2)+(e22**2)+(e33**2)+2*((e12**2)+(e13**2)+(e23**2))
        E3=e11*(e11*e11+e12*e21+e13*e31)+ &
         2*e12*(e21*e11+e22*e21+e23*e31)+ &
         2*e13*(e31*e11+e32*e21+e33*e31)+ &
           e22*(e21*e21+e22*e22+e23*e32)+ &
         2*e13*(e31*e11+e32*e21+e33*e31)+ &
           e33*(e31*e13+e32*e23+e33*e33)
        EE=1.0 - (9.0/8.0)*(E2-E3)
        EE=max(EE,0.0)
        EE=min(EE,1.0)
!        write(*,*) 'ee cycle',ee,e2,e3,e11,e12,e13,e22,e23,e33,icont,eps%n(c)
      enddo
      fss=1.0-(AA**0.5*EE**2.0)
      e11=fss*e11w+(1.0-fss)*r23*Eps % n(c)
      e22=fss*e22w+(1.0-fss)*r23*Eps % n(c)
      e33=fss*e33w+(1.0-fss)*r23*Eps % n(c)
      e12=fss*e12w
      e13=fss*e13w
      e23=fss*e23w
      CC1=CC+(AA**0.5)*(EE**2)
      CC2=0.8*(AA**0.5)
!     write(*,*) 'e-dep coeff',fss,e11,e22,e33,e12,e13,e23,cc1,cc2,EE
!
! calculate the rapid redidtribution term (used also in the second term of Gibson-Launder)
!
      VAR_211 = -CC2*(-2*(uu%n(c)*Ux(c)+uv%n(c)*Uy(c)+uw%n(c)*Uz(c))-r23*Pk(c))
      VAR_212 = -CC2*(-(uu%n(c)*Vx(c)+uv%n(c)*(Ux(c)+Vy(c))+uw%n(c)*Vx(c)+vv%n(c)*Uy(c)+vw%n(c)*Uz(c)))
      VAR_213 = -CC2*(-(uu%n(c)*Wx(c)+uw%n(c)*(Ux(c)+Wz(c))+uv%n(c)*Wy(c)+vw%n(c)*Uz(c)+ww%n(c)*Uz(c)))
      VAR_222 = -CC2*(-2*(uv%n(c)*Vx(c)+vv%n(c)*Vy(c)+vw%n(c)*Vz(c))-r23*Pk(c))
      VAR_223 = -CC2*(-(uv%n(c)*Wx(c)+uw%n(c)*Vx(c)+vv%n(c)*Wy(c)+vw%n(c)*(Vy(c)+Wz(c))+ww%n(c)*Vz(c)))
      VAR_233 = -CC2*(-2*(uw%n(c)*Wx(c)+vw%n(c)*Wy(c)+ww%n(c)*Wz(c))-r23*Pk(c))
      phi2_nn = VAR_211*n1*n1+2*VAR_212*n1*n2+2*VAR_213*n1*n3+VAR_222*n2*n2+2*VAR_223*n2*n3+VAR_233*n3*n3  
!!     write(*,*) 'VAR_2',VAR_211,VAR_212,VAR_213,VAR_222,VAR_223,VAR_233,c

!------- uu stress
    if(var == 6) then
!
      VAR_hm1 = -CC1*Eps%n(c)*a11
      VAR_hom = VAR_hm1+VAR_211
!
      VAR_wal1 = C1W*fw/Tsc(c)*(uu_nn-1.5*2.0*(uu%n(c)*n1*n1+uv%n(c)*n1*n2+uw%n(c)*n1*n3))
      VAR_wal2 = C2W*fw/Tsc(c)*(phi2_nn-1.5*2.0*(VAR_211*n1*n1+VAR_212*n1*n2+VAR_213*n1*n3))
      VAR_wall = VAR_wal1+VAR_wal2
!
      Prod = -2.0*(uv % n(c)*Uy(c) + uw % n(c)*Uz(c))  &
             -2.0*omegaY*2.0*uw%n(c) + 2.0*omegaZ*2.0*uv%n(c)   


!!      Diss_wall = uu % n(c)/Kin % n(c) * Eps % n(c) 
!!      Diss_hom  = 2.0/3.0 * Eps % n(c)

      b(c) = b(c) + (max(Prod,0.0) -  e11 + VAR_hom + VAR_wall) *volume(c) 
!
!!      b(c) = b(c) + (max(Prod,0.0) + (1.0-f22 % n(c)*f22 % n(c))*VAR_wall +&
!!             f22 % n(c)*f22 % n(c)*(VAR_hom))*volume(c)

      Aval(Adia(c)) =  Aval(Adia(c)) + (max(-Prod,0.0)/max(uu%n(c),1.0e-12)+ 2.0* Ux(c))*volume(c)
!!      Aval(Adia(c)) =  Aval(Adia(c)) + (max(-Prod,0.0)/max(uu%n(c),1.0e-12) + 2.0*Ux(c) + (1.0-f22%n(c)*f22%n(c))*&
!!                      6.0*Eps%n(c)/Kin%n(c) +&         
!!                      f22%n(c)*f22%n(c)*(Diss_hom/max(uu%n(c),1.0e-12)+g1*Eps%n(c)/(2.0*Kin%n(c))&
!!                      +g1_star*Pk(c)/(2.0*Kin%n(c))))*volume(c)

!------- vv stress
    else if(var == 7) then
!
      VAR_hm1 = -CC1*Eps%n(c)*a22
      VAR_hom = VAR_hm1+VAR_222
!
      VAR_wal1 = C1W*fw/Tsc(c)*(uu_nn-1.5*2.0*(uv%n(c)*n1*n2+vv%n(c)*n2*n2+vw%n(c)*n2*n3))
      VAR_wal2 = C2W*fw/Tsc(c)*(phi2_nn-1.5*2.0*(VAR_212*n1*n2+VAR_222*n2*n2+VAR_223*n2*n3))
      VAR_wall = VAR_wal1+VAR_wal2
!
      Prod = -2.0*(uv % n(c)*Vx(c) + vw % n(c)*Vz(c))  &
             +2.0*omegaX*2.0*vw%n(c) - 2.0*omegaZ*2.0*uw%n(c)

!!      Diss_wall = vv % n(c)/Kin % n(c) * Eps % n(c)
!!      Diss_hom  = 2.0/3.0 * Eps % n(c)

      b(c) = b(c) + (max(Prod,0.0) -  e22 + VAR_hom + VAR_wall) * volume(c)  

!!      b(c) = b(c) + (max(Prod,0.0) + (1.0-f22 % n(c)*f22 % n(c))*VAR_wall +&
!!             f22 % n(c)*f22 % n(c)*(VAR_hom))*volume(c)

      Aval(Adia(c)) =  Aval(Adia(c)) + (max(-Prod,0.0)/max(vv%n(c),1.0e-12)+2.0*Vy(c))*volume(c)
 

!!      Aval(Adia(c)) =  Aval(Adia(c)) + (max(-Prod,0.0)/max(vv%n(c),1.0e-12) + 2.0*Vy(c)+(1.0-f22%n(c)*f22%n(c))*&
!!                      6.0*Eps%n(c)/Kin%n(c) &
!!                    + f22%n(c)*f22%n(c)*(Diss_hom/max(vv%n(c),1.0e-12)+g1*Eps%n(c)/(2.0*Kin%n(c))+g1_star*Pk(c)/(2.0*Kin%n(c))))*volume(c)


!------- ww stress
    else if(var == 8) then
!
      VAR_hm1 = -CC1*Eps%n(c)*a33
      VAR_hom = VAR_hm1+VAR_233
!
      VAR_wal1 = C1W*fw/Tsc(c)*(uu_nn-1.5*2.0*(uw%n(c)*n1*n3+vw%n(c)*n2*n3+ww%n(c)*n3*n3))
      VAR_wal2 = C2W*fw/Tsc(c)*(phi2_nn-1.5*2.0*(VAR_213*n1*n3+VAR_223*n2*n3+VAR_233*n3*n3))
      VAR_wall = VAR_wal1+VAR_wal2
!

      b(c) = b(c) + (max(Prod,0.0) -  e33 + VAR_hom + VAR_wall) * volume(c)

      Aval(Adia(c)) =  Aval(Adia(c)) + (max(-Prod,0.0)/max(ww%n(c),1.0e-12)+2.0*Wz(c))*volume(c)

!!      b(c) = b(c) + (max(Prod,0.0) + (1.0-f22 % n(c)*f22 % n(c))*VAR_wall +&
!!             f22 % n(c)*f22 % n(c)*(VAR_hom))*volume(c)
!!       Aval(Adia(c)) =  Aval(Adia(c)) + (max(-Prod,0.0)/max(ww%n(c),1.0e-12)+2.0*Wz(c)+(1.0-f22%n(c)*f22%n(c))*&
!!                      6.0*Eps%n(c)/Kin%n(c)          &
!!                    + f22%n(c)*f22%n(c)*(Diss_hom/max(ww%n(c),1.0e-12)+g1*Eps%n(c)/(2.0*Kin%n(c))+g1_star*Pk(c)/(2.0*Kin%n(c))))*volume(c)

!------- uv stress
    else if(var == 9) then
!
      VAR_hm1 = -CC1*Eps%n(c)*a12
      VAR_hom = VAR_hm1+VAR_212
!
      VAR_wal1 = C1W*fw/Tsc(c)*(-1.5*(uu%n(c)*n1*n2+uv%n(c)*(n1*n1+n2*n2)+vv%n(c)*n1*n2+uw%n(c)*n2*n3+vw%n(c)*n1*n3))
      VAR_wal2 = C2W*fw/Tsc(c)*(-1.5*(VAR_211*n1*n2+VAR_212*(n1*n1+n2*n2)+VAR_222*n1*n2+VAR_213*n2*n3+VAR_223*n1*n3))
      VAR_wall = VAR_wal1+VAR_wal2
!
      Prod = -(uu % n(c)*Vx(c) + uw % n(c)*Vz(c) + &
               vv % n(c)*Uy(c) + vw % n(c)*Uz(c)) &
             +2.0*omegaX*uw%n(c) - 2.0*omegaY*vw%n(c) + 2.0*omegaZ*(vv%n(c)-uu%n(c))  

      b(c) = b(c) + (Prod  -  e12 + VAR_hom + VAR_wall) * volume(c)
      Aval(Adia(c)) =  Aval(Adia(c)) + (Vy(c)+Ux(c))*volume(c)


!------- uw stress
    else if(var == 10) then
!
      VAR_hm1 = -CC1*Eps%n(c)*a13
      VAR_hom = VAR_hm1+VAR_213
!
      VAR_wal1 = C1W*fw/Tsc(c)*(-1.5*(uu%n(c)*n1*n3+uw%n(c)*(n1*n1+n3*n3)+vw%n(c)*n1*n2+uv%n(c)*n2*n3+ww%n(c)*n1*n3))
      VAR_wal2 = C2W*fw/Tsc(c)*(-1.5*(VAR_211*n1*n3+VAR_213*(n1*n1+n3*n3)+VAR_212*n1*n3+VAR_223*n1*n2+VAR_233*n1*n3))
      VAR_wall = VAR_wal1+VAR_wal2
!
      Prod = -(uu % n(c)*Wx(c) + uv % n(c)*Wy(c)+ &
               vw % n(c)*Uy(c) + ww % n(c)*Uz(c)) & 
             -2.0*omegaX*uv%n(c)-2.0*omegaY*(ww%n(c)-uu%n(c))+2.0*omegaZ*vw%n(c)   

      b(c) = b(c) + (Prod -  e13 + VAR_hom + VAR_wall) * volume(c)
      Aval(Adia(c)) =  Aval(Adia(c)) + (Wz(c)+Ux(c))*volume(c)

!------- vw stress
    else if(var == 11) then
!
      VAR_hm1 = -CC1*Eps%n(c)*a23
      VAR_hom = VAR_hm1+VAR_223
!
      VAR_wal1 = C1W*fw/Tsc(c)*(-1.5*(uv%n(c)*n1*n3+uw%n(c)*n1*n2+vv%n(c)*n2*n3+vw%n(c)*(n2*n2+n3*n3)+ww%n(c)*n2*n3))
      VAR_wal2 = C2W*fw/Tsc(c)*(-1.5*(VAR_212*n1*n3+VAR_213*n1*n2+VAR_222*n2*n3+VAR_223*(n2*n2+n3*n3)+VAR_233*n2*n3))
      VAR_wall = VAR_wal1+VAR_wal2
!
      Prod = -(uv % n(c)*Wx(c) + vv % n(c)*Wy(c)+ &
               uw % n(c)*Vx(c) + ww % n(c)*Vz(c))  &
             -2.0*omegaX*(vw%n(c)-ww%n(c))+2.0*omegaY*uv%n(c)-2.0*omegaZ*uw%n(c)   
!
!      write(*,*) '23',cc,aa,ee,ff1,ff2  
!      write(*,*) '23',cc1,Eps%n(c),a23,VAR_hm1,VAR_hom,VAR_wal1,VAR_wal2,Prod,Wz(c),Vy(c),c
      b(c) = b(c) + (Prod -e23 + VAR_hom + VAR_wall) * volume(c)
      Aval(Adia(c)) =  Aval(Adia(c)) + (Wz(c)+Vy(c))*volume(c)
!
!------- Eps eq.
    else if(var == 13) then


      Esor = volume(c)/max(Tsc(c),1.0e-12)
!!      uiujn = max(uu%n(c)*n1*n1+uv%n(c)*n2*n1+uw%n(c)*n3*n1+&
!!                  uv%n(c)*n1*n2+vv%n(c)*n2*n2+vw%n(c)*n3*n2+&
!!                  uw%n(c)*n1*n3+vw%n(c)*n2*n3+ww%n(c)*n3*n3,1.0e-12)
!
!! PSI
!
      Diss1(c) = Ce3*VISc*Diss1(c)    
      b(c) = b(c) + (Ce1*Pk(c)+Ce3*visc*Diss1(c))*Esor 

!----- Fill in a diagonal of coefficient matrix
      Aval(Adia(c)) =  Aval(Adia(c)) + Ce2*Esor*DENc(material(c))
    end if
  end do

  if(var == 13) then
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)

!---- Calculate a values of dissipation  on wall
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER ) then
        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
          Eps%n(c2) = VISc*(uu%n(c1)+vv%n(c1)+ww%n(c1))/WallDs(c1)**2
        end if   ! end if of BC=wall
      end if    ! end if of c2<0
    end do
  end if
  deallocate (Diss1)
  RETURN
  END SUBROUTINE SourcesHJ   
