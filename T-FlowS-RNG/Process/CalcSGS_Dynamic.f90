!======================================================================!
  SUBROUTINE CalcSGS_Dynamic()
!----------------------------------------------------------------------!
!   Calculates Smagorinsky constant with dynamic procedure             !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER           :: c, j, cj, m, n 
  REAL              :: Ua, Va, Wa, UUa, VVa, WWa, UVa, UWa, VWa, DE
  REAL              :: M11a, M22a, M33a, M12a, M13a, M23a      
  REAL              :: L11, L22, L33, L12, L13, L23      
  REAL              :: M11, M22, M33, M12, M13, M23      
  REAL              :: MdotM, LdotM, Lg, Lf 
  REAL              :: MinC, MaxC, Cplus_avr, Cminus_avr
  REAL              :: fun

!--------------------------------[CVS]---------------------------------!
!  $Id: CalcSGS_Dynamic.f90,v 1.6 2008/12/10 14:10:13 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CalcSGS_Dynamic.f90,v $  
!======================================================================!
!                                                                      !
!  C is derived from:    Lij_res = Lij_mod                             !
!                                                                      !
!  res - resolved, mod - modeled                                       ! 
!                                                                      ! 
!  Lij_res = <UiUj> - <Ui><Uj>,                                        !
!                                                                      !
!  where <> denote test filter                                         !
!                                                                      !
!  Lij_mod = 2 * C * Mij, C = Csmag ** 2.0                             !
!                                                                      !
!  where Mij is:  Mij = <delta**2.0>|<Sij>|<Sij> - <delta |Sij| Sij>   !
!                                                                      ! 
!  Finaly C is :                                                       ! 
!                                                                      !
!  C = 0.5 * Lij:Mij / Mij:Mij                                         !
!                                                                      !
!  Aij : Bij = A11 * B11 + A22 * B22 + A33 * B33                       !
!              + 2.0 * A12 * B12 + 2.0 A13 * B13 + 2.0 * A23 * B23     !   
!======================================================================!
!  MinC = 1000000.0
!  MaxC = 0.0
!  Cplus_avr = 0.0
!  Cminus_avr = 0.0
!  m = 0
!  n = 0

  call Exchng(U % n)
  call Exchng(V % n)
  call Exchng(W % n)
  call Exchng(Shear)

  do c =1, NC
    Ua   = 0.0
    Va   = 0.0
    Wa   = 0.0
    DE   = 0.0

    UUa  = 0.0
    VVa  = 0.0
    WWa  = 0.0
    UVa  = 0.0
    UWa  = 0.0
    VWa  = 0.0

    M11a = 0.0
    M22a = 0.0
    M33a = 0.0
    M12a = 0.0
    M13a = 0.0
    M23a = 0.0
  
    do j = Acol(c), Acol(c + 1) - 1
      cj = Arow(j) 
      if(cj /= c) then

!--- Test velocitys
        Ua = Ua + volume(cj) * U % n(cj)
        Va = Va + volume(cj) * V % n(cj)
        Wa = Wa + volume(cj) * W % n(cj)

!--- Test stresses
        UUa = UUa + volume(cj) * U % n(cj) * U % n(cj)
        VVa = VVa + volume(cj) * V % n(cj) * V % n(cj)
        WWa = WWa + volume(cj) * W % n(cj) * W % n(cj)
        UVa = UVa + volume(cj) * U % n(cj) * V % n(cj)
        UWa = UWa + volume(cj) * U % n(cj) * W % n(cj)
        VWa = VWa + volume(cj) * V % n(cj) * W % n(cj)

!--- Test Mija
        M11a = M11a + volume(cj) * Shear(cj) * Ux(cj)
        M22a = M22a + volume(cj) * Shear(cj) * Vy(cj)
        M33a = M33a + volume(cj) * Shear(cj) * Wz(cj)
        M12a = M12a + volume(cj) * Shear(cj)                &    
               * 0.5 * ( Uy(cj) + Vx(cj) ) 
        M13a = M13a + volume(cj) * Shear(cj)                &
               * 0.5 * ( Uz(cj) + Wx(cj) )     
        M23a = M23a + volume(cj) * Shear(cj)                &
               * 0.5 * ( Vz(cj) + Wy(cj) )

!--- Test volume 
        DE = DE + volume(cj) 
      end if
    end do

!--- Now it's taking into account influence of central cell
!--- within test molecule
    
    DE = DE + volume(c)

    Ua = Ua + volume(c) * U % n(c)
    Va = Va + volume(c) * V % n(c)
    Wa = Wa + volume(c) * W % n(c)

    UUa = UUa + volume(c) * U % n(c) * U % n(c)
    VVa = VVa + volume(c) * V % n(c) * V % n(c)
    WWa = WWa + volume(c) * W % n(c) * W % n(c)
    UVa = UVa + volume(c) * U % n(c) * V % n(c)
    UWa = UWa + volume(c) * U % n(c) * W % n(c)
    VWa = VWa + volume(c) * V % n(c) * W % n(c)

    M11a = M11a + volume(c) * Shear(c) * Ux(c)
    M22a = M22a + volume(c) * Shear(c) * Vy(c)
    M33a = M33a + volume(c) * Shear(c) * Wz(c)
    M12a = M12a + volume(c) * Shear(c) * 0.5 * ( Uy(c) + Vx(c) ) 
    M13a = M13a + volume(c) * Shear(c) * 0.5 * ( Uz(c) + Wx(c) )
    M23a = M23a + volume(c) * Shear(c) * 0.5 * ( Vz(c) + Wy(c) )

    
!--- Now calculating test values
    U % filt(c) = Ua / ( DE)
    V % filt(c) = Va / ( DE)
    W % filt(c) = Wa / ( DE)

    UUf(c)      = UUa / ( DE)
    VVf(c)      = VVa / ( DE)
    WWf(c)      = WWa / ( DE)
    UVf(c)      = UVa / ( DE)
    UWf(c)      = UWa / ( DE)
    VWf(c)      = VWa / ( DE)
  
    M11f(c)     = M11a / ( DE) 
    M22f(c)     = M22a / ( DE) 
    M33f(c)     = M33a / ( DE) 
    M12f(c)     = M12a / ( DE) 
    M13f(c)     = M13a / ( DE) 
    M23f(c)     = M23a / ( DE) 
  end do

  call GraPhi(U % filt, 1, Ux,.TRUE.)    ! dU/dx
  call GraPhi(U % filt, 2, Uy,.TRUE.)    ! dU/dy
  call GraPhi(U % filt, 3, Uz,.TRUE.)    ! dU/dz
  call GraPhi(V % filt, 1, Vx,.TRUE.)    ! dV/dx
  call GraPhi(V % filt, 2, Vy,.TRUE.)    ! dV/dy
  call GraPhi(V % filt, 3, Vz,.TRUE.)    ! dV/dz
  call GraPhi(W % filt, 1, Wx,.TRUE.)    ! dW/dx
  call GraPhi(W % filt, 2, Wy,.TRUE.)    ! dW/dy
  call GraPhi(W % filt, 3, Wz,.TRUE.)    ! dW/dz

  do c=1,NC
    Lg  = volume(c)**0.33333    
    Lf  = 2.0 * Lg

    ShearTest(c) = sqrt(2.0*(Ux(c)*Ux(c) + Vy(c)*Vy(c) + Wz(c)*Wz(c) + &
             0.5*(Vz(c) + Wy(c))*(Vz(c) + Wy(c)) + &
             0.5*(Uz(c) + Wx(c))*(Uz(c) + Wx(c)) + &
             0.5*(Vx(c) + Uy(c))*(Vx(c) + Uy(c))))

    L11 = UUf(c) - U % filt(c) * U % filt(c) 
    L22 = VVf(c) - V % filt(c) * V % filt(c) 
    L33 = WWf(c) - W % filt(c) * W % filt(c) 
    L12 = UVf(c) - U % filt(c) * V % filt(c) 
    L13 = UWf(c) - U % filt(c) * W % filt(c) 
    L23 = VWf(c) - V % filt(c) * W % filt(c) 

    M11 = Lf * Lf * ShearTest(c) * Ux(c)       &
         - Lg * Lg * M11f(c) 
    M22 = Lf * Lf * ShearTest(c) * Vy(c)       &
          - Lg * Lg * M22f(c) 
    M33 = Lf * Lf * ShearTest(c) * Wz(c)       &
          - Lg * Lg * M33f(c) 
    M12 = Lf * Lf * ShearTest(c)                &
           * 0.5 * (Uy(c) + Vx(c)) - Lg * Lg * M12f(c)

    M13 = Lf * Lf * ShearTest(c)                &
          * 0.5 * (Uz(c) + Wx(c)) - Lg * Lg * M13f(c) 

    M23 = Lf * Lf * ShearTest(c)                &
          * 0.5 * (Vz(c) + Wy(c)) - Lg * Lg * M23f(c)  

    MdotM = (M11 * M11 + M22 * M22 + M33 * M33 )                        & 
            + 2.0 * M12 * M12 + 2.0 * M13 * M13 + 2.0 * M23 * M23 

    LdotM =( L11 * M11 + L22 * M22 + L33 * M33 )                       & 
            + 2.0 * L12 * M12 + 2.0 * L13 * M13 + 2.0 * L23 * M23

    Cdyn(c)  =  -0.5 * LdotM / (MdotM + tiny) 

!    MaxC = max(MaxC,Cdyn(c)) 
!    MinC = min(MinC,Cdyn(c)) 

!    if(Cdyn(c) > 0.0) then
!      Cplus_avr = Cplus_avr + Cdyn(c)
!      m = m + 1   
!    else
!      Cminus_avr = Cminus_avr + Cdyn(c)
!      n = n + 1
!    end if



    if(Cdyn(c) < 0.0) then
      Cdyn(c) = 0.0 
    else if(Cdyn(c) > 0.04) then
      Cdyn(c) = 0.04
    end if 

  end do

!  write(*,*) 'Max Cs = ', sqrt(MaxC), ' Min Cs = ', MinC/abs(MinC) * sqrt(abs(MinC))
!  write(*,*) 'Cplus_avr = ', sqrt(Cplus_avr/m), ' Cminus_avr = ', -1.0 * sqrt(abs(Cminus_avr/n)), m, n
  return
  END SUBROUTINE CalcSGS_Dynamic
