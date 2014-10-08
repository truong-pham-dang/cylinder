!======================================================================!
  SUBROUTINE SourceKinKEPSV2F()
!----------------------------------------------------------------------!

!------------------------------[Modules]-------------------------------!
!----------------------------------------------------------------------!
!   Computes the source terms in Kin transport equation,               !
!   wall shear stress (wall function approuch)                         !
!                                                                      !
!   model : k-eps-v2f                                                  !
!                                                                      !
!   Authors: Muhamed Hadziabdic and Bojan Niceno                       !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, c1, c2, s
  REAL    :: Utan, UnorSq, Unor, UtotSq, dely, Stot
  REAL    :: lf, Gblend, Ustar, Ck, yPlus, Uplus
  REAL    :: ALPHA1, Lrans, Lsgs, koeff, Lrans_integral
!--------------------------------[CVS]---------------------------------!
!  $Id: SourceKinKEPSV2F.f90,v 1.4 2008/11/19 14:59:29 IUS\mhadziabdic Exp $  
!  $UserKDSource: /home/muhamed/.CVSROOT/T-Rex/Pro/KSource.f90,v $  
!======================================================================! 

!===================!
!                   !
!     Pk terms      !
!     in domain     !
!                   !
!===================!

!----------------------------------------------------------------------!
!                                                                      !
!   Pk  is a production of turbulent kinematic energy.                 !
!                                                                      !
!   The form of Pk which is solving in this subroutine is :            !
!                                                                      !
!                                                                      !
!    Pk = 2*VISk{ [ (dU/dx)**2 + (dV/dy)**2 + (dW/dz)**2 ] +           ! 
!        0.5( dV/dx + dU/dy )**2 + 0.5( dW/dy + dV/dz )**2 +           !
!        0.5(dU/dz + dW/dx )**2 }                                      !
!                                                                      !
!   Dimension of the Pk is [ kg m /s**3 ]                              !
!   In kinetic energy equation exist two source terms which have form: !
!                                                                      !
!    /             /                                                   !
!   |             |                                                    !
!   |Pk dV  -     |DENc Eps dV                                         !
!   |             |                                                    !
!  /             /                                                     !
!                                                                      !
!----------------------------------------------------------------------!

  koeff=1.0

  if(SIMULA == HYB_ZETA) then
    do c=1,NC
      if ( xc(c)< -0.5 ) then
        lf = volume(c)**0.3333
        Lsgs  = 0.8*lf
        Lrans = 0.41*WallDs(c)
        Lrans_integral=(Kin % n(c)**1.5)/Eps % n(c)
        if ( min(Lrans/Lsgs,Lrans_integral/Lsgs/koeff)>1.0 ) then
          koeff=max(Lrans_integral/Lsgs,1.0)
        end if
      end if 
    end do

    write(*,*) 'Lrans_coef in SourceKinKEPSV2F=', koeff

    do c=1,NC
      lf = volume(c)**0.3333
      Lsgs  = 0.8*lf
      Lrans = 0.41*WallDs(c)
      Lrans_integral=(Kin % n(c)**1.5)/Eps % n(c)

      ALPHA1 = max(1.0,min(Lrans/Lsgs,Lrans_integral/Lsgs/koeff)) 
!----- Production:
      b(c) = b(c) + VISt(c) * Shear(c) * Shear(c) * volume(c)
!----------------------------------------!
!----- Dissipation:
      if(ALPHA1 < 1.05) then
        Aval(Adia(c)) = Aval(Adia(c)) +                                    &
                        DENc(material(c))*Eps%n(c)/(Kin%n(c) + tiny)*volume(c)
      else
        Aval(Adia(c)) = Aval(Adia(c)) +                                    &
                        DENc(material(c))*min(ALPHA**1.45*Eps%n(c),Kin%n(c)**1.5  &
                       /(lf*0.01))/(Kin%n(c) + tiny)*volume(c)
      end if
      Pk(c) =  VISt(c) * Shear(c) * Shear(c)
    end do
  else
    do c=1,NC
!----- Production:
      b(c) = b(c) + VISt(c) * Shear(c) * Shear(c) * volume(c)
 
!----------------------------------------!
!----- Dissipation:
      Aval(Adia(c)) = Aval(Adia(c)) +                                    &
                      DENc(material(c))*Eps%n(c)/(Kin%n(c)+tiny)*volume(c)
      Pk(c) =  VISt(c) * Shear(c) * Shear(c) 
    end do
  end if

    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
    
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
!----- Compute tangential velocity component
          UtotSq = U % n(c1) * U % n(c1) &
                 + V % n(c1) * V % n(c1) &
                 + W % n(c1) * W % n(c1)
          Unor = ( U % n(c1) * Sx(s)     &
                 + V % n(c1) * Sy(s)     &
                 + W % n(c1) * Sz(s) )   &
                 /(Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s))**0.5
          UnorSq = Unor*Unor

          if( UtotSq   >  UnorSq) then
            Utan = sqrt(UtotSq - UnorSq)
          else
            Utan = TINY
          end if
         
          Uf(c1)  = Cmu25*Kin%n(c1)**0.5
          Ynd(c1) = WallDs(c1)*Uf(c1)/VISc 

          

          TauWall(c1) =  abs(DENc(material(c1)) * kappa * Uf(c1) * Utan    &
                    / (log(Elog*Ynd(c1))))    

          Ck = sqrt(TauWall(c1))
          Ynd(c1) = WallDs(c1)*Ck/VISc

          Gblend  = 0.01*Ynd(c1)**4.0/(1.0+Ynd(c1))

            if(WallDs(c1)*Uf(c1)/VISc > 5.0) then   
              Pk(c1) = VISt(c1)*Shear(c1)*Shear(c1)*exp(-Gblend) + &
                       TauWall(c1)*0.07**0.25*Kin%n(c1)**0.5/(kappa*WallDs(c1))*exp(-1.0/Gblend)
            else
              Pk(c1) = VISt(c1)*Shear(c1)*Shear(c1)           
            end if         
          b(c1) = b(c1) + Pk(c1) * volume(c1)
          b(c1) = b(c1) - VISt(c1) * Shear(c1) * Shear(c1) * volume(c1)
        end if  ! TypeBC(c2)==WALL or WALLFL
      end if    ! c2 < 0
    end do

  RETURN
 
  END SUBROUTINE SourceKinKEPSV2F 
