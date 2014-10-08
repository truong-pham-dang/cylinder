!======================================================================!
  SUBROUTINE SourceKinKeps
!----------------------------------------------------------------------!
!   Computes the source terms in Kin transport equation.               !
!                                                                      !
!   model : k-eps                                                      !
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
  REAL    :: UtotSq, Unor, UnorSq, Utan, R
!--------------------------------[CVS]---------------------------------!
! '$Id: SourceKinKEps.f90,v 1.4 2009/06/30 12:13:57 IUS\mhadziabdic Exp $'
! '$Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/SourceKinKEps.f90,v $'
!======================================================================!


!--------------------------------------------------!
!    Compute the sources in the near wall cells    !
!--------------------------------------------------!
  if(MODE == HRe) then
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
!----- Compute nondimensional wall distance and wall-shear stress
          Ynd(c1) = sqrt(Kin % n(c1)) * Cmu25 * WallDs(c1)/VISc
          TauWall(c1) = abs(DENc(material(c1)) * kappa * sqrt(Kin % n(c1)) * Cmu25 * Utan &
                       / (log(Elog*Ynd(c1))))  
!----- Compute production in the first wall cell 
          Pk(c1) = TauWall(c1) * Cmu25 * sqrt(Kin % n(c1)) &
             / (kappa*WallDs(c1))
!----- Filling up the source term
          b(c1) = b(c1) + Pk(c1) * volume(c1) 
           
        end if  ! TypeBC(c2)==WALL or WALLFL
      end if    ! c2 < 0
    end do
!-------------------------------------------!
!    Compute the sources in the interior    !
!-------------------------------------------!
    do c=1,NC
!----- IsNearWall ensures not to put Pk twice into the near wall cells
      if(.NOT. IsNearWall(c)) then
!----- Production:
        Pk(c)= VISt(c) * Shear(c) * Shear(c)
        b(c) = b(c) + Pk(c) * volume(c)
      end if
!----- Dissipation:
      Aval(Adia(c)) = Aval(Adia(c)) + DENc(material(c))*Eps%n(c)/Kin%n(c)*volume(c)
    end do
  end if    ! end if mode = wf 

!=====================================================================
!     Jones-Launder model AND Launder - Sharma + Yap model
!=====================================================================
  if(MODE == LRe) then
    do c=1,NC
!----- Production:
      Pk(c)= VISt(c) * Shear(c) * Shear(c)
      b(c) = b(c) + Pk(c) * volume(c)
!----- Dissipation:
      Aval(Adia(c)) = Aval(Adia(c)) + DENc(material(c))*Eps%n(c)/(Kin%n(c)+tiny)*volume(c)
!----- Preparation of Kin for the boundary condition. Kin variable is temporaraly borrowed.
      Kin % n(c) = sqrt(Kin % n(c))
    end do
    call GraPhi(Kin % n,1,PHIx,.TRUE.)             ! dK/dx
    call GraPhi(Kin % n,2,PHIy,.TRUE.)             ! dK/dy
    call GraPhi(Kin % n,3,PHIz,.TRUE.)             ! dK/dz

    do c = 1, NC
!----- Turning Kin back to its real value
      Kin % n(c)    = Kin % n(c)*Kin % n(c) 
      Aval(Adia(c)) = Aval(Adia(c)) + 2.0*VISc*(PHIx(c)*PHIx(c) + &
      PHIy(c)*PHIy(c) + PHIz(c)*PHIz(c))* volume(c) / (Kin % n(c)+tiny)          
    end do
  end if

  if(SIMULA == HYB_PITM) then
    do c=1,NC
!----- Production:
      Pk(c)= VISt(c) * Shear(c) * Shear(c)
      b(c) = b(c) + Pk(c) * volume(c)
!----- Dissipation:
      Aval(Adia(c)) = Aval(Adia(c)) + DENc(material(c))*Eps%n(c)/Kin%n(c)*volume(c)
    end do
  end if

  call Exchng(Kin % n)
  RETURN
  END SUBROUTINE SourceKinKeps
