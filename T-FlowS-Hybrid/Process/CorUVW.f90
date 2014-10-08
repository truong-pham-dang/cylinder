!======================================================================!
  REAL FUNCTION CorUVW()
!----------------------------------------------------------------------!
!   Corrects the velocities, and mass fluxes on the cell faces.        !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: c, c1, c2, s, m
  REAL      :: CFLmax(256), PeMax(256)
  REAL      :: CFLs, PeS
  REAL      :: Pdrop, FluxM
!--------------------------------[CVS]---------------------------------!
!  $Id: CorUVW.f90,v 1.23 2008/12/10 14:19:45 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CorUVW.f90,v $  
!======================================================================!

!--------------------------------------------!
!     Correct velocities and fluxes with     !
!      periodic part of the pressure to      !
!      obtain divergence free velocity       !
!- - - - - - - - - - - - - - - - - - - - - - !
!     For SOLIDs, Px, Py and Pz are zero     !
!     so this loop will not correct SOLID    !
!     velocities.                            !
!--------------------------------------------!
  if(ALGOR == FRACT) then
    do c=1,NC
      U % n(c) = U % n(c) - Px(c) * volume(c) / Asave(c)
      V % n(c) = V % n(c) - Py(c) * volume(c) / Asave(c)
      W % n(c) = W % n(c) - Pz(c) * volume(c) / Asave(c)
    end do 
  else ! Algorythm is SIMPLE
    do c=1,NC
      U % n(c) = U % n(c) - Px(c) * volume(c) / Asave(c)
      V % n(c) = V % n(c) - Py(c) * volume(c) / Asave(c)
      W % n(c) = W % n(c) - Pz(c) * volume(c) / Asave(c)
    end do 
  end if

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

    if(c2  < 0) then
      if( (TypeBC(c2) == PRESSURE) ) then
        U % n(c2) = U % n(c1) 
        V % n(c2) = V % n(c1) 
        W % n(c2) = W % n(c1) 
      end if
    end if
  end do 
!----------------------------------------------------------------------!
!     Look at the following equation and you will understand why       !
!     is the matrix for pressure corrections in SIMPLE algorythm       !
!     formed from the coefficients of the velocity matrix.             !
!     Moreover, it should also be clear that pressure correction       !
!     matrix must be formed from underrelaxed velocity coefficients    !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!     Note that for FLUID-SOLID interaction, FLUX correctin is zero    !
!     because Aval(SidAij(1,s)) is also zero.                          !  
!     What will happen with parallel version ... only god knows.       !
!----------------------------------------------------------------------!
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2  > 0 .or. c2  < 0.and.TypeBC(c2) == BUFFER) then
      if(c2  > 0) then
	Flux(s)=Flux(s)+(PP % n(c2) - PP % n(c1))*Aval(SidAij(1,s))
      else 
	Flux(s)=Flux(s)+(PP % n(c2) - PP % n(c1))*Abou(c2)
      endif
    end if             !                                          !
  end do               !<---------- this is correction ---------->!

!-----------------------------------------!
!      Calculate the max mass error       !
!     with the new (corrected) fluxes     !
!-----------------------------------------!
1 do c=1,NC
    b(c) = 0.0 
  end do

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2  > 0 .or. c2  < 0 .and. TypeBC(c2) == BUFFER) then
      b(c1)=b(c1)-Flux(s)
      if(c2  > 0) b(c2)=b(c2)+Flux(s)
    else
      b(c1) = b(c1)-Flux(s)
    end if
  end do

  do c=1,NC
    b(c) = b(c) / (volume(c)*DENc(material(c)))
  end do

  errmax=0.0
  do c=1,NC
    errmax=max(errmax, abs(b(c)))
  end do
  call glomax(errmax)

!----------------------------------!
!     Calculate the CFL number     !
!       and the Peclet number      !
!----------------------------------!
  do m=1,Nmat
    CFLmax(m) = 0.0
    PeMax(m)  = 0.0
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      if( (material(c1) .eq. m) .or. (material(c2) .eq. m) ) then
        if(c2  > 0 .or. c2  < 0.and.TypeBC(c2) == BUFFER) then
          CFLs = abs( dt * Flux(s) /               &
                      ( Scoef(s) *                 &
                      (Dx(s)*Dx(s)+Dy(s)*Dy(s)+Dz(s)*Dz(s)) ) )
          PeS  = abs( Flux(s) / Scoef(s) / (VISc+TINY) )
          CFLmax(m) = max( CFLmax(m), CFLs ) 
          PeMax(m)  = max( PeMax(m),  PeS  ) 
        end if
      end if
    end do
    call glomax(CFLmax(m))
    call glomax(PeMax(m))
  end do
!->>>      
  if(Cm(1) /= 0) then     
    write(LinMon0( 19: 66), '(1PE12.3,1PE12.3,1PE12.3,1PE12.3)')  &
          U % n(Cm(1)),  V % n(Cm(1)),  W % n(Cm(1)),  P % n(Cm(1))
    if(HOT==YES) then
      write(LinMon0( 67: 78), '(1PE12.3)') T % n(Cm(1)) 
    end if
    do m=1,Nmat
      if(m .eq. 1) then
        FluxM = max(abs(FLUXx(m)),  abs(FLUXy(m)),  abs(FLUXz(m)))
        Pdrop = max(abs(PdropX(m)), abs(PdropY(m)), abs(PdropZ(m)))
        write(LinMon1( 79: 90), '(1PE12.3)') FluxM 
        write(LinMon1( 91:102), '(1PE12.3)') Pdrop    
        write(LinMon1(103:126), '(1PE12.3,1PE12.3)')  &
  	      CFLmax(m), PeMax(m)
      else if(m .eq. 2) then
        FluxM = max(abs(FLUXx(m)),  abs(FLUXy(m)),  abs(FLUXz(m)))
        Pdrop = max(abs(PdropX(m)), abs(PdropY(m)), abs(PdropZ(m)))
        write(LinMon2( 79: 90), '(1PE12.3)') FluxM 
        write(LinMon2( 91:102), '(1PE12.3)') Pdrop    
        write(LinMon2(103:126), '(1PE12.3,1PE12.3)')  &
  	      CFLmax(m), PeMax(m)
      else if(m .eq. 3) then
        FluxM = max(abs(FLUXx(m)),  abs(FLUXy(m)),  abs(FLUXz(m)))
        Pdrop = max(abs(PdropX(m)), abs(PdropY(m)), abs(PdropZ(m)))
        write(LinMon3( 79: 90), '(1PE12.3)') FluxM 
        write(LinMon3( 91:102), '(1PE12.3)') Pdrop    
        write(LinMon3(103:126), '(1PE12.3,1PE12.3)')  &
  	      CFLmax(m), PeMax(m)
      end if
    end do
    write(LineRes(  5: 16), '(1PE12.3)') errmax 
  end if 

  CorUVW = errmax ! /(velmax+TINY)

  END FUNCTION CorUVW 
