!======================================================================!
  SUBROUTINE ModOut 
!----------------------------------------------------------------------!
!   Modifies the fluxes at outflow boundaries to conserve the mass.    ! 
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: m, s, c1, c2
  REAL    :: fac(256) 
!--------------------------------[CVS]---------------------------------!
!  $Id: ModOut.f90,v 1.14 2008/12/10 14:46:52 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/ModOut.f90,v $  
!======================================================================!

!------------------------------------------!
!     Calculate the inflow mass fluxes     !
!------------------------------------------!
  do m=1,Nmat
    MassIn(m) = 0.0
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      if(c2 < 0) then
        Flux(s) = DENc(material(c1))*( U % n(c2)*Sx(s) + &
                                       V % n(c2)*Sy(s) + &
                                       W % n(c2)*Sz(s) )
        if(TypeBC(c2)  ==  INFLOW) then
	  if(material(c1) == m) MassIn(m) = MassIn(m) - Flux(s)
        endif
        if(TypeBC(c2)  ==  PRESSURE.and.Flux(s) < 0.0) then
	  if(material(c1) == m) MassIn(m) = MassIn(m) - Flux(s)
        end if
        if(TypeBC(c2)  ==  CONVECT.and.Flux(s) < 0.0) then
	  if(material(c1) == m) MassIn(m) = MassIn(m) - Flux(s)
        end if
      end if
    end do
    call glosum(MASSIN(m))
  end do                    

!-------------------------------------------!
!     Calculate the outflow mass fluxes     !
!       then correct it to satisfy the      ! 
!            overall mass balance           !
!-------------------------------------------!
  do m=1,Nmat
    MASOUT(m) = 0.0
    do s=1, NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      if(c2  < 0) then
        if(TypeBC(c2) == OUTFLOW) then
  	  U % n(c2) = U % n(c1)
	  V % n(c2) = V % n(c1)
	  W % n(c2) = W % n(c1)
	  Flux(s) = DENc(material(c1)) * ( U % n(c2)*Sx(s) +  & 
                                           V % n(c2)*Sy(s) +  &
                                           W % n(c2)*Sz(s) )
	  if(material(c1) == m) MASOUT(m) = MASOUT(m) + Flux(s)
        endif
        if(TypeBC(c2) == CONVECT.and.Flux(s) > 0.0) then
	  Flux(s) = DENc(material(c1)) * ( U % n(c2)*Sx(s) +  & 
                                           V % n(c2)*Sy(s) +  &
                                           W % n(c2)*Sz(s) )
	  if(material(c1) == m) MASOUT(m) = MASOUT(m) + Flux(s)
        endif

	Flux(s) = DENc(material(c1)) * ( U % n(c2)*Sx(s) +  & 
                                         V % n(c2)*Sy(s) +  &
                                         W % n(c2)*Sz(s) )
        if(TypeBC(c2) == PRESSURE.and.Flux(s) > 0.0) then
	  if(material(c1) == m) MASOUT(m) = MASOUT(m) + Flux(s)
        endif

      endif
    end do
    call glosum(MASOUT(m))  ! not checked
  end do

  do m=1,Nmat
    fac(m) = MASSIN(m)/(MASOUT(m)+TINY)
  end do
  do m=1,Nmat
    do s=1, NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      if(c2  < 0) then
        if(TypeBC(c2) == OUTFLOW  .or. TypeBC(c2) == CONVECT &
          .or. TypeBC(c2) == PRESSURE) then
	  if(material(c1) == m) then
  	    U % n(c2) = U % n(c2) * fac(m)
	    V % n(c2) = V % n(c2) * fac(m) 
	    W % n(c2) = W % n(c2) * fac(m)
	    Flux(s) = Flux(s)*fac(m) 
	    MASOUT(m) = MASOUT(m) + Flux(s)
	  endif
        endif
      endif
    end do
  end do

  END SUBROUTINE ModOut
