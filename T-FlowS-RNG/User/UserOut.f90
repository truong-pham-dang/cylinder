!======================================================================!
  SUBROUTINE UserOut(step)
!----------------------------------------------------------------------!
!   Calculate mass fluxes through whole domain.                        !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE par_mod
  USE les_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c1, c2, s, m, step
  REAL    :: xc1, yc1, zc1, xc2, yc2, zc2
  REAL    :: Stot
  REAL    :: AreaCross, AreaLowWall, AreaHighWall
  REAL    :: Tbulk, TlowWall, ThighWall
!--------------------------------[CVS]---------------------------------!
!  $Id: UserOut.f90,v 1.2 2002/10/30 16:30:03 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserOut.f90,v $  
!======================================================================!

!==============================!
!                              !
!   Compute bulk temperature   !
!                              !
!==============================!
    Tbulk     = 0.0
    AreaCross = 0.0
    m=1      
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      if(c2 > 0) then
        if( (material(c1)==m) .and. (material(c1) == material(c2)) ) then
          xc1=xc(c1) 
          yc1=yc(c1) 
          zc1=zc(c1) 
          xc2=xc(c1)+Dx(s) 
          yc2=yc(c1)+Dy(s) 
          zc2=zc(c1)+Dz(s)

          if((xc1 <= xp(m)) .and. (xc2 > xp(m))) then
            Tbulk = Tbulk + U % n(c1) * T % n(c1) * abs(Sx(s))
            AreaCross = AreaCross + abs(Sx(s))
          end if
 
          if((xc2 < xp(m)) .and. (xc1 >= xp(m))) then
            Tbulk = Tbulk + U % n(c2) * T % n(c2) * abs(Sx(s))
            AreaCross = AreaCross + abs(Sx(s))
          end if
        end if ! material 1&2
      else if(c2 < 0 .and. TypeBC(c2) == BUFFER) then
        if( (material(c1)==m) .and. (material(c1) == material(c2)) ) then
          xc1=xc(c1) 
          yc1=yc(c1) 
          zc1=zc(c1) 
          xc2=xc(c1)+Dx(s) 
          yc2=yc(c1)+Dy(s) 
          zc2=zc(c1)+Dz(s)

          if((xc1 <= xp(m)).and.(xc2 > xp(m))) then
            Tbulk = Tbulk + U % n(c1) * T % n(c1) * abs(Sx(s)) * 0.5
            AreaCross = AreaCross + abs(Sx(s))
          end if
 
          if((xc2 < xp(m)).and.(xc1 >= xp(m))) then
            Tbulk = Tbulk + U % n(c2) * T % n(c2) * abs(Sx(s)) * 0.5
            AreaCross = AreaCross + abs(Sx(s))
          end if
        end if ! material 1&2
      end if   ! c2 > 0
    end do
    call glosum(Tbulk)
    call glosum(AreaCross)

    Tbulk = Tbulk / (AreaCross * Ubulk(m))

!==============================!
!                              !
!   Compute wall temperature   !
!                              !
!==============================!
    AreaLowWall = 0.0
    AreaHighWall = 0.0
    TlowWall = 0.0
    ThighWall = 0.0
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
        if( (material(c1)==m) .and. (material(c1) == material(c2)) ) then
          Stot = sqrt(Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s))
          if(zc(c2) < 1.0) then
            TlowWall = TlowWall + T % n(c2) * Stot
            AreaLowWall = AreaLowWall + Stot         
          else
            ThighWall = ThighWall + T % n(c2) * Stot
            AreaHighWall = AreaHighWall + Stot      
          end if
        end if   
      end if   
    end do

    call glosum(AreaLowWall)
    call glosum(AreaHighWall)
    call glosum(TlowWall)
    call glosum(ThighWall)

    TlowWall = TlowWall / AreaLowWall
    ThighWall = ThighWall / AreaHighWall

    if(THIS < 2) &
      write(*,'(1X,I6,6F14.6,A7)') step, &
                                   Tbulk, TlowWall, ThighWall,  &
                                   AreaCross, AreaLowWall, AreaHighWall, &
                                   ' <oOOo>' 

  END SUBROUTINE UserOut
