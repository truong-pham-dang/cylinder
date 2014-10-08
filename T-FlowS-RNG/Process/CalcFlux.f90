!======================================================================!
  SUBROUTINE CalcFlux()
!----------------------------------------------------------------------!
!   Calculate mass fluxes through whole domain.                        !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c1, c2, s, m
  REAL    :: xc1, yc1, zc1, xc2, yc2, zc2
!--------------------------------[CVS]---------------------------------!
!  $Id: CalcFlux.f90,v 1.5 2008/12/05 13:07:15 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CalcFlux.f90,v $  
!======================================================================!

  do m=1,Nmat
    FLUXx(m) = 0.0
    FLUXy(m) = 0.0
    FLUXz(m) = 0.0
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

          if((xc1 <= xp(m)).and.(xc2 > xp(m))) FLUXx(m) = FLUXx(m) + Flux(s)
          if((yc1 <= yp(m)).and.(yc2 > yp(m))) FLUXy(m) = FLUXy(m) + Flux(s)
          if((zc1 <= zp(m)).and.(zc2 > zp(m))) FLUXz(m) = FLUXz(m) + Flux(s)

          if((xc2 < xp(m)).and.(xc1 >= xp(m))) FLUXx(m) = FLUXx(m) - Flux(s)
          if((yc2 < yp(m)).and.(yc1 >= yp(m))) FLUXy(m) = FLUXy(m) - Flux(s)
          if((zc2 < zp(m)).and.(zc1 >= zp(m))) FLUXz(m) = FLUXz(m) - Flux(s)
        end if ! material 1&2
      else if(c2 < 0.and.TypeBC(c2) == BUFFER) then
        if( (material(c1)==m) .and. (material(c1) == material(c2)) ) then
          xc1=xc(c1) 
          yc1=yc(c1) 
          zc1=zc(c1) 
          xc2=xc(c1)+Dx(s) 
          yc2=yc(c1)+Dy(s) 
          zc2=zc(c1)+Dz(s)

          if((xc1 <= xp(m)).and.(xc2 > xp(m))) FLUXx(m) = FLUXx(m) + .5*Flux(s)
          if((yc1 <= yp(m)).and.(yc2 > yp(m))) FLUXy(m) = FLUXy(m) + .5*Flux(s)
          if((zc1 <= zp(m)).and.(zc2 > zp(m))) FLUXz(m) = FLUXz(m) + .5*Flux(s)

          if((xc2 < xp(m)).and.(xc1 >= xp(m))) FLUXx(m) = FLUXx(m) - .5*Flux(s)
          if((yc2 < yp(m)).and.(yc1 >= yp(m))) FLUXy(m) = FLUXy(m) - .5*Flux(s)
          if((zc2 < zp(m)).and.(zc1 >= zp(m))) FLUXz(m) = FLUXz(m) - .5*Flux(s)
        end if ! material 1&2
      end if   ! c2 > 0
    end do
    call glosum(FLUXx(m))
    call glosum(FLUXy(m))
    call glosum(FLUXz(m))
    Ubulk(m) = FLUXx(m) / (AreaX(m) + TINY)
    Vbulk(m) = FLUXy(m) / (AreaY(m) + TINY)
    Wbulk(m) = FLUXz(m) / (AreaZ(m) + TINY)
  end do ! m

  END SUBROUTINE CalcFlux
