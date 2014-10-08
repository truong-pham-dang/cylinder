 SUBROUTINE SaveYplus_cyl(n0, n1) 
!----------------------------------------------------------------------!
!   						!print into file Y+ alongside the wall                 !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE les_mod
  USE pro_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Variables]------------------------------!
  INTEGER :: s, c1, c2, n
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: n0, n1
!-------------------------------[Locals]-------------------------------!
  REAL,ALLOCATABLE     :: rho(:),Uphi(:), Yplus(:) ! radial distance, Tangential velocity, Utau

	n=n1-n0

  if(n  > -1) then

    !initialization
    allocate(rho(-NbC:NC) );        rho              = 0.0 
    allocate(Uphi(-NbC:NC) );       Uphi             = 0.0 
    allocate(Yplus(-NbC:NC) );      Yplus            = 0.0 

    open(1235, FILE='Yplus.dat') ! file should consist of cells nodes

    if(this < 2 ) write(1235,'(3A16)') "x","y","Yplus"


    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)

      if(c2  < 0) then
        if( TypeBC(c2) == WALL .or. bcmark(c2) == 1 ) then !cylinders wall

            rho(c1)   = ( xc(c1) **2 + yc(c1) **2 ) **0.5
            ! vec{r}= (x,y)
            ! Urho = (vec{U}*vec{e}_{r})= (vec{U}*vec{r}) / abs(r)
            ! abs(Uphi) = (xV-yU) / abs(r)
            ! dUphi/dr=(<Uphi>-0)/(r-R_cyl), U=0 on wall
            ! abs(Utau) = sqrt(Visc*dUphi/dr) !!! Pope. pg 295
            ! y+ = Utau*(r-R_cyl)/Visc = sqrt(dUphi/dr)*(r-R_cyl)/sqrt(VISc)= &
            ! &  = (<Uphi>*(r-R_cyl)/VISc)^1/2
            Uphi(c1)     = abs( ( U % mean(c1) * yc(c1) - V % mean(c1) * xc(c1) ) / rho(c1)  )
            Yplus(c1) = sqrt( Uphi(c1) * (rho(c1) - 0.5) / VISc )
            !submiting results by writing them into file
            write(1235,'(3F16.6)') xc(c1), yc(c1), Yplus(c1)

        end if
      end if
    end do

    call Wait

    close(1235)

    deallocate(rho)
    deallocate(Uphi)
    deallocate(Yplus)

  end if !if n>1


END SUBROUTINE SaveYplus_cyl