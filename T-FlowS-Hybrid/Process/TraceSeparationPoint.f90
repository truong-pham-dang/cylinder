!======================================================================!
  SUBROUTINE TraceSeparationPoint(now)
!----------------------------------------------------------------------!
! Tracking the separation points
! on surface of cylinder in z_min - z_max range
! designed for cylinder task
! should be called in Processor.f90 in writing monitoring points block
! default name for monitoring point file was chosen as separation-points
! file must begin with z_min , z_max values on first row
!----------------------------------------------------------------------!
  USE all_mod
  USE allp_mod
  USE les_mod
  USE pro_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!  
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
	INTEGER              :: now
!------------------------------[Calling]-------------------------------!
  INTEGER              :: s, c1, c2, i, j, sp1, sp2
  REAL								 :: z_min, z_max
  REAL,ALLOCATABLE     :: rho(:), Uphi(:), gradUphi_x(:), gradUphi_y(:), &
  												LongitudinalDeviation(:)
  INTEGER, ALLOCATABLE ::	SeparationPointCandidateIndex(:)
  LOGICAL              :: FileExists
  CHARACTER            :: nameIN*80
  !--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $  
!======================================================================!
  !if(this < 2) write(*, *) ' Variable declared'

	nameIN=trim('separation-points')

	inquire(FILE=nameIN, EXIST=FileExists)
  if (FileExists) then  ! FileExists will be TRUE if the file
    !if(this < 2) write(*, *) ' # Now reading the file: ', nameIN 
    open(992, FILE=nameIN) ! file should consist of cells nodes
		read(992,*) z_min, z_max
		close(992)
		open(992, file=nameIN,POSITION='APPEND') !open file at the end of the file

		!initialization
		allocate(gradUphi_x(-NbC:NC) ); gradUphi_x = 0.0 
		allocate(gradUphi_y(-NbC:NC) ); gradUphi_y = 0.0 
		allocate(rho(-NbC:NC) );        rho        = 0.0 
		allocate(Uphi(-NbC:NC) );       Uphi       = 0.0 
		allocate(SeparationPointCandidateIndex(1:NbC)); SeparationPointCandidateIndex = 0.0
		i = 0
		j = 0
		sp1 = 0
		sp2 = 0

  else

    if(this < 2) write(*, *) ' # File ', nameIN, ' wasn`t found in current directory'
    RETURN

  end if

  call Wait

	do s=1,NS ! through all surfaces
		c1=SideC(1,s) ! outer cell, attached to boundary surface
		c2=SideC(2,s) ! inner cell, attached to boundary surface
		if(c2 < 0) then ! one of boundary point
			if(TypeBC(c2) == WALL) then ! bc = wall 
				! transition to cylindrical curvelinear coordinates rho and Phi
	      rho(c1) = ( xc(c1) **2 + yc(c1) **2 ) **0.5 + tiny

	      ! Urho = (vec{U}*vec{e}_{r})= (vec{U}*vec{r}) / abs(r)
	      ! abs(Uphi) = (xV-yU) / abs(r)
	      ! vec{r}= (x,y)
				Uphi(c1) = ( U % n(c1) * yc(c1) - V % n(c1) * xc(c1) ) / rho(c1)

				end if
		end if
	end do

!---- instantaneous derivatives
  call GraPhi(Uphi,1,gradUphi_x, .TRUE.) ! grad(U_{r})_x
  call Wait

  call GraPhi(Uphi,2,gradUphi_y, .TRUE.) ! grad(U_{r})_y
	call Wait

	!calculating grad U_{phi}_rho in attached to the boundary contition (wall) layer
	do s=1,NS ! through all surfaces
		c1=SideC(1,s) ! outer cell, attached to boundary cell
		c2=SideC(2,s) ! boundary cell
		if(c2 < 0) then ! one of boundary point
			if(TypeBC(c2) == WALL) then ! c2 = bc = wall 

	      !grad U_{phi}_rho
	      gradUphi_r(c1) = rho(c1) * ( gradUphi_x(c1) / xc(c1) + gradUphi_y(c1) / yc(c1) )
			else
				gradUphi_r(c1) = HUGE;
			end if
		end if
	end do ! end of browsing through all wall cells points

!--------BODY---------!
!---------------------!
!----- positive y ----!
!---------------------!

	!calculating grad U_{phi}_rho in attached to the boundary contition (wall) layer
  i = 0
  call Wait

	do s=1,NS ! through all surfaces
		c1=SideC(1,s) ! outer cell, attached to boundary cell
		c2=SideC(2,s) ! boundary cell
		if(c2 < 0) then ! one of boundary point
			if(TypeBC(c2) == WALL) then ! c2 = bc = wall 
	      if ( yc(c1) >= 0 .and. zc(c1) >= z_min .and. zc(c1) <= z_max) then !cutting the ring in I,II sectors !.and. abs(gradUphi_r) <= precision
	      	i = i + 1
	      	!SeparationPointCandidateIndex is the index of browsed cells in the ring
	      	SeparationPointCandidateIndex(i) = c1    ! x, y, Uphi, gradUphi_r - all that is insided in c1
		    end if

			end if
		end if
	end do ! end of browsing through all wall cells points

  call Wait
	sp1 = minloc(abs(gradUphi_r(SeparationPointCandidateIndex(1:i))),i) !return index of lowest |grad U_{phi}_rho| - separation point

	!if (sp1 == 0) ! there weren't found any points()

	write(*,*) sp1, xc(sp1), gradUphi_r(sp1), abs(gradUphi_r(SeparationPointCandidateIndex(1:i)))

  !submiting first point 	
	if(this < 2) write(992,'(A4, A2, I7, A4, I7, A2, F16.6, A2, F16.6, A5, F16.6, A11, F16.6)') "1:  ", "n=", now, &
		"sp1=",sp1,"x=", xc(sp1), "y=", yc(sp1),"Uphi=", Uphi(sp1), "gradUphi_r=", gradUphi_r(sp1)

	call Wait

!---------------------!
!----- negative y ----!
!---------------------!
  i = 0
  call Wait

	do s=1,NS ! through all surfaces
		c1=SideC(1,s) ! outer cell, attached to boundary surface
		c2=SideC(2,s) ! inner cell, attached to boundary surface
		if(c2 < 0) then ! one of boundary point
			if(TypeBC(c2) == WALL .and. yc(c1) < 0 &
			 .and. zc(c1) >= z_min .and. zc(c1) <= z_max) then ! bc = wall 
	      	i = i + 1
	      	!SeparationPointCandidateIndex is quite tricky: follow the loop to clarify it
	      	SeparationPointCandidateIndex(i) = c1    ! x, y, Uphi, gradUphi_r - all that is insided in c1
			end if
		end if
	end do ! end of browsing through all wall cells points
  
	call Wait
  
	if (this < 2) then 
		sp2 = minloc(abs(gradUphi_r(SeparationPointCandidateIndex(1:i))),i) !return index of lowest |grad U_{phi}_rho| - separation point
	end if
  
  !submiting second point
	if(this < 2) write(992,'(A4, A2, I7, A4, I7, A2, F16.6, A2, F16.6, A5, F16.6, A11, F16.6)') "1:  ", "n=", now, &
		"sp2=",sp2,"x=", xc(sp2), "y=", yc(sp2),"Uphi=", Uphi(sp2), "gradUphi_r=", gradUphi_r(sp2)

  call Wait

	close(992)

	!destructor
	deallocate(rho)
	deallocate(Uphi)
	deallocate(gradUphi_x)
	deallocate(gradUphi_y)
	deallocate(SeparationPointCandidateIndex)
	deallocate(LongitudinalDeviation)

  RETURN
  END SUBROUTINE TraceSeparationPoint
