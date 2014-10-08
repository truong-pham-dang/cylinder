!======================================================================!
  SUBROUTINE TraceSeparationPoint_new(now)
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
  INTEGER              :: s, c1, c2, i, sp
  REAL								 :: z_min, z_max, epsilon
  REAL,ALLOCATABLE     :: rho(:), Uphi(:), gradUphi_x(:), gradUphi_y(:)
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
  else

    if(this < 2) write(*, *) ' # File ', nameIN, ' wasn`t found in current directory'
    RETURN

  end if

	!if(this < 2) write(992,'(A2,7A16)') "  ","n","sp","x","y","Uphi","gradUphi_r","epsilon"

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
			end if
		end if
	end do ! end of browsing through all wall cells points

!--------BODY---------!
!---------------------!
!----- positive y ----!
!---------------------!

	!calculating grad U_{phi}_rho in attached to the boundary contition (wall) layer
  call Wait

  i = 0
  epsilon=0.001
  sp = 0

  do while (i /= 1)
	  call Wait
	  i = 0
		do s=1,NS ! through all surfaces
			c1=SideC(1,s) ! outer cell, attached to boundary cell
			c2=SideC(2,s) ! boundary cell
			if(c2 < 0) then ! one of boundary point
				if(TypeBC(c2) == WALL) then ! c2 = bc = wall 
		      if ( yc(c1) >= 0.00001 .and. zc(c1) >= z_min .and. zc(c1) <= z_max .and. abs(gradUphi_r(c1))<=epsilon) then !cutting the ring in I,II sectors !.and. abs(gradUphi_r) <= precision
						sp = c1      	
		      	i  = i+1
			    end if
				end if
			end if
		end do ! end of browsing through all wall cells points

		call Wait
		call IGlSum(i)

		if ( i > 1 ) then
			epsilon=epsilon/7
		elseif (i == 0) then
			epsilon=epsilon*2	
		end if
	end do


	!write(*,*) sp, xc(sp), gradUphi_r(sp), abs(gradUphi_r(SeparationPointCandidateIndex(1:i)))

  call Wait
	!cause if program made it's way here there is only one point among all processor is found
	call IGlSum(sp)

  !submiting first point 	
	if(this < 2) write(992,'(A2, 2I16, 5F16.6)') "1:", now, sp, xc(sp), yc(sp), Uphi(sp), gradUphi_r(sp), epsilon

!---------------------!
!----- negative y ----!
!---------------------!
  call Wait

  i = 0
  epsilon=0.001
  sp = 0

  do while (i /= 1)
	  call Wait
	  i = 0
		do s=1,NS ! through all surfaces
			c1=SideC(1,s) ! outer cell, attached to boundary cell
			c2=SideC(2,s) ! boundary cell
			if(c2 < 0) then ! one of boundary point
				if(TypeBC(c2) == WALL) then ! c2 = bc = wall 
		      if ( yc(c1) < 0.00001 .and. zc(c1) >= z_min .and. zc(c1) <= z_max .and. abs(gradUphi_r(c1))<=epsilon) then !cutting the ring in I,II sectors !.and. abs(gradUphi_r) <= precision
						sp = c1
		      	i  = i+1
			    end if
				end if
			end if
		end do ! end of browsing through all wall cells points

		call Wait
		!cause if program made it's way here there is only one point among all processor is found
		call IGlSum(i)

		if ( i > 1 ) then
			epsilon=epsilon/7
		elseif (i == 0) then
			epsilon=epsilon*2	
		end if
	end do

  call Wait

	call IGlSum(sp)
!submiting second point 	
	if(this < 2) write(992,'(A2, 2I16, 5F16.6)') "2:", now, sp, xc(sp), yc(sp), Uphi(sp), gradUphi_r(sp), epsilon

	close(992)

	!destructor
	deallocate(rho)
	deallocate(Uphi)
	deallocate(gradUphi_x)
	deallocate(gradUphi_y)
  call Wait
  
  RETURN
  END SUBROUTINE TraceSeparationPoint_new
