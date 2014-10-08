!======================================================================!
  SUBROUTINE FinPla
!----------------------------------------------------------------------!
!   Creates the file for postprocessing with Tecpot.                   !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: c1, c2, s, Ncp
  REAL      :: xc1,xc2,yc1,yc2,zc1,zc2,f1,f2,x0,y0,z0
  REAL      :: v_X, v_Y, v_U, v_V, v_W, v_P 
  CHARACTER :: namOut*80, namTem*80, dir
!--------------------------------[CVS]---------------------------------!
!  $Id: FinPla.f90,v 1.14 2008/12/10 14:36:03 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/FinPla.f90,v $  
!======================================================================!

!---- store the name
  namTem = name

1 if(THIS  < 2) then 
    write(*,*) '#-----------------------------------#'
    write(*,*) '#----- Find planes for Tecplot -----#'
    write(*,*) '#-----------------------------------#'
  end if

  if(THIS  < 2) &
    write(*,*) '# Direction (x,y or z) and coordinate, or skip'
  call ReadC(7,inp,tn,ts,te)
  call ToUppr(inp(ts(1):te(1)))
  if( inp(ts(1):te(1))  ==  'SKIP' ) then
!---- restore the name
    name = namTem
    return 
  else
    read(inp(ts(1):te(1)),'(A1)')  dir
    call ToUppr(dir)
    if(dir  ==  'X') read(inp(ts(2):te(2)),*) x0
    if(dir  ==  'Y') read(inp(ts(2):te(2)),*) y0
    if(dir  ==  'Z') read(inp(ts(2):te(2)),*) z0

    if(THIS  < 2) &
      write(*,*) '# Input the plane file name:'
    call ReadC(7,inp,tn,ts,te)
    read(inp(ts(1):te(1)), '(A80)')  namOut
    call NamFil(THIS, namOut, '.plt', len_trim('.plt'))
    open(9, FILE=namOut)
    write(*,*) '# Now creating the file:', namOut

!----- Count the values on the plane
    Ncp=0
    do s=1,NS

      c1=SideC(1,s)
      c2=SideC(2,s)

      xc1 = xc(c1)
      yc1 = yc(c1)
      zc1 = zc(c1)
      xc2 = xc(c1) + Dx(s)
      yc2 = yc(c1) + Dy(s)
      zc2 = zc(c1) + Dz(s)

      if(                                                           &
   ( (dir == 'X').and.                                              &
     ((xc1  > x0.and.xc2  < x0).or.(xc2  > x0.and.xc1  < x0)) )     &
     .or.                                                           &
   ( (dir == 'Y').and.                                              &
     ((yc1  > y0.and.yc2  < y0).or.(yc2  > y0.and.yc1  < y0)) )     &
     .or.                                                           &
   ( (dir == 'Z').and.                                              &
     ((zc1  > z0.and.zc2  < z0).or.(zc2  > z0.and.zc1  < z0)) )     &
	) Ncp = Ncp + 1

    end do

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!                                !
!      Write to Tecplot file     !
!                                !
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
    write(9,'(A13)')      'TITLE="TFlowS"'
    write(9,'(A34)')      'VARIABLES="X" "Y" "U" "V" "W" "P"'
    write(9,'(A7,I7,A8)') 'ZONE I=', Ncp, ' F=POINT'

    do s=1,NS

      c1=SideC(1,s)
      c2=SideC(2,s)

      xc1 = xc(c1)
      yc1 = yc(c1)
      zc1 = zc(c1)
      xc2 = xc(c1) + Dx(s)
      yc2 = yc(c1) + Dy(s)
      zc2 = zc(c1) + Dz(s)

      if(                                                           &
   ( (dir == 'X').and.                                              &
     ((xc1  > x0.and.xc2  < x0).or.(xc2  > x0.and.xc1  < x0)) )     &
     .or.                                                           &
   ( (dir == 'Y').and.                                              &
     ((yc1  > y0.and.yc2  < y0).or.(yc2  > y0.and.yc1  < y0)) )     &
     .or.                                                           &
   ( (dir == 'Z').and.                                              &
     ((zc1  > z0.and.zc2  < z0).or.(zc2  > z0.and.zc1  < z0)) )     &
	) then
	Ncp = Ncp + 1
	if(dir == 'X') f1 = abs(xc2-x0)/(abs(xc2-xc1))
	if(dir == 'Y') f1 = abs(yc2-y0)/(abs(yc2-yc1))
	if(dir == 'Z') f1 = abs(zc2-z0)/(abs(zc2-zc1))
	if(dir == 'X') f2 = abs(xc1-x0)/(abs(xc2-xc1))
	if(dir == 'Y') f2 = abs(yc1-y0)/(abs(yc2-yc1))
	if(dir == 'Z') f2 = abs(zc1-z0)/(abs(zc2-zc1))

	if(dir == 'X') v_X = f1*yc(c1) + f2*yc(c2)
	if(dir == 'Y') v_X = f1*xc(c1) + f2*xc(c2)
	if(dir == 'Z') v_X = f1*xc(c1) + f2*xc(c2)

	if(dir == 'X') v_Y = f1*zc(c1) + f2*zc(c2)
	if(dir == 'Y') v_Y = f1*zc(c1) + f2*zc(c2)
	if(dir == 'Z') v_Y = f1*yc(c1) + f2*yc(c2)

	v_U = f1 * U % n(c1) + f2 * U % n(c2)
	v_V = f1 * V % n(c1) + f2 * V % n(c2)
	v_W = f1 * W % n(c1) + f2 * W % n(c2)
	v_P = f1 * P % n(c1) + f2 * P % n(c2)

	write(9,*) v_X, v_Y, v_U, v_V, v_W, v_P 
      endif 

    end do    ! through sides 

    write(9,'(I8)') Ncp

    close(9)

  end if

  goto 1

  END SUBROUTINE FinPla
