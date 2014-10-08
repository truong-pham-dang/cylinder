!======================================================================!
  SUBROUTINE Calc4()
!----------------------------------------------------------------------!
!   Calculates geometrical quantities of the grid.                     !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod 
  USE gen_mod 
!----------------------------------------------------------------------!
  IMPLICIT NONE
!----------------------------------------------------------------------!
  INTERFACE
    LOGICAL FUNCTION Approx(A,B,tol)
      IMPLICIT NONE
      REAL          :: A,B
      REAL,OPTIONAL :: tol
    END FUNCTION Approx
  END INTERFACE
!------------------------------[Calling]-------------------------------!
  REAL    :: Dist       
  REAL    :: Dist2       
  REAL    :: Tol
!-------------------------------[Locals]-------------------------------!
  INTEGER                      :: c, c1, c2, n, s, ss, cc2, c_max, NNN, hh, mm
  INTEGER                      :: c11, c12, c21, c22, s1, s2,       &
                                  BouCen,                           &
                                  TypePer, Nper, NumberSides, dir, OPT
  INTEGER                      :: WallMark, rot_dir, dir_face
  REAL                         :: xt(4), yt(4), zt(4), angle_face
  REAL                         :: xs2,ys2,zs2, x_a, y_a, z_a, x_b, y_b, z_b
  REAL                         :: x_c, y_c, z_c, Det
  REAL                         :: ABi, ABj, ABk, ACi, ACj, ACk, Pi, Pj, Pk
  REAL                         :: dsc1, dsc2, PerMin, PerMax
  REAL                         :: t, SurTot, angle 
  REAL                         :: xc1, yc1, zc1, xc2, yc2, zc2 
  REAL                         :: MaxDis, TotVol, MinVol, MaxVol, Max_Coor, Min_Coor
  REAL                         :: xmin, xmax, ymin, ymax, zmin, zmax 
  DOUBLE PRECISION,ALLOCATABLE :: xnr(:), ynr(:), znr(:), xspr(:), yspr(:), zspr(:)
  DOUBLE PRECISION,ALLOCATABLE :: BCoor(:), phi_face(:)
  INTEGER,ALLOCATABLE          :: BFace(:)
!--------------------------------[CVS]---------------------------------!
  character*80 rcs1,rcs2
  data rcs1/                                                        &
  '$Id: Calc4.f90,v 1.28 2005/01/25 12:46:22 muhamed Exp $'/
  data rcs2/                                                        &
  '$Source: /home/muhamed/.CVSROOT/T-Rex/Neu2Trex/Calc4.f90,v $'/
!======================================================================!
!
!                                n3 
!                 +---------------!---------------+
!                /|              /|              /|
!               / |             / |             / |
!              /  |          n2/  |            /  |
!             +---------------!---------------+   |
!             |   |           |   |           |   |
!             |   |     o---- | s-------o     |   |      
!             |   +---- c1 ---|   !---- c2 ---|   +       
!             |  /            |  /n4          |  /
!             | /             | /             | /
!             |/              |/              |/
!             +---------------!---------------+
!                            n1
!
!   Notes:
! 
!     ! side s is oriented from cell center c1 to cell center c2     
!     ! c2 is greater then c1 inside the domain or smaller then 0
!       on the boundary
!     ! nodes are denoted with n1 - n4
!
!            c3           
!             \  4-----3
!              \/ \  . |
!              /   \  +---c2
!             /  .  \  |
!            / .     \ |
!           /.        \|
!          1-----------2
!                   |
!                   c1
!
!                                n3 
!                 +---------------!-------+         
!                /|            n2/|      /|
!               / |             !-------+ |
!              /  |            /|s|  c2 | |
!             +---------------+ | !-----| +
!             |   |           | |/n4    |/
!             |   |     c1    | !-------+
!             |   +-----------|n1 +
!             |  /            |  /
!             | /             | /
!             |/              |/ 
!             +---------------+
!                            n1
! 
!----------------------------------------------------------------------!
!   Generaly:
!
!   the equation of plane reads: A*x + B*y + C*z + D = 0
!
!   and the equation of line:  x = x0 + t*rx
!                              y = y0 + t*ry
!                              z = z0 + t*rz
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!   In our case:
!
!     line is a connection between the two cell centers:
!
!     x = xc(c1) + t*(xc(c2)-xc(c1)) = xc(c1) + t*rx
!     y = yc(c1) + t*(yc(c2)-yc(c1)) = yc(c1) + t*ry
!     z = zc(c1) + t*(zc(c2)-zc(c1)) = zc(c1) + t*rz
!    
!
!     plane is a cell face: 
!
!     Sx * x + Sy * y + Sz * z = Sx * xsp(s) + Sy * ysp(s) + Sz * zsp(s)
!  
!     and the intersection is at:
!  
!         Sx*(xsp(s)-xc(c1)) + Sy*(ysp(s)-yc(c1) + Sz*(zsp(s)-zc(c1)) 
!     t = -----------------------------------------------------------
!                           rx*Sx + ry*Sy + rz*Sz
!  
!----------------------------------------------------------------------!
 
!++++++++++++++++++++++++++++++++++++!
!     Calculate the cell centers     !
!------------------------------------!
!     => depends on: x,y,z           !
!     <= gives:      xc,yc,zc c>0    !
!++++++++++++++++++++++++++++++++++++!
    allocate(xc(-NbC:NC)); xc=0.0
    allocate(yc(-NbC:NC)); yc=0.0
    allocate(zc(-NbC:NC)); zc=0.0
    allocate(volume(NC)); volume=0.0

    do c=1,NC
      do n=1,CellN(c,0)
	xc(c) = xc(c) + x(CellN(c,n))/(1.0*CellN(c,0))
	yc(c) = yc(c) + y(CellN(c,n))/(1.0*CellN(c,0))
	zc(c) = zc(c) + z(CellN(c,n))/(1.0*CellN(c,0))
      end do
    end do

    write(*,*) 'Cell centers calculated !'

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate:                                            ! 
!        components of cell sides, cell side centers.       !
!-----------------------------------------------------------!
!     => depends on: x,y,z                                  !
!     <= gives:      Sx,Sy,Sz,xsp,yzp,zsp                   !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    allocate(Sx(NS+max(NC,NBC))); Sx=0.0
    allocate(Sy(NS+max(NC,NBC))); Sy=0.0
    allocate(Sz(NS+max(NC,NBC))); Sz=0.0
    allocate(xsp(NS+max(NC,NBC))); xsp=0.0
    allocate(ysp(NS+max(NC,NBC))); ysp=0.0
    allocate(zsp(NS+max(NC,NBC))); zsp=0.0
    allocate(Dx(NS+max(NC,NBC))); Dx=0.0
    allocate(Dy(NS+max(NC,NBC))); Dy=0.0
    allocate(Dz(NS+max(NC,NBC))); Dz=0.0

    do s=1,NS
      do n=1,SideN(s,0)    ! for quadrilateral an triangular faces
	xt(n)=x(SideN(s,n))
	yt(n)=y(SideN(s,n))
	zt(n)=z(SideN(s,n))
      end do                       

!///// cell side components

      if( SideN(s,0)  ==  4 ) then
	Sx(s)= 0.5 * ( (yt(2)-yt(1))*(zt(2)+zt(1))   &
		      +(yt(3)-yt(2))*(zt(2)+zt(3))   &
		      +(yt(4)-yt(3))*(zt(3)+zt(4))   &
		      +(yt(1)-yt(4))*(zt(4)+zt(1)) )
	Sy(s)= 0.5 * ( (zt(2)-zt(1))*(xt(2)+xt(1))   &
		      +(zt(3)-zt(2))*(xt(2)+xt(3))   &
		      +(zt(4)-zt(3))*(xt(3)+xt(4))   &
		      +(zt(1)-zt(4))*(xt(4)+xt(1)) )
	Sz(s)= 0.5 * ( (xt(2)-xt(1))*(yt(2)+yt(1))   & 
		      +(xt(3)-xt(2))*(yt(2)+yt(3))   &
		      +(xt(4)-xt(3))*(yt(3)+yt(4))   &
		      +(xt(1)-xt(4))*(yt(4)+yt(1)) )
      else if( SideN(s,0)  ==  3 ) then 
	Sx(s)= 0.5 * ( (yt(2)-yt(1))*(zt(2)+zt(1))   & 
		      +(yt(3)-yt(2))*(zt(2)+zt(3))   &
		      +(yt(1)-yt(3))*(zt(3)+zt(1)) )
	Sy(s)= 0.5 * ( (zt(2)-zt(1))*(xt(2)+xt(1))   &
		      +(zt(3)-zt(2))*(xt(2)+xt(3))   & 
		      +(zt(1)-zt(3))*(xt(3)+xt(1)) )
	Sz(s)= 0.5 * ( (xt(2)-xt(1))*(yt(2)+yt(1))   &
		      +(xt(3)-xt(2))*(yt(2)+yt(3))   & 
		      +(xt(1)-xt(3))*(yt(3)+yt(1)) )
      else
	write(*,*) 'calc4: something horrible has happened !'
	stop
      end if

!---- barycenters
      if(SideN(s,0) == 4) then  
	xsp(s) = (xt(1)+xt(2)+xt(3)+xt(4))/4.0
	ysp(s) = (yt(1)+yt(2)+yt(3)+yt(4))/4.0
	zsp(s) = (zt(1)+zt(2)+zt(3)+zt(4))/4.0
      else if(SideN(s,0) == 3) then  
	xsp(s) = (xt(1)+xt(2)+xt(3))/3.0
	ysp(s) = (yt(1)+yt(2)+yt(3))/3.0
	zsp(s) = (zt(1)+zt(2)+zt(3))/3.0
      end if 

    end do ! through sides

    write(*,*) 'Cell face components calculated !'

!+++++++++++++++++++++++++++++++++++++++++!
!     Calculate boundary cell centers     !
!-----------------------------------------!
!     => depends on: xc,yc,zc,Sx,Sy,Sz    !
!     <= gives:      xc,yc,zc for c<0     !   
!+++++++++++++++++++++++++++++++++++++++++!
  write(*,*) 'Position the boundary cell centres:'
  write(*,*) 'Type 1 for barycentric placement'
  write(*,*) 'Type 2 for orthogonal placement'
  read(*,*) BouCen 

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

    SurTot = sqrt(Sx(s)*Sx(s)+Sy(s)*Sy(s)+Sz(s)*Sz(s))

    if(c2  < 0) then
      t = (   Sx(s)*(xsp(s)-xc(c1))                               &
            + Sy(s)*(ysp(s)-yc(c1))                               &
            + Sz(s)*(zsp(s)-zc(c1)) ) / SurTot
      xc(c2) = xc(c1) + Sx(s)*t / SurTot
      yc(c2) = yc(c1) + Sy(s)*t / SurTot
      zc(c2) = zc(c1) + Sz(s)*t / SurTot
      if(BouCen == 1) then
        xc(c2) = xsp(s)
        yc(c2) = ysp(s)
        zc(c2) = zsp(s)
      end if
    endif 
  end do ! through sides

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!    Move the centers of co-planar molecules towards the walls    !
!-----------------------------------------------------------------!
!     => depends on: xc,yc,zc                                     !
!     <= gives:      xc,yc,zc                                     !
!     +  uses:       Dx,Dy,Dz                                     !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

    if(c2 > 0) then
      Dx(c1) = max( Dx(c1), abs( xc(c2)-xc(c1) ) )
      Dy(c1) = max( Dy(c1), abs( yc(c2)-yc(c1) ) )
      Dz(c1) = max( Dz(c1), abs( zc(c2)-zc(c1) ) )
      Dx(c2) = max( Dx(c2), abs( xc(c2)-xc(c1) ) )
      Dy(c2) = max( Dy(c2), abs( yc(c2)-yc(c1) ) )
      Dz(c2) = max( Dz(c2), abs( zc(c2)-zc(c1) ) )
    end if 
  end do ! through sides

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

    if(c2 < 0) then
      if( Approx(Dx(c1), 0.0, 1.e-6) ) xc(c1) = 0.75*xc(c1) + 0.25*xc(c2) 
      if( Approx(Dy(c1), 0.0, 1.e-6) ) yc(c1) = 0.75*yc(c1) + 0.25*yc(c2) 
      if( Approx(Dz(c1), 0.0, 1.e-6) ) zc(c1) = 0.75*zc(c1) + 0.25*zc(c2) 
    end if 
  end do ! through sides

  Dx = 0.0
  Dy = 0.0
  Dz = 0.0

!+++++++++++++++++++++++++++++++++++++++++++++++++!
!     Find the sides on the periodic boundary     !
!-------------------------------------------------!
!     => depends on: xc,yc,zc,Sx,Sy,Sz            !
!     <= gives:      Dx,Dy,Dz                     !
!+++++++++++++++++++++++++++++++++++++++++++++++++!
!    I am trying to do it with Gambit             !
!+++++++++++++++++++++++++++++++++++++++++++++++++!

  allocate(BCoor(NS)); BCoor=0.0
  allocate(BFace(NS)); Bface=0

!===================!
!===== Phase I =====! -> find the sides on periodic boundaries
!===================!
!  write(*,*) 'Type 1 for fast but unreliable algorithm for periodic cells search.'
!  write(*,*) 'Type 2 for slow but reliable algorithm for periodic cells search.'
!  read(*,*) OPT

  OPT = 2

2 Nper = 0 
  write(*,*) 'Type the periodic-boundary-condition marker number.'
  write(*,*) 'Type 0 if there is none !'
  write(*,*) ''
  write(*,*) '(Please note that the periodic boundaries have to be the last on the list)'
  write(*,*) 'of the boundary conditions. Their BC markers have to be larger'
  write(*,*) 'than the markers of the other boundary conditions.)'
  write(*,*) ''
  read(*,*) TypePer
  if( TypePer == 0 ) goto 1  

  if(OPT == 2) then
    write(*,*) 
    write(*,*) 'Enter the angle for the rotation of the coordiante system (in degree)'
    write(*,*) 'and the axes of rotation (1 -> x, 2 -> y, 3 -> z)'
    write(*,*) 
    write(*,*) '(if the periodic direction is not parallel to the Caresian axis (x, y and z),'
    write(*,*) 'the coordinate system has to be rotated (in 2D))'
    read(*,*) angle, rot_dir 

    angle = angle * 3.1415926 / 180.0

    write(*,*) 
    write(*,*) 'Enter the coordinates of three points that define'
    write(*,*) 'periodic plane'
    write(*,*) 
    write(*,*) 'Point 1:'
    read(*,*) x_a, y_a, z_a  !angle 
    write(*,*) 
    write(*,*) 'Point 2:'
    read(*,*) x_b, y_b, z_b  !angle 
    write(*,*) 
    write(*,*) 'Point 3:'
    read(*,*) x_c, y_c, z_c  !angle 

!    write(*,*) 
!    write(*,*) 'Enter approximative distance between the periodic faces'
!    read(*,*) PerMax


    ABi = x_b - x_a 
    ABj = y_b - y_a 
    ABk = z_b - z_a 

    ACi = x_c - x_a 
    ACj = y_c - y_a 
    ACk = z_c - z_a

    Pi = ABj*ACk - ACj*ABk 
    Pj = -ABi*ACk + ACi*ABk 
    Pk = ABi*ACj - ACi*ABj

    angle_face = angle_face * 3.1415926 / 180.0
 
    allocate(phi_face(NS)); phi_face=0.0
    allocate(xspr(NS)); xspr=0.0
    allocate(yspr(NS)); yspr=0.0
    allocate(zspr(NS)); zspr=0.0

    do s=1,NS
      c2 = SideC(2,s)
      if(c2 < 0) then
        if(BCmark(c2) == TypePer) then
          if( Approx(angle, 0.0, 1.e-6) ) then
            xspr(s) = xsp(s)  
            yspr(s) = ysp(s)  
            zspr(s) = zsp(s)  
          else 
            if(rot_dir==3) then
              xspr(s) = xsp(s)*cos(angle) + ysp(s)*sin(angle) 
              yspr(s) =-xsp(s)*sin(angle) + ysp(s)*cos(angle)
            else if(rot_dir==2) then
              xspr(s) = xsp(s)*cos(angle) + zsp(s)*sin(angle) 
              zspr(s) =-xsp(s)*sin(angle) + zsp(s)*cos(angle)
            else if(rot_dir==1) then
              yspr(s) = ysp(s)*cos(angle) + zsp(s)*sin(angle) 
              zspr(s) =-ysp(s)*sin(angle) + zsp(s)*cos(angle)
            end if
          end if
        end if
      end if
    end do
  end if

  xmin = +HUGE
  ymin = +HUGE
  zmin = +HUGE
  xmax = -HUGE
  ymax = -HUGE
  zmax = -HUGE

  write(*,*) 'Insert the periodic direction (1 -> x, 2 -> y, 3 -> z)'
  read(*,*) dir 

  BCoor=0.0
  Bface=0

  c = 0

  if(OPT == 1) then 
    do s=1,NS
      c2 = SideC(2,s)
      if(c2 < 0) then
        if(BCmark(c2) == TypePer) then
          c = c + 1
          if(dir==1) BCoor(c) = xsp(s)*1000000.0 + ysp(s)*10000.0 + zsp(s)
          if(dir==2) BCoor(c) = ysp(s)*1000000.0 + xsp(s)*10000.0 + zsp(s)
          if(dir==3) BCoor(c) = zsp(s)*1000000.0 + xsp(s)*10000.0 + ysp(s)
          BFace(c) = s
        end if
      end if
    end do
    call DISort(BCoor,BFace,c,2)
  else if(OPT == 2) then 
    c_max = 0
    do s=1,NS
      c2 = SideC(2,s)
      if(c2 < 0) then
        if(BCmark(c2) == TypePer) then
          c_max = c_max + 1
        end if 
      end if 
    end do
    Tol = 0.0001 
    
10 continue


    NNN = 0
    hh = 0 
    mm = 0 
    c = 0
  
    PerMax = -HUGE
    PerMin =  HUGE 

    do s=1,NS
      c2 = SideC(2,s)
      if(c2 < 0) then
        if(BCmark(c2) == TypePer) then
          Det = (Pi*(xsp(s)) + Pj*(ysp(s)) + Pk*(zsp(s)))/sqrt(Pi*Pi + Pj*Pj + Pk*Pk)
          PerMin = min(PerMin, Det)          
          PerMax = max(PerMax, Det)          
        end if
      end if
    end do
    
    PerMax = 0.5*(PerMax + PerMin)

    do s=1,NS
      c2 = SideC(2,s)
      if(c2 < 0) then
        if(BCmark(c2) == TypePer) then
          c = c + 1
          Det = (Pi*(xsp(s)) + Pj*(ysp(s)) + Pk*(zsp(s)))/sqrt(Pi*Pi + Pj*Pj + Pk*Pk)

          if(dir==1) then
            if((Det) < (PerMax)) then
              hh = hh + 1
              BCoor(hh) = hh
              BFace(hh) = s
              do ss=1,NS
                cc2 = SideC(2,ss)
                if(cc2 < 0) then
                  if(BCmark(cc2) == TypePer) then 
                    Det = (Pi*(xsp(ss)) + Pj*(ysp(ss)) + Pk*(zsp(ss)))/sqrt(Pi*Pi + Pj*Pj + Pk*Pk)
                    if((Det) > (PerMax)) then
                      if((abs(zsp(ss)-zsp(s))) < Tol.and.         &
                         (abs(yspr(ss)-yspr(s))) < Tol) then
                         mm = hh + c_max/2 
                         BCoor(mm) = mm
                         BFace(mm) = ss
                         NNN = NNN + 1
                      end if 
                    end if 
                  end if
                end if  
              end do
            end if          
          end if


          if(dir==2) then
            if((Det) < (PerMax)) then
              hh = hh + 1
              BCoor(hh) = hh
              BFace(hh) = s
              do ss=1,NS
                cc2 = SideC(2,ss)
                if(cc2 < 0) then
                  if(BCmark(cc2) == TypePer) then 
                    Det = (Pi*(xsp(ss)) + Pj*(ysp(ss)) + Pk*(zsp(ss)))/sqrt(Pi*Pi + Pj*Pj + Pk*Pk)
                    if((Det) > (PerMax)) then
                      if(abs((zsp(ss)-zsp(s))) < Tol.and.         &
                        abs((xspr(ss)-xspr(s))) < Tol) then
                        mm = hh + c_max/2 
                        BCoor(mm) = mm
                        BFace(mm) = ss
                        NNN = NNN + 1
                      end if 
                    end if 
                  end if
                end if  
              end do
            end if          
          end if

          if(dir==3) then
            if((Det) < (PerMax)) then
              hh = hh + 1
              BCoor(hh) = hh
              BFace(hh) = s
              do ss=1,NS
                cc2 = SideC(2,ss)
                if(cc2 < 0) then
                  if(BCmark(cc2) == TypePer) then
                    Det = (Pi*(xsp(ss)) + Pj*(ysp(ss)) + Pk*(zsp(ss)))/sqrt(Pi*Pi + Pj*Pj + Pk*Pk)
                    if((Det) > (PerMax)) then
                      if(abs((xsp(ss)-xsp(s))) < Tol.and.         &
                        abs((ysp(ss)-ysp(s))) < Tol) then
                        mm = hh + c_max/2 
                        BCoor(mm) = mm
                        BFace(mm) = ss
                        NNN = NNN + 1
                      end if 
                    end if 
                  end if
                end if  
              end do
            end if          
          end if

        end if 
      end if 
    end do

    write(*,*)'Iterating search for periodic cells: ',  &
    'Target: ', c_max/2, 'Result: ',NNN, 'Tolerance: ',Tol

    if(NNN == c_max/2) then
      continue
    else
      Tol = Tol*0.5 
      go to 10
    end if

    deallocate(phi_face)
    deallocate(xspr)
    deallocate(yspr)
    deallocate(zspr)

    call DISort(BCoor,BFace,c,2)
  end if  ! end OPT
  

  do s=1,c/2
    s1 = BFace(s)
    s2 = BFace(s+c/2)
    c11 = SideC(1,s1) ! cell 1 for side 1
    c21 = SideC(2,s1) ! cell 2 for cell 1
    c12 = SideC(1,s2) ! cell 1 for side 2
    c22 = SideC(2,s2) ! cell 2 for side 2
!>>>    write(*,*) c
!>>>    write(*,*) s1, BFace(s),     xc(c21), yc(c21), zc(c21) 
!>>>    write(*,*) s2, BFace(s+c/2), xc(c22), yc(c22), zc(c22) 
!>>>    write(*,*) '---------------------------------'
    SideC(0,s1) = s2   ! just to remember where it was coppied from
    SideC(2,s1) = c12 
    SideC(1,s2) = 0    ! c21 
    SideC(2,s2) = 0    ! c21 
  end do

  Nper = c/2
  write(*,*) 'Phase I: periodic cells: ', Nper
  go to 2

!====================!
!===== Phase II =====! -> similar to the loop in Generator
!====================!
1 continue 
  NSsh = 0
  do s=1,NS

!----- initialize
    Dx(s)=0.0
    Dy(s)=0.0
    Dz(s)=0.0

    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2   >  0) then

!----- scalar product of the side with line c1-c2 is good criteria
      if( (Sx(s) * (xc(c2)-xc(c1) )+                              &
           Sy(s) * (yc(c2)-yc(c1) )+                              &
           Sz(s) * (zc(c2)-zc(c1) ))  < 0.0 ) then

        NSsh = NSsh + 2
!----- find the coordinates of ...
        if(SideN(s,0) == 4) then
          !---- coordinates of the shadow face
          xs2=xsp(SideC(0,s))
          ys2=ysp(SideC(0,s))
          zs2=zsp(SideC(0,s))
          !---- add shadow faces
          SideN(NS+NSsh-1,0) = 4
          SideC(1,NS+NSsh-1) = c1
          SideC(2,NS+NSsh-1) = -NbC-1
          SideN(NS+NSsh-1,1) = SideN(s,1)
          SideN(NS+NSsh-1,2) = SideN(s,2)
          SideN(NS+NSsh-1,3) = SideN(s,3)
          SideN(NS+NSsh-1,4) = SideN(s,4)
          Sx(NS+NSsh-1) = Sx(s)
          Sy(NS+NSsh-1) = Sy(s)
          Sz(NS+NSsh-1) = Sz(s)
          xsp(NS+NSsh-1) = xsp(s)
          ysp(NS+NSsh-1) = ysp(s)
          zsp(NS+NSsh-1) = zsp(s)
          SideN(NS+NSsh,0) = 4
          SideC(1,NS+NSsh) = c2
          SideC(2,NS+NSsh) = -NbC-1
          SideN(NS+NSsh,1) = SideN(SideC(0,s),1) 
          SideN(NS+NSsh,2) = SideN(SideC(0,s),2)
          SideN(NS+NSsh,3) = SideN(SideC(0,s),3)
          SideN(NS+NSsh,4) = SideN(SideC(0,s),4)
          Sx(NS+NSsh) = Sx(s)
          Sy(NS+NSsh) = Sy(s)
          Sz(NS+NSsh) = Sz(s)
          xsp(NS+NSsh) = xs2
          ysp(NS+NSsh) = ys2
          zsp(NS+NSsh) = zs2
        else if(SideN(s,0) == 3) then
          !---- coordinates of the shadow face
          xs2=xsp(SideC(0,s))
          ys2=ysp(SideC(0,s))
          zs2=zsp(SideC(0,s))
          !---- add shadow faces
          SideN(NS+NSsh-1,0) = 3
          SideC(1,NS+NSsh-1) = c1
          SideC(2,NS+NSsh-1) = -NbC-1
          SideN(NS+NSsh-1,1) = SideN(s,1)
          SideN(NS+NSsh-1,2) = SideN(s,2)
          SideN(NS+NSsh-1,3) = SideN(s,3)
          Sx(NS+NSsh-1) = Sx(s)
          Sy(NS+NSsh-1) = Sy(s)
          Sz(NS+NSsh-1) = Sz(s)
          xsp(NS+NSsh-1) = xsp(s)
          ysp(NS+NSsh-1) = ysp(s)
          zsp(NS+NSsh-1) = zsp(s)
          SideN(NS+NSsh,0) = 3
          SideC(1,NS+NSsh) = c2
          SideC(2,NS+NSsh) = -NbC-1
          SideN(NS+NSsh,1) = SideN(SideC(0,s),1) 
          SideN(NS+NSsh,2) = SideN(SideC(0,s),2)
          SideN(NS+NSsh,3) = SideN(SideC(0,s),3)
          Sx(NS+NSsh) = Sx(s)
          Sy(NS+NSsh) = Sy(s)
          Sz(NS+NSsh) = Sz(s)
          xsp(NS+NSsh) = xs2
          ysp(NS+NSsh) = ys2
          zsp(NS+NSsh) = zs2
        end if

        Dx(s)=xsp(s)-xs2  !------------------------!
        Dy(s)=ysp(s)-ys2  ! later: xc2 = xc2 + Dx  !
        Dz(s)=zsp(s)-zs2  !------------------------!

!->>>     write(6,'(I5,A12,3F12.6)') s,' Dx,Dy,Dz= ', &
!->>>               Dx(s),Dy(s),Dz(s)
      endif !  S*(c2-c1) < 0.0
    end if  !  c2 > 0
  end do    !  sides
  write(*,*) 'Phase II: number of shadow faces: ', NSsh

!=====================!
!===== Phase III =====! -> find the new numbers of cell faces
!=====================!
  NumberSides = 0
  do s=1,NS+NSsh
    c1 = SideC(1,s)
    c2 = SideC(2,s)
    if(c1 > 0) then
      NumberSides = NumberSides  + 1
      NewS(s) = NumberSides 
    else
      NewS(s) = -1
    end if
  end do
  write(*,'(A21,I9,Z9)') 'Old number of sides: ', NS, NS
  write(*,'(A21,I9,Z9)') 'New number of sides: ', &
                          NumberSides-NSsh,NumberSides-NSsh
  
!====================!
!===== Phase IV =====! -> compress the sides
!====================!
  do s=1,NS+NSsh
    if(NewS(s) > 0) then
      SideC(1,NewS(s)) = SideC(1,s) 
      SideC(2,NewS(s)) = SideC(2,s)
      SideN(NewS(s),0) = SideN(s,0)
      SideN(NewS(s),1) = SideN(s,1)
      SideN(NewS(s),2) = SideN(s,2)
      SideN(NewS(s),3) = SideN(s,3)
      SideN(NewS(s),4) = SideN(s,4)
      xsp(NewS(s)) = xsp(s)
      ysp(NewS(s)) = ysp(s)
      zsp(NewS(s)) = zsp(s)
      Sx(NewS(s)) = Sx(s)
      Sy(NewS(s)) = Sy(s)
      Sz(NewS(s)) = Sz(s)
      Dx(NewS(s)) = Dx(s)
      Dy(NewS(s)) = Dy(s)
      Dz(NewS(s)) = Dz(s)
    end if
  end do 
  NS = NumberSides-NSsh

!+++++++++++++++++++++++++++++++++++++++!
!     Check the periodic boundaries     !
!+++++++++++++++++++++++++++++++++++++++!
  MaxDis = 0.0 
  do s=1,NS-NSsh
    MaxDis = max(MaxDis, (Dx(s)*Dx(s)+Dy(s)*Dy(s)+Dz(s)*Dz(s)))
  end do
  write(*,*) 'Maximal distance of periodic boundary is:', sqrt(MaxDis)

!++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate the cell volumes                   !
!--------------------------------------------------!
!     => depends on: xc,yc,zc,                     !
!                    Dx,Dy,Dz,                     !
!                    xsp, ysp, zsp                 !
!     <= gives:      volume                        !
!++++++++++++++++++++++++++++++++++++++++++++++++++!

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)   

    volume(c1) = volume(c1) + xsp(s)*Sx(s)  &
                            + ysp(s)*Sy(s)  &
                            + zsp(s)*Sz(s)
    if(c2 > 0) then
      volume(c2) = volume(c2) - (xsp(s)-Dx(s))*Sx(s)  &
                              - (ysp(s)-Dy(s))*Sy(s)  &
                              - (zsp(s)-Dz(s))*Sz(s)
    end if
  end do
  volume = volume/3.0
  c1 = 0
  MinVol =  1E+30
  MaxVol = -1E+30
  TotVol = 0.0
  do c=1,NC
    TotVol = TotVol + volume(c)
    MinVol = min(MinVol, volume(c))
    MaxVol = max(MaxVol, volume(c))
  end do
  write(*,*) 'Minimal cell volume is: ', MinVol
  write(*,*) 'Maximal cell volume is: ', MaxVol
  write(*,*) 'Total domain volume is: ', TotVol
  write(*,*) 'Cell volumes calculated !'

  if(MinVol < 0.0) then
    write(*,*) 'Negative volume occured! Another, slower, algoritham should be run !'
    write(*,*) 'Execution will be halt now! '
    stop
  end if 
 
  deallocate(BCoor)
  deallocate(Bface)
 

!++++++++++++++++++++++++++++++!
!     Calculate delta          !
!------------------------------!
!     => depends on: x,y,z     !
!     <= gives:      delta     !
!++++++++++++++++++++++++++++++!
    allocate(delta(-NbC:NC));  delta=0.0

    do c=1,NC
      delta(c)=0.0
      xmin = +HUGE
      ymin = +HUGE
      zmin = +HUGE
      xmax = -HUGE
      ymax = -HUGE
      zmax = -HUGE
      do n=1,CellN(c,0)
        xmin = min(xmin, x(CellN(c,n)))
        ymin = min(ymin, y(CellN(c,n)))
        zmin = min(zmin, z(CellN(c,n)))
        xmax = max(xmax, x(CellN(c,n)))
        ymax = max(ymax, y(CellN(c,n)))
        zmax = max(zmax, z(CellN(c,n)))
      end do
      delta(c) = xmax-xmin
      delta(c) = max(delta(c), (ymax-ymin))
      delta(c) = max(delta(c), (zmax-zmin))
    end do

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate:                                              ! 
!        distance from the cell center to the nearest wall.   !
!-------------------------------------------------------------!
!     => depends on: xc,yc,zc inside and on the boundary.     !
!     <= gives:      WallDs i                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  allocate(WallDs(-NbC:NC)); WallDs = HUGE

  write(*,*) ''
  write(*,*) 'Type the total number of wall boundary conditions:'
  write(*,*) ''
  write(*,*) '(Please note that the walls have to be the first on the list)'
  write(*,*) 'of the boundary conditions. Their BC markers have to be smaller'
  write(*,*) 'than the markers of the other boundary conditions.)'
  write(*,*) ''
  read(*,*) WallMark
 
  if(WallMark == 0) then
    WallDs = 1.0
    write(*,*) 'Distance to the wall set to 1.0 everywhere !'
  else
    do c1=1,NC
      if(mod(c1,10000) == 0) write(*,*) (100.*c1/(1.*NC)), '% complete...'
      do c2=-1,-NbC,-1
        if(BCmark(c2) <= WallMark) then
          WallDs(c1)=min(WallDs(c1),                       &
          Dist2(xc(c1),yc(c1),zc(c1),xc(c2),yc(c2),zc(c2)))
        end if
      end do
    end do

    WallDs = sqrt(WallDs)

    write(*,*) 'Distance to the wall calculated !'
  end if


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate the interpolation factors for the cell sides     !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  allocate(f(NS+max(NC,NBC))); f=0.0          

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

!----- first cell
    xc1=xc(c1)
    yc1=yc(c1)
    zc1=zc(c1)
    dsc1=dist(xc1,yc1,zc1,xsp(s), ysp(s), zsp(s))

!----- second cell (pls. check if xsi=xc on the boundary)
    xc2=xc(c2)+Dx(s)
    yc2=yc(c2)+Dy(s)
    zc2=zc(c2)+Dz(s)
    dsc2=dist(xc2,yc2,zc2,xsp(s), ysp(s), zsp(s))

!----- interpolation factor
    f(s) = dsc2 / (dsc1+dsc2)   ! not checked
  end do 

  write(*,*) 'Interpolation factors calculated !'


  return

3 write(*,*) 'Horror ! Negative volume between cells ', c1, ' and ', c2
  stop

  END SUBROUTINE Calc4
