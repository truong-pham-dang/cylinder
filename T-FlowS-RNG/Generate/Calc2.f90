!======================================================================!
  SUBROUTINE Calc2(rrun)
!----------------------------------------------------------------------!
!   Calculates geometrical quantities of the grid.                     !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  LOGICAL rrun     
!------------------------------[Calling]-------------------------------!
  REAL :: tetvol        
  REAL :: dist       
  REAL :: dist2       
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, c1, c2, m, n, s, c_1, c_2
  INTEGER :: WallMark
  REAL    :: xt(4), yt(4), zt(4)
  REAL    :: xTc, yTc, zTc       ! temporary cell coordinates 
  REAL    :: xs2,ys2,zs2
  REAL    :: dsc1, dsc2          !  for the interpolation factors
  REAL    :: t, SurTot , maxdis
  REAL    :: xc1, yc1, zc1, xc2, yc2, zc2 
  REAL    :: xmin, xmax, ymin, ymax, zmin, zmax
  INTEGER :: f4n(6,4)
  INTEGER :: f3n(4,3)
!--------------------------------[CVS]---------------------------------!
!  $Id: Calc2.f90,v 1.14 2002/10/31 11:26:46 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Generate/Calc2.f90,v $  
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
  DATA    f4n / 1, 1, 2, 4, 3, 5,                                   &
		2, 5, 6, 8, 7, 7,                                   &
		4, 6, 8, 7, 5, 8,                                   &
		3, 2, 4, 3, 1, 6  /

  DATA    f3n / 1,  1,  2,  3,                                      &
		2,  4,  4,  4,                                      &
		3,  2,  3,  1 /

!---- Without the following six lines, this procedure works for any grid
    do c=1,NC
      CellN(c,0)=8
    end do
    do s=1,NS
      SideN(s,0)=4 
    end do

!++++++++++++++++++++++++++++++++++++!
!     Calculate the cell centers     !
!------------------------------------!
!     => depends on: x,y,z           !
!     <= gives:      xc,yc,zc c>0    !
!++++++++++++++++++++++++++++++++++++!
    do c=1,NC
      xc(c)=0.0
      yc(c)=0.0
      zc(c)=0.0
      do n=1,CellN(c,0)
	xc(c) = xc(c) + x(CellN(c,n))/(1.0*CellN(c,0))
	yc(c) = yc(c) + y(CellN(c,n))/(1.0*CellN(c,0))
	zc(c) = zc(c) + z(CellN(c,n))/(1.0*CellN(c,0))
      end do
    end do

!++++++++++++++++++++++++++++++!
!     Calculate delta          !
!------------------------------!
!     => depends on: x,y,z     !
!     <= gives:      delta     !
!++++++++++++++++++++++++++++++!
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

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate:                                            ! 
!        components of cell sides, cell side centers.       !
!-----------------------------------------------------------!
!     => depends on: x,y,z                                  !
!     <= gives:      Sx,Sy,Sz,xsp,yzp,zsp                   !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
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
	write(*,*) 'calc2: something horrible has happened !'
	stop
      end if

!->>> write(*,'(3I5,3F8.4)') s,SideC(1,s),SideC(2,s),Sx(s),Sy(s),Sz(s)

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

!+++++++++++++++++++++++++++++++++++++++++!
!     Calculate boundary cell centers     !
!-----------------------------------------!
!     => depends on: xc,yc,zc,Sx,Sy,Sz    !
!     <= gives:      xc,yc,zc for c<0     !   
!+++++++++++++++++++++++++++++++++++++++++!
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
      endif 
    end do ! through sides

!->>>      do c=1,NC
!->>>        write(6,'(I7,3F8.4)') c, xc(c), yc(c), zc(c) 
!->>>      end do
!->>>      do c=-1,-NbC,-1
!->>>        write(6,'(I7,3F8.4)') c, xc(c), yc(c), zc(c)
!->>>      end do

!+++++++++++++++++++++++++++++++++++++++++++++++++!
!     Find the sides on the periodic boundary     !
!-------------------------------------------------!
!     => depends on: xc,yc,zc,Sx,Sy,Sz            !
!     <= gives:      Dx,Dy,Dz         !
!+++++++++++++++++++++++++++++++++++++++++++++++++!
    if(rrun) then
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
	  m=SideCc(s,2)

	  if(SideN(s,0) == 4) then   
            !---- coordinates of the shadow face
	    xs2=.25*(x(CellN(c2,f4n(m,1)))+x(CellN(c2,f4n(m,2)))+   &
		     x(CellN(c2,f4n(m,3)))+x(CellN(c2,f4n(m,4))) )
	    ys2=.25*(y(CellN(c2,f4n(m,1)))+y(CellN(c2,f4n(m,2)))+   &
		     y(CellN(c2,f4n(m,3)))+y(CellN(c2,f4n(m,4))) )
	    zs2=.25*(z(CellN(c2,f4n(m,1)))+z(CellN(c2,f4n(m,2)))+   &
		     z(CellN(c2,f4n(m,3)))+z(CellN(c2,f4n(m,4))) )
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
	    SideN(NS+NSsh,1) = CellN(c2,f4n(m,1)) 
	    SideN(NS+NSsh,2) = CellN(c2,f4n(m,2))
	    SideN(NS+NSsh,3) = CellN(c2,f4n(m,3))
	    SideN(NS+NSsh,4) = CellN(c2,f4n(m,4))
            Sx(NS+NSsh) = Sx(s)
            Sy(NS+NSsh) = Sy(s)
            Sz(NS+NSsh) = Sz(s)
            xsp(NS+NSsh) = xs2
            ysp(NS+NSsh) = ys2
            zsp(NS+NSsh) = zs2
 	  else if(SideN(s,0) == 3) then  
            !---- coordinates of the shadow face
	    xs2=.33333333 * (x(CellN(c2,f3n(m,1)))+                 &
		     x(CellN(c2,f3n(m,2)))+x(CellN(c2,f3n(m,3))) )
	    ys2=.33333333 * (y(CellN(c2,f4n(m,1)))+                 &
		     y(CellN(c2,f3n(m,2)))+y(CellN(c2,f3n(m,3))) )
	    zs2=.33333333 * (z(CellN(c2,f4n(m,1)))+                 &
		     z(CellN(c2,f3n(m,2)))+z(CellN(c2,f3n(m,3))) )
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
	    SideN(NS+NSsh,1) = CellN(c2,f3n(m,1)) 
	    SideN(NS+NSsh,2) = CellN(c2,f3n(m,2))
	    SideN(NS+NSsh,3) = CellN(c2,f3n(m,3))
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

!->>>      write(6,'(2I5,A12,3F12.6)') s,' Dx,Dy,Dz= ', 
!->>>     >           Dx(s),Dy(s),Dz(s)
	endif !  S*(c2-c1) < 0.0
      end if  !  c2 > 0
    end do    !  sides  
    write(*,*) 'Number of shadow faces: ', NSsh
    end if

!++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate the cell volumes                   !
!--------------------------------------------------!
!     => depends on: xc,yc,zc,                     !
!                    Dx,Dy,Dz,      !
!                    xsp, ysp, zsp                 !
!     <= gives:      volume                        !
!++++++++++++++++++++++++++++++++++++++++++++++++++!
  if(rrun) then
  do c=1,NC
    volume(c)=0.0
  end do

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)   

    do n=1,SideN(s,0)      ! for quadrilateral an triangular faces
      xt(n)=x(SideN(s,n))
      yt(n)=y(SideN(s,n))
      zt(n)=z(SideN(s,n))
    end do   

!----- first cell
    xTc=xc(c1)
    yTc=yc(c1)
    zTc=zc(c1)
    dsc1=dist(xTc,yTc,zTc,xsp(s), ysp(s), zsp(s)) 
    volume(c1)=volume(c1) + tetvol(xsp(s),ysp(s),zsp(s),            &
	       xt(1),yt(1),zt(1),xt(2),yt(2),zt(2),xTc,yTc,zTc)
    volume(c1)=volume(c1) + tetvol(xsp(s),ysp(s),zsp(s),            &
	       xt(2),yt(2),zt(2),xt(3),yt(3),zt(3),xTc,yTc,zTc)
    if(SideN(s,0) == 4) then
      volume(c1)=volume(c1) + tetvol(xsp(s),ysp(s),zsp(s),          &
		 xt(3),yt(3),zt(3),xt(4),yt(4),zt(4),xTc,yTc,zTc)
      volume(c1)=volume(c1) + tetvol(xsp(s),ysp(s),zsp(s),          &
		 xt(4),yt(4),zt(4),xt(1),yt(1),zt(1),xTc,yTc,zTc)
    else if(SideN(s,0) == 3) then
      volume(c1)=volume(c1) + tetvol(xsp(s),ysp(s),zsp(s),          &
		 xt(3),yt(3),zt(3),xt(1),yt(1),zt(1),xTc,yTc,zTc)
    end if

!----- second cell
    if(c2  > 0) then
      xTc=xc(c2)+Dx(s)
      yTc=yc(c2)+Dy(s)
      zTc=zc(c2)+Dz(s)
      dsc2=dist(xTc,yTc,zTc,xsp(s), ysp(s), zsp(s)) 
      volume(c2)=volume(c2) -tetvol(xsp(s),ysp(s),zsp(s),           &
		 xt(1),yt(1),zt(1),xt(2),yt(2),zt(2),xTc,yTc,zTc)
      volume(c2)=volume(c2) -tetvol(xsp(s),ysp(s),zsp(s),           &
		 xt(2),yt(2),zt(2),xt(3),yt(3),zt(3),xTc,yTc,zTc)
      if(SideN(s,0) == 4) then
	volume(c2)=volume(c2) -tetvol(xsp(s),ysp(s),zsp(s),         &
		   xt(3),yt(3),zt(3),xt(4),yt(4),zt(4),xTc,yTc,zTc)
	volume(c2)=volume(c2) -tetvol(xsp(s),ysp(s),zsp(s),         &
		   xt(4),yt(4),zt(4),xt(1),yt(1),zt(1),xTc,yTc,zTc)
      else if(SideN(s,0) == 3) then
	volume(c2)=volume(c2) -tetvol(xsp(s),ysp(s),zsp(s),         &
		   xt(3),yt(3),zt(3),xt(1),yt(1),zt(1),xTc,yTc,zTc)
      end if  
    else        
      dsc2=0.0
    end if

  end do
  end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate:                                              ! 
!        distance from the cell center to the nearest wall.   !
!-------------------------------------------------------------!
!     => depends on: xc,yc,zc inside and on the boundary.     !
!     <= gives:      WallDs i                              !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  WallDs = HUGE 

  write(*,'(A43)')   '#------------------------------------------'
  if(rrun) then
    write(*,'(A43)') '# Computing the distance to the walls (2/2)'           
  else            
    write(*,'(A43)') '# Computing the distance to the walls (1/2)'           
  end if 
  write(*,'(A43)')   '#------------------------------------------'
  write(*,*) 'Insert the highest BC marker which represents the wall'
  write(*,*) 'for computing the distance to the wall (0 to skip)'
  write(*,*) 'If you skip this, smoothing will not work properly'
  read(*,*) WallMark   

  if(WallMark == 0) then
    WallDs = 1.0
    write(*,*) 'Distance to the wall set to 1 everywhere !'            
  else 
    do c1=1,NC 
      do s = WallFacFst, WallFacLst      ! 1,NS
        c_1 = SideC(1,s)
        c_2 = SideC(2,s)
        if(c_2 < 0) then
          if(BCmark(c_2) <= WallMark) then
            WallDs(c1)=min(WallDs(c1), &
            dist2(xc(c1),yc(c1),zc(c1),xsp(s),ysp(s),zsp(s)))
          end if
!        else if(material(c_1) /= material(c_2) ) then
!          WallDs(c1)=min(WallDs(c1), &
!          dist2(xc(c1),yc(c1),zc(c1),xsp(s),ysp(s),zsp(s)))
        end if 
      end do
    end do

    do c=1,NC
      WallDs(c)=sqrt(WallDs(c))
    end do

    write(*,*) 'Maximal distance to the wall: ', maxval(WallDs(1:NC))
    write(*,*) 'Minimal distance to the wall: ', minval(WallDs(1:NC))
  end if

  do n=1,NN
    walln(n)=HUGE
  end do

  do c=1,NC
    do n=1,CellN(c,0)
      walln(CellN(c,n))=min(WallDs(c),walln(CellN(c,n)))
    end do
  end do

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2 < 0 .or. material(c1) /= material(c2)) then 
      do n=1,SideN(s,0)    ! for quadrilateral an triangular faces
	walln(SideN(s,n)) = 0.0
      end do
    end if
  end do 

  maxdis=0.0 
  do n=1,NN
    maxdis=max(walln(n),maxdis)
  end do

  do n=1,NN
    walln(n)=walln(n)/maxdis
!->>>	write(*,*) walln(n)
  end do

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate the interpolation factors for the cell sides     !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  if(rrun) then
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
  end if 

  END SUBROUTINE Calc2
