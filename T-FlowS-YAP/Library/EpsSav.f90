!======================================================================!
  SUBROUTINE EpsSav(xg,yg,zg,sidegDx,sidegDy,dir)
!----------------------------------------------------------------------!
! Writes: Grid in encapsulated postscript format.                      !
! ~~~~~~~                                                              !
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  REAL      :: xg(MAXN),yg(MAXN),zg(MAXN),sidegDx(MAXS),sidegDy(MAXS)
  CHARACTER :: dir
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: s, c1, c2, count, lw
  CHARACTER :: namEps*80
  REAL      :: sclf, sclp, xmax,xmin,ymax,ymin,zmax,zmin, z0, fc
  REAL      :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xin(4),yin(4) 
!--------------------------------[CVS]---------------------------------!
!  $Id: EpsSav.f90,v 1.8 2002/10/30 16:29:31 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Library/EpsSav.f90,v $            
!======================================================================!
!   Note: singnificantly improved. Keep the old version allways at     !
!         hand.                                                        !
!----------------------------------------------------------------------!

  write(*,'(A38)')       '#-------------------------------------'
  write(*,'(A9,A1,A28)') '# Making ', dir, '.eps cut through the domain'
  write(*,'(A38)')       '#-------------------------------------'

  write(6,*) 'Enter the ',dir,' coordinate for cutting or type 0 to exit: '
  call ReadC(5,inp,tn,ts,te)
  read(inp, *) z0
  if(z0 == 0) return
  write(*,*) 'Z0 = ', z0

!<<<<<<<<<<<<<<<<<<<<<<<<<!
!     create EPS file     !
!<<<<<<<<<<<<<<<<<<<<<<<<<!
  namEps = name
  namEps(len_trim(name)+1:len_trim(name)+6) = '. .eps'
  namEps(len_trim(name)+2:len_trim(name)+2) = dir
  write(6, *) 'Now creating the file:', namEps

  xmax=maxval(xg(1:NN)) 
  ymax=maxval(yg(1:NN)) 
  zmax=maxval(zg(1:NN)) 
  xmin=minval(xg(1:NN)) 
  ymin=minval(yg(1:NN)) 
  zmin=minval(zg(1:NN)) 

  sclf = 100000.0/max((xmax-xmin),(ymax-ymin))
  sclp = 0.005 

  open(9, FILE=namEps)

  write(9, '(A24)') '%!PS-Adobe-2.0 EPSF-1.2 '
  write(9, '(A24)') '%% Created by:  TFlowS %%'
  write(9, '(A15,4I7)') '%%BoundingBox: ',                          &
    int(xmin*sclf*sclp)-1,int(ymin*sclf*sclp)-1,                    &
    int(xmax*sclf*sclp)+1,int(ymax*sclf*sclp)+1 
  write(9, '(A24)') '% Abrevations:          '
  write(9, '(A24)') '/cp   {closepath}    def'
  write(9, '(A24)') '/f    {fill}         def'
  write(9, '(A24)') '/gr   {grestore}     def'
  write(9, '(A24)') '/gs   {gsave}        def'
  write(9, '(A24)') '/l    {lineto}       def'
  write(9, '(A24)') '/m    {moveto}       def'
  write(9, '(A24)') '/np   {newpath}      def'
  write(9, '(A24)') '/s    {stroke}       def'
  write(9, '(A24)') '/sc   {scale}        def'
  write(9, '(A24)') '/sg   {setgray}      def'
  write(9, '(A24)') '/slw  {setlinewidth} def'
  write(9, '(A24)') '% End of brevations.    '

  write(9,'(A4)') ' gs '
  write(9, '(A24)') '% Scale:                '
  write(9, '(2F10.6,A4)') sclp,sclp, ' sc '

  do s=1,NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)

    x1 = xg(SideN(s,1))
    y1 = yg(SideN(s,1))
    z1 = zg(SideN(s,1))

    x2 = xg(SideN(s,2))
    y2 = yg(SideN(s,2))
    z2 = zg(SideN(s,2))

    x3 = xg(SideN(s,3))
    y3 = yg(SideN(s,3))
    z3 = zg(SideN(s,3))

    x4 = xg(SideN(s,4))
    y4 = yg(SideN(s,4))
    z4 = zg(SideN(s,4))

    count=0
! 1-2
    if(z1<=z0.and.z2>=z0 .or. z1>=z0 .and. z2<=z0) then 
      count=count+1
      fc = (z2-z0)/(z2-z1+TINY)
      xin(count) = fc*x1+(1.0-fc)*x2
      yin(count) = fc*y1+(1.0-fc)*y2
    end if
! 1-3
    if(z1<=z0 .and. z3>=z0 .or. z1>=z0 .and. z3<=z0) then 
      count=count+1
      fc = (z3-z0)/(z3-z1+TINY)
      xin(count) = fc*x1+(1.0-fc)*x3
      yin(count) = fc*y1+(1.0-fc)*y3
    end if
! 2-3
    if(z2<=z0 .and. z3>=z0 .or. z2>=z0 .and. z3<=z0) then 
      count=count+1
      fc = (z3-z0)/(z3-z2+TINY)
      xin(count) = fc*x2+(1.0-fc)*x3
      yin(count) = fc*y2+(1.0-fc)*y3
    end if
!---- for quadrilateral faces
    if(SideN(s,0)==4) then
! 1-4
      if(z1<=z0 .and. z4>=z0 .or. z1>=z0 .and. z4<=z0) then 
        count=count+1
        fc = (z4-z0)/(z4-z1+TINY)
        xin(count) = fc*x1+(1.0-fc)*x4
        yin(count) = fc*y1+(1.0-fc)*y4
      end if
! 2-4
      if(z2<=z0 .and. z4>=z0 .or. z2>=z0 .and. z4<=z0) then 
        count=count+1
        fc = (z4-z0)/(z4-z2+TINY)
        xin(count) = fc*x2+(1.0-fc)*x4
        yin(count) = fc*y2+(1.0-fc)*y4
      end if
! 3-4
      if(z3<=z0 .and. z4>=z0 .or. z3>=z0 .and. z4<=z0) then 
        count=count+1
        fc = (z4-z0)/(z4-z3+TINY)
        xin(count) = fc*x3+(1.0-fc)*x4
        yin(count) = fc*y3+(1.0-fc)*y4
      end if
    end if

!----- line width 
    lw = 0 

!----- if you want to move the domain, do something like:
!----- xin = (1.0 - xin) 

    if(count == 2) then
      write(9,'(A6,I2,A5,2I8,A3,2I8,A3,A8)')                        &
		'gs np ',lw,' slw ',                                &
		int(sclf*xin(1)),int(sclf*yin(1)), ' m ',           &
		int(sclf*xin(2)),int(sclf*yin(2)), ' l ',           &
		' cp s gr' 
      if( Dx(s) /= 0.0 .or. Dy(s) /= 0.0 ) then
	write(9,'(A6,I2,A5,2I8,A3,2I8,A3,A8)')                      &
		  'gs np ',lw,' slw ',                              &
		  int(sclf*(xin(1)-sidegDx(s))),                    &
		  int(sclf*(yin(1)-sidegDy(s))), ' m ',             &
		  int(sclf*(xin(2)-sidegDx(s))),                    &
		  int(sclf*(yin(2)-sidegDy(s))), ' l ',             &
		  ' cp s gr' 
      end if
    end if

    if(count == 3) then
      write(9,'(A6,I2,A5,2I8,A3,2I8,A3,2I8,A3,A8)')                 &
		'gs np ',lw,' slw ',                                &
		int(sclf*xin(1)),int(sclf*yin(1)), ' m ',           &
		int(sclf*xin(2)),int(sclf*yin(2)), ' l ',           &
		int(sclf*xin(3)),int(sclf*yin(3)), ' l ',           &
		' cp s gr' 
    end if

    if(count == 4) then
      write(9,'(A6,I2,A5,2I8,A3,2I8,A3,2I8,A3,2I8,A3,A8)')          &
		'gs np ',lw,' slw ',                                &
		int(sclf*xin(1)),int(sclf*yin(1)), ' m ',           &
		int(sclf*xin(2)),int(sclf*yin(2)), ' l ',           &
		int(sclf*xin(3)),int(sclf*yin(3)), ' l ',           &
		int(sclf*xin(4)),int(sclf*yin(4)), ' l ',           &
		' cp s gr' 
    end if

  end do 

  write(9,'(A4)') ' gr '

  close(9)

  END SUBROUTINE EpsSav
