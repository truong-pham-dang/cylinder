!======================================================================!
  SUBROUTINE Calc1
!----------------------------------------------------------------------!
!   Calculate node coordinates inside the domain, block by block.      !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!------------------------------[Calling]-------------------------------!
  INTEGER :: IsLine
!-----------------------------[Interface]------------------------------!
  INTERFACE 
    LOGICAL FUNCTION Approx(A,B,tol)
      IMPLICIT NONE
      REAL          :: A,B 
      REAL,OPTIONAL :: tol
    END FUNCTION Approx   
  END INTERFACE 
!-------------------------------[Locals]-------------------------------!
  INTEGER :: b, bl, i, j, k, n, c, ig
  INTEGER :: l, l1, l2
  INTEGER :: is, js, ks, ie, je, ke, face 
  INTEGER :: NI, NJ, NK, CI, CJ, CK
  INTEGER :: trans(3,2)
!--------------------------------[CVS]---------------------------------!
!  $Id: Calc1.f90,v 1.10 2002/10/30 16:29:21 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Generate/Calc1.f90,v $  
!======================================================================!

  do n=1,MAXN
    x(n)=HUGE
  end do 

!++++++++++++++++++++++++++++++++++++++++!
!     calculate the node coordinates     !
!++++++++++++++++++++++++++++++++++++++++!
  NN = 0                          ! Initialize n.o.n.
  NC = 0                          ! Initialize n.o.v.
  do b=1,Nbloc
    write(*,*) 'Generating block: ', b
    NI=BlkRes(b, 1)
    NJ=BlkRes(b, 2)
    NK=BlkRes(b, 3)   

!-( 1 )-!
    n = NN + ( 1-1)*NI*NJ + ( 1-1)*NI +  1
    x(n)=xp(BlkPnt(b, 1)) 
    y(n)=yp(BlkPnt(b, 1)) 
    z(n)=zp(BlkPnt(b, 1)) 

!-( 2 )-!
    n = NN + ( 1-1)*NI*NJ + ( 1-1)*NI + NI
    x(n)=xp(BlkPnt(b, 2)) 
    y(n)=yp(BlkPnt(b, 2)) 
    z(n)=zp(BlkPnt(b, 2)) 

!-( 3 )-!
    n = NN + ( 1-1)*NI*NJ + (NJ-1)*NI +  1
    x(n)=xp(BlkPnt(b, 3)) 
    y(n)=yp(BlkPnt(b, 3)) 
    z(n)=zp(BlkPnt(b, 3)) 

!-( 4 )-!
    n = NN + ( 1-1)*NI*NJ + (NJ-1)*NI + NI
    x(n)=xp(BlkPnt(b, 4)) 
    y(n)=yp(BlkPnt(b, 4)) 
    z(n)=zp(BlkPnt(b, 4)) 

!-( 5 )-!
    n = NN + (NK-1)*NI*NJ + ( 1-1)*NI +  1
    x(n)=xp(BlkPnt(b, 5)) 
    y(n)=yp(BlkPnt(b, 5)) 
    z(n)=zp(BlkPnt(b, 5)) 

!-( 6 )-!
    n = NN + (NK-1)*NI*NJ + ( 1-1)*NI + NI
    x(n)=xp(BlkPnt(b, 6)) 
    y(n)=yp(BlkPnt(b, 6)) 
    z(n)=zp(BlkPnt(b, 6)) 

!-( 7 )-!
    n = NN + (NK-1)*NI*NJ + (NJ-1)*NI +  1
    x(n)=xp(BlkPnt(b, 7)) 
    y(n)=yp(BlkPnt(b, 7)) 
    z(n)=zp(BlkPnt(b, 7)) 

!-( 8 )-!
    n = NN + (NK-1)*NI*NJ + (NJ-1)*NI + NI
    x(n)=xp(BlkPnt(b, 8)) 
    y(n)=yp(BlkPnt(b, 8)) 
    z(n)=zp(BlkPnt(b, 8)) 

!------------------------------!
!       First on the lines     !
!     defined point by point   !
!------------------------------!
    do l=1, Nline    

      bl = IsLine(LinPnt(l,1),LinPnt(l,2),b)

      if(bl == b) then

	do n=1,8
	  if( LinPnt(l,1) == BlkPnt(b,n)  ) l1=n
	  if( LinPnt(l,2) == BlkPnt(b,n)  ) l2=n
	end do

!..... line is defined in the +i direction
	if     ( (l2-l1) == +1 ) then
	  trans(1,1)= 0
	  trans(1,2)=+1
	  trans(2,2)= 0
	  trans(3,2)= 0
	  if( (l1 == 1).or.(l1 == 5) ) then 
	    trans(2,1)=1
	  else                         
	    trans(2,1)=BlkRes(b,2)
	  endif
	  if( (l1 == 1).or.(l1 == 3) ) then
	    trans(3,1)=1
	  else                         
	    trans(3,1)=BlkRes(b,3)
	  endif

!..... line is defined in the -i direction
	else if( (l2-l1) == -1 ) then
	  trans(1,1)=BlkRes(b,1)+1   ! NI from the block + 1
	  trans(1,2)=-1
	  trans(2,2)= 0
	  trans(3,2)= 0
	  if( (l1 == 2).or.(l1 == 6) ) then
	    trans(2,1)=1
	  else                         
	    trans(2,1)=BlkRes(b,2)
	  endif
	  if( (l1 == 2).or.(l1 == 4) ) then
	    trans(3,1)=1
	  else                         
	    trans(3,1)=BlkRes(b,3)
	  endif

!..... line is defined in the +j direction
	else if( (l2-l1) == +2 ) then
	  trans(2,1)= 0
	  trans(2,2)=+1
	  trans(1,2)= 0
	  trans(3,2)= 0 
	  if( (l1 == 1).or.(l1 == 5) ) then
	    trans(1,1)=1
	  else                         
	    trans(1,1)=BlkRes(b,1)
	  endif
	  if( (l1 == 1).or.(l1 == 2) ) then
	    trans(3,1)=1
	  else                         
	    trans(3,1)=BlkRes(b,3)
	  endif

!..... line is defined in the -j direction
	else if( (l2-l1) == -2 ) then
	  trans(2,1)=BlkRes(b,2)+1   ! NJ from the block + 1
	  trans(2,2)=-1
	  trans(1,2)= 0
	  trans(3,2)= 0 
	  if( (l1 == 3).or.(l1 == 7) ) then
	    trans(1,1)=1
	  else                         
	    trans(1,1)=BlkRes(b,1)
	  endif
	  if( (l1 == 3).or.(l1 == 4) ) then
	    trans(3,1)=1
	  else                         
	    trans(3,1)=BlkRes(b,3)
	  endif

!..... line is defined in the +k direction
	else if( (l2-l1) == +4 ) then
	  trans(3,1)= 0
	  trans(3,2)=+1
	  trans(1,2)= 0
	  trans(2,2)= 0 
	  if( (l1 == 1).or.(l1 == 3) ) then
	    trans(1,1)=1
	  else                         
	    trans(1,1)=BlkRes(b,1)
	  end if
	  if( (l1 == 1).or.(l1 == 2) ) then
	    trans(2,1)=1
	  else                         
	    trans(2,1)=BlkRes(b,2)
	  endif

!..... line is defined in the -k direction
	else if( (l2-l1) == -4 ) then
	  trans(3,1)=BlkRes(b,3) + 1  ! NK from the block + 1
	  trans(3,2)=-1
	  trans(1,2)= 0
	  trans(2,2)= 0 
	  if( (l1 == 5).or.(l1 == 7) ) then
	    trans(1,1)=1
	  else                         
	    trans(1,1)=BlkRes(b,1)
	  endif
	  if( (l1 == 5).or.(l1 == 6) ) then
	    trans(2,1)=1
	  else                         
	    trans(2,1)=BlkRes(b,2)
	  endif

	endif ! l1-l2

!----- linija je zadan tocka po tocka; stara fora
	if(LinWgt(l) ==  0.0) then
	  write(6,*) 'LINIJA: ', l
	  write(6,*) 'l1= ', l1
	  write(6,*) 'l2= ', l2
	  do ig=1,LinRes(l)
	    i=trans(1,1)+trans(1,2)*ig
	    j=trans(2,1)+trans(2,2)*ig
	    k=trans(3,1)+trans(3,2)*ig

	    n = NN + (k-1)*NI*NJ + (j-1)*NI + i
	    x(n) = xl(l,ig)
	    y(n) = yl(l,ig)
	    z(n) = zl(l,ig)
	  end do

!----- linija je zadana tezinskim faktorom; nova fora
	else
	  is=trans(1,1)+trans(1,2)
	  js=trans(2,1)+trans(2,2)
	  ks=trans(3,1)+trans(3,2)
	  ie=trans(1,1)+trans(1,2)*LinRes(l)
	  je=trans(2,1)+trans(2,2)*LinRes(l)
	  ke=trans(3,1)+trans(3,2)*LinRes(l)
	  call linija( b, LinWgt(l), is, js, ks, ie, je, ke)
	endif  

      endif ! if the block contains

    end do ! for the lines

!------------------------!
!       lines...         !
!                        !
!------------------------!
     do k=1,NK,NK-1
       do j=1,NJ,NJ-1
	 call linija( b, BlkWgt(b,1), 1,j,k,NI,j,k)
       end do
     end do

     do k=1,NK,NK-1
       do i=1,NI,NI-1
	 call linija( b, BlkWgt(b,2), i,1,k,i,NJ,k)
       end do
     end do

     do j=1,NJ,NJ-1
       do i=1,NI,NI-1
	 call linija( b, BlkWgt(b,3), i,j,1,i,j,NK)
       end do
     end do

!------------------------------------------------------------!
!       surfaces...                                          !
!                                                            !
!   I think this is the propper way to calculate surfaces:   !
!   it spans the lines in the direction of higher weigh      !
!                                                            !
!------------------------------------------------------------!

!---------------------------------------------------------------------- 
! I (k=1) 
     n = (b-1)*6 + 1   ! face index
     k=1
     if( .NOT. Approx(BlfaWt(n,1),1.0 ) ) then
       do j=1,NJ
	 call linija( b, BlfaWt(n,1), 1,j,k,NI,j,k)
       end do
     else ! lines in the j direction
       do i=1,NI
	 call linija( b, BlfaWt(n,2), i,1,k,i,NJ,k)
       end do
     endif

! VI (k=NK)
     n = (b-1)*6 + 6   ! face index
     k=NK
     if( .NOT. Approx(BlfaWt(n,1),1.0 ) ) then
       do j=1,NJ
	 call linija( b, BlFaWt(n,1), 1,j,k,NI,j,k)
       end do
     else ! lines in the j direction
       do i=1,NI
	 call linija( b, BlfaWt(n,2), i,1,k,i,NJ,k)
       end do
     endif
!---------------------------------------------------------------------- 
! V (i=1)
     n = (b-1)*6 + 5   ! face index
     i=1
     if( .NOT. Approx(BlfaWt(n,3),1.0 ) ) then
       do j=1,NJ
	 call linija( b, BlFaWt(n,3), i,j,1,i,j,NK)
       end do
     else ! lines in the j direction
       do k=1,NK
	 call linija( b, BlFaWt(n,2), i,1,k,i,NJ,k)
       end do
     end if 
! III (i=NI)
     n = (b-1)*6 + 3   ! face index
     i=NI
     if( .NOT. Approx(BlfaWt(n,3),1.0 ) ) then
       do j=1,NJ
	 call linija( b, BlFaWt(n,3), i,j,1,i,j,NK)
       end do
     else ! lines in the j direction
       do k=1,NK
	 call linija( b, BlFaWt(n,2), i,1,k,i,NJ,k)
       end do
     end if 
!---------------------------------------------------------------------- 
! II (j=1)       
     n = (b-1)*6 + 2   ! face index
     j=1
     if( .NOT. Approx(BlfaWt(n,3),1.0 ) ) then
       do i=1,NI
	 call linija( b, BlFaWt(n,3), i,j,1,i,j,NK)
       end do
     else ! lines in the i direction
       do k=1,NK
	 call linija( b, BlFaWt(n,1), 1,j,k,NI,j,k)
       end do
     endif
! IV (j=NJ)       
     n = (b-1)*6 + 4   ! face index
     j=NJ
     if( .NOT. Approx(BlfaWt(n,3),1.0 ) ) then
       do i=1,NI
	 call linija( b, BlFaWt(n,3), i,j,1,i,j,NK)
       end do
     else ! lines in the i direction
       do k=1,NK
	 call linija( b, BlFaWt(n,1), 1,j,k,NI,j,k)
       end do
     endif
!---------------------------------------------------------------------- 

!------------------------!
!       volumes...       !
!                        !
!------------------------!
     if( .NOT. Approx( BlkWgt(b,3), 1.0 ) ) then
       do i=1,NI
	 do j=1,NJ
	   call linija( b, BlkWgt(b,3), i,j,1,i,j,NK)
	 end do
       end do
     else if( .NOT. Approx( BlkWgt(b,1), 1.0 ) ) then
       do k=1,NK
	 do j=1,NJ
	   call linija( b, BlkWgt(b,1), 1,j,k,NI,j,k)
	 end do
       end do
     else if( .NOT. Approx( BlkWgt(b,2), 1.0 ) ) then
       do k=1,NK
	 do i=1,NI
	   call linija( b, BlkWgt(b,2), i,1,k,i,NJ,k)
	 end do
       end do
     else

       do i=1,NI
	 do j=1,NJ
	   do k=1,NK
	     n = NN+(k-1)*NI*NJ + (j-1)*NI + i
	     call Laplac(b, i, j, k, 0.333, 0.333, 0.334,           &
				     0.333, 0.333, 0.334,           &
				     0.333, 0.333, 0.334)
	   end do
	 end do
       end do
     endif

!--------------------------------------------!
!     set the control volume nodes (CellN)  !
!     and  the control volume neighbours     !
!--------------------------------------------!
    CI = NI-1     
    CJ = NJ-1
    CK = NK-1

    do k=1,CK
      do j=1,CJ
	do i=1,CI
	  c = NC + (k-1)*CI*CJ + (j-1)*CI + i ! cell 
	  n = NN + (k-1)*NI*NJ + (j-1)*NI + i ! 1st node

!----- nodes
	  CellN(c,1)=n
	  CellN(c,2)=n+1
	  CellN(c,3)=n+NI
	  CellN(c,4)=n+NI+1
	  CellN(c,5)=CellN(c,1)+NI*NJ
	  CellN(c,6)=CellN(c,2)+NI*NJ
	  CellN(c,7)=CellN(c,3)+NI*NJ
	  CellN(c,8)=CellN(c,4)+NI*NJ

!----- neighbours
	  CellC(c,1)= c-CI*CJ
	  CellC(c,2)= c-CI 
	  CellC(c,3)= c+1
	  CellC(c,4)= c+CI
	  CellC(c,5)= c-1
	  CellC(c,6)= c+CI*CJ

!----- this value (-1) is also the default boundary marker
	  if(i == 1)  CellC(c,5)=-1
	  if(i == CI) CellC(c,3)=-1
	  if(j == 1)  CellC(c,2)=-1
	  if(j == CJ) CellC(c,4)=-1
	  if(k == 1)  CellC(c,1)=-1
	  if(k == CK) CellC(c,6)=-1

	end do
      end do
    end do

    BlkRes(b,4) = NI*NJ*NK ! is this needed ???
    BlkRes(b,5) = NN       ! old number of nodes, for fuzion 
    BlkRes(b,6) = NC       ! old number of volumes, for fuzion
    NN = NN + NI*NJ*NK
    NC = NC + CI*CJ*CK

  end do   ! through blocks 

!---------------------------------------------!
!     Insertion of the boundary condition     ! 
!           and materials information         !
!---------------------------------------------!

!----- initialize all the material markers to 1
  do c=1,NC
    material(c)=1
  end do

  do n=1,Nboun

    b  = abs(Bound(n,7))   ! block

!---- block resolution
    CI=BlkRes(b,1)-1
    CJ=BlkRes(b,2)-1
    CK=BlkRes(b,3)-1

!---- default values
    is=1
    ie=CI
    js=1
    je=CJ
    ks=1
    ke=CK

!---- boundary conditions prescribed with mnemonics
    if(BndFac(n) == 'IMIN') then
      ie=1 
      face = 5
    else if(BndFac(n) == 'IMAX') then 
      is=CI
      face = 3
    else if(BndFac(n) == 'JMIN') then 
      je=1
      face = 2
    else if(BndFac(n) == 'JMAX') then 
      js=CJ
      face = 4
    else if(BndFac(n) == 'KMIN') then 
      ke=1
      face = 1
    else if(BndFac(n) == 'KMAX') then 
      ks=CK
      face = 6
!---- boundary conditions (materials) prescribed explicitly
!     (error prone and difficult, but might be usefull)
    else   
      is = Bound(n,1)
      js = Bound(n,2)
      ks = Bound(n,3)
      ie = Bound(n,4)
      je = Bound(n,5)
      ke = Bound(n,6)
      face = 0
      if( (is == ie).and.(is ==  1) ) face=5
      if( (is == ie).and.(is == CI) ) face=3
      if( (js == je).and.(js ==  1) ) face=2
      if( (js == je).and.(js == CJ) ) face=4
      if( (ks == ke).and.(ks ==  1) ) face=1
      if( (ks == ke).and.(ks == CK) ) face=6
    end if

    do i=is,ie
      do j=js,je
	do k=ks,ke
	  c = BlkRes(b,6) + (k-1)*CI*CJ + (j-1)*CI + i   
	  if(face /= 0) then 
	    CellC(c,face)= -Bound(n,8) ! marker
	  else
	    material(c)= Bound(n,8)    ! material
	  end if
	end do
      end do
    end do

  end do  !  Nboun

  END SUBROUTINE Calc1
