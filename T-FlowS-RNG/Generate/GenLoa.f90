!======================================================================!
  SUBROUTINE GenLoa
!----------------------------------------------------------------------!
! Reads: NAME.d                                                        !
! ~~~~~~                                                               !
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
  USE par_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!------------------------------[Calling]-------------------------------!
  REAL :: TetVol   
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: b, i, l, s, fc, n, n1,n2,n3,n4
  INTEGER   :: NSchck, NNchck
  INTEGER   :: NI, NJ, NK
  INTEGER   :: dum
  CHARACTER :: namDom*80
  CHARACTER :: answer*12 
  REAL      :: xt(8), yt(8), zt(8)

  INTEGER   :: FaceN(6,4)
!--------------------------------[CVS]---------------------------------!
!  $Id: GenLoa.f90,v 1.19 2002/10/31 11:26:46 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Generate/GenLoa.f90,v $  
!======================================================================!
  DATA      FaceN / 1, 1, 2, 4, 3, 5,                               &
		    2, 5, 6, 8, 7, 7,                               &
		    4, 6, 8, 7, 5, 8,                               &
		    3, 2, 4, 3, 1, 6  /
!----------------------------------------------------------------------!

  write(6,'(A41)') '# Input problem name: (without extension)'
  call ReadC(5,inp,tn,ts,te) 
  read(inp, '(A80)')  name

  namDom = name
  namDom(len_trim(name)+1:len_trim(name)+2) = '.d'
  write(6, '(A24,A)') '# Now reading the file: ', namDom
  open(9, FILE=namDom)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Max. number of nodes (cells), boundary faces and cell faces     ! 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) MAXN, MAXB, MAXS  

!/////////////////////////!
!     Allocate memory     !
!/////////////////////////!

  write(6,'(A25)')       '# Allocating memory for: ' 
  write(6,'(A1,I8,A16)') '#', MAXN, ' nodes and cells' 
  write(6,'(A1,I8,A15)') '#', MAXB, ' boundary cells'         
  write(6,'(A1,I8,A11)') '#', MAXS, ' cell faces' 

!---- variables declared in all_mod.h90:
!---- (these are stored in .geo)
  allocate (xc(-MAXB:MAXN)); xc=0.0
  allocate (yc(-MAXB:MAXN)); yc=0.0
  allocate (zc(-MAXB:MAXN)); zc=0.0

  allocate (Sx(MAXS)); Sx=0.0
  allocate (Sy(MAXS)); Sy=0.0
  allocate (Sz(MAXS)); Sz=0.0

  allocate (Dx(MAXS)); Dx=0.0
  allocate (Dy(MAXS)); Dy=0.0
  allocate (Dz(MAXS)); Dz=0.0

  allocate (volume(-MAXB:MAXN)); volume=0.0
  allocate (delta(-MAXB:MAXN));  delta=0.0

  allocate (WallDs(MAXN)); WallDs=0.0
  allocate (f(MAXS)); f=0.0

!---- ()
  allocate (material(-MAXB:MAXN));  material=0
  allocate (SideC(0:2,MAXS)); SideC   =0

  allocate (CopyC(-MAXB:MAXN)); CopyC=0
  allocate (CopyS(2,MAXB));     CopyS=0    

  allocate (BCmark(-MAXB-1:-1)); BCmark=0;

!---- variables declared in gen_mod.h90:
  allocate (x(MAXN));     x=0 
  allocate (y(MAXN));     y=0
  allocate (z(MAXN));     z=0
  allocate (walln(MAXN)); walln=0
  allocate (xsp(MAXS));   xsp=0
  allocate (ysp(MAXS));   ysp=0
  allocate (zsp(MAXS));   zsp=0

  allocate (SideN(MAXS,0:4));       SideN =0 
  allocate (SideCc(MAXS,2));        SideCc=0 
  allocate (CellC(-MAXB:MAXN,24));  CellC =0

  allocate (NewN(-MAXB:MAXN));   NewN=0
  allocate (NewC(-MAXB:MAXN));   NewC=0
  allocate (NewS(MAXS));         NewS=0
  allocate (CelMar(-MAXB:MAXN)); CelMar=0

  allocate (CellN(MAXN,0:8));    CellN=0
  allocate (TwinN(MAXN,0:8));    TwinN=0

  allocate (NodeN2(MAXN,0:2));   NodeN2=0 
  allocate (NodeN4(MAXN,0:4));   NodeN4=0
  allocate (NodeN8(MAXN,0:8));   NodeN8=0

  allocate (level(MAXN)); level=0

!---- variables declared in pro_mod.h90:
  allocate (proces(MAXN)); proces=0
  allocate (BuSeIn(MAXS)); BuSeIn=0
  allocate (BuReIn(MAXS)); BuReIn=0
  allocate (BufPos(MAXS)); BufPos=0

  write(6,'(A26)') '# Allocation successfull !'

!////////////////////////!
!     Initialization     !
!////////////////////////!
  do i=1,120 
    do n=1,3
      BlkWgt(i,n) = 1.0
      BlFaWt(i,n) = 1.0
    end do
  end do

!>>>>>>>>>>>>>>>>>!
!     Corners     !
!>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)      NP     ! number of points
  if(NP  > MAXP) then
    write(6,*) 'ERROR MESSAGE FROM TFlowS:'
    write(6,*) 'You tried to define ', NP, ' points for'
    write(6,*) 'the domain, and the limit is: ', MAXP
    write(6,*) 'Increase the parameter MAXP in the file param.all'
    write(6,*) 'and recompile the code. Good Luck !'
    stop
  end if

  do i=1,NP
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(2):te(4)),*) xp(i), yp(i), zp(i)
  end do

!>>>>>>>>>>>>>>>>!
!     Blocks     !
!>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
   read(inp,*)      Nbloc        ! number of blocks 

  do b=1,Nbloc
    BlkPnt(b,0)=1       ! suppose it is properly oriented

    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(2):te(4)),*)  &  ! NI, NJ, NK  for a block
	 BlkRes(b, 1), BlkRes(b, 2), BlkRes(b, 3) 
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)  &            ! Block weights 
	 BlkWgt(b,1),BlkWgt(b,2),BlkWgt(b,3)
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)                       &
	 BlkPnt(b, 1), BlkPnt(b, 2),  &
	 BlkPnt(b, 3), BlkPnt(b, 4),  &
	 BlkPnt(b, 5), BlkPnt(b, 6),  &
	 BlkPnt(b, 7), BlkPnt(b, 8)

!-------------------------------!
!     Check if the block is     ! 
!       properly oriented       !
!-------------------------------!
    do n=1,8
      xt(n)=xp(BlkPnt(b, n))
      yt(n)=yp(BlkPnt(b, n))
      zt(n)=zp(BlkPnt(b, n))
!->>>     write(6,'(3F8.4)') xt(n),yt(n),zt(n)
    end do

    if(tetvol( xt(2),yt(2),zt(2), xt(5),yt(5),zt(5),  &
	       xt(3),yt(3),zt(3), xt(1),yt(1),zt(1) )  < 0) then
      BlkPnt(b,0)=-1            !  It's nor properly oriented
      call swapi(BlkPnt(b,2),BlkPnt(b,3))
      call swapi(BlkPnt(b,6),BlkPnt(b,7))
      call swapr(BlkWgt(b,1),BlkWgt(b,2))
      BlkWgt(b,1)=1.0/BlkWgt(b,1)
      BlkWgt(b,2)=1.0/BlkWgt(b,2)
      call swapi(BlkRes(b,1),BlkRes(b,2))
      write(6,*) 'Warning: Block ',b,' was not properly oriented'
    end if
  end do                 ! through blocks

!---------------------------------!
!     Set the corners of each     !
!        face of the block        !
!---------------------------------!
  do b=1,Nbloc
    do fc=1,6
      do n=1,4
	BlkFac(b, fc, n)=BlkPnt(b, FaceN(fc,n))
      end do
    end do
  end do

!>>>>>>>>>>>>>>>!
!     Lines     !   -->> Under construction
!>>>>>>>>>>>>>>>!--------------------------------!
!     Linije mogu biti zadan tocka po tocka,     ! 
!             ili samo tezinski faktor           !
!------------------------------------------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)       Nline     ! number of defined lines
  do l=1, Nline
    call ReadC(9,inp,tn,ts,te)

    read(inp(ts(1):te(3)),*) dum, LinPnt(l,1), LinPnt(l,2)

    call FinLin(LinPnt(l,1),LinPnt(l,2),LinRes(l))

    if(LinRes(l)  > MAXL) then
      write(6,*) 'ERROR MESSAGE FROM TFlowS:'
      write(6,*) 'You tried to define ', LinRes(l), ' points on'
      write(6,*) 'the line ', l, ' and the limit is: ', MAXL
      write(6,*) 'Increase the parameter MAXL in the file param.all'
      write(6,*) 'and recompile the code. Good Luck !'
      stop
    end if 

!----- zadana tocka po tocka
    if(dum  > 0) then
      do n=1,LinRes(l)
	call ReadC(9,inp,tn,ts,te)
	read(inp(ts(2):te(4)),*) xl(l,n), yl(l,n), zl(l,n)
!->>>
	write(6,*)  xl(l,n), yl(l,n), zl(l,n)
      end do
!----- zadan tezinski faktor
    else
      call ReadC(9,inp,tn,ts,te)
      read(inp,*) LinWgt(l)
    endif 

  end do

!--------------------------------!
!                                !
!                                !
!--------------------------------!
  do b=1,Nbloc
    do fc=1,6                          !  face of the block
      n = (b-1)*6 + fc                 !  surface number
      BlFaWt(n,1)=BlkWgt(b,1)
      BlFaWt(n,2)=BlkWgt(b,2)
      BlFaWt(n,3)=BlkWgt(b,3)
      BlFaLa(n)  = NO
    end do
  end do

!>>>>>>>>>>>>>>>>>>!
!     Surfaces     !  -->> Under construction
!>>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)     Nsurf     ! number of defined surfaces

  do s=1,Nsurf
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) dum, n1,n2,n3,n4
    call FinSur(n1,n2,n3,n4,b,fc)
    write(6,*) 'block: ', b, ' surf: ', fc
    n = (b-1)*6 + fc         ! surface number
    BlFaLa(n) = YES          ! perform Laplace

    call ReadC(9,inp,tn,ts,te)
    read(inp,*)  BlFaWt(n,1),BlFaWt(n,2),BlFaWt(n,3)
  end do


!??????????????????????????????????????????!
!     Is there enough allocated memory     !
!??????????????????????????????????????????!

!----- Nodes & Sides
  NNchck = 0
  NSchck = 0
  do b=1,Nbloc
    NI=BlkRes(b,1)
    NJ=BlkRes(b,2)
    NK=BlkRes(b,3)
    NNchck=NNchck + NI*NJ*NK
    NSchck=NSchck + NI*NJ*NK + 2*( (NI*NJ)+(NJ*NK)+(NI*NK) )
  end do

  if( (NSchck  > MAXS).or.(NNchck  > MAXN) ) then
    write(6,*) 'ERROR MESSAGE FROM TFlowS:'
  end if

  if( NSchck  > MAXS ) then
    write(6,*) 'The estimated number of sides is :', NSchck
    write(6,*) 'There is space available only for:', MAXS
    write(6,*) 'Increase the parameter MAXS in the input file'
    write(6,*) 'and re-run the code !'
  end if

  if( NNchck  > MAXN ) then
    write(6,*) 'The estimated number of nodes is :', NNchck
    write(6,*) 'There is space available only for:', MAXN
    write(6,*) 'Increase the parameter MAXN in the input file'
    write(6,*) 'and re-run the code !'
  end if 

  if( (NSchck  > MAXS).or.(NNchck  > MAXN) ) then
    stop
  end if

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Boundary conditions and materials     !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)     Nboun      ! number of boundary conditions 

  do n=1,Nboun
    BndFac(n)=''
    call ReadC(9,inp,tn,ts,te)
    if(tn == 7) then
      read(inp,*)  dum,                         &  
	   Bound(n,1), Bound(n,2), Bound(n,3),  &  ! is, js, ks 
	   Bound(n,4), Bound(n,5), Bound(n,6)      ! ie, je, ke
    else if(tn == 2) then
      read(inp(ts(1):te(1)),*)       dum           
      read(inp(ts(2):te(2)),'(A4)')  BndFac(n)
      call ToUppr(BndFac(n))           
    end if
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)       &  
	 Bound(n,7),  &  ! block,  
	 Bound(n,8)      ! mark         
  if( BlkPnt(Bound(n,7),0) == -1 ) then
    call swapi( Bound(n,1),Bound(n,2) )
    call swapi( Bound(n,4),Bound(n,5) )
  end if

  end do

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Periodic boundaries     !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)        Nperi      ! number of periodic boundaries

  do n=1,Nperi
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) dum,Period(n,1),Period(n,2),Period(n,3),Period(n,4) 
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)     Period(n,5),Period(n,6),Period(n,7),Period(n,8)
  end do

!>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Copy boundaries     !
!>>>>>>>>>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)        Ncopy      ! number of periodic boundaries

  do n=1,Ncopy
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) dum,Copy(n,1),Copy(n,2),Copy(n,3),Copy(n,4) 
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)     Copy(n,5),Copy(n,6),Copy(n,7),Copy(n,8)
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)     Copy(n,0)
  end do

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Refinement levels and regions     !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)        NRL ! number of refinement levels

  write(*,*) 'Number of refinement levels: ', NRL 

  do l=1,NRL
    call ReadC(9,inp,tn,ts,te)
    inp = inp(ts(2):te(2)) ! this is the only way Microsoft Fortran
                           ! compiles this part (two lines) of the code
    read(inp,*) NR(l)      ! number of regions in level n

    do n=1, NR(l)
      call ReadC(9,inp,tn,ts,te)
      read(inp(ts(3):te(3)),*) answer
      call ToUppr(answer)
      if(answer == 'RECTANGLE') then
	Fregio(l,n,0) = RECTAN
      elseif(answer == 'ELIPSOID') then
	Fregio(l,n,0) = ELIPSO 
      elseif(answer == 'PLANE') then
	Fregio(l,n,0) = PLANE
      else
	write(*,*) 'Error in input file: ', answer 
	stop
      endif 

      call ReadC(9,inp,tn,ts,te)
      read(inp,*)                                                   &
		  Fregio(l,n,1),Fregio(l,n,2),Fregio(l,n,3),        &
		  Fregio(l,n,4),Fregio(l,n,5),Fregio(l,n,6)   
    end do
  end do

!>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Smoothing regions     !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) NSR  ! number of smoothing regions 

  write(*,*) 'Number of (non)smoothing regions: ', NSR 

  do n=1, NSR
    SdirX(n) = .FALSE.
    SdirY(n) = .FALSE.
    SdirZ(n) = .FALSE.
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) Sregio(n,0)  
    if(tn == 4) then   ! smoothing in three directions
      SdirX(n) = .TRUE.
      SdirY(n) = .TRUE.
      SdirZ(n) = .TRUE.
    else if(tn == 3) then
      call ToUppr(inp(ts(2):te(2)))
      call ToUppr(inp(ts(3):te(3)))
      if( inp(ts(2):te(2))  ==  'X' ) SdirX(n) = .TRUE.
      if( inp(ts(3):te(3))  ==  'X' ) SdirX(n) = .TRUE.
      if( inp(ts(2):te(2))  ==  'Y' ) SdirY(n) = .TRUE.
      if( inp(ts(3):te(3))  ==  'Y' ) SdirY(n) = .TRUE.
      if( inp(ts(2):te(2))  ==  'Z' ) SdirZ(n) = .TRUE.
      if( inp(ts(3):te(3))  ==  'Z' ) SdirZ(n) = .TRUE.
    else if(tn == 2) then
      call ToUppr(inp(ts(2):te(2)))
      if( inp(ts(2):te(2))  ==  'X' ) SdirX(n) = .TRUE.
      if( inp(ts(2):te(2))  ==  'Y' ) SdirY(n) = .TRUE.
      if( inp(ts(2):te(2))  ==  'Z' ) SdirZ(n) = .TRUE.
    end if 

!---- read the coordinates of the (non)smoothed region
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) Siter(n), Srelax(n)
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)                                                     &
		Sregio(n,1),Sregio(n,2),Sregio(n,3),                &
		Sregio(n,4),Sregio(n,5),Sregio(n,6)   
  end do

  close(9)

  END SUBROUTINE GenLoa
