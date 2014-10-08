!======================================================================!
  SUBROUTINE GeoSav(sub, NCsub, NSsub, NBCsub, NBFsub, NCFsub)
!----------------------------------------------------------------------!
! Writes: NAME.cns, NAME.geo                                           !
! ~~~~~~~                                                              !
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: sub, NCsub, NSsub, NBCsub, NBFsub, NCFsub
!-------------------------------[Locals]-------------------------------!
  INTEGER             :: b, c, s, c1, c2, count, var, subo 
  CHARACTER           :: namOut*80
  INTEGER,ALLOCATABLE :: iwork(:,:)
  REAL,ALLOCATABLE    :: work(:)
!--------------------------------[CVS]---------------------------------!
!  $Id: GeoSav.f90,v 1.20 2002/10/31 11:26:48 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Library/GeoSav.f90,v $            
!======================================================================!
!   The files NAME.cns and NAME.geo should merge into one file in some !
!   of the future releases.                                            !
!                                                                      !
!   sub    - subdomain number                                          !
!   NCsub  - number of cells in subdomain                              !
!   NSsub  - number of sides in subdomain, but without sides on buffer !
!   NBCsub - number of physicall boundary cells in subdomain           !
!   NBFsub - number of buffer boundary faces in subdomain              !
!----------------------------------------------------------------------!

  allocate(iwork(-NbC:NS,0:2)); iwork=0
  allocate(work(NS));           work=0

!<<<<<<<<<<<<<<<<<<<<<<<<<!
!     create CNS file     !
!<<<<<<<<<<<<<<<<<<<<<<<<<!
  call NamFil( sub, namOut, '.cns', len_trim('.cns') )
  open(9, FILE=namOut,FORM='UNFORMATTED')
  write(6, *) 'Now creating the file:', namOut

!-------------------------------------------------!
!    Number of cells, boundary cells ans sides    !
!-------------------------------------------------!
  write(9) NCsub
  write(9) NBCsub+NBFsub 
  write(9) NSsub+NBFsub-NCFsub
  write(9) NSsh
  write(9) Nmat

!-----------* 
!   CELLS   * 
!-----------* 
  count=0
  do c=1,NC
    if(NewC(c) /= 0) then
      count=count+1
      iwork(count,1) = material(c)
    end if
  end do 
  write(9) (iwork(c,1), c=1,count)
!---- physicall cells
  count=0
  do c=-1,-NBC, -1
    if(NewC(c) /= 0) then
      count=count+1
      iwork(count,1) = material(c)
    end if
  end do
!---- buffer boundary cell centers
  do s=1,NBFsub
    count=count+1
    iwork(count,1) = material(BuReIn(s))
  end do
  write(9) (iwork(c,1), c=1,count)
                      
!-----------* 
!   SIDES   * 
!-----------*
  count=0

!---- NSsub physical faces
  do s=1,NS  ! OK, later chooses just sides with NewS
    if( NewS(s)  > 0  .and.  NewS(s) <= NSsub ) then
      count=count+1 
      iwork(count,0) = 0 
      iwork(count,1) = NewC(SideC(1,s))
      iwork(count,2) = NewC(SideC(2,s))
    end if
  end do 
!---- NBFsub buffer faces (copy faces here, avoid them with BufPos) 
  do s=1,NBFsub
    if(BufPos(s)  < 0) then         ! normal buffer (non-copy) 
      count=count+1 
      iwork(count,0) = BuReIn(s)    ! old cell number
      iwork(count,1) = BuSeIn(s)    ! new cell number
      iwork(count,2) = BufPos(s)    ! position in the buffer
    end if
  end do 

  write(9) (iwork(s,0), s=1,count)
  write(9) (iwork(s,1), s=1,count)
  write(9) (iwork(s,2), s=1,count)

!-----------! 
!   BOUND   !
!-----------! 
  count=0          ! count goes to negative

!---- NBCsub physical boundary cells
  do c=-1,-NbC,-1  ! OK, later chooses just cells with NewC
    if(NewC(c) /= 0) then
      count=count-1 
      ! nekad bio i: NewC(c)
      iwork(count,1) = BCmark(c)   
      iwork(count,2) = NewC(CopyC(c)) 
      if(CopyC(c) /= 0 .and. proces(CopyC(c)) /= sub) then
!	write(*,*) 'Uuups ! Sub: ', sub
	do b=1,NBFsub
	  if(BuReIn(b) == CopyC(c)) then
	    write(*,*) BufPos(b) 
	    write(*,*) xc(CopyC(c)), yc(CopyC(c)), zc(CopyC(c))  
	    iwork(count,2)=-BufPos(b) ! - sign, copy buffer
	  end if
	end do
      endif
    end if
  end do 
!---- NBFsub buffer cells
  do c=1,NBFsub
    count=count-1 
    ! nekad bio i: -NBCsub-c, 
    iwork(count,1) = BUFFER 
    iwork(count,2) = 0        ! hmm ? unused ? hmm ?
  end do 

  write(9) (iwork(c,1), c=-1,count,-1)
  write(9) (iwork(c,2), c=-1,count,-1)

!--------------!
!     COPY     !
!--------------!
  count = 0
  do s=1,Ncopy
    count = count + 1
    iwork(count,1) = CopyS(1,s) 
    iwork(count,2) = CopyS(2,s) 
  end do

  write(*,*) 'Ncopy ', count 
  write(9) count 
  write(9) (iwork(c,1), c=1,count)
  write(9) (iwork(c,2), c=1,count)

  close(9)

!<<<<<<<<<<<<<<<<<<<<<<<<<!
!     create GEO file     !
!<<<<<<<<<<<<<<<<<<<<<<<<<!
  call NamFil( sub, namOut, '.geo', len_trim('.geo') )
  open(9, FILE=namOut, FORM='UNFORMATTED')
  write(6, *) 'Now creating the file:', namOut

!---------------------------------!
!     cell center coordinates     !
!---------------------------------!
  do var=1,3
    count=0
    do c=1,NC
      if(NewC(c)  > 0) then
	count=count+1
	if(var == 1) work(count) = xc(c)
	if(var == 2) work(count) = yc(c)
	if(var == 3) work(count) = zc(c)
      end if
    end do 
    write(9) (work(c), c=1,count)
  end do

!-------------------------------!
!     boundary cell centers     !
!-------------------------------!

!---- physicall cells
  do var=1,3
    count=0
    do c=-1,-NBC, -1
      if(NewC(c) /= 0) then
	count=count+1
	if(var == 1) work(count) = xc(c)
	if(var == 2) work(count) = yc(c)
	if(var == 3) work(count) = zc(c)
      end if
    end do 
!---- buffer boundary cell centers
    do s=1,NBFsub
      count=count+1
      if(var ==  1) work(count) = xc(BuReIn(s))
      if(var ==  2) work(count) = yc(BuReIn(s))
      if(var ==  3) work(count) = zc(BuReIn(s))
    end do
    write(9) (work(c), c=1,count)
  end do

!----------------------!
!     cell volumes     !
!----------------------!
  count=0
  do c=1,NC
    if(NewC(c)  > 0) then
      count=count+1
      work(count) = volume(c)
    end if
  end do
  write(9) (work(c), c=1,count) 

!--------------------!
!     cell delta     !
!--------------------!
  count=0
  do c=1,NC
    if(NewC(c)  > 0) then
      count=count+1
      work(count) = delta(c)
    end if
  end do
  write(9) (work(c), c=1,count) 

!----------------------!
!     wall distance    !
!----------------------!
  count=0
  do c=1,NC
    if(NewC(c)  > 0) then
      count=count+1
      work(count) = WallDs(c)
    end if
  end do
  write(9) (work(c), c=1,count) 

!----------------------!
!     a1 scotti        !
!----------------------!
  count=0
  do c=1,NC
    if(NewC(c)  > 0) then
      count=count+1
      work(count) = a1(c)
    end if
  end do
  write(9) (work(c), c=1,count) 

!----------------------!
!     a2 scotti        !
!----------------------!
  count=0
  do c=1,NC
    if(NewC(c)  > 0) then
      count=count+1
      work(count) = a2(c)
    end if
  end do
  write(9) (work(c), c=1,count) 

!---------------!
!     sides     !
!---------------!

!---- from 1 to NSsub -> cell faces for which both cells are inside sub
  do var=1,10
  count=0

  do s=1,NS
    if(NewS(s)  > 0 .and. NewS(s) <= NSsub) then
      count=count+1
      if(var ==  1)  work(count) = Sx(s)
      if(var ==  2)  work(count) = Sy(s)
      if(var ==  3)  work(count) = Sz(s)
      if(var ==  4)  work(count) = Dx(s)
      if(var ==  5)  work(count) = Dy(s)
      if(var ==  6)  work(count) = Dz(s)
      if(var ==  7)  work(count) = f(s)
      if(var ==  8)  work(count) = xsp(s)
      if(var ==  9)  work(count) = ysp(s)
      if(var == 10)  work(count) = zsp(s)
    end if 
  end do

!---- from NSsub+1 to NSsub + NBFsub (think: are they in right order ?)
  do subo=1,Nsub
    do s=1,NS
      if(NewS(s)  > NSsub .and. NewS(s) <= NSsub+NBFsub) then
	c1 = SideC(1,s)
	c2 = SideC(2,s)
	if(c2  > 0) then
	  if( (proces(c1) == sub) .and. (proces(c2) == subo) ) then 
	    count=count+1
	    if(var ==  1)  work(count) = Sx(s)
	    if(var ==  2)  work(count) = Sy(s)
	    if(var ==  3)  work(count) = Sz(s)
	    if(var ==  4)  work(count) = Dx(s)
	    if(var ==  5)  work(count) = Dy(s)
	    if(var ==  6)  work(count) = Dz(s)
	    if(var ==  7)  work(count) = f(s)
            if(var ==  8)  work(count) = xsp(s)
            if(var ==  9)  work(count) = ysp(s)
            if(var == 10)  work(count) = zsp(s)
	  end if  
	  if( (proces(c2) == sub) .and. (proces(c1) == subo) ) then 
	    count=count+1
	    if(var ==  1)  work(count) = -Sx(s)
	    if(var ==  2)  work(count) = -Sy(s)
	    if(var ==  3)  work(count) = -Sz(s)
	    if(var ==  4)  work(count) = -Dx(s)
	    if(var ==  5)  work(count) = -Dy(s)
	    if(var ==  6)  work(count) = -Dz(s)
	    if(var ==  7)  work(count) = 1.0-f(s)
            if(var ==  8)  work(count) = xsp(s) - Dx(s)
            if(var ==  9)  work(count) = ysp(s) - Dy(s)
            if(var == 10)  work(count) = zsp(s) - Dz(s)
	  end if  
	end if  ! c2 > 0 
      end if    ! I think this is not really necessary 
    end do
  end do

  write(9) (work(s),s=1,count)

  end do

  close(9)

  deallocate (iwork)
  deallocate (work)

  END SUBROUTINE GeoSav
