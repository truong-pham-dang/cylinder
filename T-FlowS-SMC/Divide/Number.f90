!======================================================================!
  SUBROUTINE Number
!----------------------------------------------------------------------!
!   Number the cells in each subdomain for subsequent separate saving. !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod 
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: b, c, n, s, c1, c2
  INTEGER   :: sub,subo,NNsub,NCsub,NSsub,NBCsub,NBFsub,NCSsub,NCFsub
  CHARACTER :: namOut*80
  INTEGER,ALLOCATABLE :: SideCell(:,:)
!--------------------------------[CVS]---------------------------------!
!  $Id: Number.f90,v 1.21 2002/10/30 16:29:19 niceno Exp $'/
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Divide/Number.f90,v $'/ 
!======================================================================!
!  Each subdomain needs two buffers: a send buffer and a receive buffer.
!  A receive buffer will be stored as aditional boundary cells for each
!  subdomain. So each subdomain will have NBC physical boundary faces
!  and NBBC-NBC buffer bounndary cells. It is handy to do it that way,
!  because most of the algorythms can remain the same as they are now.
!  They won't even "know" that they use values from other processors.
!  On the other hand, a sending buffer has to be allocated in a new 
!  separate array called simply buffer(). An additional array is needed 
!  to keep track of all the indexes. That one is called buffind().
!  buffind() has stored cell numbers from it's own subdomain so that
!  later they can be copied with (well, something like that):
!  do i=1,BUFFSIZ
!    buffer(i) = U(buffind(i))
!  end do
!----------------------------------------------------------------------!

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                   !
!     Browse through subdomains     !
!                                   !
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  do sub=1,Nsub

    call NamFil(sub, namOut, '.buf', len_trim('.buf'))
    open(9, FILE=namOut)
    write(6, *) 'Now creating the file:', namOut

    write(9,'(A20)') '%%%%%%%%%%%%%%%%%%%%'
    write(9,'(A20)') '%                  %'
    write(9,'(A20)') '%  Buffer indexes  %'
    write(9,'(A20)') '%                  %'
    write(9,'(A20)') '%%%%%%%%%%%%%%%%%%%%'

!---- cells
    NCsub = 0     ! number of cells in subdomain
    do c=1,NC
      NewC(c)=0
    end do
    do c=1,NC
      if(proces(c) == sub) then
	NCsub=NCsub+1
	NewC(c)=NCsub
      end if
    end do

!---- nodes
    NNsub = 0     ! number of cells in subdomain
    do n=1,NN
      NewN(n)=0
    end do
    do c=1,NC
      if(proces(c) == sub) then
	NewN(CellN(c,1))=-1
	NewN(CellN(c,2))=-1
	NewN(CellN(c,3))=-1
	NewN(CellN(c,4))=-1
	NewN(CellN(c,5))=-1
	NewN(CellN(c,6))=-1
	NewN(CellN(c,7))=-1
	NewN(CellN(c,8))=-1
      end if
    end do
    do n=1,NN
      if(NewN(n) == -1) then
	NNsub=NNsub+1
	NewN(n)=NNsub
      end if
    end do

!---- sides & real boundary cells
    NSsub  = 0     ! number of sides in subdomain
    NBCsub = 0     ! number of real boundary cells in subdomain
    NCSsub = 0
    do s=1,NS
      NewS(s)=0
    end do
    do c=-NBC,-1
      NewC(c)=0
    end do
    do s=1,NS
      c1=SideC(1,s)  
      c2=SideC(2,s) 
      if(c2  > 0) then
	if( (proces(c1) == sub).and.(proces(c2) == sub) ) then
	  NSsub=NSsub+1
	  NewS(s)=NSsub
	end if
      else ! c2 < 0
	if( proces(c1) == sub ) then
	  NSsub =NSsub+1
	  NewS(s)=NSsub   ! new number for the side
	  NBCsub=NBCsub+1
	  NewC(c2)=-NBCsub ! new number for the boundary cell
	end if
      end if 
    end do

    do s=1,Ncopy
      c1=CopyS(1,s)
      c2=CopyS(2,s)
      if( (proces(c1) == sub).and.(proces(c2) == sub) ) then
	NCSsub=NCSsub+1
      end if
    end do

    write(*,*) 'Now saving subdomain ', sub, ' with:'
    write(*,*) NCsub, ' cells'
    write(*,*) NNsub, ' nodes' 
    write(*,*) NSsub, ' sides' 
    write(*,*) NBCsub, ' physical boundary cells' 

!------------------------!
!     Create buffers     !
!------------------------!
    NBFsub = 0
    NCFsub = 0
    write(9,'(A33)') '#--- Number of physical b. cells:'
    write(9,'(I8)')  NBCsub   
    do subo=1,Nsub
      if(subo /= sub) then
	NBBs(subo)=NBFsub+1
!----- Faces inside the domain
	do s=1,NS
	  c1=SideC(1,s)  
	  c2=SideC(2,s) 
	  if(c2  > 0) then
	    if( (proces(c1) == sub).and.(proces(c2) == subo) ) then
	      NBFsub = NBFsub+1
	      BuSeIn(NBFsub)=NewC(c1) ! Buffer Send Index 
	      BuReIn(NBFsub)=c2 ! important for coordinate
	      BufPos(NBFsub)=-NBCsub-NBFsub

	      NewS(s)=NSsub+NBFsub
	    end if
	    if( (proces(c2) == sub).and.(proces(c1) == subo) ) then
	      NBFsub = NBFsub+1
	      BuSeIn(NBFsub)=NewC(c2) ! Buffer Send Index
	      BuReIn(NBFsub)=c1 ! important for coordinate
	      BufPos(NBFsub)=-NBCsub-NBFsub

	      NewS(s)=NSsub+NBFsub
	    end if
	  end if  ! c2 > 0
	end do    ! through sides
!----- Faces on the "copy" boundary
	do s=1,Ncopy
	  c1=CopyS(1,s)  
	  c2=CopyS(2,s) 
	  if( (proces(c1) == sub).and.(proces(c2) == subo) ) then
	    NBFsub = NBFsub+1
	    NCFsub = NCFsub+1
	    BuSeIn(NBFsub)=NewC(c1) ! Buffer Send Index 
	    BuReIn(NBFsub)=c2 
	    BufPos(NBFsub)= - (-NBCsub-NBFsub) ! watch out - sign
	  end if
	  if( (proces(c2) == sub).and.(proces(c1) == subo) ) then
	    NBFsub = NBFsub+1
	    NCFsub = NCFsub+1
	    BuSeIn(NBFsub)=NewC(c2) ! Buffer Send Index
	    BuReIn(NBFsub)=c1 
	    BufPos(NBFsub)= - (-NBCsub-NBFsub) ! watch out - sign
	  end if
	end do    ! through sides
	NBBe(subo)=NBFsub

!---- write to buffer file
	write(9,'(A33)') '#===============================#' 
	write(9,'(A33)') '#   Conections with subdomain:  #' 
	write(9,'(A33)') '#===============================#' 
	write(9,'(I8)')  subo 
	write(9,'(A33)') '#--- Number of local connections:'
	write(9,'(I8)')  NBBe(subo)-NBBs(subo)+1 
	write(9,'(A40)') '#--- Local number in a buffer and index:'
	do b=NBBs(subo),NBBe(subo)
	  write(9,'(2I8)') b-NBBs(subo)+1, BuSeIn(b) 
	end do
      end if 

    end do ! for subo

    call GenSav(sub, NNsub, NCsub)
    call GeoSav(sub,        NCsub, NSsub, NBCsub, NBFsub,NCFsub)
    call TestLn(sub, NNsub, NCsub, NSsub, NBCsub, NBFsub)

    write(*,*) 'Test:'
    write(*,*) 'NNsub=',NNsub
    write(*,*) 'NCsub=',NCsub
    write(*,*) 'NSsub=',NSsub
    write(*,*) 'NBCsub=',NBCsub

    write(*,*) '====================================' 
    write(*,*) 'Subdomain ', sub
    write(*,*) 'Buffer size ', NBFsub
    do subo=1,Nsub
      if(subo /= sub) then
	write(*,*) 'Connections with ', subo ,' : ',                &
	  NBBe(subo)-NBBs(subo)+1,                                  &
	  NBCsub+NBBs(subo),                                        &
	  NBCsub+NBBe(subo) 
      end if 
    end do ! for subo
    write(*,*) '------------------------------------' 

  end do   ! through subdomains

  close(9)

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!                                                      !
!     Save the entire domain with renumbered cells     !
!                                                      !
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  do n=1,NN
    NewN(n)=n
  end do
  do c=1,NC
    NewC(c)=0
  end do
  do s=1,NS
    NewS(s)=0
  end do

  NCsub = 0     ! number of cells renumbered
  do sub=1,Nsub
    do c=1,NC
      if(proces(c) == sub) then
	NCsub=NCsub+1
	NewC(c)=NCsub
      end if
    end do
  end do

  NSsub = 0     ! number of sides renumbered
  do sub=1,Nsub
    do s=1,NS
      c1 = SideC(1,s)
      c2 = SideC(2,s)
      if(proces(c1) == sub) then
        NSsub=NSsub+1
        NewS(s)=NSsub
      end if
    end do
  end do
  write(*,*) 'Number of sides: ', NS, NSsub

  call INSort(CellN(1,0), NewC(1),NC)
  call INSort(CellN(1,1), NewC(1),NC)
  call INSort(CellN(1,2), NewC(1),NC)
  call INSort(CellN(1,3), NewC(1),NC)
  call INSort(CellN(1,4), NewC(1),NC)
  call INSort(CellN(1,5), NewC(1),NC)
  call INSort(CellN(1,6), NewC(1),NC)
  call INSort(CellN(1,7), NewC(1),NC)
  call INSort(CellN(1,8), NewC(1),NC)
  call INSort(proces(1),  NewC(1),NC)
  call INSort(material(1),NewC(1),NC)

  call INSort(SideN(1,0), NewS(1),NS)
  call INSort(SideN(1,1), NewS(1),NS)
  call INSort(SideN(1,2), NewS(1),NS)
  call INSort(SideN(1,3), NewS(1),NS)
  call INSort(SideN(1,4), NewS(1),NS)
  call RNSort(Dx(1), NewS(1), NS)  ! This is important
  call RNSort(Dy(1), NewS(1), NS)  ! for plotting the
  call RNSort(Dz(1), NewS(1), NS)  ! grid with EpsPar()
  allocate(SideCell(NS,2))
  do s=1,NS
    SideCell(s,1) = SideC(1,s)
    SideCell(s,2) = SideC(2,s)
  end do
  call INSort(SideCell(1,1), NewS(1),NS)
  call INSort(SideCell(1,2), NewS(1),NS)
  do s=1,NS
    SideC(1,s) = SideCell(s,1)
    SideC(2,s) = SideCell(s,2)
  end do
  deallocate(SideCell)

  call CouMat

  call GenSav(0, NN, NC)

  call CasSav(0, NN, NC, NS+NSsh)

  call EpsPar

  END SUBROUTINE Number
