!======================================================================!
  SUBROUTINE Fuzion
!----------------------------------------------------------------------!
!   Solve the cell connectivity after block by block grid generation   !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: i, j, n                         ! counters
  INTEGER :: b1, b2                          ! block 1 and 2
  INTEGER :: f1, f2                          ! faces of block 1 and 2
  INTEGER :: n11,n12,n13,n14,n21,n22,n23,n24 ! global node numbers
  INTEGER :: l11,l12,l13,l14,l21,l22,l23,l24 ! local  node numbers
  INTEGER :: g1, g2, g3, g4                  ! generic points
  INTEGER :: i1, j1, i2, j2, k1, k2          ! directions in blocks
  INTEGER :: ig, jg, NIG, NJG                ! generic plane 
  INTEGER :: CI1, CJ1, CK1, CI2, CJ2, CK2    ! resolution of blocks
  INTEGER :: c1, c2                          ! cells from block 1, 2
  INTEGER :: NI1, NJ1, NK1, NI2, NJ2, NK2    ! resolution of blocks
  INTEGER :: n1, n2                          !  from block 1, 2
  INTEGER :: trans1(3,3), trans2(3,3)
  INTEGER :: del                             ! number of deleted nodes
!--------------------------------[CVS]---------------------------------!
!  $Id: Fuzion.f90,v 1.9 2002/10/30 16:29:21 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Generate/Fuzion.f90,v $  
!======================================================================!

!----- initialize the NewN array
  do n=1,NN
    NewN(n)=n
  end do

  if(Nbloc == 1) return 

!----- initialize the number of deleted nodes
  del=0

  NI1=0; NJ1=0; NK1=0; NI2=0; NJ2=0; NK2=0 

!---------------------------------------------------------!
!     Search through all block and all of their faces     !
!---------------------------------------------------------!
  do b2=2,Nbloc
    do b1=1,b2-1
      do f2=1,6    ! faces of the second block
	do f1=1,6  ! faces of the first block

!----- initialize the transformation matrixes             
	  do i=1,3
	    do j=1,3
	      trans1(i,j)=0
	      trans2(i,j)=0
	    end do
	  end do

	  n11=BlkFac(b1, f1, 1)
	  n12=BlkFac(b1, f1, 2)
	  n13=BlkFac(b1, f1, 3)
	  n14=BlkFac(b1, f1, 4) 
	  n21=BlkFac(b2, f2, 1)
	  n22=BlkFac(b2, f2, 2)
	  n23=BlkFac(b2, f2, 3)
	  n24=BlkFac(b2, f2, 4)

!----- check if they are connected 
	  if( ((n11 == n21).and.(n13 == n23)) .or.  &
	      ((n11 == n24).and.(n13 == n22)) .or.  &
	      ((n11 == n23).and.(n13 == n21)) .or.  &
	      ((n11 == n22).and.(n13 == n24)) ) then

!----- definiraj genericku povrsinu (g1-g4 u biti ne trebaju)
	    g1=n11
	    g2=n12
	    g3=n13
	    g4=n14

!----- nadji lokalne cvorove (1-8) blokova 1 i 2 na generickoj povrsini
	    do n=1,8
	      if(BlkPnt(b1,n) == g1) l11=n
	      if(BlkPnt(b2,n) == g1) l21=n
	      if(BlkPnt(b1,n) == g2) l12=n
	      if(BlkPnt(b2,n) == g2) l22=n
	      if(BlkPnt(b1,n) == g3) l13=n
	      if(BlkPnt(b2,n) == g3) l23=n
	      if(BlkPnt(b1,n) == g4) l14=n
	      if(BlkPnt(b2,n) == g4) l24=n
	    end do

     	  write(6, *) 'Connecting blocks: ', b1, b2
!>>>>	  write(6, *) 'f1,f2= ', f1, f2
!>>>>	  write(6,*) n11, n12, n13, n14, n21, n22, n23, n24
!>>>>	  write(6, '(4I5)') g1,  g2,  g3,  g4
!>>>>	  write(6, '(4I5)') l11, l12, l13, l14
!>>>>	  write(6, '(4I5)') l21, l22, l23, l24

!----- direction ig, block 1
	    if((l14-l11) == +1) then
	      NIG = BlkRes(b1,1)       ! NI from block 1
	      trans1(1,2)=+1
	    elseif((l14-l11) == +2) then
	      NIG = BlkRes(b1,2)       ! NJ from block 1
	      trans1(2,2)=+1
	    elseif((l14-l11) == +4) then 
	      NIG = BlkRes(b1,3)       ! NK from block 1
	      trans1(3,2)=+1
	    elseif((l14-l11) == -1) then 
	      NIG = BlkRes(b1,1)       ! NI from block 1
	      trans1(1,1)=NIG
	      trans1(1,2)=-1
	    elseif((l14-l11) == -2) then 
	      NIG = BlkRes(b1,2)       ! NJ from block 1
	      trans1(2,1)=NIG
	      trans1(2,2)=-1
	    elseif((l14-l11) == -4) then 
	      NIG = BlkRes(b1,3)       ! NK from block 1
	      trans1(3,1)=NIG
	      trans1(3,2)=-1
	    endif

!----- direction jg, block 1 
	    if((l12-l11) == +1) then 
	      NJG = BlkRes(b1,1)       ! NI from block 1
	      trans1(1,3)=+1
	    elseif((l12-l11) == +2) then
	      NJG = BlkRes(b1,2)       ! NJ from block 1
	      trans1(2,3)=+1
	    elseif((l12-l11) == +4) then
	      NJG = BlkRes(b1,3)       ! NK from block 1
	      trans1(3,3)=+1
	    elseif((l12-l11) == -1) then
	      NJG = BlkRes(b1,1)       ! NI from block 1
	      trans1(1,1)=NJG
	      trans1(1,3)=-1
	    elseif((l12-l11) == -2) then
	      NJG = BlkRes(b1,2)       ! NJ from block 1
	      trans1(2,1)=NJG
	      trans1(2,3)=-1
	    elseif((l12-l11) == -4) then
	      NJG = BlkRes(b1,3)       ! NK from block 1
	      trans1(3,1)=NJG
	      trans1(3,3)=-1
	    endif

!----- direction ig, block 2
	    if((l24-l21) == +1) then
	      NIG = BlkRes(b2,1)       ! NI from block 2
	      trans2(1,2)=+1
	    elseif((l24-l21) == +2) then
	      NIG = BlkRes(b2,2)       ! NJ from block 2
	      trans2(2,2)=+1
	    elseif((l24-l21) == +4) then 
	      NIG = BlkRes(b2,3)       ! NK from block 2
	      trans2(3,2)=+1
	    elseif((l24-l21) == -1) then 
	      NIG = BlkRes(b2,1)       ! NI from block 2
	      trans2(1,1)=NIG
	      trans2(1,2)=-1
	    elseif((l24-l21) == -2) then 
	      NIG = BlkRes(b2,2)       ! NJ from block 2
	      trans2(2,1)=NIG
	      trans2(2,2)=-1
	    elseif((l24-l21) == -4) then 
	      NIG = BlkRes(b2,3)       ! NK from block 2
	      trans2(3,1)=NIG
	      trans2(3,2)=-1
	    endif

!----- direction jg, block 2 
	    if((l22-l21) == +1) then 
	      NJG = BlkRes(b2,1)       ! NI from block 2
	      trans2(1,3)=+1
	    elseif((l22-l21) == +2) then
	      NJG = BlkRes(b2,2)       ! NJ from block 2
	      trans2(2,3)=+1
	    elseif((l22-l21) == +4) then
	      NJG = BlkRes(b2,3)       ! NK from block 2
	      trans2(3,3)=+1
	    elseif((l22-l21) == -1) then
	      NJG = BlkRes(b2,1)       ! NI from block 2
	      trans2(1,1)=NJG
	      trans2(1,3)=-1
	    elseif((l22-l21) == -2) then
	      NJG = BlkRes(b2,2)       ! NJ from block 2
	      trans2(2,1)=NJG
	      trans2(2,3)=-1
	    elseif((l22-l21) == -4) then
	      NJG = BlkRes(b2,3)       ! NK from block 2
	      trans2(3,1)=NJG
	      trans2(3,3)=-1
	    endif

!----- set the constant directions
	    if(f1 == 1) trans1(3,1)=1
	    if(f1 == 2) trans1(2,1)=1
	    if(f1 == 3) trans1(1,1)=BlkRes(b1,1)-1
	    if(f1 == 4) trans1(2,1)=BlkRes(b1,2)-1
	    if(f1 == 5) trans1(1,1)=1
	    if(f1 == 6) trans1(3,1)=BlkRes(b1,3)-1

	    if(f2 == 1) trans2(3,1)=1
	    if(f2 == 2) trans2(2,1)=1
	    if(f2 == 3) trans2(1,1)=BlkRes(b2,1)-1
	    if(f2 == 4) trans2(2,1)=BlkRes(b2,2)-1
	    if(f2 == 5) trans2(1,1)=1
	    if(f2 == 6) trans2(3,1)=BlkRes(b2,3)-1

!>>>> ispisi to sta si dobio za provjeru                  
!>>>>     write(6, *) '   C   ig   jg'
!>>>>     do i=1,3
!>>>>       write(6, '(3I5)')  &
!>>>>	     trans1(i,1), trans1(i,2), trans1(i,3)
!>>>>     end do
!>>>>
!>>>>     write(6, *) '   C   ig   jg'
!>>>>     do i=1,3
!>>>>       write(6, '(3I5)')  &
!>>>>	     trans2(i,1), trans2(i,2), trans2(i,3)
!>>>>     end do


!----- finally conect the two blocks
	    do jg=1,NJG-1              ! through volumes only
	      do ig=1,NIG-1            ! through volumes only
		CI1=BlkRes(b1,1)-1
		CJ1=BlkRes(b1,2)-1
		CK1=BlkRes(b1,3)-1
		CI2=BlkRes(b2,1)-1
		CJ2=BlkRes(b2,2)-1
		CK2=BlkRes(b2,3)-1
		i1 = trans1(1,1) + trans1(1,2)*ig + trans1(1,3)*jg
		j1 = trans1(2,1) + trans1(2,2)*ig + trans1(2,3)*jg
		k1 = trans1(3,1) + trans1(3,2)*ig + trans1(3,3)*jg 
		i2 = trans2(1,1) + trans2(1,2)*ig + trans2(1,3)*jg
		j2 = trans2(2,1) + trans2(2,2)*ig + trans2(2,3)*jg
		k2 = trans2(3,1) + trans2(3,2)*ig + trans2(3,3)*jg
		c1 = BlkRes(b1,6)  &
		     + (k1-1)*CI1*CJ1 + (j1-1)*CI1 + i1
		c2 = BlkRes(b2,6)  &
		     + (k2-1)*CI2*CJ2 + (j2-1)*CI2 + i2
!->>>               write(6, '(2I5)') c1, c2
		CellC(c1, f1) = c2
		CellC(c2, f2) = c1
	      end do
	    end do

!----- modify the transformation matrices for nodal connection
	    if(trans1(1,1)  > 1) trans1(1,1)=trans1(1,1)+1
	    if(trans1(2,1)  > 1) trans1(2,1)=trans1(2,1)+1
	    if(trans1(3,1)  > 1) trans1(3,1)=trans1(3,1)+1
	    if(trans2(1,1)  > 1) trans2(1,1)=trans2(1,1)+1
	    if(trans2(2,1)  > 1) trans2(2,1)=trans2(2,1)+1
	    if(trans2(3,1)  > 1) trans2(3,1)=trans2(3,1)+1

!----- conect the nodes 
	    do jg=1,NJG                ! through nodes 
	      do ig=1,NIG              ! through nodes
		NI1=BlkRes(b1,1)
		NJ1=BlkRes(b1,2)
		NK1=BlkRes(b1,3)
		NI2=BlkRes(b2,1)
		NJ2=BlkRes(b2,2)
		NK2=BlkRes(b2,3)
		i1 = trans1(1,1) + trans1(1,2)*ig + trans1(1,3)*jg
		j1 = trans1(2,1) + trans1(2,2)*ig + trans1(2,3)*jg
		k1 = trans1(3,1) + trans1(3,2)*ig + trans1(3,3)*jg
		i2 = trans2(1,1) + trans2(1,2)*ig + trans2(1,3)*jg
		j2 = trans2(2,1) + trans2(2,2)*ig + trans2(2,3)*jg
		k2 = trans2(3,1) + trans2(3,2)*ig + trans2(3,3)*jg
		n1 = BlkRes(b1,5)  &
		     + (k1-1)*NI1*NJ1 + (j1-1)*NI1 + i1
		n2 = BlkRes(b2,5)  &
		     + (k2-1)*NI2*NJ2 + (j2-1)*NI2 + i2
		NewN(n2) = NewN(n1)
	      end do
	    end do

	  end if  ! are they connected ? 

	end do    ! f1
      end do      ! f2
    end do        ! b1

!----- update the node numbers
    do n=BlkRes(b2,5)+1, BlkRes(b2,5)+NI2*NJ2*NK2
      if(NewN(n) /= n) del=del+1
      if(NewN(n) == n) NewN(n)=NewN(n)-del
    end do 

  end do          ! b2 


  do n=1,NN
!->>>   write(6, '(2I8)') n, NewN(n)
    x(NewN(n))=x(n)
    y(NewN(n))=y(n)
    z(NewN(n))=z(n)
  end do

  NN=NN-del

!----- skip the merged points in the node() structure
  do i=1,NC
    do n=1,8
      CellN(i,n)=NewN(CellN(i,n))
    end do
  end do

  write(6, '(I8)') del       

  END SUBROUTINE Fuzion   
