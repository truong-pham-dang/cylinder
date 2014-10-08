!======================================================================!
  SUBROUTINE CopyBC
!----------------------------------------------------------------------!
!   Solve the cell connectivity for copy boundary conditions.          !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: i, j, n, p, Ncop                ! counters
  INTEGER :: b1, b2                          ! block 1 and 2
  INTEGER :: f1, f2                          ! faces of block 1 and 2
  INTEGER :: n11,n12,n13,n14,n21,n22,n23,n24 ! global node numbers
  INTEGER :: p11,p12,p13,p14,p21,p22,p23,p24 ! global node numbers
  INTEGER :: l11,l12,l13,l14,l21,l22,l23,l24 ! local  node numbers
  INTEGER :: i1, j1, i2, j2, k1, k2          ! directions in blocks
  INTEGER :: ig, jg, NIG, NJG                ! generic plane 
  INTEGER :: CI1, CJ1, CK1, CI2, CJ2, CK2    ! resolution of blocks
  INTEGER :: c1, c2                          ! cells from block 1, 2
  INTEGER :: trans1(3,3), trans2(3,3)
!--------------------------------[CVS]---------------------------------!
!  $Id: CopyBC.f90,v 1.7 2002/10/30 16:29:21 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Generate/CopyBC.f90,v $  
!======================================================================!

!---------------------------------------------------------!
!     Search through all block and all of their faces     !
!---------------------------------------------------------!
  Ncop = 0       

  do p=1, Ncopy    
    do b2=1,Nbloc
      do b1=1,Nbloc
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

	    p11=Copy(p, 1)
	    p12=Copy(p, 2)
	    p13=Copy(p, 3)
	    p14=Copy(p, 4) 
	    p21=Copy(p, 5)
	    p22=Copy(p, 6)
	    p23=Copy(p, 7)
	    p24=Copy(p, 8)

!----- check if they are connected 
	  if( ( ((n11 == p11).and.(n13 == p13)) .or.                &
		((n11 == p14).and.(n13 == p12)) .or.                &
		((n11 == p13).and.(n13 == p11)) .or.                &
		((n11 == p12).and.(n13 == p14)) )                   &
			     .and.                                  &
	      ( ((n21 == p21).and.(n23 == p23)) .or.                &
		((n21 == p24).and.(n23 == p22)) .or.                &
		((n21 == p23).and.(n23 == p21)) .or.                &
		((n21 == p22).and.(n23 == p24)) ) ) then

!----- nadji lokalne cvorove (1-8) blokova 1 i 2 na generickoj povrsini
	      do n=1,8
		if(BlkPnt(b1,n) == p11) l11=n
		if(BlkPnt(b1,n) == p12) l12=n
		if(BlkPnt(b1,n) == p13) l13=n
		if(BlkPnt(b1,n) == p14) l14=n
		if(BlkPnt(b2,n) == p21) l21=n
		if(BlkPnt(b2,n) == p22) l22=n
		if(BlkPnt(b2,n) == p23) l23=n
		if(BlkPnt(b2,n) == p24) l24=n
	      end do

	  write(6, *) 'Copy blocks: ', b1, b2
!>>>>	  write(6, '(2I2)') f1, f2
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
!>>>>	     write(6, *) '   C   ig   jg'
!>>>>	     do i=1,3
!>>>>	       write(6, '(3I5)')                                    &
!>>>>		     trans1(i,1), trans1(i,2), trans1(i,3)
!>>>>	     end do
!>>>>
!>>>>	     write(6, *) '   C   ig   jg'
!>>>>	     do i=1,3
!>>>>	       write(6, '(3I5)')                                    &
!>>>>		     trans2(i,1), trans2(i,2), trans2(i,3)
!>>>>	     end do


!----- finally conect the two copy boundaries
!      now, connect the cells in the interior
	      do jg=1,NJG-1              ! through volumes only
		do ig=1,NIG-1            ! through volumes only
		  Ncop = Ncop + 1
		  CI1=BlkRes(b1,1)-1
		  CJ1=BlkRes(b1,2)-1
		  CK1=BlkRes(b1,3)-1
		  CI2=BlkRes(b2,1)-1
		  CJ2=BlkRes(b2,2)-1
		  CK2=BlkRes(b2,3)-1
		  i1 = trans1(1,1)+trans1(1,2)*ig+trans1(1,3)*jg
		  j1 = trans1(2,1)+trans1(2,2)*ig+trans1(2,3)*jg
		  k1 = trans1(3,1)+trans1(3,2)*ig+trans1(3,3)*jg 
		  i2 = trans2(1,1)+trans2(1,2)*ig+trans2(1,3)*jg
		  j2 = trans2(2,1)+trans2(2,2)*ig+trans2(2,3)*jg
		  k2 = trans2(3,1)+trans2(3,2)*ig+trans2(3,3)*jg
		  c1 = BlkRes(b1,6)                                 &
		       + (k1-1)*CI1*CJ1 + (j1-1)*CI1 + i1
		  c2 = BlkRes(b2,6)                                 &
		       + (k2-1)*CI2*CJ2 + (j2-1)*CI2 + i2
!>>>>>         write(*,*) 'Copy ', c1, ' to ', c2
		  CopyC(c2) = c1 ! allway copy from c1 to c2 !
		  CopyS(1, Ncop) = c1
		  CopyS(2, Ncop) = c2
		end do
	      end do        

	    end if  ! are they connected ? 

	  end do    ! f1
	end do      ! f2
      end do        ! b1
    end do          ! b2 
  end do            ! p copies 

  Ncopy = Ncop

  END SUBROUTINE CopyBC   
