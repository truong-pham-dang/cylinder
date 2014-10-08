!======================================================================!
  SUBROUTINE Refine(lev)
!----------------------------------------------------------------------!
!   Refine the marked cells.                                           !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: lev 
!------------------------------[Calling]-------------------------------!
  LOGICAL :: IsTwin
  INTEGER :: WchNod
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, NCold, c1, c2, c3, c4, c5, c6
  INTEGER :: cr1, cr2, cr3, cr4, cr5, cr6, cr7, cr8
  INTEGER :: n, NNold, n1, n2, n3, n4, n5, n6, n7, n8
  INTEGER :: n12,n13,n24,n34,n15,n26,n37,n48,n56,n57,n68,n78
  INTEGER :: nF1, nF2, nF3, nF4, nF5, nF6, n0
  INTEGER :: del   ! number of deleted nodes 
  INTEGER :: nA, nA0, nA1, nA2, nB, nB0, nB1, nB2
!--------------------------------[CVS]---------------------------------!
!  $Id: Refine.f90,v 1.10 2002/10/30 16:29:22 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Generate/Refine.f90,v $   
!======================================================================!
!                                                                      !
!                               c6      c3                             !
!                               |      /                               !
!                         8-----|---------6                            !
!                        /  cr8 |/  cr6  /|                            !
!                       /-------+-------/ |                            !
!                      /       /       /| |                            !
!                     7-------+-------5 | |                            !
!                     |       |       | |/|                            !
!                c4---|  cr7  |  cr5  | +-----c2                       !
!                     |       |       |/| |                            !
!                     +-------+-------+ | 2                            !
!                     |      /|       | |/                             !
!                     |  cr3/ |  cr1  | /                              !
!                     |    /  |       |/                               !
!                     3---c5--+-------1                                !
!                               |                                      !
!                               c1                                     !
!                                                                      !
!----------------------------------------------------------------------!

  write(*,*) 'Refine: Number of nodes: ', NN 
  write(*,*) '        Number of cells: ', NC 

  NCold = NC 
  NNold = NN
  NN2 = 0 
  NN4 = 0
  NN8 = 0

!================================!
!                                !
!  Najprije pobroji nove celije  !
!                                !
!================================!
  do c=1,NCold
    if(CelMar(c)  ==  -1) then 
      NC = NC + 8
      CelMar(c) = NC   ! now points to cr8
    end if
  end do

  do c=1,NCold

    c1=CellC(c,1)
    c2=CellC(c,2)
    c3=CellC(c,3)
    c4=CellC(c,4)
    c5=CellC(c,5)
    c6=CellC(c,6)

    n1=CellN(c,1)
    n2=CellN(c,2)
    n3=CellN(c,3)
    n4=CellN(c,4)
    n5=CellN(c,5)
    n6=CellN(c,6)
    n7=CellN(c,7)
    n8=CellN(c,8)

!==========================!
!                          !
!     Usitnjene celije     !
!                          !
!==========================!
    if(CelMar(c)   >  0) then  ! -> samo usitnjeni 

!-------------------------!
!                         !
!     Najprije  sredi     !
!     susjedne celije     !
!                         !
!-------------------------!
      cr1=CelMar(c)-7
      cr2=CelMar(c)-6
      cr3=CelMar(c)-5
      cr4=CelMar(c)-4
      cr5=CelMar(c)-3
      cr6=CelMar(c)-2
      cr7=CelMar(c)-1
      cr8=CelMar(c)

      material(cr1) = material(c) 
      material(cr2) = material(c) 
      material(cr3) = material(c) 
      material(cr4) = material(c) 
      material(cr5) = material(c) 
      material(cr6) = material(c) 
      material(cr7) = material(c) 
      material(cr8) = material(c) 

!--------------------------------------------!
!     Interne veze ne ovise o susjedima      !
!--------------------------------------------!

!----- 6
      CellC(cr1,6) = cr5
      CellC(cr2,6) = cr6
      CellC(cr3,6) = cr7
      CellC(cr4,6) = cr8
!----- 5
      CellC(cr2,5) = cr1
      CellC(cr4,5) = cr3
      CellC(cr6,5) = cr5
      CellC(cr8,5) = cr7
!----- 4
      CellC(cr1,4) = cr3
      CellC(cr2,4) = cr4
      CellC(cr5,4) = cr7
      CellC(cr6,4) = cr8
!----- 3
      CellC(cr1,3) = cr2
      CellC(cr3,3) = cr4
      CellC(cr5,3) = cr6
      CellC(cr7,3) = cr8
!----- 2
      CellC(cr3,2) = cr1
      CellC(cr4,2) = cr2
      CellC(cr7,2) = cr5
      CellC(cr8,2) = cr6
!----- 1
      CellC(cr5,1) = cr1
      CellC(cr6,1) = cr2
      CellC(cr7,1) = cr3
      CellC(cr8,1) = cr4

!-----------------------------!
!     Level of refinement     !
!-----------------------------!
      level(cr1) = lev
      level(cr2) = lev
      level(cr3) = lev
      level(cr4) = lev
      level(cr5) = lev
      level(cr6) = lev
      level(cr7) = lev
      level(cr8) = lev

!-------------------------------------------------!         
!     Vanjske veze ovise o topologiji susjeda     !
!-------------------------------------------------!
      if(CelMar(c1)  ==  0) then  ! -> susjed 1 neusitnjen
	CellC(cr1,1) = c1
	CellC(cr2,1) = c1
	CellC(cr3,1) = c1
	CellC(cr4,1) = c1
      else                           ! -> susjed 1 usitnjen
	CellC(cr1,1) = CelMar(c1) - 8 + WchNod(c1,n1)
	CellC(cr2,1) = CelMar(c1) - 8 + WchNod(c1,n2)
	CellC(cr3,1) = CelMar(c1) - 8 + WchNod(c1,n3)
	CellC(cr4,1) = CelMar(c1) - 8 + WchNod(c1,n4)
      endif           

      if(CelMar(c2)  ==  0) then  ! -> susjed 2 neusitnjen
	CellC(cr1,2) = c2
	CellC(cr2,2) = c2
	CellC(cr5,2) = c2
	CellC(cr6,2) = c2
      else                           ! -> susjed 2 usitnjen
	CellC(cr1,2) = CelMar(c2) - 8 + WchNod(c2, n1)
	CellC(cr2,2) = CelMar(c2) - 8 + WchNod(c2, n2)
	CellC(cr5,2) = CelMar(c2) - 8 + WchNod(c2, n5)
	CellC(cr6,2) = CelMar(c2) - 8 + WchNod(c2, n6)
      endif           

      if(CelMar(c3)  ==  0) then  ! -> susjed 3 neusitnjen
	CellC(cr2,3) = c3
	CellC(cr4,3) = c3
	CellC(cr6,3) = c3
	CellC(cr8,3) = c3
      else                           ! -> susjed 3 usitnjen
	CellC(cr2,3) = CelMar(c3) - 8 + WchNod(c3, n2)
	CellC(cr4,3) = CelMar(c3) - 8 + WchNod(c3, n4)
	CellC(cr6,3) = CelMar(c3) - 8 + WchNod(c3, n6)
	CellC(cr8,3) = CelMar(c3) - 8 + WchNod(c3, n8)
      endif           

      if(CelMar(c4)  ==  0) then  ! -> susjed 4 neusitnjen
	CellC(cr3,4) = c4
	CellC(cr4,4) = c4
	CellC(cr7,4) = c4
	CellC(cr8,4) = c4
      else                           ! -> susjed 4 usitnjen
	CellC(cr3,4) = CelMar(c4) - 8 + WchNod(c4, n3)
	CellC(cr4,4) = CelMar(c4) - 8 + WchNod(c4, n4)
	CellC(cr7,4) = CelMar(c4) - 8 + WchNod(c4, n7)
	CellC(cr8,4) = CelMar(c4) - 8 + WchNod(c4, n8)
      endif           

      if(CelMar(c5)  ==  0) then  ! -> susjed 5 neusitnjen
	CellC(cr1,5) = c5
	CellC(cr3,5) = c5
	CellC(cr5,5) = c5
	CellC(cr7,5) = c5
      else                           ! -> susjed 5 usitnjen
	CellC(cr1,5) = CelMar(c5) - 8 + WchNod(c5, n1)
	CellC(cr3,5) = CelMar(c5) - 8 + WchNod(c5, n3)
	CellC(cr5,5) = CelMar(c5) - 8 + WchNod(c5, n5)
	CellC(cr7,5) = CelMar(c5) - 8 + WchNod(c5, n7)
      endif

      if(CelMar(c6)  ==  0) then  ! -> susjed 6 neusitnjen
	CellC(cr5,6) = c6
	CellC(cr6,6) = c6
	CellC(cr7,6) = c6
	CellC(cr8,6) = c6
      else                           ! -> susjed 6 usitnjen
	CellC(cr5,6) = CelMar(c6) - 8 + WchNod(c6, n5)
	CellC(cr6,6) = CelMar(c6) - 8 + WchNod(c6, n6)
	CellC(cr7,6) = CelMar(c6) - 8 + WchNod(c6, n7)
	CellC(cr8,6) = CelMar(c6) - 8 + WchNod(c6, n8)
      endif           

!-----------------------!
!                       !
!      Onda odredi      !
!     nove  cvorove     !
!                       !
!-----------------------!

!-----------------------------------!
!     Prvo cvorove na bridovima     !
!-----------------------------------!
    n12 = 0     !
    n13 = 0     !         8-----n68-----6
    n24 = 0     !        /|            /|
    n34 = 0     !      n78|          n56|
    n15 = 0     !      / n48         / n26
    n26 = 0     !     7-----n57-----5   |
    n37 = 0     !     |   |         |   |
    n48 = 0     !     |   4- - -n24-| - 2
    n56 = 0     !    n37 /         n15 / 
    n57 = 0     !     |n34          |n12
    n68 = 0     !     |/            |/
    n78 = 0     !     3-----n13-----1

!--- n12 ---!
    do n=1,NN2
      if( ( NodeN2(n,1) == n1 .and. NodeN2(n,2) == n2 ) .or.  &
	  ( NodeN2(n,1) == n2 .and. NodeN2(n,2) == n1) ) then
	n12 = NodeN2(n,0) 
      end if
    end do
    if (n12 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n12 = NN
      NodeN2(NN2,0) = n12
      NodeN2(NN2,1) = n1
      NodeN2(NN2,2) = n2
      x(n12) = 0.5 * (x(n1)+x(n2))
      y(n12) = 0.5 * (y(n1)+y(n2))
      z(n12) = 0.5 * (z(n1)+z(n2))
    end if 

!--- n13 ---!
    do n=1,NN2
      if( ( NodeN2(n,1) == n1 .and. NodeN2(n,2) == n3 ) .or.  &
	  ( NodeN2(n,1) == n3 .and. NodeN2(n,2) == n1) ) then
	n13 = NodeN2(n,0) 
      end if
    end do
    if (n13 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n13 = NN
      NodeN2(NN2,0) = n13
      NodeN2(NN2,1) = n1
      NodeN2(NN2,2) = n3
      x(n13) = 0.5 * (x(n1)+x(n3))
      y(n13) = 0.5 * (y(n1)+y(n3))
      z(n13) = 0.5 * (z(n1)+z(n3))
    end if 

!--- n24 ---!
    do n=1,NN2
      if( ( NodeN2(n,1) == n2 .and. NodeN2(n,2) == n4 ) .or.  &
	  ( NodeN2(n,1) == n4 .and. NodeN2(n,2) == n2) ) then
	n24 = NodeN2(n,0) 
      end if
    end do
    if (n24 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n24 = NN
      NodeN2(NN2,0) = n24
      NodeN2(NN2,1) = n2
      NodeN2(NN2,2) = n4
      x(n24) = 0.5 * (x(n2)+x(n4))
      y(n24) = 0.5 * (y(n2)+y(n4))
      z(n24) = 0.5 * (z(n2)+z(n4))
    end if 

!--- n34 ---!
    do n=1,NN2
      if( ( NodeN2(n,1) == n3 .and. NodeN2(n,2) == n4 ) .or.  &
	  ( NodeN2(n,1) == n4 .and. NodeN2(n,2) == n3) ) then
	n34 = NodeN2(n,0) 
      end if
    end do
    if (n34 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n34 = NN
      NodeN2(NN2,0) = n34
      NodeN2(NN2,1) = n3
      NodeN2(NN2,2) = n4
      x(n34) = 0.5 * (x(n3)+x(n4))
      y(n34) = 0.5 * (y(n3)+y(n4))
      z(n34) = 0.5 * (z(n3)+z(n4))
    end if 

!--- n15 ---!
    do n=1,NN2
      if( ( NodeN2(n,1) == n1 .and. NodeN2(n,2) == n5 ) .or.  &
	  ( NodeN2(n,1) == n5 .and. NodeN2(n,2) == n1) ) then
	n15 = NodeN2(n,0) 
      end if
    end do
    if (n15 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n15 = NN
      NodeN2(NN2,0) = n15
      NodeN2(NN2,1) = n1
      NodeN2(NN2,2) = n5
      x(n15) = 0.5 * (x(n1)+x(n5))
      y(n15) = 0.5 * (y(n1)+y(n5))
      z(n15) = 0.5 * (z(n1)+z(n5))
    end if 

!--- n26 ---!
    do n=1,NN2
      if( ( NodeN2(n,1) == n2 .and. NodeN2(n,2) == n6 ) .or.  &
	  ( NodeN2(n,1) == n6 .and. NodeN2(n,2) == n2) ) then
	n26 = NodeN2(n,0) 
      end if
    end do
    if (n26 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n26 = NN
      NodeN2(NN2,0) = n26
      NodeN2(NN2,1) = n2
      NodeN2(NN2,2) = n6
      x(n26) = 0.5 * (x(n2)+x(n6))
      y(n26) = 0.5 * (y(n2)+y(n6))
      z(n26) = 0.5 * (z(n2)+z(n6))
    end if 

!--- n37 ---!
    do n=1,NN2
      if( ( NodeN2(n,1) == n3 .and. NodeN2(n,2) == n7 ) .or.  &
	  ( NodeN2(n,1) == n7 .and. NodeN2(n,2) == n3) ) then
	n37 = NodeN2(n,0) 
      end if
    end do
    if (n37 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n37 = NN
      NodeN2(NN2,0) = n37
      NodeN2(NN2,1) = n3
      NodeN2(NN2,2) = n7
      x(n37) = 0.5 * (x(n3)+x(n7))
      y(n37) = 0.5 * (y(n3)+y(n7))
      z(n37) = 0.5 * (z(n3)+z(n7))
    end if 

!--- n48 ---!
    do n=1,NN2
      if( ( NodeN2(n,1) == n4 .and. NodeN2(n,2) == n8 ) .or.  &
	  ( NodeN2(n,1) == n8 .and. NodeN2(n,2) == n4) ) then
	n48 = NodeN2(n,0) 
      end if
    end do
    if (n48 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n48 = NN
      NodeN2(NN2,0) = n48
      NodeN2(NN2,1) = n4
      NodeN2(NN2,2) = n8
      x(n48) = 0.5 * (x(n4)+x(n8))
      y(n48) = 0.5 * (y(n4)+y(n8))
      z(n48) = 0.5 * (z(n4)+z(n8))
    end if 

!--- n56 ---!
    do n=1,NN2
      if( ( NodeN2(n,1) == n5 .and. NodeN2(n,2) == n6 ) .or.  &
	  ( NodeN2(n,1) == n6 .and. NodeN2(n,2) == n5) ) then
	n56 = NodeN2(n,0) 
      end if
    end do
    if (n56 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n56 = NN
      NodeN2(NN2,0) = n56
      NodeN2(NN2,1) = n5
      NodeN2(NN2,2) = n6
      x(n56) = 0.5 * (x(n5)+x(n6))
      y(n56) = 0.5 * (y(n5)+y(n6))
      z(n56) = 0.5 * (z(n5)+z(n6))
    end if 

!--- n57 ---!
    do n=1,NN2
      if( ( NodeN2(n,1) == n5 .and. NodeN2(n,2) == n7 ) .or.  &
	  ( NodeN2(n,1) == n7 .and. NodeN2(n,2) == n5) ) then
	n57 = NodeN2(n,0) 
      end if
    end do
    if (n57 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n57 = NN
      NodeN2(NN2,0) = n57
      NodeN2(NN2,1) = n5
      NodeN2(NN2,2) = n7
      x(n57) = 0.5 * (x(n5)+x(n7))
      y(n57) = 0.5 * (y(n5)+y(n7))
      z(n57) = 0.5 * (z(n5)+z(n7))
    end if 

!--- n68 ---!
    do n=1,NN2
      if( ( NodeN2(n,1) == n6 .and. NodeN2(n,2) == n8 ) .or.  &
	  ( NodeN2(n,1) == n8 .and. NodeN2(n,2) == n6) ) then
	n68 = NodeN2(n,0) 
      end if
    end do
    if (n68 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n68 = NN
      NodeN2(NN2,0) = n68
      NodeN2(NN2,1) = n6
      NodeN2(NN2,2) = n8
      x(n68) = 0.5 * (x(n6)+x(n8))
      y(n68) = 0.5 * (y(n6)+y(n8))
      z(n68) = 0.5 * (z(n6)+z(n8))
    end if 

!--- n78 ---!
    do n=1,NN2
      if( ( NodeN2(n,1) == n7 .and. NodeN2(n,2) == n8 ) .or.  &
	  ( NodeN2(n,1) == n8 .and. NodeN2(n,2) == n7) ) then
	n78 = NodeN2(n,0) 
      end if
    end do
    if (n78 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n78 = NN
      NodeN2(NN2,0) = n78
      NodeN2(NN2,1) = n7
      NodeN2(NN2,2) = n8
      x(n78) = 0.5 * (x(n7)+x(n8))
      y(n78) = 0.5 * (y(n7)+y(n8))
      z(n78) = 0.5 * (z(n7)+z(n8))
    end if 

!------------------------------------!
!     Onda cvorove na stranicama     !
!------------------------------------!
  nF1 = 0
  nF2 = 0
  nF3 = 0
  nF4 = 0
  nF5 = 0
  nF6 = 0

!--- nF1 ---!
    do n=1,NN4
      if( ( NodeN4(n,1) == n1 .and. NodeN4(n,4) == n4 ) .or.  &
	  ( NodeN4(n,1) == n4 .and. NodeN4(n,4) == n1 ) .or.  &
          ( NodeN4(n,1) == n2 .and. NodeN4(n,4) == n3 ) .or.  &
	  ( NodeN4(n,1) == n3 .and. NodeN4(n,4) == n2 ) ) then
	nF1 = NodeN4(n,0) 
      end if
    end do
    if (nF1 == 0) then
      NN4 = NN4 + 1
      NN  = NN  + 1
      nF1 = NN
      NodeN4(NN4,0) = nF1
      NodeN4(NN4,1) = n1
      NodeN4(NN4,2) = n2
      NodeN4(NN4,3) = n3
      NodeN4(NN4,4) = n4
      x(nF1) = 0.25 * (x(n1)+x(n2)+x(n3)+x(n4))
      y(nF1) = 0.25 * (y(n1)+y(n2)+y(n3)+y(n4))
      z(nF1) = 0.25 * (z(n1)+z(n2)+z(n3)+z(n4))
    end if 

!--- nF2 ---!
    do n=1,NN4
      if( ( NodeN4(n,1) == n1 .and. NodeN4(n,4) == n6 ) .or.  &
	  ( NodeN4(n,1) == n6 .and. NodeN4(n,4) == n1 ) .or.  &
          ( NodeN4(n,1) == n2 .and. NodeN4(n,4) == n5 ) .or.  &
	  ( NodeN4(n,1) == n5 .and. NodeN4(n,4) == n2 ) ) then
	nF2 = NodeN4(n,0) 
      end if
    end do
    if (nF2 == 0) then
      NN4 = NN4 + 1
      NN  = NN  + 1
      nF2 = NN
      NodeN4(NN4,0) = nF2
      NodeN4(NN4,1) = n1
      NodeN4(NN4,2) = n2
      NodeN4(NN4,3) = n5
      NodeN4(NN4,4) = n6
      x(nF2) = 0.25 * (x(n1)+x(n2)+x(n5)+x(n6))
      y(nF2) = 0.25 * (y(n1)+y(n2)+y(n5)+y(n6))
      z(nF2) = 0.25 * (z(n1)+z(n2)+z(n5)+z(n6))
    end if 

!--- nF3 ---!
    do n=1,NN4
      if( ( NodeN4(n,1) == n2 .and. NodeN4(n,4) == n8 ) .or.  &
	  ( NodeN4(n,1) == n8 .and. NodeN4(n,4) == n2 ) .or.  &
          ( NodeN4(n,1) == n4 .and. NodeN4(n,4) == n6 ) .or.  &
	  ( NodeN4(n,1) == n6 .and. NodeN4(n,4) == n4 ) ) then
	nF3 = NodeN4(n,0) 
      end if
    end do
    if (nF3 == 0) then
      NN4 = NN4 + 1
      NN  = NN  + 1
      nF3 = NN
      NodeN4(NN4,0) = nF3
      NodeN4(NN4,1) = n2
      NodeN4(NN4,2) = n4
      NodeN4(NN4,3) = n6
      NodeN4(NN4,4) = n8
      x(nF3) = 0.25 * (x(n2)+x(n4)+x(n6)+x(n8))
      y(nF3) = 0.25 * (y(n2)+y(n4)+y(n6)+y(n8))
      z(nF3) = 0.25 * (z(n2)+z(n4)+z(n6)+z(n8))
    end if 

!--- nF4 ---!
    do n=1,NN4
      if( ( NodeN4(n,1) == n3 .and. NodeN4(n,4) == n8 ) .or.  &
	  ( NodeN4(n,1) == n8 .and. NodeN4(n,4) == n3 ) .or.  &
          ( NodeN4(n,1) == n4 .and. NodeN4(n,4) == n7 ) .or.  &
	  ( NodeN4(n,1) == n7 .and. NodeN4(n,4) == n4 ) ) then
	nF4 = NodeN4(n,0) 
      end if
    end do
    if (nF4 == 0) then
      NN4 = NN4 + 1
      NN  = NN  + 1
      nF4 = NN
      NodeN4(NN4,0) = nF4
      NodeN4(NN4,1) = n3
      NodeN4(NN4,2) = n4
      NodeN4(NN4,3) = n7
      NodeN4(NN4,4) = n8
      x(nF4) = 0.25 * (x(n3)+x(n4)+x(n7)+x(n8))
      y(nF4) = 0.25 * (y(n3)+y(n4)+y(n7)+y(n8))
      z(nF4) = 0.25 * (z(n3)+z(n4)+z(n7)+z(n8))
    end if 

!--- nF5 ---!
    do n=1,NN4
      if( ( NodeN4(n,1) == n1 .and. NodeN4(n,4) == n7 ) .or.  &
	  ( NodeN4(n,1) == n7 .and. NodeN4(n,4) == n1 ) .or.  &
          ( NodeN4(n,1) == n3 .and. NodeN4(n,4) == n5 ) .or.  &
	  ( NodeN4(n,1) == n5 .and. NodeN4(n,4) == n3 ) ) then
	nF5 = NodeN4(n,0) 
      end if
    end do
    if (nF5 == 0) then
      NN4 = NN4 + 1
      NN  = NN  + 1
      nF5 = NN
      NodeN4(NN4,0) = nF5
      NodeN4(NN4,1) = n1
      NodeN4(NN4,2) = n3
      NodeN4(NN4,3) = n5
      NodeN4(NN4,4) = n7
      x(nF5) = 0.25 * (x(n1)+x(n3)+x(n5)+x(n7))
      y(nF5) = 0.25 * (y(n1)+y(n3)+y(n5)+y(n7))
      z(nF5) = 0.25 * (z(n1)+z(n3)+z(n5)+z(n7))
    end if 

!--- nF6 ---!
    do n=1,NN4
      if( ( NodeN4(n,1) == n5 .and. NodeN4(n,4) == n8 ) .or.  &
	  ( NodeN4(n,1) == n8 .and. NodeN4(n,4) == n5 ) .or.  &
          ( NodeN4(n,1) == n6 .and. NodeN4(n,4) == n7 ) .or.  &
	  ( NodeN4(n,1) == n7 .and. NodeN4(n,4) == n6 ) ) then
	nF6 = NodeN4(n,0) 
      end if
    end do
    if (nF6 == 0) then
      NN4 = NN4 + 1
      NN  = NN  + 1
      nF6 = NN
      NodeN4(NN4,0) = nF6
      NodeN4(NN4,1) = n5
      NodeN4(NN4,2) = n6
      NodeN4(NN4,3) = n7
      NodeN4(NN4,4) = n8
      x(nF6) = 0.25 * (x(n5)+x(n6)+x(n7)+x(n8))
      y(nF6) = 0.25 * (y(n5)+y(n6)+y(n7)+y(n8))
      z(nF6) = 0.25 * (z(n5)+z(n6)+z(n7)+z(n8))
    end if 

!---------------------------------!
!     Na kraju cvor u sredini     !
!---------------------------------!
    NN8 = NN8 + 1
    NN  = NN + 1
    n0  = NN
    NodeN8(NN8,0) = n0 
    NodeN8(NN8,1) = n1
    NodeN8(NN8,2) = n2
    NodeN8(NN8,3) = n3
    NodeN8(NN8,4) = n4
    NodeN8(NN8,5) = n5
    NodeN8(NN8,6) = n6
    NodeN8(NN8,7) = n7
    NodeN8(NN8,8) = n8
    x(n0) = 0.125*(x(n1)+x(n2)+x(n3)+x(n4)+x(n5)+x(n6)+x(n7)+x(n8))
    y(n0) = 0.125*(y(n1)+y(n2)+y(n3)+y(n4)+y(n5)+y(n6)+y(n7)+y(n8))
    z(n0) = 0.125*(z(n1)+z(n2)+z(n3)+z(n4)+z(n5)+z(n6)+z(n7)+z(n8))

!-------------------------!
!                         !
!     Postavi cvorove     !
!     novim  celijama     !
!                         !
!-------------------------!

!----- cr1 -!
    CellN(cr1,1) = n1
    CellN(cr1,2) = n12
    CellN(cr1,3) = n13
    CellN(cr1,4) = nF1
    CellN(cr1,5) = n15
    CellN(cr1,6) = nF2
    CellN(cr1,7) = nF5
    CellN(cr1,8) = n0 

!----- cr2 -!
    CellN(cr2,1) = n12
    CellN(cr2,2) = n2 
    CellN(cr2,3) = nF1
    CellN(cr2,4) = n24
    CellN(cr2,5) = nF2
    CellN(cr2,6) = n26 
    CellN(cr2,7) = n0 
    CellN(cr2,8) = nF3

!----- cr3 -!
    CellN(cr3,1) = n13
    CellN(cr3,2) = nF1
    CellN(cr3,3) = n3 
    CellN(cr3,4) = n34
    CellN(cr3,5) = nF5
    CellN(cr3,6) = n0  
    CellN(cr3,7) = n37
    CellN(cr3,8) = nF4

!----- cr4 -!
    CellN(cr4,1) = nF1
    CellN(cr4,2) = n24
    CellN(cr4,3) = n34
    CellN(cr4,4) = n4 
    CellN(cr4,5) = n0 
    CellN(cr4,6) = nF3 
    CellN(cr4,7) = nF4
    CellN(cr4,8) = n48

!----- cr5 -!
    CellN(cr5,1) = n15
    CellN(cr5,2) = nF2
    CellN(cr5,3) = nF5
    CellN(cr5,4) = n0 
    CellN(cr5,5) = n5 
    CellN(cr5,6) = n56 
    CellN(cr5,7) = n57
    CellN(cr5,8) = nF6

!----- cr6 -!
    CellN(cr6,1) = nF2
    CellN(cr6,2) = n26
    CellN(cr6,3) = n0 
    CellN(cr6,4) = nF3
    CellN(cr6,5) = n56
    CellN(cr6,6) = n6  
    CellN(cr6,7) = nF6
    CellN(cr6,8) = n68

!----- cr7 -!
    CellN(cr7,1) = nF5
    CellN(cr7,2) = n0 
    CellN(cr7,3) = n37
    CellN(cr7,4) = nF4
    CellN(cr7,5) = n57
    CellN(cr7,6) = nF6 
    CellN(cr7,7) = n7 
    CellN(cr7,8) = n78

!----- cr8 -!
    CellN(cr8,1) = n0 
    CellN(cr8,2) = nF3
    CellN(cr8,3) = nF4
    CellN(cr8,4) = n48
    CellN(cr8,5) = nF6
    CellN(cr8,6) = n68 
    CellN(cr8,7) = n78
    CellN(cr8,8) = n8 

!============================!
!                            !
!     Neusitnjene celije     !
!                            !
!============================!
    else

      if(CelMar(c1)   >  0) then  ! -> susjed 1 usitnjen
	CellC(c, 1) = CelMar(c1) - 8 + WchNod(c1,n1)
	CellC(c, 7) = CelMar(c1) - 8 + WchNod(c1,n2)
	CellC(c,13) = CelMar(c1) - 8 + WchNod(c1,n3)
	CellC(c,19) = CelMar(c1) - 8 + WchNod(c1,n4)
      endif           

      if(CelMar(c2)   >  0) then  ! -> susjed 2 usitnjen
	CellC(c, 2) = CelMar(c2) - 8 + WchNod(c2, n1)
	CellC(c, 8) = CelMar(c2) - 8 + WchNod(c2, n2)
	CellC(c,14) = CelMar(c2) - 8 + WchNod(c2, n5)
	CellC(c,20) = CelMar(c2) - 8 + WchNod(c2, n6)
      endif           

      if(CelMar(c3)   >  0) then  ! -> susjed 3 usitnjen
	CellC(c, 3) = CelMar(c3) - 8 + WchNod(c3, n2)
	CellC(c, 9) = CelMar(c3) - 8 + WchNod(c3, n4)
	CellC(c,15) = CelMar(c3) - 8 + WchNod(c3, n6)
	CellC(c,21) = CelMar(c3) - 8 + WchNod(c3, n8)
      endif           

      if(CelMar(c4)   >  0) then  ! -> susjed 4 usitnjen
	CellC(c, 4) = CelMar(c4) - 8 + WchNod(c4, n3)
	CellC(c,10) = CelMar(c4) - 8 + WchNod(c4, n4)
	CellC(c,16) = CelMar(c4) - 8 + WchNod(c4, n7)
	CellC(c,22) = CelMar(c4) - 8 + WchNod(c4, n8)
      endif           

      if(CelMar(c5)   >  0) then  ! -> susjed 5 usitnjen
	CellC(c, 5) = CelMar(c5) - 8 + WchNod(c5, n1)
	CellC(c,11) = CelMar(c5) - 8 + WchNod(c5, n3)
	CellC(c,17) = CelMar(c5) - 8 + WchNod(c5, n5)
	CellC(c,23) = CelMar(c5) - 8 + WchNod(c5, n7)
      endif

      if(CelMar(c6)   >  0) then  ! -> susjed 6 usitnjen
	CellC(c, 6) = CelMar(c6) - 8 + WchNod(c6, n5)
	CellC(c,12) = CelMar(c6) - 8 + WchNod(c6, n6)
	CellC(c,18) = CelMar(c6) - 8 + WchNod(c6, n7)
	CellC(c,24) = CelMar(c6) - 8 + WchNod(c6, n8)
      endif           

    end if   
  end do

  write(*,*) 'Number of nodes after the refinement: ', NN 
  write(*,*) 'Number of cells after the refinement: ', NC 

!==============================================!
!                                              !
!     Connect the new twins, if they exist     !
!                                              !
!==============================================!  

  do nA=1,NN2
    nA0=NodeN2(nA,0)
    nA1=NodeN2(nA,1)
    nA2=NodeN2(nA,2)

    if( (TwinN(nA1,0) /= 0).and.(TwinN(nA2,0) /= 0) ) then

      do nB=nA+1,NN2
	nB0=NodeN2(nB,0)
	nB1=NodeN2(nB,1)
	nB2=NodeN2(nB,2)

	if( (TwinN(nB1,0) /= 0).and.(TwinN(nB2,0) /= 0) ) then

	  if( (IsTwin(nA1,nB1) .and. IsTwin(nA2,nB2)) .or.          &
	      (IsTwin(nA1,nB2) .and. IsTwin(nA2,nB1))  ) then
	    if (.not. IsTwin(nA0,nB0)) then
	      TwinN(nA0,0)=TwinN(nA0,0)+1
	      TwinN(nA0,TwinN(nA0,0))=nB0
	      TwinN(nB0,0)=TwinN(nB0,0)+1
	      TwinN(nB0,TwinN(nB0,0))=nA0
	    end if
	  end if
	end if
      end do
    end if
  end do

  do nA=1,NN4
    nA0=NodeN4(nA,0)
    nA1=NodeN4(nA,1)
    nA2=NodeN4(nA,4) ! diagonal

    if( (TwinN(nA1,0) /= 0).and.(TwinN(nA2,0) /= 0) ) then

      do nB=nA+1,NN4
	nB0=NodeN4(nB,0)
	nB1=NodeN4(nB,1)
	nB2=NodeN4(nB,4) ! diagonal

	if( (TwinN(nB1,0) /= 0).and.(TwinN(nB2,0) /= 0) ) then

	  if( (IsTwin(nA1,nB1) .and. IsTwin(nA2,nB2)) .or.          &
	      (IsTwin(nA1,nB2) .and. IsTwin(nA2,nB1))  ) then
	    if (.not. IsTwin(nA0,nB0)) then
	      TwinN(nA0,0)=TwinN(nA0,0)+1
	      TwinN(nA0,TwinN(nA0,0))=nB0
	      TwinN(nB0,0)=TwinN(nB0,0)+1
	      TwinN(nB0,TwinN(nB0,0))=nA0
	    end if
	  end if
	end if
      end do
    end if
  end do

!====================================!
!                                    !
!     Delete the redundant cells     !
!                                    !
!====================================!  

!----- Initialize the new numbers for the cells
  do c=-NbC,NC
    NewN(c)=c
  end do

  del=0 
  do c=1,NC
    if(CelMar(c) /= 0) then
      NewN(c) = -1
      del = del+1
    else
      NewN(c) = c - del 
    endif 
  end do
  write(*,*) 'Deleted cells:', del

  do c=1,NC
    if(NewN(c) /= -1) then
!----- update the cell numbers. Watch out ! The numbers you are
!----- updating are old, so ...
      do n=1,24  ! n is neighbour now
	CellC( NewN(c),n ) = NewN(CellC( c,n ))
      end do

!----- update the node numbers
      do n=1,8   ! n is node now
	CellN( NewN(c),n ) = CellN( c,n ) 
      end do

      material( NewN(c) ) = material( c )  ! -> never checked !
      level( NewN(c) )    = level( c )     ! -> never checked !
    end if
  end do

  do c=NC-del+1, MAXN   ! da obrise stare podatke
    do n=1,24           ! n is neighbour now
      CellC(c,n ) = 0
    end do
  end do

  NC = NC - del    

  write(*,*) 'Number of cells after the renumeration: ', NC 

  END SUBROUTINE Refine
