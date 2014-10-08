!======================================================================!
  SUBROUTINE TopolM         
!----------------------------------------------------------------------!
! Determines the topology of the emesh.                                 !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod 
  USE gen_mod 
  USE neu_mod 
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: i, j, NbounCells
  CHARACTER :: nameOut*130
!--------------------------------[CVS]---------------------------------!
  character*80 rcs1,rcs2
  data rcs1/                                                        &
  '$Id: TopolM.f90,v 1.5 2004/06/24 14:09:45 muhamed Exp $'/ 
  data rcs2/                                                        &
  '$Source: /home/muhamed/.CVSROOT/T-Rex/Neu2Trex/TopolM.f90,v $'/
!======================================================================!

!----------------------------------------------------------------------!
!
!    *.NEU nodes are numbered as this:
!    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!        7-----------8
!       /|          /|
!      /           / |
!     /  |        /  |
!    5-----------6   |
!    |   |       |   |
!    |   3- - - -|- -4
!    |  /        |  /
!    |           | /
!    |/          |/
!    1-----------2                 
!
!    Figure 1: Numbering of cell nodes 
!
!        7-----------8             7-----------8 
!       /|          /|            /|          /| 
!      /           / |           /    (6)    / | 
!     /  |    (3) /  |          /  |        /  | 
!    5-----------6   |         5-----------6   | 
!    |(4)|       |(2)|         |   |       |   | 
!    |   3- - - -|- -4         |   3- - - -|- -4 
!    |  / (1)    |  /          |  /        |  / 
!    |           | /           |      (5)  | / 
!    |/          |/            |/          |/ 
!    1-----------2             1-----------2 
!
!    Figure 2: Numbering of directions of boundary  
!
! The directions for boundary conditions are:
!
!   direction:    nodes (not in anticlockwise order):
!   ------------------------------------------------- 
!      (1)   ->   1,2,5,6
!      (2)   ->   2,4,6,8
!      (3)   ->   3,4,7,8
!      (4)   ->   1,3,5,7
!      (5)   ->   1,2,3,4
!      (6)   ->   5,6,7,8
!
!----------------------------------------------------------------------!

!-----------------------------
!---- count the boundary cells
!-----------------------------
  NBc = 0
  NS  = 0
  do i=1,NC
    do j=1,6
      if(BCtype(i,j) /= 0) then
        NBc = NBc + 1 

!---- BCmark:
        BCmark(-NBc) = BCtype(i,j)
!---- material:
        material(-NbC) = material(i)

!---- sides
        NS  = NS  + 1
        SideC(1,NS) = i
        SideC(2,NS) = -NbC

!---- hexahedra:
        if(CellN(i,0) == 8) then
          CellN(-NBc,0) = 4 
          CellN(-NBc,1) = CellN(i,f8n(j,1))
          CellN(-NBc,2) = CellN(i,f8n(j,2))
          CellN(-NBc,3) = CellN(i,f8n(j,3))
          CellN(-NBc,4) = CellN(i,f8n(j,4))

          SideN(NS,0) = 4
          SideN(NS,1) = CellN(-NBc,1)
          SideN(NS,2) = CellN(-NBc,2)
          SideN(NS,3) = CellN(-NBc,3)
          SideN(NS,4) = CellN(-NBc,4)
        end if

!---- prisms:
        if(CellN(i,0) == 6) then
          if(j <= 3) then    ! faces (1), (2) and (3)
            CellN(-NBc,0) = 4 
            CellN(-NBc,1) = CellN(i,f6n(j,1))
            CellN(-NBc,2) = CellN(i,f6n(j,2))
            CellN(-NBc,3) = CellN(i,f6n(j,3))
            CellN(-NBc,4) = CellN(i,f6n(j,4))

            SideN(NS,0) = 4
            SideN(NS,1) = CellN(-NBc,1)
            SideN(NS,2) = CellN(-NBc,2)
            SideN(NS,3) = CellN(-NBc,3)
            SideN(NS,4) = CellN(-NBc,4)
          else if(j <= 5) then
            CellN(-NBc,0) = 3 
            CellN(-NBc,1) = CellN(i,f6n(j,1))
            CellN(-NBc,2) = CellN(i,f6n(j,2))
            CellN(-NBc,3) = CellN(i,f6n(j,3))

            SideN(NS,0) = 3
            SideN(NS,1) = CellN(-NBc,1)
            SideN(NS,2) = CellN(-NBc,2)
            SideN(NS,3) = CellN(-NBc,3)
          end if
        end if

!---- tetrahedra:
        if(CellN(i,0) == 4) then
          if(j <= 4) then
            CellN(-NBc,0) = 3 
            CellN(-NBc,1) = CellN(i,f4n(j,1))
            CellN(-NBc,2) = CellN(i,f4n(j,2))
            CellN(-NBc,3) = CellN(i,f4n(j,3))

            SideN(NS,0) = 3
            SideN(NS,1) = CellN(-NBc,1)
            SideN(NS,2) = CellN(-NBc,2)
            SideN(NS,3) = CellN(-NBc,3)
          end if
        end if

!---- pyramides:
        if(CellN(i,0) == 5) then
          if(j == 1) then    ! face (1)
            CellN(-NBc,0) = 4
            CellN(-NBc,1) = CellN(i,f5n(j,1))
            CellN(-NBc,2) = CellN(i,f5n(j,2))
            CellN(-NBc,3) = CellN(i,f5n(j,3))
            CellN(-NBc,4) = CellN(i,f5n(j,4))
 
            SideN(NS,0) = 4
            SideN(NS,1) = CellN(-NBc,1)
            SideN(NS,2) = CellN(-NBc,2)
            SideN(NS,3) = CellN(-NBc,3)
            SideN(NS,4) = CellN(-NBc,4)
          else if(j <= 5) then
            CellN(-NBc,0) = 3
            CellN(-NBc,1) = CellN(i,f5n(j,1))
            CellN(-NBc,2) = CellN(i,f5n(j,2))
            CellN(-NBc,3) = CellN(i,f5n(j,3))
 
            SideN(NS,0) = 3
            SideN(NS,1) = CellN(-NBc,1)
            SideN(NS,2) = CellN(-NBc,2)
            SideN(NS,3) = CellN(-NBc,3)
          end if
        end if

      end if
    end do 
  end do
  write(*,*) 'Boundary cells: ', NBc

  END SUBROUTINE TopolM     
