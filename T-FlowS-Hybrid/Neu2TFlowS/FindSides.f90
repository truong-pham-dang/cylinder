!======================================================================!
  SUBROUTINE FindSides
!----------------------------------------------------------------------!
!  Creates the "SideC structure"                                       !
!------------------------------[Modules]-------------------------------!
  USE all_mod 
  USE neu_mod 
  USE gen_mod 
!----------------------------------------------------------------------!
  IMPLICIT NONE
!----------------------------------------------------------------------!
  INTERFACE
    LOGICAL FUNCTION Approx(A,B,tol)
      IMPLICIT NONE
      REAL          :: A,B
      REAL,OPTIONAL :: tol
    END FUNCTION Approx
  END INTERFACE
!-------------------------------[Locals]-------------------------------!
  INTEGER             :: c, c1, c2, n1, n2, n3, n4
  INTEGER             :: Nmatch, j, MatchNodes(-1:8) 
  INTEGER             :: i1, i2, Nuber
  INTEGER             :: Broj
  INTEGER             :: fn(6,4)
  REAL   ,ALLOCATABLE :: SideCoor(:) 
  INTEGER,ALLOCATABLE :: SideCell(:), Starts(:), Ends(:) 
  REAL                :: VeryBig
!--------------------------------[CVS]---------------------------------!
  character*80 rcs1,rcs2
  data rcs1/                                                        &
  '$Id: FindSides.f90,v 1.14 2004/06/24 14:08:22 muhamed Exp $'/ 
  data rcs2/                                                        &
  '$Source: /home/muhamed/.CVSROOT/T-Rex/Neu2Trex/FindSides.f90,v $'/
!======================================================================!

  VeryBig = max(NN,NC)

  allocate(SideCoor(NC*6)); SideCoor = NN*1E+30
  allocate(SideCell(NC*6)); SideCell = 0    
  allocate(Starts(NC*6)); Starts = 0    
  allocate(Ends(NC*6));   Ends = 0    
  allocate(CellC(NC,6));  CellC = 0

!--------------------------------------------------
!---- Fill the generic coordinates with some values
!--------------------------------------------------
  do c=1,NC
    if(CellN(c,0) == 4) fn = f4n
    if(CellN(c,0) == 5) fn = f5n
    if(CellN(c,0) == 6) fn = f6n
    if(CellN(c,0) == 8) fn = f8n 
    do j=1,6
      if(BCtype(c,j) == 0) then 
        n1 = CellN(c, fn(j,1))
        n2 = CellN(c, fn(j,2))
        n3 = CellN(c, fn(j,3))
        n4 = CellN(c, fn(j,4))
        if( max(n1,n2,n3,n4)  >  0 ) then
          if(n4 > 0) then
            SideCoor((c-1)*6+j) =  VeryBig*(max(n1 , n2 , n3 , n4))   &
                                +           min(n1 , n2,  n3 , n4)
          else
            SideCoor((c-1)*6+j) =  VeryBig*(max(n1 , n2 , n3))   &
                                +           min(n1 , n2,  n3)
           end if
          SideCell((c-1)*6+j) = c 
        end if 
      end if
    end do
  end do

!-------------------------------------------------
!---- Sort the cell faces according to coordinares
!-------------------------------------------------
  call RISort(SideCoor,SideCell,NC*6,2)

!-----------------------------------------------
!---- Anotate cell faces with same coordinates
!---- (I am afraid that this might be influenced
!----      by the numerical round-off errors)  
!-----------------------------------------------
  Nuber = 1
  Starts(1) = 1
  do c=2,NC*6
    if( SideCoor(c) /= SideCoor(c-1) ) then
      Nuber = Nuber + 1
      Starts(Nuber) = c
      Ends(Nuber-1) = c-1
    end if
  end do

!==========================================
!==== Main loop to fill the SideC structure 
!==========================================
  do n3=1,Nuber
    if(Starts(n3) /= Ends(n3)) then
      do i1=Starts(n3),Ends(n3)
        do i2=i1+1,Ends(n3)
          c1 = min(SideCell(i1),SideCell(i2))
          c2 = max(SideCell(i1),SideCell(i2))
          if(c1 /= c2) then
!-----------------------------
!---- number of matching nodes
!-----------------------------
            Nmatch     = 0
            MatchNodes = 0 
            do n1=1,CellN(c1,0)
              do n2=1,CellN(c2,0)
                if(CellN(c1,n1)==CellN(c2,n2)) then
                  Nmatch = Nmatch + 1 
                  MatchNodes(n1) = 1
                end if
              end do
            end do
!----------------------
!---- general + general   
!----   c1        c2
!----------------------
            if(Nmatch > 2) then 
              if(CellN(c1,0) == 4) fn = f4n
              if(CellN(c1,0) == 5) fn = f5n
              if(CellN(c1,0) == 6) fn = f6n
              if(CellN(c1,0) == 8) fn = f8n
              do j=1,6
                if(   CellC(c1,j) == 0  .and.        & ! not set yet         
                    ( max( MatchNodes(fn(j,1)),0 ) + &
                      max( MatchNodes(fn(j,2)),0 ) + &
                      max( MatchNodes(fn(j,3)),0 ) + &
                      max( MatchNodes(fn(j,4)),0 ) == Nmatch ) ) then
                  NS = NS + 1 
                  SideC(1,NS) = c1
                  SideC(2,NS) = c2
                  SideN(NS,0) = Nmatch 
                  SideN(NS,1) = CellN(c1, fn(j,1))                
                  SideN(NS,2) = CellN(c1, fn(j,2))                
                  SideN(NS,3) = CellN(c1, fn(j,3))
                  SideN(NS,4) = CellN(c1, fn(j,4))
                  CellC(c1,j) = 1 !  -> means: set
                end if
              end do
            end if   ! Nmatch /= 2
          end if   ! c1 /= c2
        end do   ! i2
      end do   ! i1
    end if
  end do    ! do n3

  write(*,'(A30,I9,Z9)') 'FindSides: Number of sides: ', NS, NS

  END SUBROUTINE FindSides
