!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   Global variable definitions   !   Delft University of Technology   !
!         for all modules         !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!..RCS/CVS ident
! $Id: all_mod.h90,v 1.8 2002/10/31 11:26:48 niceno Exp $
! $Source: /home/muhamed/.CVSROOT/T-Rex/Library/all_mod.h90,v $
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!   Note: cell_n, parent, A_row, A_col, A_dia, side_c, side_cc, 
!         sideAij, are for all grids
!======================================================================!
MODULE all_mod

  IMPLICIT NONE

  REAL,ALLOCATABLE :: xc(:),yc(:),zc(:) 
  REAL,ALLOCATABLE :: Sx(:),Sy(:),Sz(:)
  REAL,ALLOCATABLE :: volume(:)            ! cell's volume
  REAL,ALLOCATABLE :: delta(:)             ! delta (max(dx,dy,dz))
  REAL,ALLOCATABLE :: a1(:)                ! scotti
  REAL,ALLOCATABLE :: a2(:)                ! scotti
  REAL,ALLOCATABLE :: Dx(:),Dy(:),Dz(:)
  REAL,ALLOCATABLE :: xsp(:),ysp(:),zsp(:) ! face coordinates    
  REAL,ALLOCATABLE :: WallDs(:), f(:)

  CHARACTER :: name*80
  CHARACTER :: inp*300
  INTEGER   :: tn, ts(300), te(300)

  INTEGER   :: NC, NS                    ! num. of nodes and cells 
  INTEGER   :: NbC
  INTEGER   :: MNBS
  INTEGER   :: NRL
  INTEGER   :: Ncopy
  INTEGER   :: Nmat                      ! number of materials
  LOGICAL   :: Mater(1024)               ! is the material present ?

  INTEGER,ALLOCATABLE :: material(:)     ! material markers
  INTEGER,ALLOCATABLE :: SideC(:,:)      !  c0, c1, c2

  INTEGER,ALLOCATABLE :: TypeBC(:)       ! type of boundary condition
  INTEGER,ALLOCATABLE :: bcmark(:)

  INTEGER,ALLOCATABLE :: CopyC(:)        !  might be shorter
  INTEGER,ALLOCATABLE :: CopyS(:,:)      !  similar to SideC 

  INTEGER,ALLOCATABLE :: SideC1C2(:,:)   !  similar to SideC 

  REAL, ALLOCATABLE   :: Dxsp(:,:)       !  similar to SideC 
  REAL, ALLOCATABLE   :: Dysp(:,:)       !  similar to SideC 
  REAL, ALLOCATABLE   :: Dzsp(:,:)       !  similar to SideC 

END MODULE 
