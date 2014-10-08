!+++++++++++++++++++++++++++++++++++++!
!                                     !
!     Global variable definitions     !
!       for the mesh generator        !
!                                     !
!+++++++++++++++++++++++++++++++++++++!
!..RCS/CVS ident
! $Id: gen_mod.h90,v 1.8 2000/03/22 21:10:47 bojan Exp $
! $Source: /home/muhamed/.CVSROOT/T-Rex/Library/gen_mod.h90,v $ 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
MODULE gen_mod

  USE allp_mod

  IMPLICIT NONE

  REAL,ALLOCATABLE :: x(:),  y(:),  z(:)   ! node coordinates
  REAL,ALLOCATABLE :: walln(:)             ! node distance from the wall 
  INTEGER,ALLOCATABLE :: SideN(:,:)        ! numb, n1, n2, n3, n4
  INTEGER,ALLOCATABLE :: SideCc(:,:)
						
  INTEGER,ALLOCATABLE :: CellC(:,:)        ! cell's neighbours
  INTEGER,ALLOCATABLE :: CellN(:,:)        ! cell nodes

  INTEGER,ALLOCATABLE :: TwinN(:,:)

  INTEGER,ALLOCATABLE :: NewN(:)    ! new number for the nodes and cells
  INTEGER,ALLOCATABLE :: NewC(:)    ! new number for cells
  INTEGER,ALLOCATABLE :: NewS(:)    ! new number for sides 
  INTEGER,ALLOCATABLE :: CelMar(:)  ! cell marker

  INTEGER,ALLOCATABLE :: NodeN2(:,:)    
  INTEGER,ALLOCATABLE :: NodeN4(:,:)  
  INTEGER,ALLOCATABLE :: NodeN8(:,:)    

  INTEGER,ALLOCATABLE :: level(:)   ! refinement level

  INTEGER :: MAXN, MAXB, MAXS

  INTEGER :: NR(MAXP)               ! refin. levels, refin. regions

  REAL    :: xp(MAXP), yp(MAXP), zp(MAXP)       ! point coordinates
  REAL    :: xl(MAXP,MAXL),yl(MAXP,MAXL),zl(MAXP,MAXL),LinWgt(MAXP)
  REAL    :: BlkWgt(MAXL,3), BlFaWt(MAXL,3)     ! leave this 
  REAL    :: FRegio(MAXP,MAXP,0:6)              ! levels, regions

  REAL    :: SRegio(MAXP,0:6), Srelax(MAXP)  ! levels, regions
  LOGICAL :: SdirX(MAXP), SdirY(MAXP), SdirZ(MAXP)   
  INTEGER :: Siter(MAXP)   

  INTEGER :: BlkPnt(MAXP,0:8),  & ! 0 for orientation                
	     BlkRes(MAXP,6),    & ! NI,NJ,NK,NI*NJ*NK,NNo,NVo       
	     BlkFac(MAXP,6,4),  &                                    
	     BlFaLa(MAXP),      &                                   
	     Bound(MAXP,8),     &                                  
	     Period(MAXP,8),    &                                 
	     Copy(MAXP,0:8)

  INTEGER :: LinPnt(MAXL,2), LinRes(MAXL)
  INTEGER :: Nbloc, NP, Nline, Nsurf, Nboun, Nperi
  INTEGER :: NN, NN2, NN4, NN8
  INTEGER :: NSR                   ! smoothing regions
  INTEGER :: NSsh                  ! number of shadow faces

  INTEGER :: WallFacFst, WallFacLst 

  INTEGER :: ELIPSO, RECTAN, PLANE,YES,NO

  CHARACTER*4 :: BndFac(MAXP)

END MODULE
