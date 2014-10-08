!======================================================================!
  SUBROUTINE Smooth
!----------------------------------------------------------------------!
!   Smooths the grid lines by a Laplacian-like algorythm.              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER             :: c, n, i,j,k,m
  REAL                :: xTnew, yTnew, zTnew    ! temporary new values
  REAL                :: xmax, ymax, zmax, xmin, ymin, zmin 
  REAL                :: x1,y1,z1,x8,y8,z8 
  INTEGER             :: regio
  REAL,ALLOCATABLE    :: xnew(:), ynew(:), znew(:) ! new values
  INTEGER,ALLOCATABLE :: NodeN(:,:)
!--------------------------------[CVS]---------------------------------!
!  $Id: Smooth.f90,v 1.8 2002/10/30 16:29:22 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Generate/Smooth.f90,v $  
!======================================================================!

!---- allocate memory for additional arrays
  allocate(xnew(MAXN)); xnew=0
  allocate(ynew(MAXN)); ynew=0
  allocate(znew(MAXN)); znew=0
  allocate(NodeN(MAXN,0:40)); NodeN=0

  write(*,*) 'Now smoothing the cells. This may take a while !' 

  do n=1,Nn
    NodeN(n,0) = 0
  end do 

  xmax=-HUGE
  ymax=-HUGE
  zmax=-HUGE
  xmin=+HUGE
  ymin=+HUGE
  zmin=+HUGE
  do n=1,Nn
    xmax=max(x(n),xmax) 
    ymax=max(y(n),ymax) 
    zmax=max(z(n),zmax) 
    xmin=min(x(n),xmin) 
    ymin=min(y(n),ymin) 
    zmin=min(z(n),zmin) 
  end do 

!---------------------------!
!     Connect the nodes     !
!---------------------------!
  do c=1,Nc            ! through cells
    do i=1,8           ! through nodes of a cell 
      n = CellN(c,i)   ! first cell
      do j=1,8         ! through nodes of a cell 
	m = CellN(c,j) ! second cell 
	if(n  /=  m) then 
	  do k=1,NodeN(n,0)
	    if(NodeN(n,k) == m) goto 10            
	  end do
	  NodeN(n,0)=NodeN(n,0)+1
	  NodeN(n,NodeN(n,0)) = m
	end if
10      end do
    end do
  end do 

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Browse through regions     !
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  do regio=1,NSR
    if( ( .not. SdirX(regio) ) .and.  &
	( .not. SdirY(regio) ) .and.  &
	( .not. SdirZ(regio) ) ) then
      do n=1,Nn
	x1=SRegio(regio,1)
	y1=SRegio(regio,2)
	z1=SRegio(regio,3)
	x8=SRegio(regio,4)
	y8=SRegio(regio,5)
	z8=SRegio(regio,6)
	if( (x1 <= x(n)) .and. (x(n) <= x8) .and.  &
	    (y1 <= y(n)) .and. (y(n) <= y8) .and.  &
	    (z1 <= z(n)) .and. (z(n) <= z8) ) then
	  NodeN(n,0) = 0
	endif
      end do
    end if
  end do

!->>> do n=1,Nn
!->>>   write(*,*) '=>', n, NodeN(n,0)
!->>>   write(*,*) (NodeN(n,i), i=1,NodeN(n,0) )
!->>> end do

!-------------------------!
!     Smooth the grid     !
!-------------------------!
  do regio=1,NSR
    write(*,*) 'Now smoothing region ',regio,' with:',              &
		Siter(regio), ' iterations.'

    do j=1,Siter(regio)         

!---- calculate new coordinates using the old values (x(),y(),z())
      do n=1,Nn
	if(NodeN(n,0)   >  0) then
	  xTnew=0.0
	  yTnew=0.0
	  zTnew=0.0
	  do i=1,NodeN(n,0)
	    xTnew = xTnew + x(NodeN(n,i))
	    yTnew = yTnew + y(NodeN(n,i))
	    zTnew = zTnew + z(NodeN(n,i))
	  end do
	  xTnew = xTnew / (1.0*NodeN(n,0)) 
	  yTnew = yTnew / (1.0*NodeN(n,0)) 
	  zTnew = zTnew / (1.0*NodeN(n,0)) 
	  if(x(n)  > 0.001*xmin .and. x(n)  < 0.999*xmax)           &
	  xnew(n) = (1.0-Srelax(regio)*walln(n))*x(n)               &
		  +      Srelax(regio)*walln(n) *xTnew
	  if(y(n)  > 0.001*ymin .and. y(n)  < 0.999*ymax)           &
	  ynew(n) = (1.0-Srelax(regio)*walln(n))*y(n)               &
		  +      Srelax(regio)*walln(n) *yTnew
	  if(z(n)  > 0.001*zmin .and. z(n)  < 0.999*zmax)           &
	  znew(n) = (1.0-Srelax(regio)*walln(n))*z(n)               &
		  +      Srelax(regio)*walln(n)* zTnew
	end if
      end do

!---- update coordinates
      do n=1,Nn
	if(NodeN(n,0)   >  0) then

	  x1=SRegio(regio,1)
	  y1=SRegio(regio,2)
	  z1=SRegio(regio,3)
	  x8=SRegio(regio,4)
	  y8=SRegio(regio,5)
	  z8=SRegio(regio,6)

	  if( (x1 <= x(n)) .and. (x(n) <= x8) .and.                 &
	      (y1 <= y(n)) .and. (y(n) <= y8) .and.                 &
	      (z1 <= z(n)) .and. (z(n) <= z8) ) then

	    if(SdirX(regio)) then
	      if(x(n)  > 0.001*xmin .and. x(n)  < 0.999*xmax)       &
	      x(n)=xnew(n)
	    end if

	    if(SdirY(regio)) then
	      if(y(n)  > 0.001*ymin .and. y(n)  < 0.999*ymax)       &
	      y(n)=ynew(n)
	    end if 

	    if(SdirZ(regio)) then
	      if(z(n)  > 0.001*zmin .and. z(n)  < 0.999*zmax)       &
	      z(n)=znew(n)
	    end if 

	  end if  ! if the point belongs to region
	end if    ! if the point is not excluded from smoothing
      end do      ! through nodes
    end do
  end do

  deallocate(xnew)
  deallocate(ynew)
  deallocate(znew)
  deallocate(NodeN)

  END SUBROUTINE Smooth
