!======================================================================!
  PROGRAM Linija
!----------------------------------------------------------------------!
!   Places the nodes on the line defined with local block position     !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------[Calling]-------------------------------!
  REAL     atanh
!-------------------------------[Locals]-------------------------------!
  INTEGER  N, i, node, case
  REAL     x0, delx, t, dt, ddt, pr, xi
  REAL     x(0:1000)
  REAL     w
!--------------------------------[CVS]---------------------------------!
  character*80 rcs1,rcs2
  data rcs1/                                                        &
  '$Id: Stretch.f90,v 1.1 2000/09/07 07:39:08 niceno Exp $'/
  data rcs2/                                                        &
  '$Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Test/Stretch.f90,v $'/
!======================================================================!

 character*20 file1 

  x0=0.0
 
  write(*,*) 'Distance: '
  read(*,*)  delx 
  write(*,*) 'Weigth: '
  read(*,*)  w
  write(*,*) 'Number of cells: '
  read(*,*)  N

  x(0) = x0

!-----------------------------!
!     Linear distribution     !
!-----------------------------!
  if(w  > 0.0) then
    ddt = ( 2.0*(1.0-w) ) / ( 1.0*N*(1.0*N-1.0)*(1.0+w) )
    t=0.0
    do i=1,N
      dt=1.0/(1.0*N)+(1.0*i-0.5*(1.0*N+1)) * ddt
      t=t+dt
      x(i) = x0 + t*delx
    end do

!---------------------------------!
!     Hyperbolic distribution     !
!---------------------------------!
  else
    case = 0
    if     ((w  >  -0.5).and.(w <=  -0.25)) then
      pr = 1.0 - abs(0.5 - abs(w))    
      case = 1
    else if((w >=  -0.75).and.(w  <  -0.5)) then
      pr = 1.0 - abs(0.5 - abs(w))            
      case = 2
    else
      pr = -w
      case = 3 
    endif

    do i=1,N
      if(case == 1) xi = -1.0*(1.0*i)/(1.0*N)
      if(case == 2) xi =  1.0 - 1.0*(1.0*i)/(1.0*N)
      if(case == 3) xi = -1.0 + 2.0*(1.0*i)/(1.0*N)
      if    (case == 1) then
        x(i) = x0 - (tanh(xi*atanh(pr))/pr)*delx
      elseif(case == 2) then
        x(i) = x0 + delx - (tanh(xi*atanh(pr))/pr)*delx
      elseif(case == 3) then
        x(i) = x0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delx
      endif 
    end do

  endif

!------------------------------------!
!     Print out the distribution     !
!------------------------------------!
  write(*,*) ' number:   node:                  cell:'
  do i=N,1,-1
    write(*,'(I6,F12.6,A10)') i+1, x(i),              ' +-------+'
    write(*,'(I6,A23,F12.6)') i,  '|   o   | ', 0.5*(x(i)+x(i-1))
  end do
  write(*,'(I6,F12.6,A10)') i+1, x(0),              ' +-------+'
  write(*,*) 'wall:    //////////////////////////////'

  file1 = 'fileout'

  open(3,file=file1) 
  
  do i = 1, N
  
   write(3,*) 0.5*(x(i)+x(i-1)) 

  end do

  close (3)
 
  END PROGRAM 

!======================================================================!
  REAL FUNCTION Atanh(x)
!----------------------------------------------------------------------!
!   Calculates inverse of hyperbolic tangens.                          !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  REAL   x
!--------------------------------[CVS]---------------------------------!
  character*80 rcs1,rcs2
  data rcs1/                                                        &
  '$Id: Stretch.f90,v 1.1 2000/09/07 07:39:08 niceno Exp $'/
  data rcs2/                                                        &
  '$Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Test/Stretch.f90,v $'/
!======================================================================!
  if(x  > 1.0) then
    write(*,*) 'Error message from atanh: bad argument'
    stop
  end if 

  atanh=log( sqrt( (1.0+x)/(1.0-x) ) )

  END FUNCTION
