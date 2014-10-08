!======================================================================!
  SUBROUTINE Probe2D 
!----------------------------------------------------------------------!
! Finds coordinates of all the planes for the channel flow.            !
! It assumes that homogeneous directions of the flow are x and y.      !
!----------------------------------------------------------------------!
  USE all_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------[Calling]-------------------------------! 
  INTERFACE
    LOGICAL FUNCTION Approx(A,B,tol)
      REAL           :: A,B
      REAL, OPTIONAL :: tol
    END FUNCTION Approx
  END INTERFACE
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: Nprob, p, c
  REAL      :: yp(20000), zp(20000)
  CHARACTER :: namPro*80
  CHARACTER :: answer*80
!--------------------------------[CVS]---------------------------------!
!  $Id: Probe2D.f90,v 1.4 2002/10/30 16:29:33 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Library/Probe2D.f90,v $  
!======================================================================!

  write(*,*) '==============================='
  write(*,*) ' Looking for homogeneous plane '
  write(*,*) '-------------------------------'
  write(*,*) 'Insert homogeneous direction (xy,yz,zx or skip)'
  read(*,*) answer
  call touppr(answer)
  if(answer=='SKIP') return

  NProb = 0
  zp=0.0
  yp=0.0

  do c=1,NC

!---- try to find the cell among the probes
    do p=1,Nprob
      if(answer=='YZ') then
        if( Approx(yc(c), yp(p)) .and. &  
            Approx(zc(c), zp(p)) ) go to 1
      else if(answer=='ZX') then
        if( Approx(xc(c), yp(p)) .and. &  
            Approx(zc(c), zp(p)) ) go to 1
      else if(answer=='XY') then
        if( Approx(xc(c), yp(p)) .and. &  
            Approx(yc(c), zp(p)) ) go to 1
      end if
    end do 

!---- couldn't find a cell among the probes, add a new one
    Nprob = Nprob+1
    if(answer=='YZ') then
      yp(Nprob)=yc(c)
      zp(Nprob)=zc(c)
    else if(answer=='ZX') then
      yp(Nprob)=xc(c)
      zp(Nprob)=zc(c)
    else if(answer=='XY') then
      yp(Nprob)=xc(c)
      zp(Nprob)=yc(c)
    end if 

    if(Nprob == 20000) then
      write(*,*) 'Probe 2D: Not a 2D problem.'
      return
    end if
1 end do

!<<<<<<<<<<<<<<<<<<<<<<<<!
!     create 2D file     !
!<<<<<<<<<<<<<<<<<<<<<<<<!
  namPro = name
  namPro(len_trim(name)+1:len_trim(name)+3) = '.2D'
  write(6, *) 'Now creating the file:', namPro
  open(9, FILE=namPro)

!---- write the number of probes 
  write(9,'(I8)') Nprob

!---- write the probe coordinates out
  do p=1,Nprob
    write(9,'(I8,1PE17.8,1PE17.8)') p, yp(p), zp(p)
  end do

  close(9)

  END SUBROUTINE Probe2D
