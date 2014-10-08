!======================================================================!
  SUBROUTINE Probe1D
!----------------------------------------------------------------------!
! This program finds the coordinate of cell-centers in non-homogeneous
! direction and write them in file name.1D
!----------------------------------------------------------------------!
  USE all_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  LOGICAL :: isit
!------------------------------[Calling]-------------------------------! 
  INTERFACE
    LOGICAL FUNCTION Approx(A,B,tol)
      REAL           :: A,B
      REAL, OPTIONAL :: tol
    END FUNCTION Approx
  END INTERFACE
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: Nprob, p, c
  REAL      :: zp(1000)
  CHARACTER :: namPro*80
  CHARACTER :: answer*80
!--------------------------------[CVS]---------------------------------!
!  $Id: Probe1D.f90,v 1.5 2002/10/30 16:29:33 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Library/Probe1D.f90,v $    
!======================================================================!

  write(*,*) '========================================'
  write(*,*) ' Looking for non-homogeneous directions '
  write(*,*) '----------------------------------------'
  write(*,*) 'Insert non-homogeneous direction (x,y,z or skip)'
  read(*,*) answer
  call touppr(answer)
  if(answer=='SKIP') return

  NProb = 0
  zp=0.0

  do c=1,NC

!---- try to find the cell among the probes
    do p=1,Nprob
      if(answer == 'X') then
        if( Approx(xc(c), zp(p)) ) go to 1
      else if(answer == 'Y') then
        if( Approx(yc(c), zp(p)) ) go to 1
      else if(answer == 'Z') then
        if( Approx(zc(c), zp(p)) ) go to 1
      end if
    end do 

!---- couldn't find a cell among the probes, add a new one
    Nprob = Nprob+1
    if(answer=='X') zp(Nprob)=xc(c)
    if(answer=='Y') zp(Nprob)=yc(c)
    if(answer=='Z') zp(Nprob)=zc(c)

    if(Nprob == 1000) then
      write(*,*) 'Probe 1D: Not a 1D (channel flow) problem.'
      isit = .false.
      return
    end if
1 end do

  isit = .true.

!<<<<<<<<<<<<<<<<<<<<<<<<!
!     create 1D file     !
!<<<<<<<<<<<<<<<<<<<<<<<<!
  namPro = name
  namPro(len_trim(name)+1:len_trim(name)+4) = '.1Dc'
  write(6, *) 'Now creating the file:', namPro
  open(9, FILE=namPro)

!---- write the number of probes 
  write(9,'(I8)') Nprob

!---- write the probe coordinates out
  do p=1,Nprob
    write(9,'(I8,1PE17.8)') p, zp(p)
  end do

  close(9)

  END SUBROUTINE Probe1D
