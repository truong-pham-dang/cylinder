!======================================================================!
  SUBROUTINE Probe1D_nodes_gen
!----------------------------------------------------------------------!
! This program finds the coordinate of nodes in non-homogeneous  
! direction and write them in file name.1D
!----------------------------------------------------------------------!
  USE all_mod
  USE gen_mod
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
  INTEGER   :: Nprob, p, c, n
  REAL      :: N_p(10000)
  CHARACTER :: namPro*80
  CHARACTER :: answer*80
!--------------------------------[CVS]---------------------------------!
!  $Id: Probe1D.f90,v 1.5 2002/10/30 16:29:33 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Library/Probe1D.f90,v $    
!======================================================================!

  write(*,*) '==========================================='
  write(*,*) ' Creating 1D file with the node '
  write(*,*) ' coordinates in non-homogeneous directions '
  write(*,*) '-------------------------------------------'
  write(*,*) 'Insert non-homogeneous direction (x, y, z, Rx, Ry, Rz or skip)'
  read(*,*) answer
  call touppr(answer)
  if(answer=='SKIP') return
 
  NProb = 0
  N_p   = 0.0

  do n =  1, NN 
!---- try to find the cell among the probes
    do p=1,Nprob
        if(answer == 'X') then
          if( Approx(x(n), N_p(p)) ) go to 1
        else if(answer == 'Y') then
          if( Approx(y(n), N_p(p)) ) go to 1
        else if(answer == 'Z') then
          if( Approx(z(n), N_p(p)) ) go to 1
        else if(answer == 'RX') then
          if( Approx( (z(n)**2.0 + y(n)**2.0)**0.5, N_p(p)) ) go to 1
        else if(answer == 'RY') then
          if( Approx( (x(n)**2.0 + z(n)**2.0)**0.5, N_p(p)) ) go to 1
        else if(answer == 'RZ') then
          if( Approx( (x(n)*x(n) + y(n)*y(n))**0.5, N_p(p)) ) go to 1
        end if
    end do 
  
!---- couldn't find a cell among the probes, add a new one
    Nprob = Nprob+1
    if(answer=='X') N_p(Nprob)= x(n)
    if(answer=='Y') N_p(Nprob)= y(n)
    if(answer=='Z') N_p(Nprob)= z(n)

    if(answer=='RX') N_p(Nprob)= (z(n)**2.0 + y(n)**2.0)**0.5
    if(answer=='RY') N_p(Nprob)= (x(n)**2.0 + z(n)**2.0)**0.5
    if(answer=='RZ') N_p(Nprob)= (x(n)*x(n) + y(n)*y(n))**0.5

    if(Nprob == 10000) then
      write(*,*) 'Probe 1D: Not a 1D (channel flow) problem.'
      isit = .false.
      return
    end if
1  end do

  isit = .true.

!<<<<<<<<<<<<<<<<<<<<<<<<!
!     create 1D file     !
!<<<<<<<<<<<<<<<<<<<<<<<<!
  namPro = name
  namPro(len_trim(name)+1:len_trim(name)+3) = '.1D'
  write(6, *) 'Now creating the file:', namPro
  open(9, FILE=namPro)
!---- write the number of probes 
  write(9,'(I8)') Nprob

!  call SSORT (N_p, 1, Nprob*3, 1)
!  call  RISort(N_p, Nprob, Nprob*6,2)

  call SORT2(N_p, Nprob*2, Nprob)

!---- write the probe coordinates out
  do p=1, Nprob
    write(9,'(I8,1E17.8)') p, N_p(p)
  end do

  close(9)

  END SUBROUTINE Probe1D_nodes_gen
