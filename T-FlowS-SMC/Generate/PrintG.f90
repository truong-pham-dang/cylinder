!======================================================================!
  SUBROUTINE PrintG 
!----------------------------------------------------------------------!
!   Prints some statistical data about the grid on the standard output !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: i, j, k, numb, nonz, stencw
!--------------------------------[CVS]---------------------------------!
!  $Id: PrintG.f90,v 1.7 2002/10/30 16:29:22 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Generate/PrintG.f90,v $  
!======================================================================!

  write(*,*) '==============='
  write(*,*) 'Grid statistics'
  write(*,*) '==============='
  write(*,*) '  number of nodes         :', NN
  write(*,*) '  number of cells         :', NC
  write(*,*) '  number of sides         :', NS
  write(*,*) '  number of boundary cells:', NbC
  write(*,*) '----------------------------------'

!----- Find the number of non zero entries
  nonz=0
  do i = 1,NC
    stencw=1            ! prije je bilo 0
    do j=1,24
      if( CellC(i,j)  > 0 ) stencw=stencw+1
    end do
    nonz = nonz + stencw
  end do

  write(*,*) '  number of non zero matrix entries:', nonz
  write(*,*) '  average stencil size:', REAL(nonz)/REAL(NC)
  write(*,*) '  max number of nodes and cells:',   MAXN
  write(*,*) '  max number of boundary cells:',    MAXB
  write(*,*) '----------------------------------'

!----- Neighbours
  do j=1,24
    numb=0
    do i=1,NC
      stencw=0
      do k=1,24
	if( CellC(i,k)  > 0 ) stencw=stencw+1
      end do
      if(stencw  ==  j) numb=numb+1
    end do
    if(numb /= 0) then
    write(*,*) '  number of cells with ',j, ' neighbours: ',numb
    endif
  end do

!----- Twins
  do j=1,8
    numb=0
    do i=1,NN
      if(TwinN(i,0)  ==  j) numb=numb+1
    end do
    if(numb /= 0) then
    write(*,*) '  number of nodes with ',j, ' twins: ',numb
    endif 
  end do

  END SUBROUTINE PrintG
