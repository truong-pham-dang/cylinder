!======================================================================!
  SUBROUTINE GenSav(sub, NNsub, NCsub)
!----------------------------------------------------------------------!
! Writes: NAME.gmv, NAME.faces.gmv, NAME.shadow.gmv                    !
! ~~~~~~~                                                              ! 
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: sub, NNsub, NCsub, NmaterBC
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: c,  c1,  c2,  n, s
  CHARACTER :: namOut*80
!--------------------------------[CVS]---------------------------------!
!  $Id: GenSav.f90,v 1.24 2002/10/30 16:29:33 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Library/GenSav.f90,v $             
!======================================================================!
!   See also: number                                                   !
!----------------------------------------------------------------------!

!<<<<<<<<<<<<<<<<<<<<<<<<<!
!                         !
!     create GMV file     !
!                         !
!<<<<<<<<<<<<<<<<<<<<<<<<<!
  call NamFil(sub, namOut, '.gmv', len_trim('.gmv'))
  open(9, FILE=namOut)
  write(6, *) 'Now creating the file:', namOut

!---------------!
!     start     !
!---------------!
  write(9,'(A14)') 'gmvinput ascii'

!---------------!
!     nodes     !
!---------------!
  write(9,*) 'nodes', NNsub

  do n=1,NN
    if(NewN(n) /= 0) write(9, '(1PE14.7)') x(n)
  end do
  do n=1,NN
    if(NewN(n) /= 0) write(9, '(1PE14.7)') y(n)
  end do
  do n=1,NN
    if(NewN(n) /= 0) write(9, '(1PE14.7)') z(n)
  end do

!----------------------!
!     cell section     !
!----------------------!
  write(9,*) 'cells', NCsub
  do c=1,NC
    if(NewC(c) /= 0) then
      if(CellN(c,0) == 8) then
        write(9,*) 'hex 8'
        write(9,'(8I9)')                        &
  	  NewN(CellN(c,1)), NewN(CellN(c,2)),   &
	  NewN(CellN(c,4)), NewN(CellN(c,3)),   &
	  NewN(CellN(c,5)), NewN(CellN(c,6)),   &
	  NewN(CellN(c,8)), NewN(CellN(c,7))
      else if(CellN(c,0) == 6) then
        write(9,*) 'prism 6'
        write(9,'(6I9)')                        &
  	  NewN(CellN(c,1)), NewN(CellN(c,2)),   &
	  NewN(CellN(c,3)), NewN(CellN(c,4)),   &
	  NewN(CellN(c,5)), NewN(CellN(c,6))
      else if(CellN(c,0) == 4) then
        write(9,*) 'tet 4'
        write(9,'(4I9)')                        &
  	  NewN(CellN(c,1)), NewN(CellN(c,2)),   &
	  NewN(CellN(c,3)), NewN(CellN(c,4))
      else if(CellN(c,0) == 5) then
        write(9,*) 'pyramid 5'
        write(9,'(5I9)')                        &
  	  NewN(CellN(c,5)), NewN(CellN(c,1)),   &
	  NewN(CellN(c,2)), NewN(CellN(c,4)),   &
	  NewN(CellN(c,3))
      else
        write(*,*) 'Unsupported cell type ', CellN(c,0), ' nodes.'
        write(*,*) 'Exiting'
        stop 
      end if 
    end if
  end do  

!->>>  write(9,*) 'variables'
!->>>  write(9,*) 'wallds 0'
!->>>  do c=1,NC
!->>>    write(9,*) WallDs(c)
!->>>  end do
!->>>  write(9,*) 'endvars'
 
!--------------------------!
!     material section     !
!--------------------------!
  write(9,'(A10,2I5)') 'materials', Nmat, 0
  do n=1,1024
    if(Mater(n)) write(9,*) n 
  end do	
  do c=1,NC
    if(NewC(c) /= 0) then
      write(9,*) material(c)
    end if
  end do	

  write(9,'(A6)') 'endgmv'            !  end the GMV file
  close(9)

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!                                            !
!     create boundary condition GMV file     !
!                                            !
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  if(sub /= 0) return

  call NamFil(sub, namOut, '.faces.gmv', len_trim('.faces.gmv'))
  open(9, FILE=namOut)
  write(6, *) 'Now creating the file:', namOut

!---------------!
!     start     !
!---------------!
  write(9,'(A14)') 'gmvinput ascii'

!---------------!
!     nodes     !
!---------------!
  write(9,*) 'nodes', NNsub

  do n=1,NN
    write(9, '(1PE14.7)') x(n)
  end do
  do n=1,NN
    write(9, '(1PE14.7)') y(n)
  end do
  do n=1,NN
    write(9, '(1PE14.7)') z(n)
  end do

!----------------------!
!     cell section     !
!----------------------!

!---- count the cell faces on the periodic boundaries
  write(9,*) 'cells', NS
  do s=1,NS
    if(SideN(s,0) == 4) then
      write(9,*) 'quad 4'
      write(9,'(4I9)')           &
        SideN(s,1), SideN(s,2),  &
        SideN(s,3), SideN(s,4)
    else if(SideN(s,0) == 3) then
      write(9,*) 'tri 3'
      write(9,'(3I9)')           &
        SideN(s,1), SideN(s,2),  &
        SideN(s,3)
    else
      write(*,*) 'Unsupported cell type ', CellN(c,0), ' nodes.'
      write(*,*) 'Exiting'
      stop 
    end if
  end do  

!--------------------------!
!     material section     !
!--------------------------!
  NmaterBC = 0
  do s=1,NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)
    if( c2 < 0 ) then 
      NmaterBC=max(NmaterBC,bcmark(c2))
    end if
  end do  

  write(9,*) 'materials', NmaterBC + 1, 0
  do n = 1, NmaterBC + 1
    write(9,*) n 
  end do	

  do s=1,NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)
    if( c2 < 0 ) then 
      write(9,*) bcmark(c2) 
    else 
      write(9,*) NmaterBC+1 
    end if
  end do

  write(9,'(A6)') 'endgmv'            !  end the GMV file
  close(9)

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!                            !
!     create shadow file     !
!                            !
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  if(sub /= 0) return

  call NamFil(sub, namOut, '.shadow.gmv', len_trim('.shadow.gmv'))
  open(9, FILE=namOut)
  write(6, *) 'Now creating the file:', namOut

  do s=NS+1,NS+NSsh
    write(9,*) SideN(s,0) 
    if(SideN(s,0)==3) then
      write(9,*) SideN(s,1), SideN(s,2), SideN(s,3),             &
                 SideC(1,s), SideC(2,s) 
    else if(SideN(s,0)==4) then
      write(9,*) SideN(s,1), SideN(s,2), SideN(s,3), SideN(s,4), &
                 SideC(1,s), SideC(2,s) 
    end if
  end do  

  close(9)

  END SUBROUTINE GenSav
