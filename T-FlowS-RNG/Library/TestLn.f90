!======================================================================!
  SUBROUTINE TestLn(sub, NNsub, NCsub, NSsub, NBCsub, NBFsub) 
!----------------------------------------------------------------------!
!   Creates the file "test.link.gmv" to check the cell connections.    !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod 
  USE gen_mod 
  USE par_mod 
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: sub, NNsub, NCsub, NSsub, NBCsub, NBFsub
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: n,c,c1,c2,s 
  INTEGER   :: NSsubNonPer, NSsubPer
  CHARACTER :: namOut*80
!--------------------------------[CVS]---------------------------------!
!  $Id: TestLn.f90,v 1.12 2002/10/30 16:29:33 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Library/TestLn.f90,v $  
!======================================================================!

!<<<<<<<<<<<<<<<<<<<<<<<<<!
!     create GMV file     !
!<<<<<<<<<<<<<<<<<<<<<<<<<!--------------------------------------------!
!     Big changes have taken place in this section. Links between      !
!     the computational cells have been introduced as aditional        !
!     cells of general type. Cell centers are introduced as aditional  !
!     nodes with numbers NN+1 to NN+NC. Material of this links is      !
!     different than from the cells, so that they can be visualised    !
!     more easily in GMV. I still have to check this changes with      !
!     divide4 program.                                                 !
!----------------------------------------------------------------------!
  namOut = name         

  call NamFil(sub, namOut, '.ln.gmv', len_trim('.ln.gmv'))
  open(9, FILE=namOut)
  write(6, *) 'Now creating the file:', namOut

!---------------!
!     nodes     !
!---------------!
  write(9,'(A14)') 'gmvinput ascii'          !  start of GMV file
  write(9,*) 'nodes', NNsub + NCsub + NBCsub + NBFsub

!---- x
  do n=1,NN
    if(NewN(n) /= 0) write(9, '(1PE14.7)') x(n)
  end do
  do c=1,NC
    if(NewC(c)  > 0) write(9, '(1PE14.7)') xc(c)
  end do 
  do c=-1,-NBC,-1
    if(NewC(c) /= 0) write(9, '(1PE14.7)') xc(c)
  end do 
  do c=1,NBFsub
    write(9, '(1PE14.7)') xc(BuReIn(c))
  end do

!---- y
  do n=1,NN
    if(NewN(n) /= 0) write(9, '(1PE14.7)') y(n)
  end do
  do c=1,NC
    if(NewC(c)  > 0) write(9, '(1PE14.7)') yc(c)
  end do 
  do c=-1,-NBC,-1
    if(NewC(c) /= 0) write(9, '(1PE14.7)') yc(c)
  end do 
  do c=1,NBFsub
    write(9, '(1PE14.7)') yc(BuReIn(c))
  end do

!---- z
  do n=1,NN
    if(NewN(n) /= 0) write(9, '(1PE14.7)') z(n)
  end do
  do c=1,NC
    if(NewC(c)  > 0) write(9, '(1PE14.7)') zc(c)
  end do 
  do c=-1,-NBC,-1
    if(NewC(c) /= 0) write(9, '(1PE14.7)') zc(c)
  end do 
  do c=1,NBFsub
    write(9, '(1PE14.7)') zc(BuReIn(c))
  end do

!----------------------!
!     cell section     !
!----------------------!

!---- regular (ordinary) cells
  write(9,*) 'cells', NCsub + NSsub + NBFsub ! + NBFsub
  do c=1,NC
    if(NewC(c)  > 0) then
      if(CellN(c,0) == 8) then
        write(9,*) 'hex 8'
        write(9,*)                                                  &
	      NewN(CellN(c,1)), NewN(CellN(c,2)),                   &
	      NewN(CellN(c,4)), NewN(CellN(c,3)),                   &
	      NewN(CellN(c,5)), NewN(CellN(c,6)),                   &
	      NewN(CellN(c,8)), NewN(CellN(c,7))
      else if(CellN(c,0) == 6) then
        write(9,*) 'prism 6'
        write(9,*)                                                  &
	      NewN(CellN(c,1)), NewN(CellN(c,2)),                   &
	      NewN(CellN(c,3)), NewN(CellN(c,4)),                   &
	      NewN(CellN(c,5)), NewN(CellN(c,6))
      else if(CellN(c,0) == 4) then
        write(9,*) 'tet 4'
        write(9,*)                                                  &
	      NewN(CellN(c,1)), NewN(CellN(c,2)),                   &
	      NewN(CellN(c,3)), NewN(CellN(c,4))
      else if(CellN(c,0) == 5) then
        write(9,*) 'pyramid 5'
        write(9,*)                                                  &
              NewN(CellN(c,5)), NewN(CellN(c,1)),                   &
              NewN(CellN(c,2)), NewN(CellN(c,4)),                   &
              NewN(CellN(c,3))      
      else
        write(*,*) 'Unsupported cell type with ', CellN(c,0), 'nodes. '
        write(*,*) 'Exiting !'
        stop
      end if
    end if 
  end do  

!---- physical links; non-periodic
  NSsubNonPer = 0 
  do s=1, NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)

    if( (NewS(s) > 0) .and. (NewS(s) <= NSsub) ) then

      if( (Sx(s) * (xc(c2)-xc(c1) )+               &
           Sy(s) * (yc(c2)-yc(c1) )+               &
           Sz(s) * (zc(c2)-zc(c1) ))  > 0.0 ) then 

        NSsubNonPer = NSsubNonPer + 1

        c1 = NewC(SideC(1,s))
        c2 = NewC(SideC(2,s))
        if( c2  > 0 ) then
	  write(9,*) 'general 1'
	  write(9,*) '  2'
	  write(9,*) NNsub+c1, NNsub+c2
        else
	  write(9,*) 'general 1'
	  write(9,*) '  2'
	  write(9,*) NNsub+c1, NNsub+NCsub-c2
        end if
      end if

    end if
  end do  

!---- physical links; periodic
  NSsubPer    = 0 
  do s=1, NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)

    if( (NewS(s) > 0) .and. (NewS(s) <= NSsub) ) then

      if( (Sx(s) * (xc(c2)-xc(c1) )+               &
           Sy(s) * (yc(c2)-yc(c1) )+               &
           Sz(s) * (zc(c2)-zc(c1) ))  < 0.0 ) then 

        NSsubPer = NSsubPer + 1

        c1 = NewC(SideC(1,s))
        c2 = NewC(SideC(2,s))
        if( c2  > 0 ) then
	  write(9,*) 'general 1'
	  write(9,*) '  2'
	  write(9,*) NNsub+c1, NNsub+c2
        else
	  write(9,*) 'general 1'
	  write(9,*) '  2'
	  write(9,*) NNsub+c1, NNsub+NCsub-c2
        end if
      end if

    end if
  end do  

  write(*,*) 'Non-periodic links:', NSsubNonPer
  write(*,*) 'Periodic links    :', NSsubPer
  write(*,*) 'Number of sides   :', NSsub

!---- interprocessor links
  do c=1,NBFsub
    c1 = BuSeIn(c) 
    write(9,*) 'general 1'
    write(9,*) '  2'
    write(9,*) NNsub+c1, NNsub+NCsub+NBCsub+c
  end do  

!--------------------------!
!     material section     !
!--------------------------!
  write(9,*) 'material  4  0'
  write(9,*) 'cells'
  write(9,*) 'links'
  write(9,*) 'periodic'
  write(9,*) 'buffers'

  do c=1, NCsub
    write(9,*) 1
  end do
  do s=1, NSsubNonPer
    write(9,*) 2
  end do
  do s=1, NSsubPer
    write(9,*) 3
  end do
  do s=1,NBFsub
    write(9,*) 4
  end do

  write(9,'(A6)') 'endgmv'            !  end the GMV file
  close(9)

  END SUBROUTINE TestLn
