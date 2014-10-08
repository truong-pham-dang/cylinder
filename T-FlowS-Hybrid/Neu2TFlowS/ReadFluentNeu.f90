!======================================================================!
  SUBROUTINE ReadFluentNeu() 
!----------------------------------------------------------------------!
! Reads the Fluents (Gambits) neutral file format.                     !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod 
  USE gen_mod 
  USE neu_mod 
  USE par_mod 
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  CHARACTER           :: Line*300
  CHARACTER           :: nameIn*130
  INTEGER             :: i, j, dum1, dum2, dum3
  INTEGER,ALLOCATABLE :: temp(:)
  INTEGER             :: c, dir, type
!--------------------------------[CVS]---------------------------------!
  character*80 rcs1,rcs2
  data rcs1/                                                        &
  '$Id: ReadFluentNeu.f90,v 1.14 2005/01/25 12:45:44 muhamed Exp $'/ 
  data rcs2/                                                        &
  '$Source: /home/muhamed/.CVSROOT/T-Rex/Neu2Trex/ReadFluentNeu.f90,v $'/
!======================================================================!

  write(*,*) 'Enter the Fluent''s (*.NEU) file name (without ext.):'
  read(*,*) name

  nameIn = name
  nameIn(len_trim(name)+1:len_trim(name)+4) = '.neu'

  open(9,FILE=nameIn)
  write(*,*) 'Now reading the file: ', nameIn

!---- skip first 6 lines
  do i=1,6
    call Read9Ln(Line,tn,ts,te)
  end do 

!---- read the line which contains usefull information  
  call Read9Ln(Line,tn,ts,te)
  read(Line(ts(1):te(1)),*) NN  
  read(Line(ts(2):te(2)),*) NC
  read(Line(ts(3):te(3)),*) NBloc
  read(Line(ts(4):te(4)),*) NBS

  write(*,*) 'Total number of nodes:  ', NN
  write(*,*) 'Total number of cells:  ', NC
  write(*,*) 'Total number of blocks: ', NBloc
  write(*,*) 'Total number of boundary sections: ', NBS

!---- count the boundary cells
  NbC = 0
  do 
    call Read9Ln(Line,tn,ts,te)
    if( Line(ts(1):te(1)) == 'BOUNDARY' ) then
      do j=1,NBS
        if(j>1) call Read9Ln(Line,tn,ts,te) ! BOUNDARY CONDITIONS
        call Read9Ln(Line,tn,ts,te)
        read(Line(ts(3):te(3)),*) dum1  
        NBc = NBc + dum1 
        do i=1,dum1
          read(9,*) c, dum2, dir
        end do
        call Read9Ln(Line,tn,ts,te)         ! ENDOFSECTION
      end do
      write(*,*) 'Total number of boundary cells: ', NbC
      go to 1
    end if
  end do 

1 rewind(9)

!---- skip first 7 lines
  do i=1,7
    call Read9Ln(Line,tn,ts,te)
  end do 

!---- allocate memory ==> ATTENTION: NO CHECKING
  allocate(x(NN)); x=0
  allocate(y(NN)); y=0
  allocate(z(NN)); z=0

  allocate(material(-NbC:NC));     material=0 
  allocate(BCtype(NC,6));          BCtype=0
  allocate(BCmark(-NbC-1:-1));     BCmark=0
  allocate(CellN(-NbC-1:NC,-1:8)); CellN=0
  allocate(SideC(0:2,NC*5));       SideC=0    
  allocate(SideN(NC*5,0:4));       SideN=0
  Ncopy = 1000000 
  allocate(CopyC(-NbC:-1));  CopyC=0
  allocate(CopyS(2,Ncopy));  CopyS=0

  allocate(NewN(NN));        NewN=0  
  allocate(NewC(-NbC-1:NC)); NewC=0  
  allocate(NewS(NC*5));      NewS=0  

  allocate(proces(0:NC)); proces=0

  allocate(temp(NC)); temp=0

!--- skip one line 
  call Read9Ln(Line,tn,ts,te)

!%%%%%%%%%%%%%%%%%%%!
! NODAL COORDINATES !
!%%%%%%%%%%%%%%%%%%%!
!---- read the nodal coordinates
  call Read9Ln(Line,tn,ts,te)          ! NODAL COORDINATES
  do i=1,NN
    call Read9Ln(Line,tn,ts,te)
    read(Line(ts(2):te(2)),*) x(i)  
    read(Line(ts(3):te(3)),*) y(i)
    read(Line(ts(4):te(4)),*) z(i) 
  end do
  call Read9Ln(Line,tn,ts,te)          ! ENDOFSECTION

!%%%%%%%%%%%%%%%%!
! ELEMENTS/CELLS !
!%%%%%%%%%%%%%%%%!
!---- read nodes of each cell
  call Read9Ln(Line,tn,ts,te)          ! ELEMENTS/CELLS
  do i=1,NC
    read(9,'(I8,1X,I2,1X,I2,1X,7I8:/(15X,7I8:))') dum1, dum2, CellN(i,0), (CellN(i,j), j=1,CellN(i,0))
  end do
  call Read9Ln(Line,tn,ts,te)          ! ENDOFSECTION
!---- read block data
  do j=1,NBloc
    call Read9Ln(Line,tn,ts,te)        ! ELEMENT GROUP
    call Read9Ln(Line,tn,ts,te)
    read(Line(ts(4):te(4)),'(I10)') dum1  
    call Read9Ln(Line,tn,ts,te)        ! block*
    call Read9Ln(Line,tn,ts,te)        ! 0
    read(9,'(10I8)') (temp(i), i=1,dum1)
    do i=1,dum1
      material(temp(i)) = j
    end do
    call Read9Ln(Line,tn,ts,te)        ! ENDOFSECTION
  end do

!%%%%%%%%%%%%%%%%%%%%%!
! BOUNDARY CONDITIONS !
!%%%%%%%%%%%%%%%%%%%%%!
  write(*,*) 'NBS=', NBS 
  do j=1,NBS
    call Read9Ln(Line,tn,ts,te)        ! BOUNDARY CONDITIONS
    call Read9Ln(Line,tn,ts,te)
    if( Line(ts(1):te(1)) == 'one')      type =  1
    if( Line(ts(1):te(1)) == 'two')      type =  2
    if( Line(ts(1):te(1)) == 'three')    type =  3
    if( Line(ts(1):te(1)) == 'four')     type =  4
    if( Line(ts(1):te(1)) == 'five')     type =  5
    if( Line(ts(1):te(1)) == 'six')      type =  6
    if( Line(ts(1):te(1)) == 'seven')    type =  7
    if( Line(ts(1):te(1)) == 'eight')    type =  8
    if( Line(ts(1):te(1)) == 'nine')     type =  9
    if( Line(ts(1):te(1)) == 'ten')      type = 10
    if( Line(ts(1):te(1)) == 'oneone')   type = 11
    if( Line(ts(1):te(1)) == 'onetwo')   type = 12
    if( Line(ts(1):te(1)) == 'onethree') type = 13
    if( Line(ts(1):te(1)) == 'onefour')  type = 14
    if( Line(ts(1):te(1)) == 'onefive')  type = 15
    if( Line(ts(1):te(1)) == 'onesix')   type = 16
    if( Line(ts(1):te(1)) == 'oneseven') type = 17
    if( Line(ts(1):te(1)) == 'oneeight') type = 18
    if( Line(ts(1):te(1)) == 'onenine')  type = 19
    write(*,*)  Line(ts(1):te(1)), type 
    read(Line(ts(3):te(3)),'(I8)') dum1  
    do i=1,dum1
      read(9,*) c, dum2, dir
      BCtype(c,dir) = type 
    end do
    call Read9Ln(Line,tn,ts,te)        ! ENDOFSECTION
  end do

  close(9)

  END SUBROUTINE ReadFluentNeu
