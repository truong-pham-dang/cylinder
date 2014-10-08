!======================================================================!
  SUBROUTINE CnsLoa
!----------------------------------------------------------------------!
! Reads:  NAME.cns                                                     !
! ~~~~~~                                                               !
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE sol_mod
  USE les_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER      :: c, s, it
  CHARACTER*80 :: nameIn
  CHARACTER*8  :: answer
!--------------------------------[CVS]---------------------------------!
!  $Id: CnsLoa.f90,v 1.7 2009/06/30 12:01:16 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CnsLoa.f90,v $  
!======================================================================!

  if(this < 2) write(*,*) '# Input problem name:'
  call ReadC(7,inp,tn,ts,te)  
  read(inp(ts(1):te(1)), '(A80)')  name

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!       Read the file with the      !
!     connections between cells     !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  call NamFil(THIS, nameIn, '.cns', len_trim('.cns'))  

  open(9, FILE=nameIn,FORM='UNFORMATTED')
  if(this < 2) write(*,*) '# Now reading the file:'
  if(this < 2) write(*,*) '  ', nameIn

!///// number of cells, boundary cells and sides
  read(9) NC                                 
  read(9) NbC
  read(9) NS
  read(9) c   ! NSsh is not used in Processor. 
  read(9) Nmat

!///// cell materials
  allocate (material(-NbC:NC))
  read(9) (material(c), c=1,NC)        
  read(9) (material(c), c=-1,-NBC,-1) 

!///// sides
  allocate (SideC(0:2,NS))
  read(9) (SideC(0,s), s=1,NS)
  read(9) (SideC(1,s), s=1,NS)
  read(9) (SideC(2,s), s=1,NS)

!///// boundary cells
  allocate (TypeBC(-NbC:-1)); TypeBC=0 
  allocate (bcmark(-NbC:-1))
  read(9) (bcmark(c), c=-1,-NbC, -1) 

!///// boundary copy cells
  allocate (CopyC(-NbC:-1))
  read(9) (CopyC(c), c=-1,-NbC, -1)   

  close(9)

!----- Type of the problem
  if(this  < 2) write(*,*) '# Type of problem: '
  call ReadC(7,inp,tn,ts,te)
  do it=1,tn
    read(inp(ts(it):te(it)),'(A8)')  answer
    call ToUppr(answer)
    if(answer == 'CHANNEL') then
      CHANNEL = YES
      if(this  < 2) write(*,*) '  CHANNEL -> Channel flow'
    else if(answer == 'PIPE') then
      PIPE = YES
      if(this  < 2) write(*,*) '  PIPE -> Pipe flow'
    else if(answer == 'JET') then
      JET = YES
      if(this  < 2) write(*,*) '  JET -> Impinging jet flow'
    else if(answer == 'TEST') then
      TEST = YES
      if(this  < 2) write(*,*) '  TEST -> ???'
    else if(answer == 'OTHER') then
      OTHER = YES
      if(this  < 2) write(*,*) '  OTHER -> All the other problems'
    else if(answer == 'HOT') then
      HOT = YES
      if(this  < 2) write(*,*) '  HOT -> Problems with temperature'
    else if(answer == 'XHOM') then
      XHOM = YES
      if(this  < 2) write(*,*) '  XHOM -> Homogeneous direction'
    else if(answer == 'YHOM') then
      YHOM = YES
      if(this  < 2) write(*,*) '  YHOM -> Homogeneous direction'
    else if(answer == 'ZHOM') then
      ZHOM = YES
      if(this  < 2) write(*,*) '  ZHOM -> Homogeneous direction'
    else if(answer == 'TGV') then
      TGV = YES
      if(this  < 2) write(*,*) '  TGV -> Taylor-Green Vortex test case'
    else if(answer == 'ROT') then
      ROT = YES
      if(this  < 2) write(*,*) '  ROT -> Rotation'
    else if(answer == 'BUDG') then
      BUDG = YES
      if(this  < 2) write(*,*) '  BUDG -> Test Budgets'
    else
      write(*,*) 'Error in input ! Exiting'
      stop
    endif
  end do

  END SUBROUTINE CnsLoa
