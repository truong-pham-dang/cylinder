!======================================================================!
  PROGRAM Neu2TFlowS 
!------------------------------[Modules]-------------------------------!
  USE all_mod 
  USE gen_mod 
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, n, s
!--------------------------------[CVS]---------------------------------!
  character*80 rcs1,rcs2
  data rcs1/                                                        &
  '$Id: Neu2TFlowS.f90,v 1.15 2004/06/24 14:09:07 muhamed Exp $'/ 
  data rcs2/                                                        &
  '$Source: /home/muhamed/.CVSROOT/TFlowS/Neu2TFlowS/Neu2TFlowS.f90,v $'/
!======================================================================!

  call logo

!---- Test the precision
  open(90,FORM='UNFORMATTED',FILE='Neu2FlowS.real');
  write(90) 3.1451592
  close(90)  

  call ReadFluentNeu
  call TopolM
  call FindSides
  call Calc4
  
  call Connect3

  do n=1,NN
    NewN(n) = n 
  end do  
  do c=-NbC,NC
    NewC(c) = c 
  end do  
  do s=1,NS 
    NewS(s) = s
  end do  

!---- caount all the materials
  call CouMat

  call GenSav(0, NN, NC, NS, NbC)
  call GeoSav(0, NC, NS, NBC, 0, 0) 
  call TestLn(0, NN, NC, NS, NbC, 0)

!---- create output for Fluent
  NewC(-NBC-1) = -NBC-1
  call CasSav(0, NN, NC, NS+NSsh, NBC) ! Save grid for postprocessing
                                       ! with Fluent

!---- create 1D file (used for channel or pipe flow) 
  call Probe1D_nodes

!---- make eps figures
  write(*,*) 'Making three .eps cuts through the domain.'
  call EpsSav(y,z,x,Dy,Dz,'x')
  call EpsSav(z,x,y,Dz,Dx,'y')
  call EpsSav(x,y,z,Dx,Dy,'z')
 
  write(*,*) 'Making a 3D shaded .eps figure of the domain.'
  call EpsWho(NSsh)  ! Draw the domain with shadows

  END PROGRAM Neu2TFlowS 
