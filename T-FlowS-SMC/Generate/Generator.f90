!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                    ________  ______                                  !
!                   |        ||      \                                 !
!                   `--.  .--'|  ,-.  \________  ___                   !
!                      |  |___|  |_/  /  _  \  \/  /                   !
!                      |  |___|      /  _____\    /                    !
!                      |  |   |  |\  \  \____/    \                    !
!                      |__|   |__| \__\_____/__/\__\                   !
!                                                                      !
!                                                                      !
!           BLOCK-STRUCTURED 3D HEXAHEDRAL MESH GENERATOR              !
!                                +                                     !
!                  UNSTRUCTURED CELL REFINEMENT                        !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
!                                     Bojan Niceno                     !
!                                     Delft University of Technology   !
!                                     Section Heat Transfer            !
!                                     niceno@duttwta.wt.tn.tudelft.nl  !
!                                                                      !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                                                                      !
!  Block structure:                                                    !
!                                                                      !
!                               b6                                     !
!                               |                                      !
!                         8-----|---------6                            !
!                        /|     |        /|                            !
!                       / |     +   b3  / |                            !
!                      /  |        /   /  |                            !
!                     7---------------5   |                            !
!                     |   |     |/    |   |                            !
!              b4---- | +-------o-----| +-------b2                     !
!                     |   |    /|     |   |                            !
!                     |   4---/-|-----|---2                            !
!                     |  /   /  |     |  /                             !
!                     | /   b5  +     | /                              !
!                     |/              |/                               !
!                     3---------------1                                !
!                               |                                      !
!                               b1                                     !
!  Faces are defined as:
!
!     I   :
!     II  :
!     III :
!     IV  :
!     V   :
!     VI  :
!
!  Local coordinate directions:                                        !
!                                                                      !
!     i: 1 -> 2                                                        !
!     j: 1 -> 3                                                        !
!     k: 1 -> 5                                                        !
!                                                                      !
!  Notes:                                                              !
!                                                                      !
!    - can't handle domains with less then 3x3x3 cells properly        !
!    - nodes of a cell (and block) are deonoted with numbers 1 - 8     !
!    - neighbouring cells are denoted with c1 - c6                     !
!    - local coordinare directions (for blocks) are defined with:      !
!                                                                      !
!======================================================================!
  PROGRAM Generator
!----------------------------------------------------------------------!
!   Block structured mesh generation and unstructured cell refinement. *
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c,s,n
!--------------------------------[CVS]---------------------------------!
!  $Id: Generator.f90,v 1.20 2002/10/30 16:29:21 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Generate/Generator.f90,v $  
!======================================================================!

!---- Test the precision
  open(90,FORM='UNFORMATTED',FILE='Generator.real'); 
  write(90) 3.1451592
  close(90)

  call logo

  call IniGen 

  call GenLoa
  call Calc1  
  call Fuzion
  call PeriBC
  call CopyBC

  call TopSys(.FALSE.)  ! trial run 
  call Calc2(.FALSE.)
  call Smooth 

  call Mark

  call TopSys(.TRUE.) ! real run
  call Calc2(.TRUE.)

!---- prepare for saving
  do n=1,NN
    NewN(n)=n
  end do
  do c=-NBC,NC
    NewC(c)=c
  end do
  do s=1,NS
    NewS(s)=s
  end do

!---- count the materials in the grid
  call CouMat

!---- save the grid
  call GenSav(0, NN, NC)            ! Save grid for postprocessing
  call GeoSav(0, NC, NS, NBC, 0, 0) ! Saved data for processing

  call TestLn(0, NN, NC, NS, NbC, 0)

!---- save the 1D probe (good for the channel flow)
  call Probe1D_nodes_gen

!---- save the 2D probe (good for the channel flow)
  call Probe2D

!---- create output for Fluent
  NewC(-NBC-1) = -NBC-1
  call CasSav(0, NN, NC, NS+NSsh) ! Save grid for postprocessing
                                  ! with Fluent
!---- make eps figures
  call EpsSav(y,z,x,Dy,Dz,'x') 
  call EpsSav(z,x,y,Dz,Dx,'y') 
  call EpsSav(x,y,z,Dx,Dy,'z') 

  call EpsWho(NSsh)  ! Draw the domain with shadows

!---- write something on the screen
  call PrintG

  END PROGRAM Generator 
