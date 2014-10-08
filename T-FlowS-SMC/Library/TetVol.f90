!======================================================================!
  REAL FUNCTION TetVol(xA,yA,zA,xB,yB,zB,xC,yC,zC,xD,yD,zD)
!----------------------------------------------------------------------!
!   Returns the volume of tethraedra spanned with A, B, C and D.       !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  REAL :: xA,yA,zA,xB,yB,zB,xC,yC,zC,xD,yD,zD
!--------------------------------[CVS]---------------------------------!
!  $Id: TetVol.f90,v 1.4 2002/10/30 16:29:33 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Library/TetVol.f90,v $  
!======================================================================!
!                                                                      !
!   The order of nodes (A,B,C and D) DOES matters.                     !
!                                                                      !
!                D-----C                                               !
!               / \  . |                                               !
!              /   \   |                                               !
!             /  .  \  |    I am not 100% sure that the figure is OK   !
!            / .     \ |                                               !
!           /.        \|                                               !
!          A-----------B                                               !
!                                                                      !
!----------------------------------------------------------------------!

  TetVol=( ((yB-yA)*(zC-zA)-(yC-yA)*(zB-zA))*(xD-xA) +              &
	   ((xC-xA)*(zB-zA)-(xB-xA)*(zC-zA))*(yD-yA) +              &
	   ((xB-xA)*(yC-yA)-(xC-xA)*(yB-yA))*(zD-zA) )/6.0

  END FUNCTION TetVol
