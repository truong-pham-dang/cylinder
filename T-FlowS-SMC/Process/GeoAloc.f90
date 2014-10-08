!======================================================================!
  SUBROUTINE GeoAloc
!----------------------------------------------------------------------!
! Alocates memory for geometrical quantities.                          !
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE sol_mod
  USE les_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!--------------------------------[CVS]---------------------------------!
!  $Id: GeoAloc.f90,v 1.6 2008/11/19 14:50:07 IUS\mhadziabdic Exp $
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/GeoAloc.f90,v $
!======================================================================!

!----- variables defined in all.h90:
  allocate (xc(-NbC:NC)); xc=0.0       
  allocate (yc(-NbC:NC)); yc=0.0   
  allocate (zc(-NbC:NC)); zc=0.0  
  allocate (Sx(NS)); Sx=0.0       
  allocate (Sy(NS)); Sy=0.0       
  allocate (Sz(NS)); Sz=0.0  
  allocate (Dx(NS)); Dx=0.0       
  allocate (Dy(NS)); Dy=0.0       
  allocate (Dz(NS)); Dz=0.0  
  allocate (xsp(NS)); xsp=0.0       
  allocate (ysp(NS)); ysp=0.0       
  allocate (zsp(NS)); zsp=0.0  
  allocate (volume(-NbC:NC)); volume=0.  
  allocate (delta(-NbC:NC));  delta=0.  
  allocate (WallDs(-NbC:NC)); WallDs=0.       
  allocate (f(NS));  f=0.0  
  allocate (fF(NS)); fF=0.0  

  allocate (a1(-NbC:NC));  a1=0. !IZMJENA 
  allocate (a2(-NbC:NC));  a2=0.  
!----- variables defined in sol.h90:
  allocate (D(NC)); D=0
  allocate (p1(-NbC:NC)); p1=0
  allocate (p2(-NbC:NC)); p2=0
  allocate (q1(NC)); q1=0
  allocate (q2(NC)); q2=0
  allocate (r2(NC)); r2=0
  allocate (u1(NC)); u1=0
  allocate (u2(NC)); u2=0
  allocate (v1(NC)); v1=0
  allocate (v2(NC)); v2=0

!---- variables defined in pro_mod.h90:
  allocate (Acol(NC+1));     Acol=0
  allocate (Adia(NC));       Adia=0
  allocate (Asave(-NbC:NC)); Asave=0
  allocate (Abou(-NbC:-1));  Abou=0
  allocate (b(NC));          b=0

  allocate (Scoef(NS)); Scoef=0.
  allocate (SidAij(2,NS)); SidAij=0

  allocate (xp(Nmat));    xp   =0.0
  allocate (yp(Nmat));    yp   =0.0
  allocate (zp(Nmat));    zp   =0.0
  allocate (AreaX(Nmat)); AreaX=0.0
  allocate (AreaY(Nmat)); AreaY=0.0
  allocate (AreaZ(Nmat)); AreaZ=0.0

!---- variables defined in par_mod.h90:
  allocate (BufInd(-NbC:-1)); BufInd=0

!??????????????????????????????????????????!
!     Is there enough allocated memory     !
!??????????????????????????????????????????!
! Do something !  

  END SUBROUTINE GeoAloc
