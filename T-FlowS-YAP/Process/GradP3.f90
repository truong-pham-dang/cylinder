!======================================================================!
  SUBROUTINE GradP3(PHI, PHI_x, PHI_y, PHI_z)
!----------------------------------------------------------------------!
! Calculates gradient of generic variable PHI. PHI may stand either    !
! for pressure (P) or pressure corrections (PP). This procedure also   !
! handles different materials.                                         ! 
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  REAL :: PHI(-NbC:NC),                                             &
          PHI_x(-NbC:NC),                                           &
          PHI_y(-NbC:NC),                                           &
          PHI_z(-NbC:NC)
!-------------------------------[Locals]-------------------------------!
  INTEGER :: s, c, c1, c2
  REAL    :: PHI_s, xs, ys, zs 
!--------------------------------[CVS]---------------------------------!
!  $Id: GradP3.f90,v 1.8 2008/12/10 14:40:18 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/GradP3.f90,v $      
!======================================================================!
						 
  call Exchng(PHI)

  Ps = 0.0

  do c=1,NC
    PHI_x(c)=0.0
    PHI_y(c)=0.0
    PHI_z(c)=0.0
  end do

!----------------------------------------------------------------!
!     First step: without any wall influence, except outflow     !
!----------------------------------------------------------------!
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2 > 0                            .or. &
       c2 < 0 .and. TypeBC(c2) == BUFFER .or. &
       c2 < 0 .and. TypeBC(c2) == OUTFLOW) then  
      if( StateMat(material(c1))==FLUID .and. &
          StateMat(material(c2))==FLUID ) then  
        PHI_s = f(s)*PHI(c1)+(1.0-f(s))*PHI(c2)  
        PHI_x(c1) = PHI_x(c1) + PHI_s * Sx(s)
        PHI_y(c1) = PHI_y(c1) + PHI_s * Sy(s)
        PHI_z(c1) = PHI_z(c1) + PHI_s * Sz(s)
        PHI_x(c2) = PHI_x(c2) - PHI_s * Sx(s)
        PHI_y(c2) = PHI_y(c2) - PHI_s * Sy(s)
        PHI_z(c2) = PHI_z(c2) - PHI_s * Sz(s)
      end if
    end if
  end do

  do c=1,NC
    if(StateMat(material(c))==FLUID) then
      PHI_x(c)=PHI_x(c)/volume(c)
      PHI_y(c)=PHI_y(c)/volume(c)
      PHI_z(c)=PHI_z(c)/volume(c)
    end if
  end do

!----------------------------------------------------------------!
!     Second step: extrapolate to boundaries, except outflow     !
!----------------------------------------------------------------!
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2 < 0               .and. &
       TypeBC(c2) /= BUFFER .and. &
       TypeBC(c2) /= OUTFLOW) then  
      if(StateMat(material(c1))==FLUID) then
        Ps(s) = (   PHI(c1)                     +     &
                    PHI_x(c1) * (xc(c2)-xc(c1)) +     &
                    PHI_y(c1) * (yc(c2)-yc(c1)) +     &
                    PHI_z(c1) * (zc(c2)-zc(c1))   ) / &
              ( 1.0 - (  Sx(s) * (xc(c2)-xc(c1))      & 
                       + Sy(s) * (yc(c2)-yc(c1))      &
                       + Sz(s) * (zc(c2)-zc(c1))  ) / volume(c1)  )
      end if
    end if

!---- Handle two materials
    if(c2 > 0 .or. c2 < 0 .and. TypeBC(c2) == BUFFER) then  
      if( StateMat(material(c1))==FLUID .and. &
          StateMat(material(c2))==SOLID ) then  
        xs = xsp(s) 
        ys = ysp(s) 
        zs = zsp(s) 
        Ps(s) = (   PHI(c1)                 +     &
                    PHI_x(c1) * (xs-xc(c1)) +     &
                    PHI_y(c1) * (ys-yc(c1)) +     &
                    PHI_z(c1) * (zs-zc(c1))   ) / &
              ( 1.0 - (  Sx(s) * (xs-xc(c1))      & 
                       + Sy(s) * (ys-yc(c1))      &
                       + Sz(s) * (zs-zc(c1))  ) / volume(c1)  )
      end if
      if( StateMat(material(c1))==SOLID .and. &
          StateMat(material(c2))==FLUID ) then  
        xs = xsp(s) 
        ys = ysp(s) 
        zs = zsp(s) 
        Ps(s) = (   PHI(c2)                 +     &
                    PHI_x(c2) * (xs-xc(c2)) +     &
                    PHI_y(c2) * (ys-yc(c2)) +     &
                    PHI_z(c2) * (zs-zc(c2))   ) / &
              ( 1.0 + (  Sx(s) * (xs-xc(c2))      & 
                       + Sy(s) * (ys-yc(c2))      &
                       + Sz(s) * (zs-zc(c2))  ) / volume(c2)  )
      end if
    end if ! c2 < 0
  end do

!-------------------------------------------------!
!     Third step: compute the final gradients     !
!-------------------------------------------------!
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2 < 0               .and. &
       TypeBC(c2) /= BUFFER .and. &
       TypeBC(c2) /= OUTFLOW) then  
      PHI_x(c1) = PHI_x(c1) + Ps(s) * Sx(s)/volume(c1)
      PHI_y(c1) = PHI_y(c1) + Ps(s) * Sy(s)/volume(c1)
      PHI_z(c1) = PHI_z(c1) + Ps(s) * Sz(s)/volume(c1)
    end if

!---- Handle two materials
    if(c2 > 0 .or. c2 < 0 .and. TypeBC(c2) == BUFFER) then  
      if( StateMat(material(c1))==FLUID .and. &
          StateMat(material(c2))==SOLID ) then  
        PHI_x(c1) = PHI_x(c1) + Ps(s) * Sx(s)/volume(c1)
        PHI_y(c1) = PHI_y(c1) + Ps(s) * Sy(s)/volume(c1)
        PHI_z(c1) = PHI_z(c1) + Ps(s) * Sz(s)/volume(c1)
      end if 
      if( StateMat(material(c1))==SOLID .and. &
          StateMat(material(c2))==FLUID ) then  
        PHI_x(c2) = PHI_x(c2) - Ps(s) * Sx(s)/volume(c2)
        PHI_y(c2) = PHI_y(c2) - Ps(s) * Sy(s)/volume(c2)
        PHI_z(c2) = PHI_z(c2) - Ps(s) * Sz(s)/volume(c2)
      end if 
    end if  ! c2 < 0
  end do

  call Exchng(PHI_x)
  call Exchng(PHI_y)
  call Exchng(PHI_z)

  END SUBROUTINE GradP3
