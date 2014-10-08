!======================================================================!
  SUBROUTINE GradP(PHI, PHI_x, PHI_y, PHI_z)
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
  INTEGER :: s, c, c1, c2, iter
!--------------------------------[CVS]---------------------------------!
!  $Id: GradP.f90,v 1.15 2008/12/10 14:40:03 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/GradP.f90,v $      
!======================================================================!
						 
  call Exchng(PHI)

  Ps = 0.0

  do c=1,NC
    PHI_x(c)=0.0
    PHI_y(c)=0.0
    PHI_z(c)=0.0
  end do

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
      if(TypeBC(c2) /= PRESSURE) then  !IZMJENA 1
        PHI(c2) = PHI(c1)
      end if
    end if  
  end do

  call GraPhi(PHI,1,PHI_x,.TRUE.) ! dP/dx
  call GraPhi(PHI,2,PHI_y,.TRUE.) ! dP/dy
  call GraPhi(PHI,3,PHI_z,.TRUE.) ! dP/dz    

  do iter=1,1 

    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
        if(TypeBC(c2) /= PRESSURE) then                  ! IZMJENA
          PHI(c2) = PHI(c1) + 1.2 * ( PHI_x(c1) * (xc(c2)-xc(c1)) &  
                            +         PHI_y(c1) * (yc(c2)-yc(c1)) &  
                            +         PHI_z(c1) * (zc(c2)-zc(c1)) )   
        end if  
      end if  
    end do

    call GraPhi(PHI,1,PHI_x,.TRUE.) ! dP/dx
    call GraPhi(PHI,2,PHI_y,.TRUE.) ! dP/dy
    call GraPhi(PHI,3,PHI_z,.TRUE.) ! dP/dz 

  end do

  END SUBROUTINE GradP
