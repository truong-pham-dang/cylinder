!======================================================================!
  SUBROUTINE GraPhiCor(PHI, PHI_x, PHI_y, PHI_z)
!----------------------------------------------------------------------!
! Calculates gradient in the cells adjencent to material interface     !
! boundaries. Assumes that "tentative" gradients are just calculated   !
! and stored in "PHI_x", "PHI_y" and "PHI_z" arrays.                   !
!                                                                      !
! It also assumes that the gradients "PHI_x", "PHI_y" and "PHI_z"      !
! are fresh in buffers.                                                !
!                                                                      !
! This entire procedure is for two materials.                          !
!                                                                      !
! It is not desiged, and probably won't work in the following          !
! situations:                                                          !
!                                                                      !
!    +---+---+---+                                                     !
!    | F | F | F |                                                     !
!    +---+---+---+                                                     !
!    | F | S | S |  ->  The cell in the middle will probably not work  !
!    +---+---+---+                                                     !
!    | F | S | S |                                                     !
!    +---+---+---+                                                     !
!                                                                      !
! Further, it will probably not work on periodic boundaries.           !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE sol_mod  ! needed for p1 and p2 arrays
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  REAL :: PHI(-NbC:NC),    &
          PHI_x(-NbC:NC),  &
          PHI_y(-NbC:NC),  &
          PHI_z(-NbC:NC)
!-------------------------------[Locals]-------------------------------!
  INTEGER :: s, c, c1, c2
  REAL    :: DPHI1, DPHI2, Dxc1, Dyc1, Dzc1, Dxc2, Dyc2, Dzc2 
  REAL    :: f1, f2, PHIs
!--------------------------------[CVS]---------------------------------!
!  $Id: GraPhiCor.f90,v 1.9 2008/12/10 14:39:46 IUS\mhadziabdic Exp $
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/GraPhiCor.f90,v $ 
!======================================================================!

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

!---- Take care of material interfaces          
    if( StateMat(material(c1))==FLUID .and. &  
        StateMat(material(c2))==SOLID       &  
        .or.                                &  
        StateMat(material(c1))==SOLID .and. &  
        StateMat(material(c2))==FLUID ) then   

      Dxc1 = xsp(s)-xc(c1)                     
      Dyc1 = ysp(s)-yc(c1)                     
      Dzc1 = zsp(s)-zc(c1)                     
      Dxc2 = xsp(s)-xc(c2)                     
      Dyc2 = ysp(s)-yc(c2)                     
      Dzc2 = zsp(s)-zc(c2)                     

!---- Missing parts of the gradient vector
      p1(c1) = CONc(material(c1)) *  &
           ( (G(1,c1)*Dxc1+G(4,c1)*Dyc1+G(5,c1)*Dzc1) * Sx(s) + &
             (G(4,c1)*Dxc1+G(2,c1)*Dyc1+G(6,c1)*Dzc1) * Sy(s) + & 
             (G(5,c1)*Dxc1+G(6,c1)*Dyc1+G(3,c1)*Dzc1) * Sz(s) )
      if(c2 > 0) then               
        p2(c2) = CONc(material(c2)) *  &
              ( (G(1,c2)*Dxc2+G(4,c2)*Dyc2+G(5,c2)*Dzc2) * Sx(s) + &
                (G(4,c2)*Dxc2+G(2,c2)*Dyc2+G(6,c2)*Dzc2) * Sy(s) + &
                (G(5,c2)*Dxc2+G(6,c2)*Dyc2+G(3,c2)*Dzc2) * Sz(s) )
      else if(TypeBC(c2) == BUFFER) then ! prepare to exchange
        p2(c1) = -p1(c1)
      end if
    end if    
  end do
 
  call Exchng(p2)

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

!---- Take care of material interfaces          
    if( StateMat(material(c1))==FLUID .and. &  
        StateMat(material(c2))==SOLID       &  
        .or.                                &  
        StateMat(material(c1))==SOLID .and. &  
        StateMat(material(c2))==FLUID ) then   

!---- Flux from cell 1 towards the material interface
      f1 = CONc(material(c1)) *  &
           (PHI_x(c1)*Sx(s) + PHI_y(c1)*Sy(s) + PHI_z(c1)*Sz(s))

!---- Flux from cell 2 towards the material interface
      f2 = CONc(material(c2)) *  &
           (PHI_x(c2)*Sx(s) + PHI_y(c2)*Sy(s) + PHI_z(c2)*Sz(s))

!---- The two fluxes (q1 and q2) should be the same
      PHIs = (f2 - f1) / (p1(c1) - p2(c2) + TINY)
      PHIside(s) = PHIs

      Dxc1 = xsp(s)-xc(c1)                     
      Dyc1 = ysp(s)-yc(c1)                     
      Dzc1 = zsp(s)-zc(c1)                     
      Dxc2 = xsp(s)-xc(c2)                     
      Dyc2 = ysp(s)-yc(c2)                     
      Dzc2 = zsp(s)-zc(c2)                     

!---- Now update the gradients
      PHI_x(c1)=PHI_x(c1)+PHIs*(G(1,c1)*Dxc1+G(4,c1)*Dyc1+G(5,c1)*Dzc1)
      PHI_y(c1)=PHI_y(c1)+PHIs*(G(4,c1)*Dxc1+G(2,c1)*Dyc1+G(6,c1)*Dzc1) 
      PHI_z(c1)=PHI_z(c1)+PHIs*(G(5,c1)*Dxc1+G(6,c1)*Dyc1+G(3,c1)*Dzc1)

      if(c2 > 0) then
       PHI_x(c2)=PHI_x(c2)+PHIs*(G(1,c2)*Dxc2+G(4,c2)*Dyc2+G(5,c2)*Dzc2)
       PHI_y(c2)=PHI_y(c2)+PHIs*(G(4,c2)*Dxc2+G(2,c2)*Dyc2+G(6,c2)*Dzc2)
       PHI_z(c2)=PHI_z(c2)+PHIs*(G(5,c2)*Dxc2+G(6,c2)*Dyc2+G(3,c2)*Dzc2)
      end if

    end if    
  end do

  call Exchng(PHI_x)
  call Exchng(PHI_y)
  call Exchng(PHI_z)

  END SUBROUTINE GraPhiCor
