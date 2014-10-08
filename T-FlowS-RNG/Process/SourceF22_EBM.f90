!=====================================================================!
  SUBROUTINE SourceF22_EBM()
!---------------------------------------------------------------------!
!  Calculate source terms in eliptic relaxation  equation             !
!  for vi2 and imposing  boundary condition for f22                   !
!---------------------------------------------------------------------!
!------------------------------[Moules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER s, c, c1, c2, j
  REAL    Sor11,  f22hg
!=============================================
  REAL    A0
!--------------------------------[CVS]---------------------------------!
!  $Id: SourceF22KEPSV2F.f90,v 1.4 2008/11/19 15:00:41 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/SourceF22KEPSV2F.f90,v $    
!----------------------------------------------------------------------!
!                                                                      !
!  The form of source terms are :                                      ! 
!                                                                      !
!     /                                                                !
!    |                                                                 !  
!    | f22hg*dV ; f22hg - f22hg homogenious is placed in a source      !
!    |                     coefficients b(c)                           !
!   /                                                                  !
!      f22hg = (1.0 - Cv_1)*(vi2(c)/Kin(c) - 2.0/3.0)/Tsc(c)     &     !
!              + 2.0*Cv_2*Pk(c)/(3.0*Kin(c))                           !        
!                                                                      ! 
!     /                                                                ! 
!    |                                                                 !
!    | f22*dV ;this term is placed in a diagonal of coefficient matrice!  
!    |                                                                 !
!   /                                                                  !
!                                                                      !
!                                                                      ! 
!  Dimensions of certain variables                                     !
!                                                                      !
!     Tsc            [s]                                               !     
!     Kin            [m^2/s^2]                                         !
!     Eps            [m^3/s^2]                                         !
!     vi2            [m^2/s^2]                                         !
!     f22            [-]                                               !
!     Lsc            [m]                                               !
!======================================================================!
   call Scale()

!---- izvorni clan f22hg
   do c = 1,NC
     f22hg = 1.0
     Sor11 = volume(c)/Lsc(c)**2
     Aval(Adia(c)) = Aval(Adia(c)) + Sor11     
     b(c) = b(c) + f22hg*volume(c)/Lsc(c)**2
   end do
!---- izvorni clan
   do s=1,NS
     c1=SideC(1,s)
     c2=SideC(2,s)
     if(c2 < 0 .and. TypeBC(c2) /= BUFFER ) then
       if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then

           f22 % n(c2) = 0.0

       end if   ! end if of BC=wall
     end if    ! end if of c2<0
   end do

 RETURN
 END SUBROUTINE  SourceF22_EBM
