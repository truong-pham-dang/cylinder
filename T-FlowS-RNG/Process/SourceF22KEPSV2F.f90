!=====================================================================!
  SUBROUTINE SourceF22KEPSV2F()
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
 if(SIMULA == ZETA.or.SIMULA==HYB_ZETA) then 
   do c = 1,NC
     f22hg = (1.0 - Cf_1 - 0.65*Pk(c)/(Eps % n(c) + tiny) )*(v_2 % n(c) - 2.0/3.0)/(Tsc(c) + tiny)   &
             + 0.0085 * Pk(c)/(Kin % n(c) + tiny)
     b(c) = b(c) + f22hg*volume(c)/(Lsc(c)**2.0 + tiny) 
   end do
!---- izvorni clan

   do c = 1, NC
     Sor11 = volume(c)/(Lsc(c)**2 + tiny)
     Aval(Adia(c)) = Aval(Adia(c)) + Sor11     
   end do

!  Imposing boundary condition for f22 on the wall
   do s=1,NS
     c1=SideC(1,s)
     c2=SideC(2,s)
     if(c2 < 0 .and. TypeBC(c2) /= BUFFER ) then
       if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then


         f22 % n(c2) = -2.0*VISc*v_2 % n(c1)/                     &
                        WallDs(c1)**2

!----   Fill in a source coefficients
!---------------------------------------------------------------------  
!----   Linearization of the near wall terms helps to get more  
!----   stable solution, especially for small wall distance.
!---------------------------------------------------------------------
         A0 = Scoef(s) 
         b(c1) = b(c1) + A0*f22%n(c2)
       end if   ! end if of BC=wall
     end if    ! end if of c2<0
   end do
 else if(SIMULA == K_EPS_VV) then
   do c = 1,NC
     f22hg = (1.0 - Cf_1)*(v_2 % n(c)/(Kin % n(c)+tiny) - 2.0/3.0)/(Tsc(c)+tiny)   &
                + Cf_2*Pk(c)/(Kin % n(c)+tiny)
     b(c) = b(c) + f22hg*volume(c)/(Lsc(c)+tiny)**2.0
     Sor11 = volume(c)/Lsc(c)**2.0
     Aval(Adia(c)) = Aval(Adia(c)) + Sor11
   end do

!  Imposing boundary condition for f22 on the wall
   do s=1,NS
     c1=SideC(1,s)
     c2=SideC(2,s)
     if(c2 < 0 .and. TypeBC(c2) /= BUFFER ) then
       if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
         f22 % n(c2) = -20.0*VISc**2*v_2 % n(c1)/                     &
                       (WallDs(c1)**4*Eps % n(c2))
!----   Fill in a source coefficients
!---------------------------------------------------------------------  
!----   Linearization of the near wall terms helps to get more  
!----   stable solution, especially for small wall distance.
!---------------------------------------------------------------------
         A0 = Scoef(s)
         b(c1) = b(c1) + A0*f22%n(c2)
       end if   ! end if of BC=wall
     end if    ! end if of c2<0
   end do
 end if 
 RETURN
 END SUBROUTINE  SourceF22KEPSV2F
