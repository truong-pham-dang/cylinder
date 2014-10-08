!======================================================================!
  SUBROUTINE UserPerturb2(fac,n, Dom)
!----------------------------------------------------------------------!
!   Perturbs the flow field for any flow.                              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE par_mod
  USE les_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------[Calling]-------------------------------!
  INTEGER :: n, Dom
  REAL    :: fac
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, seed(1), q
  REAL    :: randn
  integer, allocatable :: a(:) 
!--------------------------------[CVS]---------------------------------!
!  $Id: UserPerturb2.f90,v 1.11 2002/10/30 16:30:03 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserPerturb2.f90,v $  
!======================================================================!

  call random_seed(size=q)
  allocate(a(q)) 
  a = This*100 + n
  call random_seed(PUT = a)    ! Set user seed

!----------------------------------!
!      add fluctuating values      !
!----------------------------------!
  do c=1,NC
    if(material(c) == Dom) then
    call random_number(randn)

!---- 10 % of the maximum values
    U % n(c)  = U % n(c) + .1*fac*(0.5-randn) & 
              * max(abs(U % n(c)), abs(V % n(c)), abs(W % n(c))) 
    U % o(c)  = U % n(c)
    U % oo(c) = U % n(c)

!---- 10 % of the maximum values
    V % n(c)  = V % n(c)   + .1*fac*(0.5-randn) &
              * max(abs(U % n(c)), abs(V % n(c)), abs(W % n(c))) 
    V % o(c)  = V % n(c)
    V % oo(c) = V % n(c) 

!---- 10 % of the maximum values
    W % n(c)  = W % n(c)   + .1*fac*(0.5-randn) &
              * max(abs(U % n(c)), abs(V % n(c)), abs(W % n(c))) 
    W % o(c)  = W % n(c)
    W % oo(c) = W % n(c)
    end if
  end do


  END SUBROUTINE UserPerturb2
