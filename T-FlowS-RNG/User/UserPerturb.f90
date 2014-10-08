!======================================================================!
  SUBROUTINE UserPerturb(fac,n)
!----------------------------------------------------------------------!
!   Perturbs the flow field for the channel flow.                      !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: n
  REAL    :: fac
!------------------------------[Calling]-------------------------------!
  REAL    :: Uta, UserURms, UserVRms, UserWRms
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, seed(1)
  REAL    :: randn
!--------------------------------[CVS]---------------------------------!
!  $Id: UserPerturb.f90,v 1.10 2002/10/30 16:30:03 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserPerturb.f90,v $  
!======================================================================!
!   See also: Uplus, URms, VRms, WRms, perturb                         !
!----------------------------------------------------------------------!

!  seed(1) = This*100 + n
!  call random_seed(PUT = seed(1:1))    ! Set user seed

  do c=1,NC

!----------------------------------!
!      add fluctuating values      !
!----------------------------------!
    call random_number(randn)
    Uta = sqrt(   Utau(material(c))**2  &
                + Vtau(material(c))**2  &
                + Wtau(material(c))**2 )

    U % n(c)  = U % n(c)   + 2.0*(0.5-randn) *                  &
	     fac * Uta * UserURms(1.0-abs(zc(c)))
    call random_number(randn)
    U % o(c)  = U % o(c)  + 2.0*(0.5-randn) *                   &
	     fac * Uta * UserURms(1.0-abs(zc(c)))
    call random_number(randn)
    U % oo(c) = U % oo(c) + 2.0*(0.5-randn) *                   &
	     fac * Uta * UserURms(1.0-abs(zc(c)))

    call random_number(randn)
    V % n(c)  = V % n(c)   + 2.0*(0.5-randn) *                  &
	     fac * Uta * UserVRms(1.0-abs(zc(c)))
    call random_number(randn)
    V % o(c)  = V % o(c)  + 2.0*(0.5-randn) *                   &
	     fac * Uta * UserVRms(1.0-abs(zc(c)))
    call random_number(randn)
    V % oo(c) = V % oo(c) + 2.0*(0.5-randn) *                   &
	     fac * Uta * UserVRms(1.0-abs(zc(c)))

    call random_number(randn)
    W % n(c)  = W % n(c)   + 2.0*(0.5-randn) *                  &
	     fac * Uta * UserWRms(1.0-abs(zc(c)))
    call random_number(randn)
    W % o(c)  = W % o(c)  + 2.0*(0.5-randn) *                   &
	     fac * Uta * UserWRms(1.0-abs(zc(c)))
    call random_number(randn)
    W % oo(c) = W % oo(c) + 2.0*(0.5-randn) *                   &
	     fac * Uta * UserWRms(1.0-abs(zc(c)))

!->>> write(*,*) U(c), V(c), W(c)

    write(LinMon0(127:138),'(A6)') ' <<-- '           

  end do

  END SUBROUTINE UserPerturb
