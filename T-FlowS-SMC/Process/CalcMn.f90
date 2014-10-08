!======================================================================!
  SUBROUTINE CalcMn(n0, n1,restart)   
!----------------------------------------------------------------------!
!   Calculates time averaged velocity and velocity fluctuations.       !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE les_mod
  USE pro_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: n0, n1
  LOGICAL :: restart
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, n
!--------------------------------[CVS]---------------------------------!
!  $Id: CalcMn.f90,v 1.17 2008/12/05 13:22:06 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CalcMn.f90,v $   
!======================================================================!
 

  n=n1-n0

  if(n  > -1) then

  do c=1,NC
!-----------------------!
!      mean values      !
!-----------------------!

    U   % mean(c) = ( U  % mean(c) * (1.*n) + U % n(c) ) / ( 1. *(n+1) )
    V   % mean(c) = ( V  % mean(c) * (1.*n) + V % n(c) ) / ( 1. *(n+1) )
    W   % mean(c) = ( W  % mean(c) * (1.*n) + W % n(c) ) / ( 1. *(n+1) )
    P   % mean(c) = ( P  % mean(c) * (1.*n) + P  % n(c) ) / ( 1. *(n+1) )

    Ls  % mean(c) = ( Ls % mean(c) * (1.*n) + ( ( Eps % n(c) / ( Shear(c) / ( 2**0.5 ) ) ** 3 )**0.5 ) ) / ( 1. *(n+1) )
    Eps % mean(c) = ( Eps% mean(c) * (1.*n) + Eps % n(c) ) / ( 1. *(n+1) )
    Kin % mean(c) = ( Kin% mean(c) * (1.*n) + Kin % n(c) ) / ( 1. *(n+1) )
    VISt_mean(c) = ( VISt_mean(c) * (1.*n) + VISt(c) ) / ( 1. *(n+1) )

    uu  % mean(c) = ( uu % mean(c) * (1.*n) + uu % n(c) ) / ( 1. *(n+1) )
    vv  % mean(c) = ( vv % mean(c) * (1.*n) + vv % n(c) ) / ( 1. *(n+1) )
    ww  % mean(c) = ( ww % mean(c) * (1.*n) + ww % n(c) ) / ( 1. *(n+1) )
    uv  % mean(c) = ( uv % mean(c) * (1.*n) + uv % n(c) ) / ( 1. *(n+1) )
    uw  % mean(c) = ( uw % mean(c) * (1.*n) + uw % n(c) ) / ( 1. *(n+1) )
    vw  % mean(c) = ( vw % mean(c) * (1.*n) + vw % n(c) ) / ( 1. *(n+1) )

  end do
  end if ! n > -1
!----------------------------------------------------------------------------------------------------------------------!
!!    Tracking local minimum of Fpy to determine the phase of Uy
!----------------------------------------------------------------------------------------------------------------------!

call Wait

END SUBROUTINE CalcMn
