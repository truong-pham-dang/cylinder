!======================================================================!
  SUBROUTINE ConvScheme(Phis, s,                            &
                        Phi,                                &
                        dPhidi, dPhidj, dPhidk, Di, Dj, Dk, &
                        blenda) 
!----------------------------------------------------------------------!
! Computes the value at the cell face using different convective       !
! schemes. In this subroutine I try to follow the nomenclature from    !
! Basara's and Przulj's AIAA paper.                                    !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  REAL          :: Phis
  INTEGER       :: s
  REAL          :: Phi(-NbC:NC)
  REAL          :: dPhidi(-NbC:NC), dPhidj(-NbC:NC), dPhidk(-NbC:NC)
  REAL          :: Di(NS),          Dj(NS),          Dk(NS)
  INTEGER       :: blenda
!-------------------------------[Locals]-------------------------------!
  INTEGER       :: c1, c2, C, D
  REAL          :: fj ! flow oriented interpolation factor
  REAL          :: gD, gU, alfa, beta1, beta2 
  REAL          :: Phij, PhiU, PhiUstar, rj, sign, GammaC, Beta
!--------------------------------[CVS]---------------------------------!
!  $Id: ConvScheme.f90,v 1.6 2009/06/30 12:06:21 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/ConvScheme.f90,v $  
!======================================================================!
!
!               Flux > 0
!   +---------+---------+---------+---------+
!   |         |         |         |         |
!   |         |   c1    |   c2    |         |
!   |    o    |    o  ==s=>  o    |    o    |----> xi
!   |    U    |    C    |    D    |         |
!   |         |         |         |         |
!   +---------+---------+---------+---------+   
!
!
!               Flux < 0
!   +---------+---------+---------+---------+
!   |         |         |         |         |
!   |         |   c1    |   c2    |         |
!   |    o    |    o  <=s==  o    |    o    |----> xi
!   |         |    D    |    C    |    U    |
!   |         |         |         |         |
!   +---------+---------+---------+---------+   
!
!----------------------------------------------------------------------!

  c1 = SideC(1,s)
  c2 = SideC(2,s)

  if(Flux(s) > 0.0) then ! goes from c1 to c2
    fj   = 1.0 - f(s)
    C    = c1
    D    = c2
    sign = +1.0
  else ! Flux(s) < 0.0   ! goes from c2 to c1
    fj = f(s)
    C    = c2
    D    = c1
    sign = -1.0
  end if

  if(Flux(s) > 0.0) then
    PhiUstar = Phi(D) - 2.0 * ( dPhidi(C)*Di(s) &
                               +dPhidj(C)*Dj(s) &
                               +dPhidk(C)*Dk(s) )
  else
    PhiUstar = Phi(D) + 2.0 * ( dPhidi(C)*Di(s) &
                               +dPhidj(C)*Dj(s) &
                               +dPhidk(C)*Dk(s) )
  end if

  PhiU = max( PhiMin(C), min(PhiUstar,PhiMax(C)) )

  rj = ( Phi(C) - PhiU ) / ( Phi(D)-Phi(C) + 1.0e-16 )

  gD = 0.5 * fj * (1.0+fj)
  gU = 0.5 * fj * (1.0-fj)

  if(blenda == CDS) then
    Phij = fj
  else if(blenda == QUICK) then
    rj = ( Phi(C) - PhiU ) / ( Phi(D)-Phi(C) + 1.0e-12 )
    alfa = 0.0
    Phij = (gD - alfa) + (gU + alfa) * rj
  else if(blenda == LUDS) then
    alfa = 0.5 * fj * (1+fj)
    Phij = (gD - alfa) + (gU + alfa) * rj
  else if(blenda == MINMOD) then
    Phij = fj * max(0.0, min(rj,1.0))
  else if(blenda == SMART) then
    beta1 = 3.0
    beta2 = 1.0
    Phij = max( 0.0, min( (beta1-1.0)*rj, gD+gU*rj, beta2 ) )
  else if(blenda == AVL_SMART) then
    beta1 = 1.0 + fj*(2.0+fj) 
    beta2 = fj*(2.0-fj) 
    Phij = max( 0.0, min( (beta1-1.0)*rj, gD+gU*rj, beta2 ) )
  else if(blenda == SUPERBEE) then
    Phij = 0.5 * max( 0.0, min( 2.0*rj,1.0 ), min( rj,2.0 ) )
  else if(blenda == YES) then
    return
  end if

  PhiS = Phi(C) + Phij * sign * (Phi(c2)-Phi(c1))

  if(blenda == GAMMA) then
    Beta = 0.1

    if(Flux(s) > 0.0) then
      PhiUstar = 1.0 - (Phi(D) - Phi(C))/(2.0 * ( dPhidi(C)*Di(s) &
                                 +dPhidj(C)*Dj(s) &
                                 +dPhidk(C)*Dk(s)))
    else
      PhiUstar = 1.0 + (Phi(D) - Phi(C))/(2.0 * ( dPhidi(C)*Di(s) &
                                 +dPhidj(C)*Dj(s) &
                                 +dPhidk(C)*Dk(s)))
    end if

    GammaC = PhiUstar/Beta

    if(PhiUstar < Beta.and.PhiUstar > 0.0) then
      PhiS = (1.0 - GammaC*(1.0 - f(s)))*Phi(C) + GammaC*(1.0 - f(s))*Phi(D)
    else if(PhiUstar < 1.0.and.PhiUstar >= Beta) then
       PhiS = f(s)*Phi(C) + (1.0 - f(s))*Phi(D)
    else
      PhiS = Phi(C)
    end if
  end if 

  RETURN 

  END SUBROUTINE ConvScheme
