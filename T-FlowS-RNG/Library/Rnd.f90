!======================================================================!
  REAL FUNCTION Rnd (r)
!----------------------------------------------------------------------!
!   Computes random double precision number in the range 0.0-1.0.      !
!----------------------------------------------------------------------!
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: r
!--------------------------------[CVS]---------------------------------!
!  $Id: Rnd.f90,v 1.4 2002/10/30 16:29:33 niceno Exp $
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Library/Rnd.f90,v $
!======================================================================!
!   Input parameter r must be 0, otherwise it was not working.         !
!   One of authors of TFlowS (B. Niceno) found it somewhere on         !
!   the web and never really understood how and why it works.          !
!----------------------------------------------------------------------!
  data ia1, ia0, ia1ma0 /1536, 1029, 507/
  data ic /1731/
  data ix1, ix0 /0, 0/
!----------------------------------------------------------------------!
  if (r  < 0.) go to 1
  if (r  > 0.) go to 2

  iy0 = ia0*ix0
  iy1 = ia1*ix1 + ia1ma0*(ix0-ix1) + iy0
  iy0 = iy0 + ic
  ix0 = mod (iy0, 2048)
  iy1 = iy1 + (iy0-ix0)/2048
  ix1 = mod (iy1, 2048)

1 rnd = ix1*2048 + ix0
  rnd = rnd / 4194304.
  return

2 ix1 = mod(r,1)*4194304. + 0.5
  ix0 = mod (ix1, 2048)
  ix1 = (ix1-ix0)/2048
  go to 1

  END FUNCTION Rnd
