!======================================================================!
  SUBROUTINE UserPotent(PHI,                   &
                        XPHI,                  &
                        dPHIdx,  dPHIdy,   dPHIdz)
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Parameters]---------------------------!
  TYPE(Unknown) :: PHI
  REAL XPHI(NC)
  REAL dPHIdx(-NbC:NC), dPHIdy(-NbC:NC), dPHIdz(-NbC:NC)
!-------------------------------[Locals]-------------------------------!
  INTEGER n,c,s,c1,c2,niter,miter
  REAL    A0,A12,A21,error
  REAL    CONeff, FUex1, FUim1, PHIxS1, PHIyS1, PHIzS1
  REAL            FUex2, FUim2, PHIxS2, PHIyS2, PHIzS2
  REAL    Stot, PHIs
!--------------------------------[CVS]---------------------------------!
!  $Id: UserPotent.f90,v 1.2 2002/10/30 16:30:03 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserPotent.f90,v $  
!======================================================================!

  do n=1,Acol(NC+1) ! to je broj nonzero + 1
    Aval(n) = 0.0
  end do

  do c=1,NC
    b(c)=0.0
  end do

  do c=-NbC,-1
    Abou(c)=0.0
  end do

  call GraPhi(PHI % n, 1, dPHIdx, .TRUE.) ! dp/dx
  call GraPhi(PHI % n, 2, dPHIdy, .TRUE.) ! dp/dy
  call GraPhi(PHI % n, 3, dPHIdz, .TRUE.) ! dp/dz

  do c=1,NC
    U % n(c) = -dPHIdx(c)
    V % n(c) = -dPHIdy(c)
    W % n(c) = -dPHIdz(c)
  end do

!==================!
!                  !
!     Difusion     !
!                  !
!==================!

!----- Set XPHI back to zero 
  do c=1,NC
    XPHI(c) = 0.0  
  end do

!--------------------------------!
!     Spatial Discretization     !
!--------------------------------!

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

!----- gradients on the cell face 
    if(c2  > 0 .or. c2  < 0.and.TypeBC(c2) == BUFFER) then
      PHIxS1 = dPHIdx(c1) 
      PHIyS1 = dPHIdy(c1) 
      PHIzS1 = dPHIdz(c1) 
      PHIxS2 = dPHIdx(c2) 
      PHIyS2 = dPHIdy(c2) 
      PHIzS2 = dPHIdz(c2) 
      CONeff = 1.
    else
      PHIxS1 = dPHIdx(c1) 
      PHIyS1 = dPHIdy(c1) 
      PHIzS1 = dPHIdz(c1) 
      PHIxS2 = PHIxS1 
      PHIyS2 = PHIyS1 
      PHIzS2 = PHIzS1 
      CONeff = 1. 
    endif

!---- total (exact) diffusive flux
    FUex1 = CONeff*(PHIxS1*Sx(s)+PHIyS1*Sy(s)+PHIzS1*Sz(s))
    FUex2 = CONeff*(PHIxS2*Sx(s)+PHIyS2*Sy(s)+PHIzS2*Sz(s))

!---- implicit diffusive flux
    FUim1 = CONeff*Scoef(s)*     &
            (  PHIxS1*Dx(s)      &
             + PHIyS1*Dy(s)      &
             + PHIzS1*Dz(s) )
    FUim2 = CONeff*Scoef(s)*     &
            (  PHIxS2*Dx(s)      &
             + PHIyS2*Dy(s)      &
             + PHIzS2*Dz(s) )

!---- cross diffusion part
    XPHI(c1) = XPHI(c1) + FUex1 - FUim1 
    if(c2.gt.0) then
      XPHI(c2) = XPHI(c2) - FUex2 + FUim2 
    end if 

!----- calculate the coefficients for the sysytem matrix
    A12 = CONeff*Scoef(s)  
    A21 = CONeff*Scoef(s)  

!----- fill the system matrix
    if(c2.gt.0) then
      Aval(Adia(c1)) = Aval(Adia(c1)) + A12
      Aval(Adia(c2)) = Aval(Adia(c2)) + A21
      Aval(SidAij(1,s)) = Aval(SidAij(1,s)) - A12
      Aval(SidAij(2,s)) = Aval(SidAij(2,s)) - A21
    else if(c2.lt.0) then
      if( (TypeBC(c2).eq.INFLOW) .or. (TypeBC(c2).eq.OUTFLOW) ) then        
        Aval(Adia(c1)) = Aval(Adia(c1)) + A12
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !     Warning: this is case-dependent     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        if(TypeBC(c2).eq.INFLOW) then 
          b(c1) = b(c1) + 2.0 * A12
          PHI % n(c2) = 2.0     
        end if
        if(TypeBC(c2).eq.OUTFLOW) then 
          b(c1) = b(c1) + 1.0 * A12
          PHI % n(c2) = 1.0 
        end if
      else if( TypeBC(c2) == BUFFER ) then
        Aval(Adia(c1)) = Aval(Adia(c1)) + A12
        Abou(c2) = - A12  ! cool parallel stuff
      endif
    end if
  end do  ! through sides

  do c=1,NC
    b(c)  = b(c) + XPHI(c)
  end do

!===================================!
!                                   !
!     Solve the equations for PHI   !
!                                   !    
!===================================!
  miter=210 

  niter=miter
  call bicg(NC, Nbc, NONZERO, Aval,Acol,Arow,Adia,Abou,   &
	    PHI % n, b, PREC,                                 &
	    niter, PHI % STol, res(1), error)

  write(*,*) res(1)
  write(*,*) niter       
 
!--------------------------------------!
!    Handle non inflow and outflow     !
!--------------------------------------!
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2 < 0) then
      if( (TypeBC(c2).ne.INFLOW)  .and. &
          (TypeBC(c2).ne.OUTFLOW) .and. & 
          (TypeBC(c2).ne.BUFFER)  ) then 
        PHI % n(c2) = PHI % n(c1)
      end if
    end if ! c2 < 0 
  end do

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2 < 0) then
      if( (TypeBC(c2).eq.OUTFLOW) ) then   
        U % n(c2) = U % n(c1)
        V % n(c2) = V % n(c1)
        W % n(c2) = W % n(c1)
      end if
    end if ! c2 < 0 
  end do
  do c=1 ,NC
    U % o(c)  = U % n(c)
    U % oo(c) = U % n(c)
    V % o(c)  = V % n(c)
    V % oo(c) = V % n(c)
    W % o(c)  = W % n(c)
    W % oo(c) = W % n(c)
  end do

  call Exchng(PHI % n)

  END SUBROUTINE
