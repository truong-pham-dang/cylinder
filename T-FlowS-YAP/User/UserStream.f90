!======================================================================!
  SUBROUTINE UserStream(PHI,                   &
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
  INTEGER :: n,c,s,c1,c2,niter,miter
  REAL    :: A12,A21,error
  REAL    :: CONeff, FUex1, FUim1, PHIxS1, PHIyS1, PHIzS1
  REAL    ::         FUex2, FUim2, PHIxS2, PHIyS2, PHIzS2
  REAL    :: H, half
!--------------------------------[CVS]---------------------------------!
!  $Id: UserStream.f90,v 1.5 2002/10/30 16:30:03 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserStream.f90,v $  
!======================================================================!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     Warning: this is case-dependent     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!  half = 0.15 ! splits the domain in lower and upper half
!  H    = 0.19 ! characteristic height of the domain

  half = 3.0 ! splits the domain in lower and upper half
  H    = 2.0 ! characteristic height of the domain

  write(*,'(A)') '==> In UserStream'
  write(*,'(A,F8.3)') 'H   =', H
  write(*,'(A,F8.3)') 'half=', half

  do ini = 1,10 

    call GraPhi(PHI % n, 1, Px, .TRUE.) ! dp/dx
    call GraPhi(PHI % n, 2, Py, .TRUE.) ! dp/dy
    call GraPhi(PHI % n, 3, Pz, .TRUE.) ! dp/dz

    do n=1,Acol(NC+1) ! to je broj nonzero + 1
      Aval(n) = 0.0
    end do

    do c=1,NC
      b(c)=0.0
    end do

    do c=-NbC,-1
      Abou(c)=0.0
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
        if( (TypeBC(c2).eq.WALL)   .or. & 
            (TypeBC(c2).eq.WALLFL) .or. & 
            (TypeBC(c2).eq.SYMMETRY).or. &         
            (TypeBC(c2).eq.CONVECT) ) then        
          Aval(Adia(c1)) = Aval(Adia(c1)) + A12
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !     Warning: this is case-dependent     !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
          if(zc(c2) < half) then 
!          if(zc(c2) > half) then 
!            b(c1) = b(c1) + Wbulk(1) * H * A12
!            PHI % n(c2) = Wbulk(1) * H
!            write(*,*) Wbulk(1)      
            b(c1) = b(c1) + 1.0 * H * A12
            PHI % n(c2) = 1.0 * H
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
  
    call GraPhi(U % n,3,PHIz, .TRUE.) ! dU/dz
    call GraPhi(W % n,1,PHIx, .TRUE.) ! dW/dx 

    do c=1,NC
      b(c)  = b(c) + volume(c)*(PHIx(c) - PHIz(c)) 
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

    write(*,*) res(5)
    write(*,*) niter       
 
!----------------------------------!
!    Handle inflow and outflow     !
!----------------------------------!
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      if(c2 < 0) then
        if( (TypeBC(c2).eq.INFLOW) .or. (TypeBC(c2).eq.OUTFLOW) ) then 
          PHI % n(c2) = PHI % n(c1)
        end if
      end if ! c2 < 0 
    end do

    call Exchng(PHI % n)

  end do ! ini

  END SUBROUTINE
