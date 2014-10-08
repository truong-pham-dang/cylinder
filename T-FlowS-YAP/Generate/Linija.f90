!======================================================================!
  SUBROUTINE Linija(b, w, is,js,ks,ie,je,ke)
!----------------------------------------------------------------------!
!   Places the nodes on the line defined with local block position     !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: b, is,js,ks,ie,je,ke
  REAL    :: w
!------------------------------[Calling]-------------------------------!
  REAL :: atanh
!-------------------------------[Locals]-------------------------------!
  INTEGER :: N, NI, NJ, NK, i, j, k, node, case
  REAL    :: x0, y0, z0, delx, dely, delz, t, dt, ddt, pr, xi
!--------------------------------[CVS]---------------------------------!
!  $Id: Linija.f90,v 1.7 2002/10/30 16:29:22 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Generate/Linija.f90,v $  
!======================================================================!

  NI=BlkRes(b,1)
  NJ=BlkRes(b,2)
  NK=BlkRes(b,3)   

  x0=x(NN+(ks-1)*NI*NJ+(js-1)*NI+is)
  y0=y(NN+(ks-1)*NI*NJ+(js-1)*NI+is)
  z0=z(NN+(ks-1)*NI*NJ+(js-1)*NI+is)
  delx=x(NN+(ke-1)*NI*NJ+(je-1)*NI+ie)-x0
  dely=y(NN+(ke-1)*NI*NJ+(je-1)*NI+ie)-y0
  delz=z(NN+(ke-1)*NI*NJ+(je-1)*NI+ie)-z0

  N=max( (ie-is), (je-js),  (ke-ks) )

!-----------------------------!
!     Linear distribution     !
!-----------------------------!
  if(w  > 0.0) then
    ddt = ( 2.0*(1.0-w) ) / ( 1.0*N*(1.0*N-1.0)*(1.0+w) )
    t=0.0
    do i=is,ie
      do j=js,je
	do k=ks,ke
	  if( ie /= is ) then
	    dt=1.0/(1.0*N)+(1.0*i-0.5*(1.0*N+1)) * ddt
	    t=t+dt
	    node = NN + (k-1)*NI*NJ + (j-1)*NI + i+1
	    if( (i  < ie).and.(x(node) == HUGE) ) then 
	      x(node) = x0 + t*delx
	      y(node) = y0 + t*dely
	      z(node) = z0 + t*delz
	    endif
	  end if 
	  if( je /= js ) then
	    dt=1.0/(1.0*N)+(1.0*j-0.5*(1.0*N+1)) * ddt
	    t=t+dt
	    node = NN + (k-1)*NI*NJ + (j-0)*NI + i 
	    if( (j  < je).and.(x(node) == HUGE) ) then 
	      x(node) = x0 + t*delx
	      y(node) = y0 + t*dely
	      z(node) = z0 + t*delz
	    endif
	  end if 
	  if( ke /= ks ) then
	    dt=1.0/(1.0*N)+(1.0*k-0.5*(1.0*N+1)) * ddt
	    t=t+dt
	    node = NN + (k-0)*NI*NJ + (j-1)*NI + i 
	    if( (k  < ke).and.(x(node) == HUGE) ) then 
	      x(node) = x0 + t*delx
	      y(node) = y0 + t*dely
	      z(node) = z0 + t*delz
	    endif
	  end if 
	end do
      end do
    end do

!---------------------------------!
!     Hyperbolic distribution     !
!---------------------------------!
  else
    case = 0
    if     ((w  >  -0.5).and.(w <=  -0.25)) then
      pr = 1.0 - abs(0.5 - abs(w))    
      case = 1
    else if((w >=  -0.75).and.(w  <  -0.5)) then
      pr = 1.0 - abs(0.5 - abs(w))            
      case = 2
    else
      pr = -w
      case = 3 
    endif

    do i=is,ie
      do j=js,je
	do k=ks,ke
	  if( ie /= is ) then
    if(case == 1) xi = -1.0*(1.0*i)/(1.0*N)
    if(case == 2) xi =  1.0 - 1.0*(1.0*i)/(1.0*N)
    if(case == 3) xi = -1.0 + 2.0*(1.0*i)/(1.0*N)
	    node = NN + (k-1)*NI*NJ + (j-1)*NI + i+1
	    if( (i  < ie).and.(x(node) == HUGE) ) then 
    if    (case == 1) then
      x(node) = x0 - (tanh(xi*atanh(pr))/pr)*delx
      y(node) = y0 - (tanh(xi*atanh(pr))/pr)*dely
      z(node) = z0 - (tanh(xi*atanh(pr))/pr)*delz
    elseif(case == 2) then
      x(node) = x0 + delx - (tanh(xi*atanh(pr))/pr)*delx
      y(node) = y0 + dely - (tanh(xi*atanh(pr))/pr)*dely
      z(node) = z0 + delz - (tanh(xi*atanh(pr))/pr)*delz
    elseif(case == 3) then
      x(node) = x0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delx
      y(node) = y0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*dely
      z(node) = z0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delz
    endif 
	    endif
	  end if 
	  if( je /= js ) then
    if(case == 1) xi = -1.0*(1.0*j)/(1.0*N)
    if(case == 2) xi =  1.0 - 1.0*(1.0*j)/(1.0*N)
    if(case == 3) xi = -1.0 + 2.0*(1.0*j)/(1.0*N)
	    node = NN + (k-1)*NI*NJ + (j-0)*NI + i 
	    if( (j  < je).and.(x(node) == HUGE) ) then 
    if    (case == 1) then
      x(node) = x0 - (tanh(xi*atanh(pr))/pr)*delx
      y(node) = y0 - (tanh(xi*atanh(pr))/pr)*dely
      z(node) = z0 - (tanh(xi*atanh(pr))/pr)*delz
    elseif(case == 2) then
      x(node) = x0 + delx - (tanh(xi*atanh(pr))/pr)*delx
      y(node) = y0 + dely - (tanh(xi*atanh(pr))/pr)*dely
      z(node) = z0 + delz - (tanh(xi*atanh(pr))/pr)*delz
    elseif(case == 3) then
      x(node) = x0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delx
      y(node) = y0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*dely
      z(node) = z0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delz
    endif 
	    endif
	  end if 
	  if( ke /= ks ) then
    if(case == 1) xi = -1.0*(1.0*k)/(1.0*N)
    if(case == 2) xi =  1.0 - 1.0*(1.0*k)/(1.0*N)
    if(case == 3) xi = -1.0 + 2.0*(1.0*k)/(1.0*N)
	    node = NN + (k-0)*NI*NJ + (j-1)*NI + i 
	    if( (k  < ke).and.(x(node) == HUGE) ) then 
    if    (case == 1) then
      x(node) = x0 - (tanh(xi*atanh(pr))/pr)*delx
      y(node) = y0 - (tanh(xi*atanh(pr))/pr)*dely
      z(node) = z0 - (tanh(xi*atanh(pr))/pr)*delz
    elseif(case == 2) then
      x(node) = x0 + delx - (tanh(xi*atanh(pr))/pr)*delx
      y(node) = y0 + dely - (tanh(xi*atanh(pr))/pr)*dely
      z(node) = z0 + delz - (tanh(xi*atanh(pr))/pr)*delz
    elseif(case == 3) then
      x(node) = x0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delx
      y(node) = y0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*dely
      z(node) = z0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delz
    endif 
	    endif
	  end if 
	end do
      end do
    end do

  endif

  END SUBROUTINE Linija
