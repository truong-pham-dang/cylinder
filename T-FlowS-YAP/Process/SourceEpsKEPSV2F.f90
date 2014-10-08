!======================================================================!
  SUBROUTINE SourceEpsKEPSV2F()
!-----------------------------------------------------------------------!
! Purpose:                                                              !
! Calculates source terms in equation of dissipation of turbulent energy! 
!  and imposes  a boundary condition                                    ! 
! model : k-eps-v2f                                                     !  
! Authors: Muhamed Hadziabdic and Bojan Niceno                          !
!-----------------------------------------------------------------------!
!------------------------------[Modules]--------------------------------!
  USE all_mod
  USE pro_mod
  USE rans_mod
  USE les_mod
  !use ieee_arithmetic
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, s, c1, c2,j	 
  REAL    :: Esor, Ce_11, Gblend, fp, fa, Rey, Ret, Ck
  REAL    :: Utan, UnorSq, Unor, UtotSq, dely, Stot ,F_YAP
!--------------------------------[CVS]---------------------------------!
!  $Id: SourceEpsKEPSV2F.f90,v 1.4 2008/11/19 15:03:34 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/SourceEpsKEPSV2F.f90,v $  
!======================================================================!


!----------------------------------------------------------------------!
!    In dissipation of turbulent kinetic energy equation exist two     !
!    source terms which have form:                                     !
!                                                                      !
!    /                                                                 !
!   |                                                                  !
!   |((Cv_e1*PkE - Cv_11 DENc*Eps)/Tsc)*dV                             !
!   |                                                                  !
!  /                                                                   !
!                                                                      !
!  First, positive , source term is solved and added to source         !
!  coefficient b(c) on right hand side.                                !
!  Second, negative, source term is added to main diagonal left hand   !
!  side coefficient matrix in order to increase stability of solver    !
!  It is nessesary to calculate coefficient Cv_11 using Kin, Cv_e2, vi2!
!  and coefficient A1                                                  !
!----------------------------------------------------------------------!      

  call Scale()

  if(SIMULA == ZETA.or.SIMULA==HYB_ZETA) then
    !call GraPhi((Kin % n**1.5) / (Eps % n + tiny),1,grad_L_x, .TRUE.) ! grad(L)_x
    !call GraPhi((Kin % n**1.5) / (Eps % n + tiny),2,grad_L_y, .TRUE.) ! grad(L)_y
    !call GraPhi((Kin % n**1.5) / (Eps % n + tiny),3,grad_L_z, .TRUE.) ! grad(L)_z

call wait
    Eps % mean(-1) = 0.0

    do c = 1,NC 

      F_YAP = (abs(Kin % n(c)) ** 1.5)/(Eps % n(c) + tiny)/(2.44*abs(WallDs(c)))
      ! to make it out of here     
      Eps % mean(-1) = max(F_YAP,Eps % mean(-1))
 
      Ce_11 = Ce1*(1.0 + alpha*(1.0/(v_2%n(c)+tiny) ))

      !F_YAP=(((grad_L_x(c)**2+grad_L_y(c)**2+grad_L_z(c)**2)**0.5)/2.44)      ! F_YAP. ref. 155 pg.M.T.in.E.a.E. Hanjalic
      Esor = volume(c)/(Tsc(c)+tiny)

      b(c) = b(c) + Ce_11*Pk(c)*Esor +&
      max(0.00,0.83*(Eps % n(c)**2) / abs(Kin % n(c)+ tiny) * (F_YAP - 1)*(F_YAP**2))*volume(c)
 
!----- Fill in a diagonal of coefficient matrix 
      Aval(Adia(c)) =  Aval(Adia(c)) + Ce2*Esor*DENc(material(c))
    end do                   
  end if
!----  Imposing a boundary condition on wall for Eps 

    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER ) then
        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
          UtotSq = U % n(c1) * U % n(c1) &
                 + V % n(c1) * V % n(c1) &
                 + W % n(c1) * W % n(c1) 
          Unor = ( U % n(c1) * Sx(s)     &
                 + V % n(c1) * Sy(s)     &
                 + W % n(c1) * Sz(s) )   &
                 /(Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s))**0.5  
          UnorSq = Unor*Unor

          if( UtotSq   >  UnorSq) then
            Utan = sqrt(UtotSq - UnorSq)
          else
            Utan = TINY
          end if

          Uf(c1)  = Cmu**0.25*Kin%n(c1)**0.5
          Ynd(c1) = WallDs(c1)*Uf(c1)/VISc
          TauWall(c1) =  abs(abs(DENc(material(c1)) * kappa * Uf(c1) * Utan    &
                    / (log(Elog*Ynd(c1)))))
      
          Ck = sqrt(TauWall(c1))
          Ynd(c1) = WallDs(c1)*Ck/VISc

          fa = min(Ck**3.0/(kappa*WallDs(c1)*Pk(c1)),1.0)
          if(WallDs(c1)*Ck/VISc>5.0) then
            Eps%n(c1) = (1.0-fa)*2.0*VISc*Kin%n(c1)/WallDs(c1)**2 + &
            fa*0.07**0.75*Kin%n(c1)**1.5/(kappa*WallDs(c1))
!----- Adjusting a coefficient in order to get a fixed
!----- value of Eps in near wall calls
            do j=Acol(c1), Acol(c1+1) -1
              Aval(j) = 0.0
            end do
            Aval(Adia(c1)) = 1.0 
            b(c1)       = Eps % n(c1)
          else
            Eps%n(c2) = 2.0*VISc*Kin%n(c1)/WallDs(c1)**2         
          end if  
       end if
      end if
    end do  

  if(SIMULA == K_EPS_VV) then
    do c = 1,NC
      Esor = volume(c)/Tsc(c)
      Ce_11 = Ce1*(1.0 + alpha*(Kin%n(c)/(v_2%n(c)) + tiny)**0.5)
      b(c) = b(c) + Ce_11*Pk(c)*Esor

!----- Fill in a diagonal of coefficient matrix
      Aval(Adia(c)) =  Aval(Adia(c)) + Ce2*Esor*DENc(material(c))
    end do
  end if 

  RETURN
  END SUBROUTINE SourceEpsKEPSV2F 
