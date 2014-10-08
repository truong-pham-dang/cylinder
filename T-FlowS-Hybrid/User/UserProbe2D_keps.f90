!======================================================================!
  SUBROUTINE UserProbe2D(namAut) 
!----------------------------------------------------------------------!
! Reads the ".1D" file created by the "Generator" and  creates cut     !
! lines which position is defined in *.cmn file.It works for all models!
! implemented by now. For LES it creates additional files which        !
! contains data of average variables in cut lines as well as           !
! instantenious.                                                       !
! It is necessary to defined exact values od x,y coordinates in*.cmn   !
!----------------------------------------------------------------------!
  USE all_mod
  USE allp_mod
  USE les_mod
  USE pro_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------[Calling]-------------------------------!
  INTERFACE
    LOGICAL FUNCTION Approx(A,B,tol)
      REAL           :: A,B
      REAL, OPTIONAL :: tol
    END FUNCTION Approx
  END INTERFACE 
!-------------------------------[Locals]-------------------------------!
  INTEGER             :: pl, c, dummy, K, Npoints
  CHARACTER           :: namCoo*80, answer*80
  REAL                :: u_tau, Re_tau, tau_w,u_tau_i,tau_w_i, Uref, WD
  REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:),         & 
                                 uup(:), vvp(:), wwp(:),         &
                                 uvp(:), uwp(:), vwp(:),         &
                                 Umean(:), Vmean(:), Wmean(:),   &
                                 Kin_mp(:), Eps_mp(:), v_2_mp(:),&
                                 f22_mp(:), Utau_a(:), Utau_i(:),&
                                 VISt_mp(:),VIS_mp(:), T_mp(:),  &
                                 Pk_mp(:), PHI(:)
  INTEGER,ALLOCATABLE :: Np(:), order(:) 
  CHARACTER, OPTIONAL :: namAut*(*)            
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe2D.f90,v 1.6 2003/01/16 13:35:26 muhamed Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe2D.f90,v $  
!======================================================================!
  if(PRESENT(namAut)) then
    write(*,*) namAut
!---- save the name
    namCut = namAut
  else
    if(this  < 2)  & 
    call ReadC(7,inp,tn,ts,te)  
    read(inp(ts(1):te(1)),'(A80)') namCut
    answer=namCut
    call wait
    read(inp(ts(2):te(2)),*) x_o 
    read(inp(ts(3):te(3)),*) y_o
    read(inp(ts(4):te(4)),*) Href
  end if
  call wait
  k = 0
  do c=1, NC
    if( Approx( xc(c),x_o,1.0e-4 ) ) then
      if( Approx( yc(c),y_o,1.0e-4)  ) then
        k = k + 1
      end if
    end if
  end do
  Npoints = k  
  write(*,*)' ovdje i probe2d'
  allocate(z_p(k));      z_p     =0.0
  allocate(Np(k));       Np      =0 
  allocate(Ump(k));      Ump     =0.0
  allocate(Vmp(k));      Vmp     =0.0
  allocate(Wmp(k));      Wmp     =0.0
  allocate(Kin_mp(k));   Kin_mp  =0.0
  allocate(Eps_mp(k));   Eps_mp  =0.0        
  allocate(v_2_mp(k));   v_2_mp  =0.0
  allocate(f22_mp(k));   f22_mp  =0.0
  allocate(uup(k))   ;   uup     =0.0
  allocate(vvp(k))   ;   vvp     =0.0
  allocate(wwp(k))   ;   wwp     =0.0
  allocate(uvp(k))   ;   uvp     =0.0
  allocate(uwp(k))   ;   uwp     =0.0
  allocate(vwp(k))   ;   vwp     =0.0  
  allocate(Utau_a(k));   Utau_a  =0.0  
  allocate(Umean(k));    Umean   =0.0
  allocate(Vmean(k));    Vmean   =0.0
  allocate(Wmean(k));    Wmean   =0.0
  allocate(Utau_i(k));   Utau_i  =0.0 
  allocate(VISt_mp(k));  VISt_mp =0.0 
  allocate(VIS_mp(k));   VIS_mp  =0.0 
  allocate(order(k));    order   =0 
  allocate(T_mp(k));     T_mp    =0.0 
  allocate(Pk_mp(k));    Pk_mp    =0.0 

!+++++++++++++++++++++++++++++!
!   Finding a cut lines cells !
!   values                    !
!+++++++++++++++++++++++++++++!
  k = 0
  do c = 1, NC
    if( Approx( xc(c),x_o,1.0e-4 ) ) then
      if( Approx( yc(c),y_o,1.0e-4)  ) then
          Np(k)    =  Np(k) + 1
          k        =  k + 1
          order(k) =  k
          z_p(k)   =  zc(c)
      end if
    end if
  end do 
   
  if(K == 0) then
    write(*,*) 'Did not find any point in define cut line position'
    return
  else
    call RISort(z_p, order, Npoints, 2)
  end if
  k = 0
  do pl=1,Npoints
  do c = 1, NC
    if( Approx( xc(c),x_o,1.0e-4 ) ) then
      if( Approx( yc(c),y_o,1.0e-4)  ) then
        if( Approx( zc(c),z_p(pl),1.0e-5)  ) then
          k        =  k + 1
          Ump(k)   =  U % n(c)
          Vmp(k)   =  V % n(c)
          Wmp(k)   =  W % n(c)
          if(pl == 1) WD  =  WallDs(c)
          VISt_mp(k) = VISt(c)
!          uvp(k)    =  VISt_r(c)
          uvp(k)    =  VISt(c)*Uz(c)
!IZMJENA
          if(HOT == YES) T_mp(k) = T % n(c)
          if(SIMULA==SPA_ALL) then
            VIS_mp(k) = VIS % n(c)
          end if
          if(SIMULA==K_EPS.or.SIMULA==LES) then
            Kin_mp(k) = Kin % n(c)
            Eps_mp(k) = Eps % n(c)  
          end if
          if(SIMULA==K_EPS_VV) then
            Kin_mp(k) = Kin % n(c)
            Eps_mp(k) = Eps % n(c)
            v_2_mp(k) = v_2 % n(c)
            f22_mp(k) = f22 % n(c)
            Pk_mp(k)  = Pk(c)
            Utau_a(k) = Ufmean(c)
          end if   
!          if(SIMULA==WOLF) then
!            u_tau  = (VISc*Ump(1)/WD)**0.5
!            Re_tau = u_tau*Href/VISc
!          end if 
!          if(SIMULA==LES) then
!            Umean(k)  =  U % mean(c)
!            Vmean(k)  =  V % mean(c)
!            Wmean(k)  =  W % mean(c)
!            uup(k)    =  uu % n(c) - U % mean(c)*U % mean(c)
!            vvp(k)    =  vv % n(c) - V % mean(c)*V % mean(c)   
!            wwp(k)    =  ww % n(c) - W % mean(c)*W % mean(c) 
!            uvp(k)    =  uv % n(c) - U % mean(c)*V % mean(c)  
!            uwp(k)    =  uw % n(c) - U % mean(c)*W % mean(c)  
!            vwp(k)    =  vw % n(c) - V % mean(c)*W % mean(c) 
!            Kin_mp(k) =  0.5*(uup(k) + vvp(k) + wwp(k))
!          end if 
          end if
        end if
      end if
    end do
    end do
    write(*,*)'finished sorting'
!+++++++++++++++++++++++++++++++++++!
!   Calculating a friction velocity !
!+++++++++++++++++++++++++++++++++++!
  

!  if(SIMULA==LES) then
!    tau_w  = VISc * 0.5*(Umean(1)+Umean(Npoints))/(WD)
!    tau_w_i= VISc * 0.5*(Ump(1)+Ump(Npoints))/(WD)
!    u_tau  = sqrt(abs(tau_w))
!    u_tau_i= sqrt(abs(tau_w_i))
!    Re_tau = u_tau*Href/VISc
!  end if 

  if(SIMULA==K_EPS.or.SIMULA==LES) then
    u_tau  = Kin_mp(1)**0.5*Cmu**0.25 
    Re_tau = u_tau*Href/VISc
    if(MODE == J_L.or.MODE == S_L_Y.or.MODE==NAG) then
      u_tau  = (VISc*Ump(1)/WD)**0.5 
      Re_tau = u_tau*Href/VISc
    end if
  end if 

!  if(SIMULA==WOLF) then
!    u_tau  = (VISc*Ump(1)/WD)**0.5
!    Re_tau = u_tau*Href/VISc
!  end if

  if(SIMULA==SPA_ALL.or.SIMULA==DES_SPA) then
    u_tau  = (VISc*Ump(1)/WD)**0.5
    Re_tau = u_tau*Href/VISc
  end if 

  if(SIMULA==K_EPS_VV) then
    u_tau  = (VISc*abs(Ump(1))/WD)**0.5
    Re_tau = u_tau*Href/VISc
    if(MODE == WF) then
      u_tau  = Utau_a(1) 
      Re_tau = u_tau*Href/VISc
    end if
  end if 

!<<<<<<<<<<<<<<<<<<<<<<<<<!
!     creating  files     !
!<<<<<<<<<<<<<<<<<<<<<<<<<!
 
    namCut(len_trim(namCut)+1:len_trim(namCut)+3) = 'CL.'
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!---- write the values at the probes
    if(SIMULA==K_EPS.or.SIMULA==LES) then                         
      open(9, FILE=namCut)
      write(9,'(A1,I8)') '#', Npoints
      write(9,'(A80)') '# 1-z  2-U  3-V  4-W  5-Kin  6-Eps  7-uv'  
      do pl=1,Npoints 
        write(9,'(8F16.8)')                      &
                     z_p(pl),                    &
                     Ump(pl), Vmp(pl), Wmp(pl),  &
                     Kin_mp(pl), Eps_mp(pl),     &
                     uvp(pl), VISt_mp(pl)                    
      end do
      close(9)
    end if
       
    if(SIMULA==SPA_ALL.or.SIMULA==DES_SPA) then                         
      open(9, FILE=namCut)
      write(9,'(A1,I8)') '#', Npoints
      write(9,'(A80)') '# 1-z  2-U  3-V  4-W  '  
      do pl=1,Npoints 
        write(9,'(7F12.6)')                      &
                     z_p(pl),                    &
                     Ump(pl), Vmp(pl), Wmp(pl), VISt_mp(pl), VIS_mp(pl), uvp(pl)
      end do
      close(9)
    end if

!    if(SIMULA==WOLF) then
!      open(9, FILE=namCut)
!      write(9,'(A1,I8)') '#', Npoints
!      write(9,'(A80)') '# 1-z  2-U  3-V  4-W  '
!      do pl=1,Npoints
!        write(9,'(7F12.6)')                      &
!                     z_p(pl),                    &
!                     Ump(pl), Vmp(pl), Wmp(pl), Kin_mp(pl), &
!                     VISt_mp(pl), uvp(pl)
!      end do
!      close(9)
!    end if
       
    if(SIMULA==K_EPS_VV) then                      
      open(9, FILE=namCut)
      if(name /= 'Backstep') then
      write(9,'(A1,I8)') '#', Npoints
      write(9,'(A80)') '# 1-z  2-U  3-V  4-W  5-Kin  6-Eps  7-uv  8-vv  9-f22'  
      if(HOT == YES) then
        do pl=1,Npoints 
          write(9,'(12F16.7)')                      &
                       z_p(pl),                    &
                       Ump(pl), Vmp(pl), Wmp(pl),T_mp(pl),  &
                       Kin_mp(pl), Eps_mp(pl),     &  
                       uvp(pl), v_2_mp(pl),        &
                       f22_mp(pl),   &
                       VISt_mp(pl), Pk_mp(pl)
        end do
      else
        do pl=1,Npoints
          write(9,'(11F16.7)')                      &
                       z_p(pl),                    &
                       Ump(pl), Vmp(pl), Wmp(pl),  &
                       Kin_mp(pl), Eps_mp(pl),     &  
                       uvp(pl), v_2_mp(pl),        &
                       f22_mp(pl),   &
                       VISt_mp(pl), Pk_mp(pl)
        end do
        close(9)
      end if 
      end if
      if(name == 'Backstep') then
        open(9, FILE=namCut)
        Uref = 11.3
        Uref = 2.0
        write(9,'(A1,I8)') '#', Npoints
        write(9,'(A80)') '# 1-z  2-U  3-V  4-W  5-Kin  6-Eps  7-uv  8-vv  9-f22'  
        if(HOT == YES) then
          do pl=1,Npoints 
            write(9,'(11F16.7)')                      &
                         z_p(pl)/0.038,                    &
                         Ump(pl)/Uref, Vmp(pl)/Uref, Wmp(pl)/Uref,T_mp(pl),   &
                         Kin_mp(pl)/Uref**2.0, Eps_mp(pl)*0.038/Uref**3.0,     &  
                         uvp(pl)/Uref**2.0, v_2_mp(pl)/Uref**2.0,        &
                         f22_mp(pl)*0.038/Uref**2.0
          end do
        else  
          do pl=1,Npoints
            write(9,'(10F16.7)')                      &
                         z_p(pl)/0.038,                    &
                         Ump(pl)/Uref, Vmp(pl)/Uref, Wmp(pl)/Uref,  &
                         Kin_mp(pl)/Uref**2.0, Eps_mp(pl)*0.038/Uref**3.0,     &  
                         uvp(pl)/Uref**2.0, v_2_mp(pl)/Uref**2.0,        &
                         f22_mp(pl)*0.038/Uref**2.0
          end do
        end if
      end if
      close(9)
    end if

!    if(SIMULA==LES) then
!      open(9, FILE=namCut)
!      write(9,'(A1,I8)') '#', Npoints
!      write(9,'(A80)') '# 1-z  2-Umean  3-Vmean  4-Wmean  5-U   6-V   7-W  8-uu  9-vv 10-ww 11-uv 12-uw  13-vw'
!      do pl=1,Npoints 
!        write(9,'(13F12.6)')                           &
!                     z_p(pl),                          &
!                     Umean(pl), Vmean(pl), Wmean(pl),  &
!                     Ump(pl), Vmp(pl), Wmp(pl),        &
!                     uup(pl), vvp(pl), wwp(pl),        &
!                     uvp(pl), uwp(pl), vwp(pl)
!      end do
!      close(9)
!    end if

!    open(10,file = 'data_norm.dat')
!    if(TEST==YES) then
!      write(10,'(A1,I8)') '#', Npoints
!      write(10,*) '# variables'
!      do pl=1,Npoints 
!        write(10,'(9F16.7)')                      &
!                     z_p(pl)/0.038,                    &
!                     Ump(pl)/2.0, Vmp(pl)/2.0, Wmp(pl)/2.0,  &
!                     Kin_mp(pl)/4.0, Eps_mp(pl)*0.038/8.0,     &  
!                     uvp(pl)/4.0, v_2_mp(pl)/4.0,        &
!                     f22_mp(pl)*0.038/2.0  
!      end do
!    end if
!    close(10)
        
!+++++++++++++++++++++++++++++++++++
!   Calculate normalise values
!+++++++++++++++++++++++++++++++++++

  do pl=1,Npoints
            if(SIMULA==K_EPS.or.SIMULA==LES) then
              z_p(pl)   =  (z_p(pl))*u_tau/VISc
              Ump(pl)   =  Ump(pl)/u_tau
              Vmp(pl)   =  Vmp(pl)/u_tau
              Wmp(pl)   =  Wmp(pl)/u_tau      
              Kin_mp(pl) = Kin_mp(pl)/u_tau**2.0
              Eps_mp(pl) = Eps_mp(pl)*VISc/u_tau**4.0
              uvp(pl)    = uvp(pl)/u_tau**2.0
              VISt_mp(pl) = VISt_mp(pl)/VISc
            end if
            if(SIMULA==SPA_ALL.or.SIMULA==DES_SPA) then
              z_p(pl)   =  (z_p(pl))*u_tau/VISc
              Ump(pl)   =  Ump(pl)/u_tau
              Vmp(pl)   =  Vmp(pl)/u_tau
              Wmp(pl)   =  Wmp(pl)/u_tau      
              VISt_mp(pl) = VISt_mp(pl)/VISc
              VIS_mp(pl) = VIS_mp(pl)/VISc
              uvp(pl)    = uvp(pl)/u_tau**2.0
            end if
            if(SIMULA==K_EPS_VV) then
              z_p(pl)   =  z_p(pl)*u_tau/VISc
              Ump(pl)   =  Ump(pl)/u_tau
              Vmp(pl)   =  Vmp(pl)/u_tau
              Wmp(pl)   =  Wmp(pl)/u_tau      
              Kin_mp(pl) = Kin_mp(pl)/u_tau**2.0
              Eps_mp(pl) = Eps_mp(pl)*VISc/u_tau**4.0
              Pk_mp(pl)  = Pk_mp(pl)*VISc/u_tau**4.0
              uvp(pl)    = uvp(pl)/u_tau**2.0 
              v_2_mp(pl) = v_2_mp(pl)/u_tau**2.0
              f22_mp(pl) =  f22_mp(pl)*VISc/u_tau**2.0 
              VISt_mp(pl) = VISt_mp(pl)/VISc
            end if
!            if(SIMULA==WOLF) then
!              z_p(pl)   =  (z_p(pl))*u_tau/VISc
!              Ump(pl)   =  Ump(pl)/u_tau
!              Vmp(pl)   =  Vmp(pl)/u_tau
!              Wmp(pl)   =  Wmp(pl)/u_tau
!              Kin_mp(pl) = Kin_mp(pl)/u_tau**2.0
!              VISt_mp(pl) = VISt_mp(pl)/VISc
!              uvp(pl)    = uvp(pl)/u_tau**2.0
!            end if
!            if(SIMULA==LES) then
!              z_p(pl)   =  (z_p(pl)+1.0)*u_tau/VISc
!              Ump(pl)   =  Ump(pl)/u_tau
!              Vmp(pl)   =  Vmp(pl)/u_tau
!              Wmp(pl)   =  Wmp(pl)/u_tau       
!              Umean(pl)  =  Umean(pl)/u_tau
!              Vmean(pl)  =  Vmean(pl)/u_tau 
!              Wmean(pl)  =  Wmean(pl)/u_tau 
!              uup(pl)    =  uup(pl)/u_tau**2.0 
!              vvp(pl)    =  vvp(pl)/u_tau**2.0  
!              wwp(pl)    =  wwp(pl)/u_tau**2.0  
!              uvp(pl)    =  uvp(pl)/u_tau**2.0  
!              uwp(pl)    =  uwp(pl)/u_tau**2.0  
!              vwp(pl)    =  vwp(pl)/u_tau**2.0  
!              Kin_mp(pl) =  0.5*(uup(pl)+vvp(pl)+wwp(pl))/u_tau**2.0  
!            end if
      end do

!<<<<<<<<<<<<<<<<<<<<<<<<<!
!     create 2Drplus file !
!<<<<<<<<<<<<<<<<<<<<<<<<<!
  if(This < 2) then
    namCut(len_trim(namCut)+1:len_trim(namCut)+7) = 'dat'
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!---- write the values at the probes
    if(SIMULA==K_EPS.or.SIMULA==LES) then
      open(9, FILE=namCut)
      write(9,'(A1,I8,A10,F10.7)') '#', Npoints, 'u_tau=',u_tau 
      write(9,'(A40)')'# 1-z  2-U  3-V  4-W  5-Kin  6-Eps  7-uv'
      do pl=1,Npoints
        write(9,'(8E18.6)')                      &
                     z_p(pl),                    &
                     Ump(pl), Vmp(pl), Wmp(pl),  &
                     Kin_mp(pl), Eps_mp(pl),     &
                     uvp(pl), VISt_mp(pl)
      end do
      close(9)
    end if
 
    if(SIMULA==SPA_ALL.or.SIMULA==DES_SPA) then
      open(9, FILE=namCut)
      write(9,'(A1,I8,A10,F8.4)') '#', Npoints, 'u_tau=',u_tau 
      write(9,'(A80)') '# 1-z  2-U  3-V  4-W  '
      do pl=1,Npoints
        write(9,'(7F12.6)')                      &
                     z_p(pl),                    &
                     Ump(pl), Vmp(pl), Wmp(pl), VISt_mp(pl),VIS_mp(pl), uvp(pl)
      end do
      close(9)
    end if
!    if(SIMULA==WOLF) then
!      open(9, FILE=namCut)
!      write(9,'(A1,I8,A10,F8.4)') '#', Npoints, 'u_tau=',u_tau
!      write(9,'(A80)') '# 1-z  2-U  3-V  4-W  '
!      do pl=1,Npoints
!        write(9,'(7F12.6)')                      &
!                     z_p(pl),                    &
!                     Ump(pl), Vmp(pl), Wmp(pl), Kin_mp(pl)  &
!                   , VISt_mp(pl), uvp(pl)
!      end do
!      close(9)
!    end if
    if(SIMULA==K_EPS_VV) then
      open(9, FILE=namCut)
      write(9,*) '#', Npoints, 'u_tau=',u_tau, 'Ret=', Re_tau 
      write(9,'(A80)') '# 1-z  2-U  3-V  4-W  5-Kin  6-Eps  7-uv  8-vv  9-f22'
      if(HOT == YES) then
        do pl=1,Npoints
          write(9,'(12F16.8)')                      &
                       z_p(pl),                    &
                       Ump(pl), Vmp(pl), Wmp(pl), T_mp(pl), &
                       Kin_mp(pl), Eps_mp(pl),     &
                       uvp(pl), v_2_mp(pl),        &
                       f22_mp(pl), VISt_mp(pl), Pk_mp(pl)
        end do
      else
        do pl=1,Npoints 
          write(9,'(11F16.8)')                      &
                       z_p(pl),                    &
                       Ump(pl), Vmp(pl), Wmp(pl),  &
                       Kin_mp(pl), Eps_mp(pl),     &
                       uvp(pl), v_2_mp(pl),        &
                       f22_mp(pl), VISt_mp(pl), Pk_mp(pl)
        end do
      end if
      close(9)
    end if
 
!    if(SIMULA==LES) then
!      open(9, FILE=namCut)
!      write(9,'(A1,I8,A10,F8.4)') '#', Npoints, 'u_tau=',u_tau 
!      write(9,'(A80)') '#  1-z  2-Umean  3-Vmean  4-Wmean  5-U   6-V   7-W  8-uu  9-vv 10-ww 11-uv 12-uw  13-vw'
!      do pl=1,Npoints
!        write(9,'(13F12.6)')                       &
!                     z_p(pl),                          &
!                     Umean(pl), Vmean(pl), Wmean(pl),  &
!                     Ump(pl), Vmp(pl), Wmp(pl),        &
!                     uup(pl), vvp(pl), wwp(pl),        &
!                     uvp(pl), uwp(pl), vwp(pl)
!      end do
!      close(9)
!    end if 
  end if
  namCut = name
 
  deallocate(Np)
  deallocate(Ump)
  deallocate(Vmp)
  deallocate(Wmp)
  deallocate(Kin_mp)
  deallocate(Eps_mp)
  deallocate(v_2_mp)
  deallocate(f22_mp)
  deallocate(uup)
  deallocate(vvp)
  deallocate(wwp)
  deallocate(uvp)
  deallocate(uwp)
  deallocate(vwp)
  deallocate(Utau_a)
  deallocate(Umean)
  deallocate(Vmean)
  deallocate(Wmean)
  deallocate(Utau_i)   
  deallocate(VISt_mp)   
  deallocate(VIS_mp)   
  deallocate(Pk_mp)   

  END SUBROUTINE UserProbe2D
