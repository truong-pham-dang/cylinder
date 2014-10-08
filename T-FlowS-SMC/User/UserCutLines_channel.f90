!======================================================================!
  SUBROUTINE UserCutLines_channel(y) 
!----------------------------------------------------------------------!
! This program reads name.1D file created by NEU or GEN and averages    
! the results in the homogeneous directions. 
! The results are writen in files name_res.dat and name_res_plus.dat
!----------------------------------------------------------------------!
    USE all_mod
    USE allp_mod
    USE les_mod
    USE pro_mod
    USE par_mod
    USE rans_mod
!----------------------------------------------------------------------!
    IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
    REAL :: y(-NbC:NC)
    REAL :: Ufric, Wall_near 
!------------------------------[Calling]-------------------------------!
    INTERFACE
      LOGICAL FUNCTION Approx(A,B,tol)
        REAL           :: A,B
        REAL, OPTIONAL :: tol
      END FUNCTION Approx
    END INTERFACE 
!-------------------------------[Locals]-------------------------------!
    INTEGER             :: Nprob, pl, c, i, count, s, c1, c2
    CHARACTER           :: namCoo*80, namPro*80, answer*80, namRes*80
    CHARACTER           :: namRes_plus*80
    REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:), ind(:),       & 
                                 VIStp(:), var_2(:), var_3(:), &
                                 Wall_p(:), Ufric_p(:), Kinp(:), Epsp(:), &
                                 v2p(:), fp(:) 
    INTEGER,ALLOCATABLE :: Np(:), Ncount(:)
    REAL                :: R, Urad_mean, Utan_mean, dummy, Lscale, Twall
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $  
!======================================================================!

    namPro = name

    namRes = name
    namRes(len_trim(name)+1:len_trim(name)+8) = '_res.dat'
    namRes_plus = name
    namRes_plus(len_trim(name)+1:len_trim(name)+13) = '_res_plus.dat'
!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
    namCoo = name
    namCoo(len_trim(name)+1:len_trim(name)+3) = '.1D'
    if(this < 2)  write(6, *) '# Now reading the file:', namCoo
    open(9, FILE=namCoo)

!---- write the number of searching intervals 
    read(9,*) Nprob
    allocate(z_p(Nprob*2))
    allocate(ind(Nprob*2))

!---- read the intervals positions
    do pl=1,Nprob
      read(9,*) ind(pl), z_p(pl) 
    end do
    close(9)

    call GraPhi(U % n, 1, Ux,.TRUE.)    ! dU/dx
    call GraPhi(U % n, 2, Uy,.TRUE.)    ! dU/dy
    call GraPhi(U % n, 3, Uz,.TRUE.)    ! dU/dz
    call GraPhi(V % n, 1, Vx,.TRUE.)    ! dV/dx
    call GraPhi(V % n, 2, Vy,.TRUE.)    ! dV/dy
    call GraPhi(V % n, 3, Vz,.TRUE.)    ! dV/dz
    call GraPhi(W % n, 1, Wx,.TRUE.)    ! dW/dx
    call GraPhi(W % n, 2, Wy,.TRUE.)    ! dW/dy
    call GraPhi(W % n, 3, Wz,.TRUE.)    ! dW/dz


    call SSORT (z_p, ind, Nprob, 0)
    write(*,*) Nprob
    
    call CalcShear(U % n, V % n, W % n, Shear)

    allocate(Np(Nprob));    Np=0 
    allocate(Wall_p(Nprob));Wall_p=0.0
    allocate(Ump(Nprob));   Ump=0.0
    allocate(Vmp(Nprob));   Vmp=0.0
    allocate(Wmp(Nprob));   Wmp=0.0
    allocate(uup(Nprob));   uup=0.0
    allocate(vvp(Nprob));   vvp=0.0
    allocate(wwp(Nprob));   wwp=0.0
    allocate(uvp(Nprob));   uvp=0.0
    allocate(uwp(Nprob));   uwp=0.0
    allocate(vwp(Nprob));   vwp=0.0
    allocate(Kinp(Nprob));  Kinp=0.0
    allocate(Epsp(Nprob));  Epsp=0.0
    allocate(v2p(Nprob));   v2p=0.0
    allocate(fp(Nprob));    fp=0.0
    allocate(Ksgsp(Nprob)); Ksgsp=0.0
    allocate(VIStp(Nprob));  VIStp=0.0
    allocate(Ufric_p(Nprob)); Ufric_p = 0.0

    allocate(Ncount(Nprob)); Ncount=0
    count = 0
    if(HOT==YES) then
      allocate(Tmp(Nprob));   Tmp=0.0
      allocate(TTp(Nprob));   TTp=0.0
      allocate(uTp(Nprob));   uTp=0.0
      allocate(vTp(Nprob));   vTp=0.0
      allocate(wTp(Nprob));   wTp=0.0
    end if  

    Lscale = 0.0
!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
    do i = 1, Nprob-1
      do c=1, NC
        Lscale = max(Lscale,WallDs(c))
        if(y(c) > (z_p(i+1)) .and. y(c) < (z_p(i))) then
          if(ROT == YES) then
            Wall_p(i) = Wall_p(i) + y(c) 
          else
            Wall_p(i) = Wall_p(i) + WallDs(c) 
          end if 
          if(SIMULA==K_EPS_VV) then
            Ump(i)   = Ump(i) + U % n(c)
            Vmp(i)   = Vmp(i) + V % n(c)
            Wmp(i)   = Wmp(i) + W % n(c)
            Kinp(i)  = Kinp(i) + Kin % n(c) 
            Epsp(i)  = Epsp(i) + Eps % n(c) 
            v2p(i)   = v2p(i) + v_2 % n(c) 
            fp(i)    = fp(i) + f22 % n(c) 
            uup(i)   = uup(i) + 2.0*VISt(c)*Ux(c) 
            vvp(i)   = vvp(i) + 2.0*VISt(c)*Vy(c)
            wwp(i)   = wwp(i) + 2.0*VISt(c)*Wz(c)
            uvp(i)   = uvp(i) + VISt(c)*(Uy(c) + Vx(c))
            uwp(i)   = uwp(i) + VISt(c)*(Uz(c) + Wx(c))
            vwp(i)   = vwp(i) + VISt(c)!*(Vz(c) + Wy(c))
            if(IsNearWall(c)) then
              if(Ynd(c)>5.0) then      
                Ufric_p(i) = Ufric_p(i) + sqrt(TauWall(c))
              else
                Ufric_p(i) = Ufric_p(i) + (VISc * (U % n(c)**2 + V % n(c)**2 + &
                W% n(c)**2)**0.5/WallDs(c))**0.5      
              end if
            end if 
          else if(SIMULA==EBM) then
            Ump(i)   = Ump(i) + U % n(c)
            Vmp(i)   = Vmp(i) + V % n(c)
            Wmp(i)   = Wmp(i) + W % n(c)
            Kinp(i)  = Kinp(i) + Kin % n(c) 
            Epsp(i)  = Epsp(i) + Eps % n(c) 
            fp(i)    = fp(i) + 0.0 
            uup(i)   = uup(i) + uu%n(c) 
            vvp(i)   = vvp(i) + vv%n(c) 
            wwp(i)   = wwp(i) + ww%n(c)
            uvp(i)   = uvp(i) + uv%n(c)
            uwp(i)   = uwp(i) + uw%n(c)
            vwp(i)   = vwp(i) + vw%n(c)
            if(IsNearWall(c)) then
              Ufric_p(i) = Ufric_p(i) + (VISc * (U % n(c)**2 + V % n(c)**2 + W % n(c)**2)**0.5/WallDs(c))**0.5
            end if 
          else if(SIMULA==ZETA) then
            Ump(i)   = Ump(i) + U % n(c)
            Vmp(i)   = Vmp(i) + V % n(c)
            Wmp(i)   = Wmp(i) + W % n(c)
            Kinp(i)  = Kinp(i) + Kin % n(c) 
            Epsp(i)  = Epsp(i) + Eps % n(c) 
            v2p(i)   = v2p(i) + v_2 % n(c) 
            fp(i)    = fp(i) + f22 % n(c) 
            uup(i)   = uup(i) + 2.0*VISt(c)*Ux(c) 
            vvp(i)   = vvp(i) + 2.0*VISt(c)*Vy(c)
            wwp(i)   = wwp(i) + 2.0*VISt(c)*Wz(c)
            uvp(i)   = uvp(i) + VISt(c)*(Uy(c) + Vx(c))
            uwp(i)   = uwp(i) + VISt(c)*(Uz(c) + Wx(c))
            vwp(i)   = vwp(i) + VISt(c)!*(Vz(c) + Wy(c))
            if(IsNearWall(c)) then
              if(Ynd(c)>5.0) then      
                Ufric_p(i) = Ufric_p(i) + sqrt(TauWall(c))
              else
                Ufric_p(i) = Ufric_p(i) + (VISc * (U % n(c)**2 + V % n(c)**2 + &
                W% n(c)**2)**0.5/WallDs(c))**0.5      
              end if
            end if 
          else if(SIMULA==K_EPS) then
            Ump(i)   = Ump(i) + U % n(c)
            Vmp(i)   = Vmp(i) + V % n(c)
            Wmp(i)   = Wmp(i) + W % n(c)
            Kinp(i)  = Kinp(i) + Kin % n(c) 
            Epsp(i)  = Epsp(i) + Eps % n(c) 
            uup(i)   = uup(i) + 2.0*VISt(c)*Ux(c) 
            vvp(i)   = vvp(i) + 2.0*VISt(c)*Vy(c)
            wwp(i)   = wwp(i) + 2.0*VISt(c)*Wz(c)
            uvp(i)   = uvp(i) + VISt(c)*(Uy(c) + Vx(c))
            uwp(i)   = uwp(i) + VISt(c)*(Uz(c) + Wx(c))
            vwp(i)   = vwp(i) + VISt(c)*(Vz(c) + Wy(c))
            if(MODE == HRe) then
              if(IsNearWall(c)) then
                Ufric_p(i) = Ufric_p(i) + Cmu**0.25*Kin % n(c)**0.5
              end if  
            else
              if(IsNearWall(c)) then
                Ufric_p(i) = Ufric_p(i) + (VISc * (U % n(c)**2 + V % n(c)**2 + W % n(c)**2)**0.5/WallDs(c))**0.5
              end if  
            end if
          else if(SIMULA==SPA_ALL) then
            Ump(i)   = Ump(i) + U % n(c)
            Vmp(i)   = Vmp(i) + V % n(c)
            Wmp(i)   = Wmp(i) + W % n(c)
            VIStp(i) = VIStp(i) + VISt(c)
            uup(i)   = uup(i) + 2.0*VISt(c)*Ux(c) 
            vvp(i)   = vvp(i) + 2.0*VISt(c)*Vy(c)
            wwp(i)   = wwp(i) + 2.0*VISt(c)*Wz(c)
            uvp(i)   = uvp(i) + VISt(c)*(Uy(c) + Vx(c))
            uwp(i)   = uwp(i) + VISt(c)*(Uz(c) + Wx(c))
            vwp(i)   = vwp(i) + VISt(c)*(Vz(c) + Wy(c))
            if(IsNearWall(c)) then
              Ufric_p(i) = Ufric_p(i) + (VISc * (U % n(c)**2 + V % n(c)**2 + W % n(c)**2)**0.5/WallDs(c))**0.5
            end if 
          else if(SIMULA==LES.or.SIMULA==DES_SPA.or.SIMULA == DNS) then
            Ump(i)   = Ump(i) + U % mean(c)
            Vmp(i)   = Vmp(i) + V % mean(c)
            Wmp(i)   = Wmp(i) + W % mean(c)
            uup(i)   = uup(i) + (uu % mean(c)- U % mean(c) * U % mean(c))
            vvp(i)   = vvp(i) + (vv % mean(c)- V % mean(c) * V % mean(c))
            wwp(i)   = wwp(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
            uvp(i)   = uvp(i) + (uv % mean(c)- U % mean(c) * V % mean(c))
            uwp(i)   = uwp(i) + (uw % mean(c)- U % mean(c) * W % mean(c))
            vwp(i)   = vwp(i) + (vw % mean(c)- V % mean(c) * W % mean(c))
            if(IsNearWall(c)) then
              Ufric_p(i) = Ufric_p(i) + (VISc * (U % mean(c)**2 + V % mean(c)**2 + W % mean(c)**2)**0.5/WallDs(c))**0.5
            end if 
          else if(SIMULA==HYB_ZETA) then
            Ump(i)   = Ump(i) + U % mean(c)
            Vmp(i)   = Vmp(i) + V % mean(c)
            Wmp(i)   = Wmp(i) + W % mean(c)
            uup(i)   = uup(i) + (uu % mean(c)- U % mean(c) * U % mean(c))
            vvp(i)   = vvp(i) + (vv % mean(c)- V % mean(c) * V % mean(c))
            wwp(i)   = wwp(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
            uvp(i)   = uvp(i) + (uv % mean(c)- U % mean(c) * V % mean(c))
            uwp(i)   = uwp(i) + (uw % mean(c)- U % mean(c) * W % mean(c))
            vwp(i)   = vwp(i) + (vw % mean(c)- V % mean(c) * W % mean(c))
            Kinp(i)  = Kinp(i) + Kin % n(c) 
            Epsp(i)  = Epsp(i) + Eps % n(c) 
            v2p(i)   = v2p(i) + v_2 % n(c) 
            fp(i)    = fp(i) + f22 % n(c) 
            if(IsNearWall(c)) then
              Ufric_p(i) = Ufric_p(i) + (VISc * (U % mean(c)**2 + V % mean(c)**2 + W % mean(c)**2)**0.5/WallDs(c))**0.5
            end if 
          else if(SIMULA==HYB_PITM) then
            Ump(i)   = Ump(i) + U % mean(c)
            Vmp(i)   = Vmp(i) + V % mean(c)
            Wmp(i)   = Wmp(i) + W % mean(c)
            uup(i)   = uup(i) + (uu % mean(c)- U % mean(c) * U % mean(c))
            vvp(i)   = vvp(i) + (vv % mean(c)- V % mean(c) * V % mean(c))
            wwp(i)   = wwp(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
            uvp(i)   = uvp(i) + (uv % mean(c)- U % mean(c) * V % mean(c))
            uwp(i)   = uwp(i) + (uw % mean(c)- U % mean(c) * W % mean(c))
            vwp(i)   = vwp(i) + (vw % mean(c)- V % mean(c) * W % mean(c))
            Kinp(i)  = Kinp(i) + Kin % n(c) 
            Epsp(i)  = Epsp(i) + Eps % n(c) 
            fp(i)    = fp(i) + 1.5 + 0.4/(1.0 + 3.0*(0.41*WallDs(c)/volume(c)**0.33333)**0.5 ) 
            if(IsNearWall(c)) then
              Ufric_p(i) = Ufric_p(i) + (VISc * (U % mean(c)**2 + V % mean(c)**2 + W % mean(c)**2)**0.5/WallDs(c))**0.5
            end if 
          end if 
 
          If(SIMULA==LES.or.SIMULA==DES_SPA.or.SIMULA==HYB_PITM) then
            VIStp(i) = VIStp(i) + VISt(c)
          end if
          If(SIMULA==HYB_ZETA) then
            VIStp(i) = VIStp(i) + VISt_eff(c)
          end if

          if(HOT==YES) then
            if(SIMULA == LES.or.SIMULA == DES_SPA.or.SIMULA == DNS &
               .or.SIMULA==HYB_PITM.or.SIMULA==HYB_ZETA) then
              Tmp(i)   = Tmp(i) + T % mean(c)
              TTp(i)   = TTp(i) + (TT % mean(c) - T % mean(c) * T % mean(c))
              uTp(i)   = uTp(i) + (uT % mean(c) - u % mean(c) * T % mean(c))
              vTp(i)   = vTp(i) + (vT % mean(c) - v % mean(c) * T % mean(c))
              wTp(i)   = wTp(i) + (wT % mean(c) - w % mean(c) * T % mean(c))
            else
              Tmp(i)   = Tmp(i) + T % n(c)
            end if
          end if
          Ncount(i) = Ncount(i) + 1
        end if
      end do 
    end do 
   if(HOT==YES) then
     do s=1,NS
       c1=SideC(1,s)
       c2=SideC(2,s)
       if(c2  < 0) then
         if( TypeBC(c2) ==  WALL.or. TypeBC(c2) ==  WALLFL) then
           Twall = T % n(c2) 
         end if
       end if
     end do
   end if

!---- average over all processors
    do pl=1, Nprob-1
      call IGlSum(Ncount(pl))

      call GloSum(Wall_p(pl))

      call GloSum(Ump(pl))
      call GloSum(Vmp(pl))
      call GloSum(Wmp(pl))

      call GloSum(Kinp(pl))
      call GloSum(Epsp(pl))
      call GloSum(fp(pl))
      call GloSum(v2p(pl))

      call GloSum(uup(pl))
      call GloSum(vvp(pl))
      call GloSum(wwp(pl))

      call GloSum(uvp(pl))
      call GloSum(uwp(pl))
      call GloSum(vwp(pl))
      call GloSum(VIStp(pl))
      call GloSum(Ufric_p(pl))

      count =  count + Ncount(pl) 

      if(HOT==YES) then
        call GloSum(Tmp(pl))
        call GloSum(TTp(pl))
        call GloSum(uTp(pl))
        call GloSum(vTp(pl))
        call GloSum(wTp(pl))
      end if
    end do

    call wait

    do i = 1, Nprob-1
      if(Ncount(i) /= 0) then
        Wall_p(i) =  Wall_p(i)/Ncount(i)
        Ump(i)    =  Ump(i)/Ncount(i)
        Vmp(i)    =  Vmp(i)/Ncount(i)
        Wmp(i)    =  Wmp(i)/Ncount(i)
        uup(i)    =  uup(i)/Ncount(i)
        vvp(i)    =  vvp(i)/Ncount(i)
        wwp(i)    =  wwp(i)/Ncount(i)
        uvp(i)    =  uvp(i)/Ncount(i)
        uwp(i)    =  uwp(i)/Ncount(i)
        vwp(i)    =  vwp(i)/Ncount(i)
        Kinp(i)   =  Kinp(i)/Ncount(i)
        Epsp(i)   =  Epsp(i)/Ncount(i)
        v2p(i)    =  v2p(i)/Ncount(i)
        fp(i)     =  fp(i)/Ncount(i)
        VIStp(i)  =  VIStp(i)/Ncount(i)
        Ufric_p(i)=  Ufric_p(i)/Ncount(i)
        if(HOT == YES) then
          Tmp(i)    =  Tmp(i)/Ncount(i)
          TTp(i)    =  TTp(i)/Ncount(i)
          uTp(i)    =  uTp(i)/Ncount(i)
          vTp(i)    =  vTp(i)/Ncount(i)
          wTp(i)    =  wTp(i)/Ncount(i)
        end if
      end if
    end do 


    if(SIMULA == LES.or.SIMULA==DES_SPA.or.SIMULA==DNS) then
      do i = 1, (Nprob-1)/2
        Ump(i)   =  (Ump(i)+Ump(Nprob-i))/2.0
        Vmp(i)   =  (Vmp(i)+Vmp(Nprob-i))/2.0 
        Wmp(i)   =  (Wmp(i)+Wmp(Nprob-i))/2.0
        uup(i)   =  (uup(i)+uup(Nprob-i))/2.0
        vvp(i)   =  (vvp(i)+vvp(Nprob-i))/2.0
        wwp(i)   =  (wwp(i)+wwp(Nprob-i))/2.0
        uvp(i)   =  (uvp(i)+abs(uvp(Nprob-i)))/2.0
        uwp(i)   =  (uwp(i)+abs(uwp(Nprob-i)))/2.0
        vwp(i)   =  (vwp(i)+vwp(Nprob-i))/2.0
        Ufric_p(i) = (Ufric_p(i) + Ufric_p(Nprob-i))/2
      end do

!      Nprob = (Nprob-1)/2

      Ufric = Ufric_p(1)
    end if

    Ufric = 0.0
    do i = 1, Nprob
      Ufric = max(Ufric_p(i),Ufric)
    end do
10  continue

    if(Ufric == 0.0) then
      if(this < 2) write(*,*) 'Friction velocity is zero in UserCutLines_channel.f90 !'
      return
    end if

    open(3,FILE=namRes)
    write(3,'(A1,2(A10, F13.5))') '#', 'Utau = ', Ufric, 'Re_tau = ', Ufric/VISc
    if(SIMULA == DNS) then
      write(3,'(A1,2X,A80)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8:uv, 9:uw, 10:vw, 11:Kin' 
      do i = 1, Nprob
        if(Ncount(i) /= 0) then
          write(3,'(11E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i),   &
          uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
          vwp(i), 0.5*(uup(i) + vvp(i) + wwp(i))
        end if
      end do 
    else if(SIMULA == LES.or.SIMULA == DES_SPA) then
      if(HOT == YES) then 
        write(3,'(A1,2X,A150)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:T, 6:uu, 7:vv, 8:ww, &
        9:uv, 10:uw, 11:vw, 12: t, 13: uT, 14: vT, 15: wT, 16:Kin, 17: VISt' 
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(3,'(17E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i), Tmp(i), &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
            vwp(i), TTp(i), uTp(i), vTp(i), wTp(i), 0.5*(uup(i) + vvp(i) + wwp(i)), VIStp(i)
          end if
        end do 
      else
        write(3,'(A1,2X,A88)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8:uv, 9:uw, 10:vw, 11:Kin, 12: VISt' 
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(3,'(12E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i),   &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
            vwp(i), 0.5*(uup(i) + vvp(i) + wwp(i)), VIStp(i)
          end if
        end do 
      end if 
    else if(SIMULA == HYB_ZETA) then
      if(HOT==YES) then
        write(3,'(A1,2X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5: T, 6:Kin, 7:Eps, 8:v2, 9:f22, &
        10:uu, 11:vv, 12:ww, 13:uv, 14:uw, 15:vw, 16:t, 17:uT, 18:vT, 19:wT, 20:Kin_res, 21:VISt' 
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(3,'(21E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i), Tmp(i), &
            Kinp(i), Epsp(i), v2p(i), fp(i),      &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
            vwp(i), TTp(i), uTp(i), vTp(i), wTp(i), &
            0.5*(uup(i) + vvp(i) + wwp(i)), VIStp(i)
          end if
        end do 
      else 
        write(3,'(A1,2X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:Kin, 6:Eps, 7:v2, 8:f22, &
        9:uu, 10:vv, 11:ww, 12:uv, 13:uw, 14:vw, 15:Kin_res, 16: VISt' 
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(3,'(16E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i),   &
            Kinp(i), Epsp(i), v2p(i), fp(i),      &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
            vwp(i), 0.5*(uup(i) + vvp(i) + wwp(i)), VIStp(i)
          end if
        end do 
      end if
    else if(SIMULA == HYB_PITM) then
      if(HOT==YES) then
        write(3,'(A1,2X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5: T, 6:Kin, 7:Eps,  &
        10:uu, 11:vv, 12:ww, 13:uv, 14:uw, 15:vw, 16:t, 17:uT, 18:vT, 19:wT, 20:Kin_res, 21:VISt' 
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(3,'(21E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i), Tmp(i), &
            Kinp(i), Epsp(i),      &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
            vwp(i), TTp(i), uTp(i), vTp(i), wTp(i), &
            0.5*(uup(i) + vvp(i) + wwp(i)), VIStp(i)
          end if
        end do 
      else 
        write(3,'(A1,2X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:Kin, 6:Eps, &
        9:uu, 10:vv, 11:ww, 12:uv, 13:uw, 14:vw, 15:Kin_res, 16: VISt' 
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(3,'(16E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i),   &
            Kinp(i), Epsp(i),  &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
            vwp(i), 0.5*(uup(i) + vvp(i) + wwp(i)), VIStp(i)
          end if
        end do 
      end if
    else if(SIMULA == ZETA.or.SIMULA==EBM) then
      if(HOT==YES) then
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8: uv, &
        9: uw, 10: vw, 11:Kin, 12:Eps, 13:zeta, 14:f' 
                                
      do i = Nprob, 1, -1 
          if(Ncount(i) /= 0) then
            write(3,'(15E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i), Tmp(i),  &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i), vwp(i), &
            Kinp(i), Epsp(i), &
            v2p(i), fp(i)
          end if
        end do 
      else
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8: uv, &
                                9: uw, 10: vw, 11:Kin, 12:Eps, 13:zeta, 14:f' 
        do i = Nprob, 1, -1 
          if(Ncount(i) /= 0) then
            write(3,'(14E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i),   &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i), vwp(i), Kinp(i), Epsp(i), &
            v2p(i), fp(i)
          end if
        end do 
      end if
    else if(SIMULA == K_EPS_VV) then
      if(HOT==YES) then
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8: uv, &
        9: uw, 10: vw, 11:Kin, 12:Eps, 13:zeta, 14:f' 
        do i = Nprob, 1, -1
          if(Ncount(i) /= 0) then
            write(3,'(15E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i), Tmp(i),  &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i), vwp(i), &
            Kinp(i), Epsp(i), &
            v2p(i), fp(i)
          end if
        end do 
      else
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8: uv, &
                                9: uw, 10: vw, 11:Kin, 12:Eps, 13:zeta, 14:f' 
        do i = 1, Nprob  !, 1, -1
          if(Ncount(i) /= 0) then
            write(3,'(14E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i),   &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i), vwp(i), Kinp(i), Epsp(i), &
            v2p(i), fp(i)
          end if
        end do 
      end if
    else if(SIMULA == K_EPS) then
      if(HOT==YES) then 
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8: uv, &
                                9: uw, 10: vw, 11:Kin, 12:Eps, 13: VISt' 
        do i = Nprob, 1, -1
          if(Ncount(i) /= 0) then
            write(3,'(14E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i), Tmp(i),  &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i), vwp(i), Kinp(i), Epsp(i), VIStp(i)
          end if
        end do 
      else
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8: uv, &
                                9: uw, 10: vw, 11:Kin, 12:Eps, 13: VISt' 
        do i = Nprob, 1, -1
          if(Ncount(i) /= 0) then
            write(3,'(13E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i), &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i), vwp(i), Kinp(i), Epsp(i), VIStp(i)
          end if
        end do 
      end if
    else if(SIMULA == SPA_ALL) then
      if(HOT==YES) then
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:T, 6:uu, 7:vv, 8:ww, 9: uv, &
                                10: uw, 11: vw, 12:VISt' 
        do i = Nprob, 1, -1
          if(Ncount(i) /= 0) then
            write(3,'(12E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i), Tmp(i),  &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i), vwp(i), VIStp(i)
          end if
        end do 
      else
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8: uv, &
                                9: uw, 10: vw, 11:VISt' 
        do i = Nprob, 1, -1
          if(Ncount(i) /= 0) then
            write(3,'(11E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i),   &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i), vwp(i), VIStp(i)
          end if
        end do 
      end if
    end if

    close(3)

    open(3,FILE=namRes_plus)
    write(3,'(A1,2(A10, F13.5))') '#', 'Utau = ', Ufric, 'Re_tau = ', Ufric/VISc
    if(SIMULA == DNS) then
      write(3,'(A1,2X,A80)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8:uv, 9:uw, 10:vw, 11:Kin' 
      do i = 1, Nprob
        if(Ncount(i) /= 0) then
          write(3,'(11E15.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric,   &
          uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2,      &
          vwp(i)/Ufric**2, 0.5*(uup(i) + vvp(i) + wwp(i))/Ufric**2
        end if
      end do 
    else if(SIMULA == LES.or.SIMULA == DES_SPA) then
      if(HOT == YES) then 
        write(3,'(A1,2X,A110)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:T, 6:uu, 7:vv, 8:ww, &
                               9:uv, 10:uw, 11:vw, 12: t, 13: uT, 14: vT, 15: wT, 16:Kin, 17: VISt' 
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(3,'(17E15.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric, Tmp(i), &
            uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2,      &
            vwp(i)/Ufric**2, TTp(i), uTp(i), vTp(i), wTp(i), &
            0.5*(uup(i) + vvp(i) + wwp(i))/Ufric**2, VIStp(i)/VISc
          end if
        end do 
      else
        write(3,'(A1,2X,A88)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8:uv, 9:uw, 10:vw, 11:Kin, 12: VISt' 
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(3,'(12E15.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric,   &
            uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2,      &
            vwp(i)/Ufric**2, 0.5*(uup(i) + vvp(i) + wwp(i))/Ufric**2, VIStp(i)/VISc
          end if
        end do 
      end if 
    else if(SIMULA == HYB_ZETA) then
      if(HOT==YES) then
        write(3,'(A1,2X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5: T, 6:Kin, 7:Eps, 8:v2, 9:f22, &
                                     10:uu, 11:vv, 12:ww, 13:uv, 14:uw, 15:vw, 16:t, 17:uT, 18:vT, 19:wT, 20:Kin_res, 21:VISt' 
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(3,'(21E15.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric, Tmp(i), &
            Kinp(i)/Ufric, Epsp(i)*VISc/Ufric**4.0, v2p(i), fp(i)*VISc/Ufric**2.0, &
            uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2,      &
            vwp(i)/Ufric**2, TTp(i), uTp(i), vTp(i), wTp(i), &
            0.5*(uup(i) + vvp(i) + wwp(i))/Ufric**2, VIStp(i)/VISc
          end if
        end do 
      else 
        write(3,'(A1,2X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:Kin, 6:Eps, 7:v2, 8:f22, &
                                     9:uu, 10:vv, 11:ww, 12:uv, 13:uw, 14:vw, 15:Kin_res, 16: VISt' 
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(3,'(16E15.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric,   &
            Kinp(i)/Ufric**2, Epsp(i)*VISc/Ufric**4.0, v2p(i), fp(i)*VISc/Ufric**2.0, &
            uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2,      &
            vwp(i)/Ufric**2, 0.5*(uup(i) + vvp(i) + wwp(i))/Ufric**2, VIStp(i)/VISc
          end if
        end do 
      end if
    else if(SIMULA == HYB_PITM) then
      if(HOT==YES) then
        write(3,'(A1,2X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5: T, 6:Kin, 7:Eps,  &
                                     10:uu, 11:vv, 12:ww, 13:uv, 14:uw, 15:vw, 16:t, 17:uT, 18:vT, 19:wT, 20:Kin_res, 21:VISt' 
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(3,'(21E15.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric, Tmp(i), &
            Kinp(i)/Ufric**2, Epsp(i)*VISc/Ufric**4.0,      &
            uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2,      &
            vwp(i)/Ufric**2, TTp(i), uTp(i), vTp(i), wTp(i), &
            0.5*(uup(i) + vvp(i) + wwp(i))/Ufric**2, VIStp(i)/VISc
          end if
        end do 
      else 
        write(3,'(A1,2X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:Kin, 6:Eps, &
                                     9:uu, 10:vv, 11:ww, 12:uv, 13:uw, 14:vw, 15:Kin_res, 16: VISt' 
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(3,'(16E15.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric,   &
            Kinp(i)/Ufric**2, Epsp(i)*VISc/Ufric**4.0,  &
            uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2,      &
            vwp(i)/Ufric**2, 0.5*(uup(i) + vvp(i) + wwp(i))/Ufric**2, VIStp(i)/VISc
          end if
        end do 
      end if
    else if(SIMULA == ZETA.or.SIMULA==EBM) then
      if(HOT==YES) then
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8: uv, &
                                9: uw, 10: vw, 11:Kin, 12:Eps, 13:zeta, 14:f' 
        do i = Nprob, 1, -1 
          if(Ncount(i) /= 0) then
            write(3,'(15E15.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric, (Twall - Tmp(i))*Ufric/0.1,  &
            uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2, vwp(i)/Ufric**2, &
            Kinp(i)/Ufric**2, Epsp(i)*VISc/Ufric**4.0, &
            v2p(i), fp(i)
          end if
        end do 
      else
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8: uv, &
                                9: uw, 10: vw, 11:Kin, 12:Eps, 13:zeta, 14:f' 
        do i = Nprob, 1, -1 
          if(Ncount(i) /= 0) then
            write(3,'(14E15.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric,   &
            uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2, vwp(i)/Ufric**2,&
            Kinp(i)/Ufric**2, Epsp(i)*VISc/Ufric**4.0, &
            v2p(i), fp(i)
          end if
        end do 
      end if
    else if(SIMULA == K_EPS_VV) then
      if(HOT==YES) then
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8: uv, &
                                9: uw, 10: vw, 11:Kin, 12:Eps, 13:zeta, 14:f' 
        do i = Nprob, 1, -1
          if(Ncount(i) /= 0) then
            write(3,'(15E15.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric, Tmp(i),  &
            uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2, vwp(i)/Ufric**2, &
            Kinp(i)/Ufric**2, Epsp(i)*VISc/Ufric**4.0, &
            v2p(i)/Ufric**2, fp(i)*VISc/Ufric**2.0
          end if
        end do 
      else
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8: uv, &
                                9: uw, 10: vw, 11:Kin, 12:Eps, 13:zeta, 14:f' 
        do i = Nprob, 1, -1
          if(Ncount(i) /= 0) then
            write(3,'(14E15.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric,   &
            uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2, &
            vwp(i)/Ufric**2, Kinp(i)/Ufric**2, Epsp(i)*VISc/Ufric**4.0, &
            v2p(i)/Ufric**2, fp(i)*VISc/Ufric**2.0
          end if
        end do 
      end if
    else if(SIMULA == K_EPS) then
      if(HOT==YES) then 
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8: uv, &
                                9: uw, 10: vw, 11:Kin, 12:Eps, 13: VISt' 
        do i = Nprob, 1, -1
          if(Ncount(i) /= 0) then
            write(3,'(14E15.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric, Tmp(i),  &
            uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2, &
            vwp(i)/Ufric**2, Kinp(i)/Ufric**2, Epsp(i)*VISc/Ufric**4.0, VIStp(i)/VISc
          end if
        end do 
      else
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8: uv, &
                                9: uw, 10: vw, 11:Kin, 12:Eps, 13: VISt' 
        do i = Nprob, 1, -1
          if(Ncount(i) /= 0) then
            write(3,'(13E15.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric, &
            uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2, &
            vwp(i)/Ufric**2, Kinp(i)/Ufric**2, Epsp(i)*VISc/Ufric**4.0, VIStp(i)/VISc
          end if
        end do 
      end if
    else if(SIMULA == SPA_ALL) then
      if(HOT==YES) then
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:T, 6:uu, 7:vv, 8:ww, 9: uv, &
                                10: uw, 11: vw, 12:VISt' 
        do i = Nprob, 1, -1
          if(Ncount(i) /= 0) then
            write(3,'(12E15.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric, Tmp(i),  &
            uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2, vwp(i)/Ufric**2, VIStp(i)/VISc
          end if
        end do 
      else
        write(3,'(A1,1X,A120)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8: uv, &
                                9: uw, 10: vw, 11:VISt' 
        do i = Nprob, 1, -1
          if(Ncount(i) /= 0) then
            write(3,'(11E15.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric,   &
            uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2, vwp(i)/Ufric**2, VIStp(i)/VISc
          end if
        end do 
      end if
    end if

    close(3)

    deallocate(Np)
    deallocate(z_p)
    deallocate(Ump)
    deallocate(Vmp)
    deallocate(Wmp)
    deallocate(uup)
    deallocate(vvp)
    deallocate(wwp)
    deallocate(uvp)
    deallocate(uwp)
    deallocate(vwp)
    deallocate(Kinp)
    deallocate(Epsp)
    deallocate(v2p)
    deallocate(fp)
    deallocate(VIStp)
    deallocate(Ksgsp)
    if(HOT==YES) then
      deallocate(Tmp)
      deallocate(TTp)
      deallocate(uTp)
      deallocate(vTp)
      deallocate(wTp)
    end if
  END SUBROUTINE UserCutLines_channel
