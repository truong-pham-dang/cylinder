!======================================================================!
  SUBROUTINE UserProbe1D_pipe() 
!----------------------------------------------------------------------!
! Reads the "pipe_rad.1D" file created by yourself and averages the    !
! results in the homogeneous azimutal direction and axial direction.   !
! The values of Umean, Vmean, Wmean, uu, vv, ww, uv, uw and vw are     !
! writen  into file "pipe_mean.dat" and "pipe_mean_plus.dat".          !
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
    REAL :: Rad_2, Ufric 
!------------------------------[Calling]-------------------------------!
    INTERFACE
      LOGICAL FUNCTION Approx(A,B,tol)
        REAL           :: A,B
        REAL, OPTIONAL :: tol
      END FUNCTION Approx
    END INTERFACE 
!-------------------------------[Locals]-------------------------------!
    INTEGER             :: Nprob, pl, c, i, count
    CHARACTER           :: namCoo*80, namPro*80, answer*80
    REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:),               & 
                                 var_1(:), var_2(:), var_3(:), Rad_1(:), Rad_mp(:)
    INTEGER,ALLOCATABLE :: Np(:), Ncount(:)
    REAL                :: R, Urad_mean, Utan_mean, dummy
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $  
!======================================================================!

    namPro = name

!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
    namCoo = name
    namCoo(len_trim(name)+1:len_trim(name)+3) = '.1D'
    write(6, *) '# Now reading the file:', namCoo
    open(9, FILE=namCoo)

!---- write the number of probes 
    read(9,*) Nprob
    allocate(z_p(Nprob*2))

!---- write the probe coordinates out
    do pl=1,Nprob
      read(9,*) dummy, z_p(pl) 
    end do
    close(9)
    call CalcShear(U % n, V % n, W % n, Shear)

    allocate(Np(Nprob));    Np=0 
    allocate(Ump(Nprob));   Ump=0.0
    allocate(Vmp(Nprob));   Vmp=0.0
    allocate(Wmp(Nprob));   Wmp=0.0
    allocate(uup(Nprob));   uup=0.0
    allocate(vvp(Nprob));   vvp=0.0
    allocate(wwp(Nprob));   wwp=0.0
    allocate(uvp(Nprob));   uvp=0.0
    allocate(uwp(Nprob));   uwp=0.0
    allocate(vwp(Nprob));   vwp=0.0
    allocate(Ksgsp(Nprob)); Ksgsp=0.0
    allocate(Rad_mp(Nprob));  Rad_mp=0.0
    allocate(var_1(Nprob));  var_1=0.0

    allocate(Rad_1(Nprob));  Rad_1=0.0
    allocate(Ncount(Nprob)); Ncount=0
    count = 0

    if(HOT==YES) then
      allocate(Tmp(Nprob));   Tmp=0.0
      allocate(TTp(Nprob));   TTp=0.0
      allocate(uTp(Nprob));   uTp=0.0
      allocate(vTp(Nprob));   vTp=0.0
      allocate(wTp(Nprob));   wTp=0.0
    end if  

!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
    do i = 1, Nprob
      Rad_1(i) = 1.0 - abs(z_p(i)) 
    end do 

    do i = 2, Nprob
      do c=1, NC
        Rad_2 = WallDs(c)   
        if(Rad_2 < Rad_1(i-1) .and. Rad_2 > Rad_1(i)) then
          R           = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 
          Urad_mean   = (U % n(c) * xc(c) / R  + V % n(c) * yc(c) / R)
          Utan_mean   = (-U % n(c) * yc(c) / R  + V % n(c) * xc(c) / R) 
 
          Ump(i)   = Ump(i) + U % n(c)
          Vmp(i)   = Vmp(i) + V % n(c)
          Wmp(i)   = Wmp(i) + W % n(c)
          if(SIMULA==K_EPS_VV) then
            uup(i)   = uup(i) + Kin % n(c) 
            vvp(i)   = vvp(i) + Eps % n(c) 
            wwp(i)   = wwp(i) + v_2 % n(c) 
            uvp(i)   = uvp(i) + f22 % n(c) 
            uwp(i)   = uwp(i) + 0.0 
          else if(SIMULA==ZETA) then
            uup(i)   = uup(i) + Kin % n(c) 
            vvp(i)   = vvp(i) + Eps % n(c) 
            wwp(i)   = wwp(i) + v_2 % n(c) 
            uvp(i)   = uvp(i) + f22 % n(c) 
            uwp(i)   = uwp(i) + 0.0 
          else if(SIMULA==K_EPS) then
            uup(i)   = uup(i) + Kin % n(c) 
            vvp(i)   = vvp(i) + Eps % n(c) 
            wwp(i)   = wwp(i) + 0.0
            uvp(i)   = uvp(i) + 0.0 
            uwp(i)   = uwp(i) + 0.0 
          else if(SIMULA==LES.or.SIMULA==DES_SPA) then
            uup(i)   = uup(i) + (uu % mean(c)- Urad_mean * Urad_mean)
            vvp(i)   = vvp(i) + (vv % mean(c)- Utan_mean * Utan_mean)
            wwp(i)   = wwp(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
            uvp(i)   = uvp(i) + (uv % mean(c)- Urad_mean * Utan_mean )
            uwp(i)   = uwp(i) + (uw % mean(c)- Urad_mean * W % mean(c))
            vwp(i)   = uwp(i) + (vw % mean(c)- Utan_mean * W % mean(c))
          end if 
 
          if(HOT==YES) then
            Tmp(i)   = Tmp(i) + T % mean(c)
            if(SIMULA == LES.or.SIMULA == DES_SPA) then
              TTp(i)   = TTp(i) + (TT % mean(c) - T % mean(c) * T % mean(c))
              uTp(i)   = uTp(i) + (uT % mean(c) - u % mean(c) * T % mean(c))
              vTp(i)   = vTp(i) + (vT % mean(c) - v % mean(c) * T % mean(c))
              wTp(i)   = wTp(i) + (wT % mean(c) - w % mean(c) * T % mean(c))
            end if
          end if
       
          Rad_mp(i) = Rad_mp(i) + WallDs(c) !(xc(c)*xc(c) + yc(c)*yc(c))**0.5
          Ncount(i) = Ncount(i) + 1
        end if
      end do 
    end do 

!---- average over all processors
    do pl=2, Nprob
      call IGlSum(Ncount(pl))

      call GloSum(Ump(pl))
      call GloSum(Vmp(pl))
      call GloSum(Wmp(pl))
      call GloSum(Rad_mp(pl))

      call GloSum(uup(pl))
      call GloSum(vvp(pl))
      call GloSum(wwp(pl))

      call GloSum(uvp(pl))
      call GloSum(uwp(pl))
      call GloSum(vwp(pl))
      call GloSum(var_1(pl))

      count =  count + Ncount(pl) 

      if(HOT==YES) then
        call GloSum(Tmp(pl))
        call GloSum(TTp(pl))
        call GloSum(uTp(pl))
        call GloSum(vTp(pl))
        call GloSum(wTp(pl))
      end if
    end do

    open(3,FILE='pipe_mean.dat')

    do i = 2, Nprob
      if(Ncount(i) /= 0) then
        Ump(i)    =  Ump(i)/Ncount(i)
        Vmp(i)    =  Vmp(i)/Ncount(i)
        Wmp(i)    =  Wmp(i)/Ncount(i)
        Ump(i)    =  Ump(i)/Ncount(i)
        Vmp(i)    =  Vmp(i)/Ncount(i)
        uup(i)    =  uup(i)/Ncount(i)
        vvp(i)    =  vvp(i)/Ncount(i)
        wwp(i)    =  wwp(i)/Ncount(i)
        uvp(i)    =  uvp(i)/Ncount(i)
        uwp(i)    =  uwp(i)/Ncount(i)
        var_1(i)  =  var_1(i)/Ncount(i)
        Rad_mp(i) =  Rad_mp(i)/Ncount(i)
      end if
    end do 

    Ufric = (VISc*abs(Wmp(Nprob))/(Rad_mp(Nprob)))**0.5 

10  continue
 
    open(3,FILE='pipe_mean.dat')
    write(3,'(A1,2(A10, F13.5))') '#', 'Utau = ', Ufric, 'Re_tau = ', Ufric/VISc
    if(SIMULA == LES.or.SIMULA == DES_SPA) then
      write(3,'(A1,2X,A80)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8:uv, 9:uw, 10:vw, 11:Kin' 
      do i = Nprob, 2, -1
        if(Ncount(i) /= 0) then
          write(3,'(11E15.7)') Rad_mp(i), abs(Ump(i)), abs(Vmp(i)), abs(Wmp(i)),   &
          uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
          vwp(i), (uup(i) + vvp(i) + wwp(i))
        end if
      end do 
    else if(SIMULA == ZETA) then
      write(3,'(A1,1X,A60)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:Kin, 6:Eps, 7:zeta, 8:f' 
      do i = Nprob, 2, -1
        if(Ncount(i) /= 0) then
          write(3,'(9E15.7)') Rad_mp(i), abs(Ump(i)), abs(Vmp(i)), abs(Wmp(i)),   &
          uup(i), vvp(i), wwp(i), uvp(i), uwp(i)
        end if
      end do 
    end if
    close(3)

    open(3,FILE='pipe_mean_plus.dat')
    write(3,'(A1,2(A10, F13.5))') '#', 'Utau = ', Ufric, 'Re_tau = ', Ufric/VISc
    if(SIMULA == LES.or.SIMULA == DES_SPA) then
      write(3,'(A1,2X,A80)') '#', '1:Xrad+, 2:U+, 3:V+, 4:W+, 5:uu+, 6:vv+, 7:ww+, 8:uv+, 9:uw+, 10:vw+, 11:Kin+' 
      do i = Nprob, 2, -1
        if(Ncount(i) /= 0) then
          write(3,'(11E15.7)') Rad_mp(i)*Ufric/VISc, abs(Ump(i))/Ufric, abs(Vmp(i))/Ufric, abs(Wmp(i))/Ufric,   &
          uup(i)/Ufric**2.0, vvp(i)/Ufric**2.0, wwp(i)/Ufric**2.0, uvp(i)/Ufric**2.0, uwp(i)/Ufric**2.0,      &
          vwp(i)/Ufric**2.0, (uup(i) + vvp(i) + wwp(i))/Ufric**2.0
        end if
      end do 
    else if(SIMULA == ZETA) then
      write(3,'(A1,1X,A60)') '#', '1:Xrad+, 2:U+, 3:V+, 4:W+, 5:Kin+, 6:Eps+, 7:zeta+, 8:f+' 
      do i = Nprob, 2, -1
        if(Ncount(i) /= 0) then
          write(3,'(9E15.7)') (Rad_mp(i))*Ufric/VISc, abs(Ump(i))/Ufric, abs(Vmp(i))/Ufric, abs(Wmp(i))/Ufric,   &
          uup(i)/Ufric**2.0, vvp(i)*VISc/Ufric**4.0, wwp(i), uvp(i)*VISc/Ufric**2.0, uwp(i)/Ufric**2.0
        end if
      end do 
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
    deallocate(var_1)
    deallocate(Ksgsp)
    if(HOT==YES) then
      deallocate(Tmp)
      deallocate(TTp)
      deallocate(uTp)
      deallocate(vTp)
      deallocate(wTp)
    end if
  END SUBROUTINE UserProbe1D_pipe
