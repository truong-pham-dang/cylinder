!======================================================================!
  SUBROUTINE UserProbe1D_pipe_LES() 
!----------------------------------------------------------------------!
! Reads the ".1D" file created by the "Generator" and averages the     !
! results in the planes defined by coordinates in it. Then averages    !
! the values of Umean, Vmean, Wmean, uu, vv, ww, uv, uw and vw and     !
! writes them into file ".1Dr".                                        !
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
  INTEGER             :: Nprob, pl, c, i, count, Ncount_wall
  CHARACTER           :: namCoo*80, namPro*80, answer*80
  REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                         uup(:), vvp(:), wwp(:), &
                         uvp(:), uwp(:), vwp(:), &
                         Tmp(:), TTp(:),         &
                         uTp(:), vTp(:), wTp(:), &
                         Ksgsp(:),               & 
                         var_1(:), var_2(:), var_3(:), Rad_1(:), Rad_mp(:), &
                         var_10(:), var_11(:), var_12(:), var_13(:)
  
  INTEGER,ALLOCATABLE :: Np(:), Ncount(:)
  REAL                :: R, Urad_mean, Utan_mean, NF, vol1
  REAL                :: Tfric, Twall, Twall_p, FFF, dummy, tan1
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $  
!======================================================================!
!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
  write(6, *) '# Now reading the file: pipe_rad.1D ' 
  open(9, FILE='pipe_rad.1D')

!---- write the number of probes 
  read(9,*) Nprob
  allocate(z_p(Nprob))

!---- write the probe coordinates out
  do pl=1,Nprob
    read(9,*) dummy, z_p(pl) 
  end do
  close(9)

  call CalcShear(U % n, V % n, W % n, Shear)

  allocate(Np(Nprob));    Np=0 
  allocate(Ncount(Nprob));    Ncount=0 
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
  allocate(var_10(Nprob));  var_10=0.0
  allocate(var_11(Nprob));  var_11=0.0
  allocate(var_12(Nprob));  var_12=0.0
  allocate(var_13(Nprob));  var_13=0.0
  allocate(Rad_1(Nprob));  Rad_1=0.0

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
      Rad_1(i) = abs(z_p(i)) 
    end do 
    Ncount_wall = 0

    if(HOT==YES) then
    Twall_p     = 0.0
    do c = -Nbc,1
      if(TypeBC(c) == WALL) then  
        Twall_p = Twall_p + T % mean(c)
        Ncount_wall = Ncount_wall + 1
      end if
    end do
    call GloSum(Twall_p)
    call IGlSum(Ncount_wall)

    Twall =  Twall_p / Ncount_wall
    end if

    do i = 1, Nprob-1
      do c=1,NC
        Rad_2 = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 + tiny
        tan1 = yc(c)/xc(c)
        if(Rad_2 < Rad_1(i) .and. Rad_2 >= Rad_1(i+1)) then
          if(tan1 > 0.36397.and.tan1 < 2.7474) then
          R           = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 + tiny
          Urad_mean   = (U % mean(c) * xc(c) / R  + V % mean(c) * yc(c) / R)
          Utan_mean   = (-U % mean(c) * yc(c) / R  + V % mean(c) * xc(c) / R) 
 
          Ump(i)   = Ump(i) + U % mean(c)
          Vmp(i)   = Vmp(i) + V % mean(c)
          Wmp(i)   = Wmp(i) + W % mean(c)
          uup(i)   = uup(i) + (uu % mean(c)- Urad_mean * Urad_mean)
          vvp(i)   = vvp(i) + (vv % mean(c)- Utan_mean * Utan_mean)
          wwp(i)   = wwp(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
          uvp(i)   = uvp(i) + (uv % mean(c)- Urad_mean * Utan_mean )
          uwp(i)   = uwp(i) + (uw % mean(c)- Urad_mean * W % mean(c))


          if(HOT==YES) then
            Tmp(i)   = Tmp(i) + T % mean(c)
            TTp(i)   = TTp(i) + (TT % mean(c) - T % mean(c) * T % mean(c))
            uTp(i)   = uTp(i) + (uT % mean(c) - Urad_mean * T % mean(c))
            vTp(i)   = vTp(i) + (vT % mean(c) - Utan_mean * T % mean(c))
            wTp(i)   = wTp(i) + (wT % mean(c) - w % mean(c) * T % mean(c))
          end if
       
          Rad_mp(i) = Rad_mp(i) + (xc(c)*xc(c) + yc(c)*yc(c))**0.5
          Ncount(i) = Ncount(i) + 1
        end if
        end if
      end do 
    end do 
!---- average over all processors
  do pl=1, Nprob
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

    call GloSum(var_10(pl))
    call GloSum(var_11(pl))
    call GloSum(var_12(pl))
    call GloSum(var_13(pl))

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
    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        Wmp(i)    =  Wmp(i)/Ncount(i)
        Ump(i)    =  Ump(i)/Ncount(i)
        Vmp(i)    =  Vmp(i)/Ncount(i)
        uup(i)    =  uup(i)/Ncount(i)
        vvp(i)    =  vvp(i)/Ncount(i)
        wwp(i)    =  wwp(i)/Ncount(i)
        uvp(i)    =  uvp(i)/Ncount(i)
        uwp(i)    =  uwp(i)/Ncount(i)
        var_1(i)  =  var_1(i)/Ncount(i)
        var_10(i)  =  var_10(i)/Ncount(i)
        var_11(i)  =  var_11(i)/Ncount(i)
        var_12(i)  =  var_12(i)/Ncount(i)
        var_13(i)  =  var_13(i)/Ncount(i)
        Rad_mp(i)  =  Rad_mp(i)/Ncount(i)

        write(3,'(7E15.7,I8)') 1.0 - Rad_mp(i), Wmp(i) , uup(i), vvp(i), wwp(i), uvp(i), uwp(i),&
                                Ncount(i)
      end if
    end do 
    close(3)

    write(*,*) Wmp(1), 1.0-Rad_mp(1) 
      Ufric = (VISc*abs(Wmp(1))/(1.0-Rad_mp(1)))**0.5 
      if(this < 2)write(*,*) 'mean Utau = ', Ufric,         &
       'ukupno nasao celija: ', count

    open(3,FILE='pipe_mean_plus.dat')
    write(3,*)'# Utau = ', Ufric, 'Re_tau = ', Ufric/VISc
    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        write(3,'(6E15.7)') (1.0 - Rad_mp(i))*Ufric/VISc, abs(Wmp(i))/Ufric,  &
         (abs(uup(i))/Ufric**2.0)**0.5, (abs(vvp(i))/Ufric**2.0)**0.5,  &
         (abs(wwp(i))/Ufric**2.0)**0.5, uwp(i)/Ufric**2.0
      end if
    end do 
    close(3)

  deallocate(Np)
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
  deallocate(var_10)
  deallocate(var_11)
  deallocate(var_12)
  deallocate(var_13)
  deallocate(Ksgsp)
  deallocate(Rad_1)

  if(HOT==YES) then
    deallocate(Tmp)
    deallocate(TTp)
    deallocate(uTp)
    deallocate(vTp)
    deallocate(wTp)
  end if
  END SUBROUTINE UserProbe1D_pipe_LES
