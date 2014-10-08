!======================================================================!
  SUBROUTINE UserProbe1D_Nusselt_jet()
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
  INTEGER             :: Nprob, pl, c, dummy, i, count, c1, c2, s
  CHARACTER           :: namCoo*80, namPro*80, answer*80
  REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 Ump1(:), Vmp1(:), Wmp1(:), &
                                 Ump2(:), Vmp2(:), Wmp2(:), &
                                 Ump3(:), Vmp3(:), Wmp3(:), &
                                 Ump4(:), Vmp4(:), Wmp4(:), &
                                 Ump5(:), Vmp5(:), Wmp5(:), &
                                 Ump6(:), Vmp6(:), Wmp6(:), &
                                 Ump7(:), Vmp7(:), Wmp7(:), &
                                 uup1(:), vvp1(:), wwp1(:), &
                                 uup2(:), vvp2(:), wwp2(:), &
                                 uup3(:), vvp3(:), wwp3(:), &
                                 uup4(:), vvp4(:), wwp4(:), &
                                 uup5(:), vvp5(:), wwp5(:), &
                                 uup6(:), vvp6(:), wwp6(:), &
                                 uup7(:), vvp7(:), wwp7(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 uwp1(:), uwp2(:), uwp3(:), &
                                 uwp4(:), uwp5(:), uwp6(:), &
                                 uwp7(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:),               & 
                                 var_1(:), var_2(:), var_3(:), var_4(:), Rad_1(:), Rad_mp(:)
                                  
  INTEGER,ALLOCATABLE :: Np(:), Ncount(:), Ncount2(:)
  REAL                :: R, Urad_mean, Utan_mean
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $  
!======================================================================!


!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
  write(6, *) '# Now reading the file: pipe_rad_coordinate.1D ' 
  open(9, FILE='pipe_rad_coordinate.1D')

!---- write the number of probes 
  read(9,*) Nprob
  allocate(z_p(Nprob))

!---- write the probe coordinates out
  do pl=1,Nprob
    read(9,*) dummy, z_p(pl) 
  end do
  close(9)
  if(this < 2) write(*,*)'finished reading'
  allocate(Np(Nprob));     Np=0 
  allocate(Ump(Nprob));    Ump=0.0
  allocate(Vmp(Nprob));    Vmp=0.0
  allocate(Wmp(Nprob));    Wmp=0.0
  allocate(Ump1(Nprob));   Ump1=0.0
  allocate(Vmp1(Nprob));   Vmp1=0.0
  allocate(Wmp1(Nprob));   Wmp1=0.0
  allocate(Ump2(Nprob));   Ump2=0.0
  allocate(Vmp2(Nprob));   Vmp2=0.0
  allocate(Wmp2(Nprob));   Wmp2=0.0
  allocate(Ump3(Nprob));   Ump3=0.0
  allocate(Vmp3(Nprob));   Vmp3=0.0
  allocate(Wmp3(Nprob));   Wmp3=0.0
  allocate(Ump4(Nprob));   Ump4=0.0
  allocate(Vmp4(Nprob));   Vmp4=0.0
  allocate(Wmp4(Nprob));   Wmp4=0.0
  allocate(Ump5(Nprob));   Ump5=0.0
  allocate(Vmp5(Nprob));   Vmp5=0.0
  allocate(Wmp5(Nprob));   Wmp5=0.0
  allocate(Ump6(Nprob));   Ump6=0.0
  allocate(Vmp6(Nprob));   Vmp6=0.0
  allocate(Wmp6(Nprob));   Wmp6=0.0
  allocate(Ump7(Nprob));   Ump7=0.0
  allocate(Vmp7(Nprob));   Vmp7=0.0
  allocate(Wmp7(Nprob));   Wmp7=0.0
  allocate(uup1(Nprob));   uup1=0.0
  allocate(vvp1(Nprob));   vvp1=0.0
  allocate(wwp1(Nprob));   wwp1=0.0
  allocate(uup2(Nprob));   uup2=0.0
  allocate(vvp2(Nprob));   vvp2=0.0
  allocate(wwp2(Nprob));   wwp2=0.0
  allocate(uup3(Nprob));   uup3=0.0
  allocate(vvp3(Nprob));   vvp3=0.0
  allocate(wwp3(Nprob));   wwp3=0.0
  allocate(uup4(Nprob));   uup4=0.0
  allocate(vvp4(Nprob));   vvp4=0.0
  allocate(wwp4(Nprob));   wwp4=0.0
  allocate(uup5(Nprob));   uup5=0.0
  allocate(vvp5(Nprob));   vvp5=0.0
  allocate(wwp5(Nprob));   wwp5=0.0
  allocate(uup6(Nprob));   uup6=0.0
  allocate(vvp6(Nprob));   vvp6=0.0
  allocate(wwp6(Nprob));   wwp6=0.0
  allocate(uup7(Nprob));   uup7=0.0
  allocate(vvp7(Nprob));   vvp7=0.0
  allocate(wwp7(Nprob));   wwp7=0.0
  allocate(uwp1(Nprob));   uwp1=0.0
  allocate(uwp2(Nprob));   uwp2=0.0
  allocate(uwp3(Nprob));   uwp3=0.0
  allocate(uwp4(Nprob));   uwp4=0.0
  allocate(uwp5(Nprob));   uwp5=0.0
  allocate(uwp6(Nprob));   uwp6=0.0
  allocate(uwp7(Nprob));   uwp7=0.0
  allocate(uvp(Nprob));    uvp=0.0
  allocate(uwp(Nprob));    uwp=0.0
  allocate(vwp(Nprob));    vwp=0.0
  allocate(Rad_mp(Nprob));  Rad_mp=0.0
  allocate(var_1(Nprob));  var_1=0.0
  allocate(var_2(Nprob));  var_2=0.0
  allocate(var_3(Nprob));  var_3=0.0


  allocate(Rad_1(Nprob));  Rad_1=0.0
  allocate(Ncount(Nprob)); Ncount=0
  allocate(Ncount2(Nprob)); Ncount2=0
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

    do i = 1, Nprob-1
      do s=1,NS
        c1=SideC(1,s)
        c2=SideC(2,s)
        if(c2 < 0) then
          if(TypeBC(c2).eq.WALLFL) then
            Rad_2 = (xc(c1)*xc(c1) + yc(c1)*yc(c1))**0.5 + tiny
            if(Rad_2 < Rad_1(i) .and. Rad_2 > Rad_1(i+1).and.zc(c1) < 0.1) then
              R           = (xc(c1)*xc(c1) + yc(c1)*yc(c1))**0.5 + tiny 
              Ump(i)   = Ump(i) + (U % mean(c1) * U % mean(c1)            &
              + V % mean(c1) * V % mean(c1) + W % mean(c1)*W % mean(c1))**0.5 
              Wmp(i)   = Wmp(i) + W % mean(c1)
              var_1(i) = var_1(i) + zc(c1)
              var_3(i) = var_3(i) +  (U % n(c1) * U % n(c1) + V % n(c1) * V % n(c1) + W % n(c1)*W % n(c1))**0.5
              Rad_mp(i)= Rad_mp(i) + (xc(c1)*xc(c1) + yc(c1)*yc(c1))**0.5
              Tmp(i)   = Tmp(i) + T % mean(c2)
              var_2(i) = var_2(i) + T % n(c2)
              Ncount(i)= Ncount(i) + 1
            end if
          end if
        end if
      end do

      do c = 1, NC
        Rad_2 = (xc(c)*xc(c) + yc(c)*yc(c))**0.5  + tiny
        if(Rad_2 < Rad_1(i) .and. Rad_2 > Rad_1(i+1)) then
            R           = (xc(c)*xc(c) + yc(c)*yc(c))**0.5
            Urad_mean   = (U % mean(c) * xc(c) / R  + V % mean(c) * yc(c) / R)
            Utan_mean   = (-U % mean(c) * yc(c) / R  + V % mean(c) * xc(c) / R)

          if(zc(c) < 0.0284.and.zc(c) > 0.0236)then
            Ump1(i)   = Ump1(i) + (U % mean(c) * U % mean(c)            &
                       + V % mean(c) * V % mean(c) + W % mean(c)*W % mean(c))**0.5 
            Vmp1(i)   = Vmp1(i) + Urad_mean 
            Wmp1(i)   = Wmp1(i) + W % mean(c)
            uup1(i)   = uup1(i) + (uu % mean(c)- Urad_mean * Urad_mean)
            vvp1(i)   = vvp1(i) + (vv % mean(c)- Utan_mean * Utan_mean)
            wwp1(i)   = wwp1(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
            uwp1(i)   = uwp1(i) + (uw % mean(c)- Urad_mean * W % mean(c))
            Ncount2(i)= Ncount2(i) + 1
          end if
          if(zc(c) < 0.109.and.zc(c) > 0.0956)then
            Ump2(i)   = Ump2(i) + (U % mean(c) * U % mean(c)            &
                       + V % mean(c) * V % mean(c) + W % mean(c)*W % mean(c))**0.5 
            Vmp2(i)   = Vmp2(i) + Urad_mean 
            Wmp2(i)   = Wmp2(i) + W % mean(c)
            uup2(i)   = uup2(i) + (uu % mean(c)- Urad_mean * Urad_mean)
            vvp2(i)   = vvp2(i) + (vv % mean(c)- Utan_mean * Utan_mean)
            wwp2(i)   = wwp2(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
            uwp2(i)   = uwp2(i) + (uw % mean(c)- Urad_mean * W % mean(c))
          end if
          if(zc(c) < 0.5517.and.zc(c) > 0.4893)then
            Ump3(i)   = Ump3(i) + (U % mean(c) * U % mean(c)            &
                       + V % mean(c) * V % mean(c) + W % mean(c)*W % mean(c))**0.5 
            Vmp3(i)   = Vmp3(i) + Urad_mean 
            Wmp3(i)   = Wmp3(i) + W % mean(c)
            uup3(i)   = uup3(i) + (uu % mean(c)- Urad_mean * Urad_mean)
            vvp3(i)   = vvp3(i) + (vv % mean(c)- Utan_mean * Utan_mean)
            wwp3(i)   = wwp3(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
            uwp3(i)   = uwp3(i) + (uw % mean(c)- Urad_mean * W % mean(c))
          end if
          if(zc(c) < 1.0.and.zc(c) > 0.8883)then
            Ump4(i)   = Ump4(i) + (U % mean(c) * U % mean(c)            &
                      + V % mean(c) * V % mean(c) + W % mean(c)*W % mean(c))**0.5 
            Vmp4(i)   = Vmp4(i) + Urad_mean 
            Wmp4(i)   = Wmp4(i) + W % mean(c)
            uup4(i)   = uup4(i) + (uu % mean(c)- Urad_mean * Urad_mean)
            vvp4(i)   = vvp4(i) + (vv % mean(c)- Utan_mean * Utan_mean)
            wwp4(i)   = wwp4(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
            uwp4(i)   = uwp4(i) + (uw % mean(c)- Urad_mean * W % mean(c))
          end if
          if(zc(c) < 2.0.and.zc(c) > 0.1875000E+01)then
            Ump5(i)   = Ump5(i) + (U % mean(c) * U % mean(c)            &
                       + V % mean(c) * V % mean(c) + W % mean(c)*W % mean(c))**0.5 
            Vmp5(i)   = Vmp5(i) + Urad_mean 
            Wmp5(i)   = Wmp5(i) + W % mean(c)
            uup5(i)   = uup5(i) + (uu % mean(c)- Urad_mean * Urad_mean)
            vvp5(i)   = vvp5(i) + (vv % mean(c)- Utan_mean * Utan_mean)
            wwp5(i)   = wwp5(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
            uwp5(i)   = uwp5(i) + (uw % mean(c)- Urad_mean * W % mean(c))
          end if
          if(zc(c) < 3.0.and.zc(c) > 2.875)then
            Ump6(i)   = Ump6(i) + (U % mean(c) * U % mean(c)            &
                       + V % mean(c) * V % mean(c) + W % mean(c)*W % mean(c))**0.5 
            Vmp6(i)   = Vmp6(i) + Urad_mean 
            Wmp6(i)   = Wmp6(i) + W % mean(c)
            uup6(i)   = uup6(i) + (uu % mean(c)- Urad_mean * Urad_mean)
            vvp6(i)   = vvp6(i) + (vv % mean(c)- Utan_mean * Utan_mean)
            wwp6(i)   = wwp6(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
            uwp6(i)   = uwp6(i) + (uw % mean(c)- Urad_mean * W % mean(c))
          end if
          if(zc(c) < 4.0.and.zc(c) > 3.875)then
            Ump7(i)   = Ump7(i) + (U % mean(c) * U % mean(c)            &
                       + V % mean(c) * V % mean(c) + W % mean(c)*W % mean(c))**0.5 
            Vmp7(i)   = Vmp7(i) + Urad_mean 
            Wmp7(i)   = Wmp7(i) + W % mean(c)
            uup7(i)   = uup7(i) + (uu % mean(c)- Urad_mean * Urad_mean)
            vvp7(i)   = vvp7(i) + (vv % mean(c)- Utan_mean * Utan_mean)
            wwp7(i)   = wwp7(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
            uwp7(i)   = uwp7(i) + (uw % mean(c)- Urad_mean * W % mean(c))
          end if
        end if
      end do 

    end do 
!---- average over all processors
  do pl=1, Nprob
    call IGlSum(Ncount(pl))
    call IGlSum(Ncount2(pl))

    call GloSum(Ump(pl))
    call GloSum(Wmp(pl))

    call GloSum(Ump1(pl))
    call GloSum(Vmp1(pl))
    call GloSum(Wmp1(pl))

    call GloSum(Ump2(pl))
    call GloSum(Vmp2(pl))
    call GloSum(Wmp2(pl))

    call GloSum(Ump3(pl))
    call GloSum(Vmp3(pl))
    call GloSum(Wmp3(pl))

    call GloSum(Ump4(pl))
    call GloSum(Vmp4(pl))
    call GloSum(Wmp4(pl))

    call GloSum(Ump5(pl))
    call GloSum(Vmp5(pl))
    call GloSum(Wmp5(pl))

    call GloSum(Ump6(pl))
    call GloSum(Vmp6(pl))
    call GloSum(Wmp6(pl))

    call GloSum(Ump7(pl))
    call GloSum(Vmp7(pl))
    call GloSum(Wmp7(pl))

    call GloSum(uup1(pl))
    call GloSum(vvp1(pl))
    call GloSum(wwp1(pl))

    call GloSum(uup2(pl))
    call GloSum(vvp2(pl))
    call GloSum(wwp2(pl))

    call GloSum(uup3(pl))
    call GloSum(vvp3(pl))
    call GloSum(wwp3(pl))

    call GloSum(uup4(pl))
    call GloSum(vvp4(pl))
    call GloSum(wwp4(pl))

    call GloSum(uup5(pl))
    call GloSum(vvp5(pl))
    call GloSum(wwp5(pl))

    call GloSum(uup6(pl))
    call GloSum(vvp6(pl))
    call GloSum(wwp6(pl))

    call GloSum(uup7(pl))
    call GloSum(vvp7(pl))
    call GloSum(wwp7(pl))

    call GloSum(uwp1(pl))
    call GloSum(uwp2(pl))
    call GloSum(uwp3(pl))
    call GloSum(uwp4(pl))
    call GloSum(uwp5(pl))
    call GloSum(uwp6(pl))
    call GloSum(uwp7(pl))

    call GloSum(var_1(pl))
    call GloSum(var_2(pl))
    call GloSum(var_3(pl))

    call GloSum(Rad_mp(pl))
    call GloSum(Tmp(pl))

    count =  count + Ncount(pl) 
  end do

    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        Wmp(i)    =  Wmp(i)/Ncount(i)
        Ump(i)    =  Ump(i)/Ncount(i)
        Wmp1(i)    =  Wmp1(i)/Ncount2(i)
        Vmp1(i)    =  Vmp1(i)/Ncount2(i)
        Ump1(i)    =  Ump1(i)/Ncount2(i)
        Wmp2(i)    =  Wmp2(i)/Ncount2(i)
        Vmp2(i)    =  Vmp2(i)/Ncount2(i)
        Ump2(i)    =  Ump2(i)/Ncount2(i)
        Wmp3(i)    =  Wmp3(i)/Ncount2(i)
        Vmp3(i)    =  Vmp3(i)/Ncount2(i)
        Ump3(i)    =  Ump3(i)/Ncount2(i)
        Wmp4(i)    =  Wmp4(i)/Ncount2(i)
        Vmp4(i)    =  Vmp4(i)/Ncount2(i)
        Ump4(i)    =  Ump4(i)/Ncount2(i)
        Wmp5(i)    =  Wmp5(i)/Ncount2(i)
        Vmp5(i)    =  Vmp5(i)/Ncount2(i)
        Ump5(i)    =  Ump5(i)/Ncount2(i)
        Wmp6(i)    =  Wmp6(i)/Ncount2(i)
        Vmp6(i)    =  Vmp6(i)/Ncount2(i)
        Ump6(i)    =  Ump6(i)/Ncount2(i)
        Wmp7(i)    =  Wmp7(i)/Ncount2(i)
        Vmp7(i)    =  Vmp7(i)/Ncount2(i)
        Ump7(i)    =  Ump7(i)/Ncount2(i)
        Tmp(i)    =  Tmp(i)/Ncount(i)
        var_1(i)  =  var_1(i)/Ncount(i)
        var_2(i)  =  var_2(i)/Ncount(i)
        var_3(i)  =  var_3(i)/Ncount(i)
        Rad_mp(i) =  Rad_mp(i)/Ncount(i)
        uup1(i) =  uup1(i)/Ncount2(i)
        vvp1(i) =  vvp1(i)/Ncount2(i)
        wwp1(i) =  wwp1(i)/Ncount2(i)
        uup2(i) =  uup2(i)/Ncount2(i)
        vvp2(i) =  vvp2(i)/Ncount2(i)
        wwp2(i) =  wwp2(i)/Ncount2(i)
        uup3(i) =  uup3(i)/Ncount2(i)
        vvp3(i) =  vvp3(i)/Ncount2(i)
        wwp3(i) =  wwp3(i)/Ncount2(i)
        uup4(i) =  uup4(i)/Ncount2(i)
        vvp4(i) =  vvp4(i)/Ncount2(i)
        wwp4(i) =  wwp4(i)/Ncount2(i)
        uup5(i) =  uup5(i)/Ncount2(i)
        vvp5(i) =  vvp5(i)/Ncount2(i)
        wwp5(i) =  wwp5(i)/Ncount2(i)
        uup6(i) =  uup6(i)/Ncount2(i)
        vvp6(i) =  vvp6(i)/Ncount2(i)
        wwp6(i) =  wwp6(i)/Ncount2(i)
        uup7(i) =  uup7(i)/Ncount2(i)
        vvp7(i) =  vvp7(i)/Ncount2(i)
        wwp7(i) =  wwp7(i)/Ncount2(i)
        uwp1(i) =  uwp1(i)/Ncount2(i)
        uwp2(i) =  uwp2(i)/Ncount2(i)
        uwp3(i) =  uwp3(i)/Ncount2(i)
        uwp4(i) =  uwp4(i)/Ncount2(i)
        uwp5(i) =  uwp5(i)/Ncount2(i)
        uwp6(i) =  uwp6(i)/Ncount2(i)
        uwp7(i) =  uwp7(i)/Ncount2(i)
      end if
    end do 
    call wait

    open(3,FILE='Nusselt_Utau.dat')
    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        write(3,'(5E15.7,I6)') Rad_mp(i)/2.0, 0.2/(1.4e-4*(Tmp(i)-20.0)), &
                            0.2/(1.4e-4*(var_2(i)-20.0)),              &   
                            (abs(Ump(i)*0.0001/var_1(i)))**0.5,        &
                            (abs(var_3(i)*0.0001/var_1(i)))**0.5, &
                            Ncount(i)       
      end if
    end do 
    close(3)

    open(3,FILE='Kinetic_energy.dat')
    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        write(3,'(6E15.7,1X,I6)') Rad_mp(i)/2.0, & 
                            0.5*(uup1(i) + vvp1(i) + wwp1(i)), &     
                            0.5*(uup2(i) + vvp2(i) + wwp2(i)), &     
                            0.5*(uup3(i) + vvp3(i) + wwp3(i)), &     
                            0.5*(uup4(i) + vvp4(i) + wwp4(i)), &     
                            0.5*(uup4(i) + vvp5(i) + wwp5(i)), Ncount2(i)     
      end if
    end do 
    close(3)

    open(3,FILE='Res_position_0.025.dat')
    do i = 1, Nprob
      if(Ncount2(i) /= 0) then
        write(3,'(9E15.7,1X,I6)') Rad_mp(i)/2.0, & 
                            Ump1(i), Vmp1(i), Wmp1(i), & 
                            uup1(i), vvp1(i), wwp1(i), uwp1(i), &     
                            0.5*(uup1(i) + vvp1(i) + wwp1(i)), &     
                            Ncount2(i)     
      end if
    end do 
    close(3)

    open(3,FILE='Res_position_0.050.dat')
    do i = 1, Nprob
      if(Ncount2(i) /= 0) then
        write(3,'(9E15.7,1X,I6)') Rad_mp(i)/2.0, & 
                            Ump2(i), Vmp2(i), Wmp2(i), & 
                            uup2(i), vvp2(i), wwp2(i), uwp2(i), &     
                            0.5*(uup2(i) + vvp2(i) + wwp2(i)), &     
                            Ncount2(i)     
      end if
    end do 
    close(3)

    open(3,FILE='Res_position_0.100.dat')
    do i = 1, Nprob
      if(Ncount2(i) /= 0) then
        write(3,'(9E15.7,1X,I6)') Rad_mp(i)/2.0, & 
                            Ump3(i), Vmp3(i), Wmp3(i), & 
                            uup3(i), vvp3(i), wwp3(i), uwp3(i), &     
                            0.5*(uup3(i) + vvp3(i) + wwp3(i)), &     
                            Ncount2(i)     
      end if
    end do 
    close(3)

    open(3,FILE='Res_position_1.000.dat')
    do i = 1, Nprob
      if(Ncount2(i) /= 0) then
        write(3,'(9E15.7,1X,I6)') Rad_mp(i)/2.0, & 
                            Ump4(i), Vmp4(i), Wmp4(i), & 
                            uup4(i), vvp4(i), wwp4(i), uwp4(i), &     
                            0.5*(uup4(i) + vvp4(i) + wwp4(i)), &     
                            Ncount2(i)     
      end if
    end do 
    close(3)

    open(3,FILE='Res_position_2.000.dat')
    do i = 1, Nprob
      if(Ncount2(i) /= 0) then
        write(3,'(9E15.7,1X,I6)') Rad_mp(i)/2.0, & 
                            Ump5(i), Vmp5(i), Wmp5(i), & 
                            uup5(i), vvp5(i), wwp5(i), uwp5(i), &     
                            0.5*(uup5(i) + vvp5(i) + wwp5(i)), &     
                            Ncount2(i)     
      end if
    end do 
    close(3)

    open(3,FILE='Res_position_3.000.dat')
    do i = 1, Nprob
      if(Ncount2(i) /= 0) then
        write(3,'(9E15.7,1X,I6)') Rad_mp(i)/2.0, & 
                            Ump6(i), Vmp6(i), Wmp6(i), & 
                            uup6(i), vvp6(i), wwp6(i), uwp6(i), &     
                            0.5*(uup6(i) + vvp6(i) + wwp6(i)), &     
                            Ncount2(i)     
      end if
    end do 
    close(3)

    open(3,FILE='Res_position_4.000.dat')
    do i = 1, Nprob
      if(Ncount2(i) /= 0) then
        write(3,'(9E15.7,1X,I6)') Rad_mp(i)/2.0, & 
                            Ump7(i), Vmp7(i), Wmp7(i), & 
                            uup7(i), vvp7(i), wwp7(i), uwp7(i), &     
                            0.5*(uup7(i) + vvp7(i) + wwp7(i)), &     
                            Ncount2(i)     
      end if
    end do 
    close(3)



  deallocate(Np)
  deallocate(Ump)
  deallocate(Vmp)
  deallocate(Wmp)

  deallocate(Ump1)
  deallocate(Vmp1)
  deallocate(Wmp1)

  deallocate(Ump2)
  deallocate(Vmp2)
  deallocate(Wmp2)

  deallocate(Ump3)
  deallocate(Vmp3)
  deallocate(Wmp3)

  deallocate(Ump4)
  deallocate(Vmp4)
  deallocate(Wmp4)

  deallocate(Ump5)
  deallocate(Vmp5)
  deallocate(Wmp5)

  deallocate(Ump6)
  deallocate(Vmp6)
  deallocate(Wmp6)

  deallocate(Ump7)
  deallocate(Vmp7)
  deallocate(Wmp7)

  deallocate(uup1)
  deallocate(vvp1)
  deallocate(wwp1)

  deallocate(uup2)
  deallocate(vvp2)
  deallocate(wwp2)

  deallocate(uup3)
  deallocate(vvp3)
  deallocate(wwp3)

  deallocate(uup4)
  deallocate(vvp4)
  deallocate(wwp4)

  deallocate(uup5)
  deallocate(vvp5)
  deallocate(wwp5)

  deallocate(uup6)
  deallocate(vvp6)
  deallocate(wwp6)

  deallocate(uup7)
  deallocate(vvp7)
  deallocate(wwp7)

  deallocate(uvp)
  deallocate(uwp)
  deallocate(vwp)

  deallocate(uwp1)
  deallocate(uwp2)
  deallocate(uwp3)
  deallocate(uwp4)
  deallocate(uwp5)
  deallocate(uwp6)
  deallocate(uwp7)

  deallocate(var_1)
  deallocate(var_2)
  deallocate(var_3)
  deallocate(Rad_mp)
  deallocate(Rad_1)
  deallocate(Ncount)
  deallocate(Ncount2)

  if(HOT==YES) then
    deallocate(Tmp)
    deallocate(TTp)
    deallocate(uTp)
    deallocate(vTp)
    deallocate(wTp)
  end if
  END SUBROUTINE UserProbe1D_Nusselt_jet
