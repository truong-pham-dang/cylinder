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
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:),               & 
                                 var_1(:), var_2(:), var_3(:), var_4(:), Rad_1(:), Rad_mp(:)
  INTEGER,ALLOCATABLE :: Np(:), Ncount(:)
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
  allocate(Rad_mp(Nprob));  Rad_mp=0.0
  allocate(var_1(Nprob));  var_1=0.0
  allocate(var_2(Nprob));  var_2=0.0
  allocate(var_3(Nprob));  var_3=0.0
  allocate(var_4(Nprob));  var_4=0.0

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
      Rad_1(i) = abs(z_p(i)) 
    end do 

    do i = 1, Nprob-1
      do s=1,NS
        c1=SideC(1,s)
        c2=SideC(2,s)
        if(c2 < 0) then
          if(TypeBC(c2).eq.WALLFL) then
            Rad_2 = (xc(c1)*xc(c1) + yc(c1)*yc(c1))**0.5
            if(Rad_2 < Rad_1(i) .and. Rad_2 > Rad_1(i+1).and.zc(c1) < 0.1) then
              R           = (xc(c1)*xc(c1) + yc(c1)*yc(c1))**0.5 
              Ump(i)   = Ump(i) + (U % mean(c1) * U % mean(c1) + V % mean(c1) * V % mean(c1) + W % mean(c1)*W % mean(c1))**0.5 
              Wmp(i)   = Wmp(i) + W % mean(c1)
              var_1(i) = var_1(i) + zc(c1)
              var_3(i) = var_3(i) +  (U % n(c1) * U % n(c1) + V % n(c1) * V % n(c1) + W % n(c1)*W % n(c1))**0.5
              Rad_mp(i)= Rad_mp(i) + (xc(c1)*xc(c1) + yc(c1)*yc(c1))**0.5
              Tmp(i)   = Tmp(i) + T % mean(c2)
              var_2(i) = var_2(i) + T % n(c2)
              var_4(i) = Px(c1) * xc(c1)/R + Py(c1) * yc(c1)/R  
              Ncount(i)= Ncount(i) + 1
            end if
          end if
        end if
      end do
    end do 
!---- average over all processors
  do pl=1, Nprob
    call IGlSum(Ncount(pl))

    call GloSum(Ump(pl))
    call GloSum(Wmp(pl))
    call GloSum(var_1(pl))
    call GloSum(var_2(pl))
    call GloSum(var_3(pl))
    call GloSum(var_4(pl))
    call GloSum(Rad_mp(pl))
    call GloSum(Tmp(pl))

    count =  count + Ncount(pl) 
  end do

    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        Wmp(i)    =  Wmp(i)/Ncount(i)
        Ump(i)    =  Ump(i)/Ncount(i)
        Tmp(i)    =  Tmp(i)/Ncount(i)
        var_1(i)  =  var_1(i)/Ncount(i)
        var_2(i)  =  var_2(i)/Ncount(i)
        var_3(i)  =  var_3(i)/Ncount(i)
        var_4(i)  =  var_4(i)/Ncount(i)
        Rad_mp(i) =  Rad_mp(i)/Ncount(i)
      end if
    end do 
    call wait

    open(3,FILE='Nusselt_Utau.dat')
    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        write(3,'(6E15.7,I6)') Rad_mp(i)/2.0, 0.2/(1.4e-4*(Tmp(i)-20.0)), &
                            0.2/(1.4e-4*(var_2(i)-20.0)),              &   
                            (abs(Ump(i)*0.0001/var_1(i)))**0.5,        &
                            (abs(var_3(i)*0.0001/var_1(i)))**0.5, var_4(i), Ncount(i)       
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
  deallocate(var_2)
  deallocate(var_3)
  deallocate(var_4)
  deallocate(Rad_mp)
  deallocate(Rad_1)
  deallocate(Ncount)

  if(HOT==YES) then
    deallocate(Tmp)
    deallocate(TTp)
    deallocate(uTp)
    deallocate(vTp)
    deallocate(wTp)
  end if
  write(*,*)'Finished Nusselt_jet'
  END SUBROUTINE UserProbe1D_Nusselt_jet
