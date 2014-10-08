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
  INTEGER             :: Nprob, pl, c, dummy, i, count, k
  CHARACTER           :: namCoo*80, namPro*80, answer*80
  REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:),               & 
                                 var_1(:), var_2(:), var_3(:), Rad_1(:), Rad_mp(:)
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

  allocate(Np(Nprob));    Np=0 
  allocate(Ump(Nprob));   Ump=0.0
  allocate(Wmp(Nprob));   Wmp=0.0
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
      k = Nprob -i +1
      Rad_1(i) = abs(z_p(k)) 
    end do 

    do i = 1, Nprob-1
      do c=1,NC
        if(zc(c) < 9.2 .and. zc(c) > 8.8) then
          Rad_2 = (xc(c)*xc(c) + yc(c)*yc(c))**0.5
          if(Rad_2 > Rad_1(i) .and. Rad_2 < Rad_1(i+1)) then
            R           = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 
            Urad_mean   = (U % mean(c) * xc(c) / R  + V % mean(c) * yc(c) / R)
            Utan_mean   = (-U % mean(c) * yc(c) / R  + V % mean(c) * xc(c) / R) 
 
            Ump(i)   = Ump(i) + Urad_mean 
            Wmp(i)   = Wmp(i) + W % n(c)

            var_1(i) = var_1(i) + W % n(c)
            if(HOT==YES) then
              Tmp(i)   = Tmp(i) + T % mean(c)
              TTp(i)   = TTp(i) + (TT % mean(c) - T % mean(c) * T % mean(c))
              uTp(i)   = uTp(i) + (uT % mean(c) - u % mean(c) * T % mean(c))
              vTp(i)   = vTp(i) + (vT % mean(c) - v % mean(c) * T % mean(c))
              wTp(i)   = wTp(i) + (wT % mean(c) - w % mean(c) * T % mean(c))
            end if
       
            Rad_mp(i) = Rad_mp(i) + R 
            Ncount(i) = Ncount(i) + 1
          end if
        end if
      end do 
    end do 
!---- average over all processors
  do pl=1, Nprob
    call IGlSum(Ncount(pl))

    call GloSum(Ump(pl))
    call GloSum(Wmp(pl))
    call GloSum(Rad_mp(pl))


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
  call wait
  write(*,*) Ncount(1), Ncount(70)
10 continue
 
    open(3,FILE='Nusselt_number_Utau.dat')
    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        write(3,'(2E15.7)') Rad_mp(i)/2.0, abs(Wmp(i))
      end if
    end do 
    close(3)

  deallocate(Np)
  deallocate(Ump)
  deallocate(Wmp)
  deallocate(var_1)
  if(HOT==YES) then
    deallocate(Tmp)
    deallocate(TTp)
    deallocate(uTp)
    deallocate(vTp)
    deallocate(wTp)
  end if
  END SUBROUTINE UserProbe1D_Nusselt_jet
