!======================================================================!
  SUBROUTINE UserProbe1D(y) 
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
!------------------------------[Calling]-------------------------------!
  INTERFACE
    LOGICAL FUNCTION Approx(A,B,tol)
      REAL           :: A,B
      REAL, OPTIONAL :: tol
    END FUNCTION Approx
  END INTERFACE 
!-------------------------------[Locals]-------------------------------!
  INTEGER             :: Nprob, pl, c, dummy
  CHARACTER           :: namCoo*80, namPro*80, answer*80
  REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:), VISt_m(:), VISt_m1(:) ,VISt_m2(:),&
                                 Kin_mp(:), Eps_mp(:), v_2_mp(:), f22_mp(:), &
                                 var_1(:), var_2(:), var_3(:)
  INTEGER,ALLOCATABLE :: Np(:)
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $  
!======================================================================!

!  if(this  < 2)  & 
!    write(*,*) '# Input probe file name [skip cancels]:'
!  call ReadC(7,inp,tn,ts,te)  
!  read(inp(ts(1):te(1)),'(A80)') namPro
!  answer=namPro
!  call ToUppr(answer) 
!  if(answer == 'SKIP') return 

!  call wait 

  
  namPro = name
!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
  namCoo = name
  namCoo(len_trim(name)+1:len_trim(name)+3) = '.1D'
  write(6, *) '# Now reading the file:', namCoo
  open(9, FILE=namCoo)

!---- write the number of probes 
  read(9,'(I8)') Nprob
  allocate(z_p(Nprob))

!---- write the probe coordinates out
  do pl=1,Nprob
    read(9,'(I8,1PE17.8)') dummy, z_p(pl)
  end do

  close(9)

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
  allocate(VISt_m(Nprob)); VISt_m=0.0
  allocate(VISt_m1(Nprob)); VISt_m1=0.0
  allocate(VISt_m2(Nprob)); VISt_m2=0.0
  allocate(Kin_mp(Nprob)); Kin_mp=0.0
  allocate(Eps_mp(Nprob)); Eps_mp=0.0
  allocate(v_2_mp(Nprob)); v_2_mp=0.0
  allocate(f22_mp(Nprob)); f22_mp=0.0
  allocate(var_1(Nprob));  var_1=0.0
  allocate(var_2(Nprob));  var_2=0.0
  allocate(var_3(Nprob));  var_3=0.0
  if(HOT==YES) then
    allocate(Tmp(Nprob));   Tmp=0.0
    allocate(TTp(Nprob));   TTp=0.0
    allocate(uTp(Nprob));   uTp=0.0
    allocate(vTp(Nprob));   vTp=0.0
    allocate(wTp(Nprob));   wTp=0.0
  end if  

  call CalcShear(U % n, V % n, W % n, Shear)
  call GraPhi(U % n, 1, Ux,.TRUE.)    ! dU/dx
  call GraPhi(U % n, 2, Uy,.TRUE.)    ! dU/dy
  call GraPhi(U % mean, 3, Uz,.TRUE.)    ! dU/dz
  call GraPhi(V % n, 1, Vx,.TRUE.)    ! dV/dx
  call GraPhi(V % n, 2, Vy,.TRUE.)    ! dV/dy
  call GraPhi(V % n, 3, Vz,.TRUE.)    ! dV/dz
  call GraPhi(W % mean, 1, Wx,.TRUE.)    ! dW/dx
  call GraPhi(W % n, 2, Wy,.TRUE.)    ! dW/dy
  call GraPhi(W % mean, 3, Wz,.TRUE.)    ! dW/dz
!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
  do pl=1,Nprob
    do c=1,NC
      if( Approx( y(c),z_p(pl) ) ) then
        Np(pl)    = Np(pl) + 1
        Ump(pl)   = Ump(pl) + U % mean(c)
        Vmp(pl)   = Vmp(pl) + V % mean(c)
        Wmp(pl)   = Wmp(pl) + W % mean(c)
        uup(pl)   = uup(pl) + (uu % mean(c)-U % mean(c) * U % mean(c))
        vvp(pl)   = vvp(pl) + (vv % mean(c)-V % mean(c) * V % mean(c))
        wwp(pl)   = wwp(pl) + (ww % mean(c)-W % mean(c) * W % mean(c))
        uvp(pl)   = uvp(pl) + (uv % mean(c)-U % mean(c) * V % mean(c)) 
        uwp(pl)   = uwp(pl) + (uw % mean(c)-U % mean(c) * W % mean(c)) 
        vwp(pl)   = vwp(pl) + (vw % mean(c)-V % mean(c) * W % mean(c))
        VISt_m(pl)= VISt_m(pl) + VISt_mean(c)
        VISt_m1(pl)= VISt_m1(pl) + VISt(c) 

        if(HOT==YES) then
          Tmp(pl)   = Tmp(pl) + T % mean(c) 
          TTp(pl)   = TTp(pl) + (TT % mean(c) - T % mean(c) * T % mean(c))
          uTp(pl)   = uTp(pl) + (uT % mean(c) - u % mean(c) * T % mean(c))
          vTp(pl)   = vTp(pl) + (vT % mean(c) - v % mean(c) * T % mean(c))
          wTp(pl)   = wTp(pl) + (wT % mean(c) - w % mean(c) * T % mean(c))
        end if  
      end if
    end do
  end do

!---- average over all processors
  do pl=1,Nprob
    call IGlSum(Np(pl))

    call GloSum(Ump(pl))
    call GloSum(Vmp(pl))
    call GloSum(Wmp(pl))

    call GloSum(uup(pl))
    call GloSum(vvp(pl))
    call GloSum(wwp(pl))

    call GloSum(uvp(pl))
    call GloSum(uwp(pl))
    call GloSum(vwp(pl))

    call GloSum(VISt_m(pl))
    call GloSum(VISt_m1(pl))
    call GloSum(VISt_m2(pl))

    call GloSum(Ksgsp(pl))

    if(HOT==YES) then
      call GloSum(Tmp(pl))
      call GloSum(TTp(pl))
      call GloSum(uTp(pl))
      call GloSum(vTp(pl))
      call GloSum(wTp(pl))
    end if
  end do

  do pl=1,Nprob
    Ump(pl)   = Ump(pl) / Np(pl) 
    Vmp(pl)   = Vmp(pl) / Np(pl)
    Wmp(pl)   = Wmp(pl) / Np(pl)
    uup(pl)   = uup(pl) / Np(pl)
    vvp(pl)   = vvp(pl) / Np(pl)
    wwp(pl)   = wwp(pl) / Np(pl)
    uvp(pl)   = uvp(pl) / Np(pl)
    uwp(pl)   = uwp(pl) / Np(pl)
    vwp(pl)   = vwp(pl) / Np(pl)
    Ksgsp(pl) = Ksgsp(pl) / Np(pl)
    VISt_m(pl) = VISt_m(pl) / Np(pl)
    VISt_m1(pl) = VISt_m1(pl) / Np(pl)
    VISt_m2(pl) = VISt_m2(pl) / Np(pl)
    if(HOT==YES) then
      Tmp(pl) = Tmp(pl) / Np(pl)
      TTp(pl) = TTp(pl) / Np(pl)
      uTp(pl) = uTp(pl) / Np(pl)
      vTp(pl) = vTp(pl) / Np(pl)
      wTp(pl) = wTp(pl) / Np(pl)
    end if
  end do

!<<<<<<<<<<<<<<<<<<<<<<<<<!
!     create 1Dr file     !
!<<<<<<<<<<<<<<<<<<<<<<<<<!
  if(This < 2) then
    namPro(len_trim(namPro)+1:len_trim(namPro)+4) = '.1Dr'
    write(6, *) '# Now creating the file:', namPro
    open(9, FILE=namPro)

!---- write the number of probes
    write(9,'(A1,I8)') '#', Nprob
    write(9,'(A80)')  '# 1-z  2-Umean  3-Vmean  4-Wmean  5-uu  6-vv  7-ww  8-uv  9-uw  10-vw  11-Ksgs'
    if(HOT==YES) then
      write(9,'(A38)') '# 12-Tmean  13-TT  14-ut  15-vt  16-wt'
    end if

!---- write the values at the probes
    do pl=1,Nprob
      if(HOT==YES) then
        write(9,'(16E17.7,2I8)')                    &
                       z_p(pl),                     &
                       Ump(pl), Vmp(pl), Wmp(pl),   &
                       uup(pl), vvp(pl), wwp(pl),   &
                       uvp(pl), uwp(pl), vwp(pl),   &
                       Ksgsp(pl), Tmp(pl), TTp(pl), & 
                       uTp(pl), vTp(pl), wTp(pl),   & 
                       pl, Np(pl)
      else 
        write(9,'(12E14.6)')                   &
                       z_p(pl),                    &
                       Ump(pl), Vmp(pl), Wmp(pl),  &
                       uup(pl), vvp(pl), wwp(pl),  &
                       uvp(pl), uwp(pl), vwp(pl),  &
                       VISt_m(pl), VISt_m1(pl)
      end if
    end do

    close(9)
  end if


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
  deallocate(Ksgsp)
  deallocate(VISt_m)
  deallocate(VISt_m1)
  deallocate(VISt_m2)
  deallocate(Kin_mp)
  deallocate(Eps_mp)
  deallocate(v_2_mp)
  deallocate(f22_mp)
  deallocate(var_1)
  deallocate(var_2)
  deallocate(var_3)
  if(HOT==YES) then
    deallocate(Tmp)
    deallocate(TTp)
    deallocate(uTp)
    deallocate(vTp)
    deallocate(wTp)
  end if
  END SUBROUTINE UserProbe1D
