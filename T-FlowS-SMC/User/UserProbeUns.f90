!======================================================================!
  SUBROUTINE UserProbeUns
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
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER   Nprob, pl, s, c1, c2, dummy
  CHARACTER namCoo*80, namPro*80, answer*80
  REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                         uup(:), vvp(:), wwp(:),         &
                         uvp(:), uwp(:), vwp(:),         &
                         Ksgsp(:)
  REAL                :: f1, f2, zc1, zc2
  REAL                :: Um_p, Vm_p, Wm_p, &
                         uu_p, vv_p, ww_p, &
                         uv_p, uw_p, vw_p, &
                         Ksgs_p  
  INTEGER,ALLOCATABLE :: Np(:)
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbeUns.f90,v 1.12 2002/10/30 16:30:03 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbeUns.f90,v $  
!======================================================================!

  if(this  < 2)  & 
    write(*,*) '# Input probe file name [skip cancels]:'
  call ReadC(7,inp,tn,ts,te)  
  read(inp(ts(1):te(1)),'(A80)') namPro
  answer=namPro
  call ToUppr(answer) 
  if(answer == 'SKIP') return 

  call wait 

  call Exchng(U % mean)
  call Exchng(V % mean)
  call Exchng(W % mean)
  call Exchng(uu % mean)
  call Exchng(vv % mean)
  call Exchng(ww % mean)
  call Exchng(uv % mean)
  call Exchng(uw % mean)
  call Exchng(vw % mean)
  call Exchng(Ksgs)

!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
  namCoo = name
  namCoo(len_trim(name)+1:len_trim(name)+3) = '.1D'
  write(6, *) 'Now reading the file:', namCoo
  open(9, FILE=namCoo)

!---- write the number of probes 
  read(9,'(I8)') Nprob
  allocate(z_p(Nprob))

!---- write the probe coordinates out
  do pl=1,Nprob
    read(9,*) dummy, z_p(pl)
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

!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
  do pl=1,Nprob

    do s=1,NS

      c1=SideC(1,s)
      c2=SideC(2,s)

      if( (c2>0)  .or.  ((c2<0) .and. (TypeBC(c2)==BUFFER)) ) then

        zc1 = zc(c1)
        zc2 = zc(c1) + Dz(s)

        if( (zc1 > z_p(pl) .and. zc2 < z_p(pl)) .or. & 
            (zc1 < z_p(pl) .and. zc2 > z_p(pl)) ) then
	  Np(pl) = Np(pl) + 1
	  f1 = abs(zc2-z_p(pl))/(abs(zc2-zc1))
	  f2 = abs(zc1-z_p(pl))/(abs(zc2-zc1))

	  Um_p   = f1 *  U % mean(c1) + f2 *  U % mean(c2)
	  Vm_p   = f1 *  V % mean(c1) + f2 *  V % mean(c2)
	  Wm_p   = f1 *  W % mean(c1) + f2 *  W % mean(c2)
          uu_p   = f1 * uu % mean(c1) + f2 * uu % mean(c2)
          vv_p   = f1 * vv % mean(c1) + f2 * vv % mean(c2)
          ww_p   = f1 * ww % mean(c1) + f2 * ww % mean(c2)
          uv_p   = f1 * uv % mean(c1) + f2 * uv % mean(c2)
          uw_p   = f1 * uw % mean(c1) + f2 * uw % mean(c2)
          vw_p   = f1 * vw % mean(c1) + f2 * vw % mean(c2)
          Ksgs_p = f1 *Ksgs(c1) + f2*Ksgs(c2)

   	  Ump(pl)   = Ump(pl)   + Um_p
	  Vmp(pl)   = Vmp(pl)   + Vm_p
	  Wmp(pl)   = Wmp(pl)   + Wm_p
          uup(pl)   = uup(pl)   + uu_p - Um_p * Um_p
          vvp(pl)   = vvp(pl)   + vv_p - Vm_p * Vm_p
          wwp(pl)   = wwp(pl)   + ww_p - Wm_p * Wm_p
          uvp(pl)   = uvp(pl)   + uv_p - Um_p * Vm_p
          uwp(pl)   = uwp(pl)   + uw_p - Um_p * Wm_p
          vwp(pl)   = vwp(pl)   + vw_p - Vm_p * Wm_p
	  Ksgsp(pl) = Ksgsp(pl) + Ksgs_p
        end if

      endif 

    end do    ! through sides 

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

    call GloSum(Ksgsp(pl))
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
  end do

!<<<<<<<<<<<<<<<<<<<<<<<<<!
!     create 1Dr file     !
!<<<<<<<<<<<<<<<<<<<<<<<<<!
  if(This < 2) then
    namPro(len_trim(namPro)+1:len_trim(namPro)+4) = '.1Dr'
    write(6, *) 'Now creating the file:', namPro
    open(9, FILE=namPro)

!---- write the number of probes
    write(9,'(A1,I8)') '#', Nprob
    write(9,'(A80)') '# 1-z  2-Umean  3-Vmean  4-Wmean  5-uu  6-vv  7-ww  8-uv  9-uw  10-vw  11-Ksgs'

!---- write the values at the probes
    do pl=1,Nprob
      write(9,'(11E12.4,2I8)')                   &
                     z_p(pl),                    &
                     Ump(pl), Vmp(pl), Wmp(pl),  &
                     uup(pl), vvp(pl), wwp(pl),  &
                     uvp(pl), uwp(pl), vwp(pl),  &
                     Ksgsp(pl),                  &  
                     pl, Np(pl)
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

  END SUBROUTINE UserProbeUns
