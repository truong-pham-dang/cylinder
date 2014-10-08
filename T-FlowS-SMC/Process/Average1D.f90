!======================================================================!
  SUBROUTINE Average1D(y)   
!----------------------------------------------------------------------!
! Reads the ".1D" file created by the "Generator" and averages the     !
! results in the points defined by coordinates in it and then          !
! distributes teh result back to cells.                                !
!----------------------------------------------------------------------!
  USE all_mod
  USE allp_mod
  USE les_mod
  USE pro_mod
  USE par_mod
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
                                 Ksgsp(:)
  INTEGER,ALLOCATABLE :: Np(:), AveragePoint(:)
!--------------------------------[CVS]---------------------------------!
!  $Id: Average1D.f90,v 1.5 2009/06/30 11:40:17 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/Average1D.f90,v $  
!======================================================================!

  if(this  < 2)  & 
    write(*,*) '# Input probe file name [skip cancels]:'
  call ReadC(7,inp,tn,ts,te)  
  read(inp(ts(1):te(1)),'(A80)') namPro
  answer=namPro
  call ToUppr(answer) 
  if(answer == 'SKIP') return 

  call wait 

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
    read(9,'(I8,1PE17.8,1PE17.8)') dummy, z_p(pl)
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
  if(HOT==YES) then
    allocate(Tmp(Nprob));   Tmp=0.0
    allocate(TTp(Nprob));   TTp=0.0
    allocate(uTp(Nprob));   uTp=0.0
    allocate(vTp(Nprob));   vTp=0.0
    allocate(wTp(Nprob));   wTp=0.0
  end if  
  allocate(AveragePoint(NC)); AveragePoint = 0

!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
  do pl=1,Nprob
    do c=1,NC
      if( Approx( y(c),z_p(pl) ) ) then  
        AveragePoint(c) = pl
        Np(pl)    = Np(pl) + 1
        Ump(pl)   = Ump(pl) + U % mean(c)
        Vmp(pl)   = Vmp(pl) + V % mean(c)
        Wmp(pl)   = Wmp(pl) + W % mean(c)
        uup(pl)   = uup(pl) + uu % mean(c)
        vvp(pl)   = vvp(pl) + vv % mean(c)
        wwp(pl)   = wwp(pl) + ww % mean(c)
        uvp(pl)   = uvp(pl) + uv % mean(c)
        uwp(pl)   = uwp(pl) + uw % mean(c)
        vwp(pl)   = vwp(pl) + vw % mean(c)
        Ksgsp(pl) = Ksgsp(pl) + Ksgs(c) 
        if(HOT==YES) then
          Tmp(pl)   = Tmp(pl) + T % mean(c) 
          TTp(pl)   = TTp(pl) + TT % mean(c)
          uTp(pl)   = uTp(pl) + uT % mean(c)
          vTp(pl)   = vTp(pl) + vT % mean(c)
          wTp(pl)   = wTp(pl) + wT % mean(c)
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
    if(HOT==YES) then
      Tmp(pl) = Tmp(pl) / Np(pl)
      TTp(pl) = TTp(pl) / Np(pl)
      uTp(pl) = uTp(pl) / Np(pl)
      vTp(pl) = vTp(pl) / Np(pl)
      wTp(pl) = wTp(pl) / Np(pl)
    end if
  end do

!--------------------------------!
!     distribute the results     !
!--------------------------------!
  do c=1,NC
    U % mean(c)  = Ump(AveragePoint(c))
    V % mean(c)  = Vmp(AveragePoint(c))
    W % mean(c)  = Wmp(AveragePoint(c))
    uu % mean(c) = uup(AveragePoint(c))
    vv % mean(c) = vvp(AveragePoint(c))
    ww % mean(c) = wwp(AveragePoint(c))
    uv % mean(c) = uvp(AveragePoint(c))
    uw % mean(c) = uwp(AveragePoint(c))
    vw % mean(c) = vwp(AveragePoint(c))
    Ksgs(c)      = Ksgsp(AveragePoint(c))
    if(HOT==YES) then
      T % mean(c)  = Tmp(AveragePoint(c))
      TT % mean(c) = TTp(AveragePoint(c))
      uT % mean(c) = uTp(AveragePoint(c))
      vT % mean(c) = vTp(AveragePoint(c))
      wT % mean(c) = wTp(AveragePoint(c))
    end if  
  end do

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
  if(HOT==YES) then
    deallocate(Tmp)
    deallocate(TTp)
    deallocate(uTp)
    deallocate(vTp)
    deallocate(wTp)
  end if
  deallocate(AveragePoint)

  END SUBROUTINE Average1D  
