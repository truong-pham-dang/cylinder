!======================================================================!
  SUBROUTINE UserHill() 
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
  INTEGER             :: Nprob, pl, c, i, count, k, c1, c2, s
  CHARACTER           :: namCoo*80, namPro*80, answer*80, line_name*11, namOut*80
  REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:)
  INTEGER,ALLOCATABLE :: Np(:), Ncount(:)
  REAL                :: R, dummy, R1, R2
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $  
!======================================================================!



  do k = 0, 5

!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
  if(k == 0) then  
    if(this < 2) write(6, *) '# Now reading the file: Hill_y_0.2.1D ' 
    open(9, FILE='Hill_y_0.2.1D')

!---- write the number of probes 
    read(9,*) Nprob
    allocate(z_p(Nprob))

!---- write the probe coordinates out
    do pl=1,Nprob
      read(9,*) dummy, z_p(pl)
    end do
    close(9)
  else if(k == 1) then
    if(this < 2) write(6, *) '# Now reading the file: Hill_y_0.6.1D ' 
    open(9, FILE='Hill_y_0.6.1D')

!---- write the number of probes 
    read(9,*) Nprob
    allocate(z_p(Nprob))

!---- write the probe coordinates out
    do pl=1,Nprob
      read(9,*) dummy, z_p(pl)
    end do
    close(9)
  else if(k == 2) then
    if(this < 2) write(6, *) '# Now reading the file: Hill_y_1.0.1D ' 
    open(9, FILE='Hill_y_1.0.1D')

!---- write the number of probes 
    read(9,*) Nprob
    allocate(z_p(Nprob))

!---- write the probe coordinates out
    do pl=1,Nprob
      read(9,*) dummy, z_p(pl)
    end do
    close(9)
  else if(k == 3) then
    if(this < 2) write(6, *) '# Now reading the file: Hill_y_1.4.1D ' 
    open(9, FILE='Hill_y_1.4.1D')

!---- write the number of probes 
    read(9,*) Nprob
    allocate(z_p(Nprob))

!---- write the probe coordinates out
    do pl=1,Nprob
      read(9,*) dummy, z_p(pl)
    end do
    close(9)
  else if(k == 4) then
    if(this < 2) write(6, *) '# Now reading the file: Hill_y_1.8.1D ' 
    open(9, FILE='Hill_y_1.8.1D')

!---- write the number of probes 
    read(9,*) Nprob
    allocate(z_p(Nprob))

!---- write the probe coordinates out
    do pl=1,Nprob
      read(9,*) dummy, z_p(pl)
    end do
    close(9)
  else if(k == 5) then
    if(this < 2) write(6, *) '# Now reading the file: Hill_y_2.0.1D ' 
    open(9, FILE='Hill_y_2.0.1D')

!---- write the number of probes 
    read(9,*) Nprob
    allocate(z_p(Nprob))

!---- write the probe coordinates out
    do pl=1,Nprob
      read(9,*) dummy, z_p(pl)
    end do
    close(9)
  end if
  
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

    if(k == 0) then
      R1 = 0.16885  
      R2 = 0.210227
    else if(k == 1) then
      R1 = 0.582824
      R2 = 0.624227
    else if(k == 2) then
      R1 = 0.957342
      R2 = 1.042677
    else if(k == 3) then
      R1 = 1.375784
      R2 = 1.417197
    else if(k == 4) then
      R1 = 1.789784
      R2 = 1.831177
    else if(k == 5) then
      R1 = 1.957214
      R2 = 2.000000
    end if  

    

    do i = 1, Nprob-1
      do c=1,NC
        if( xc(c) > R1 .and. xc(c) < R2) then
          if(yc(c) > z_p(i) .and. yc(c) < z_p(i+1)) then
            Ump(i)   = Ump(i) + U % mean(c) 
            Vmp(i)   = Vmp(i) + V % mean(c)
            Wmp(i)   = Wmp(i) + W % mean(c)
            uuP(i)   = uup(i) + (uu % mean(c)- U % mean(c) * U % mean(c))
            vvp(i)   = vvp(i) + (vv % mean(c)- V % mean(c) * V % mean(c))
            wwp(i)   = wwp(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
            uvp(i)   = uvp(i) + (uv % mean(c)- U % mean(c) * V % mean(c))
            uwp(i)   = uwp(i) + (uw % mean(c)- U % mean(c) * W % mean(c))
            vwp(i)   = vwp(i) + (vw % mean(c)- V % mean(c) * W % mean(c))

            if(HOT==YES) then
              Tmp(i)   = Tmp(i) + T % mean(c)
              TTp(i)   = TTp(i) + (TT % mean(c) - T % mean(c) * T % mean(c))
              uTp(i)   = uTp(i) + (uT % mean(c) - U % mean(c) * T % mean(c))
              vTp(i)   = vTp(i) + (vT % mean(c) - V % mean(c) * T % mean(c))
              wTp(i)   = wTp(i) + (wT % mean(c) - w % mean(c) * T % mean(c))
            end if
       
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

      call GloSum(uup(pl))
      call GloSum(vvp(pl))
      call GloSum(wwp(pl))

      call GloSum(uvp(pl))
      call GloSum(uwp(pl))
      call GloSum(vwp(pl))

      count =  count + Ncount(pl) 

      if(HOT==YES) then
        call GloSum(Tmp(pl))
        call GloSum(TTp(pl))
        call GloSum(uTp(pl))
        call GloSum(vTp(pl))
        call GloSum(wTp(pl))
      end if
    end do

    line_name  = 'Hill_line_x'
    write(line_name(11:11),'(I1)') k

    do i = 1, Nprob
      Ncount(i) = max(Ncount(i),1) 
    end do

    open(3,FILE=line_name)
    do i = 1, Nprob-1
      Wmp(i)    =  Wmp(i)/Ncount(i)
      Ump(i)    =  Ump(i)/Ncount(i)
      Vmp(i)    =  Vmp(i)/Ncount(i)
      uup(i)    =  uup(i)/Ncount(i)
      vvp(i)    =  vvp(i)/Ncount(i)
      wwp(i)    =  wwp(i)/Ncount(i)
      uvp(i)    =  uvp(i)/Ncount(i)
      uwp(i)    =  uwp(i)/Ncount(i)
      vwp(i)    =  vwp(i)/Ncount(i)
      Tmp(i)    =  Tmp(i)/Ncount(i)
      TTp(i)    =  TTp(i)/Ncount(i)
      uTp(i)    =  uTp(i)/Ncount(i)
      vTp(i)    =  vTp(i)/Ncount(i)
      wTp(i)    =  wTp(i)/Ncount(i)

      write(3,'(15E15.7,1X,I6)') (z_p(i)+z_p(i+1))/1.0, Ump(i), Vmp(i), Wmp(i) , &
                                 uup(i), vvp(i), wwp(i), uvp(i), uwp(i), vwp(i), &
                                 Tmp(i), TTp(i), uTp(i), vTp(i), wTp(i),         &
                                 Ncount(i) 
    end do 
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
    deallocate(Ncount)
    if(HOT==YES) then
      deallocate(Tmp)
      deallocate(TTp)
      deallocate(uTp)
      deallocate(vTp)
      deallocate(wTp)
    end if
  end do   !end number of radius

  if(this < 2) write(*,*) 'Finished with UserHill '

  END SUBROUTINE UserHill
