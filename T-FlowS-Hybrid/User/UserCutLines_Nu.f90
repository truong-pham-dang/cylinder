!======================================================================!
  SUBROUTINE UserCutLines_Nu()     
!----------------------------------------------------------------------!
! This program reads name.1D file created by NEU or GEN and averages    
! the results in the homogeneous direction directions.            
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
    REAL :: Ufric, Wall_near
!------------------------------[Calling]-------------------------------!
    INTERFACE
      LOGICAL FUNCTION Approx(A,B,tol)
        REAL           :: A,B
        REAL, OPTIONAL :: tol
      END FUNCTION Approx
    END INTERFACE
!-------------------------------[Locals]-------------------------------!
    INTEGER             :: Nprob, pl, c, i, count, kk, s, c1, c2
    CHARACTER           :: namCoo*80, namPro*80, answer*80, namRes*80
    CHARACTER           :: namRes_plus*80
    REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), &
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:), ind(:),       &
                                 var_1(:), var_2(:), Rad_mp(:),     &
                                 var_3(:), Wall_p(:), Rad_1(:), &
                                 Ufric_p(:)
    INTEGER,ALLOCATABLE :: Np(:), Ncount(:), Ncount2(:)
    REAL                :: R, Urad_mean, Utan_mean, dummy, Lscale, R_max, Rad_2
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exp $
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $
!======================================================================!

!    namPro = name

!    namRes = name
!    namRes(len_trim(name)+1:len_trim(name)+8) = '_res.dat'
!    namRes_plus = name
!    namRes_plus(len_trim(name)+1:len_trim(name)+13) = '_res_plus.dat'
!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
    namCoo = name
1    namCoo(len_trim(name)+1:len_trim(name)+3) = '.1D'
    if(this < 2) write(6, *) '# Now reading the file:', namCoo
    open(9, FILE='rad_coordinate.dat')

!---- write the number of searching intervals 
    read(9,*) Nprob
    if(Nprob == 0) then
      if(this < 2) write(*,*) "You have to create an ascii file with cell-faces coordinates"
      if(this < 2) write(*,*) "named rad_coordinate.dat, like this:"
      if(this < 2) write(*,*) "10  ! number of cells + 1"
      if(this < 2) write(*,*) "0.0"
      if(this < 2) write(*,*) "0.1"
      if(this < 2) write(*,*) "0.2"
      if(this < 2) write(*,*) "... "
      return
    end if
    allocate(z_p(Nprob*2))
    allocate(ind(Nprob*2))

!---- read the intervals positions
    do pl=1,Nprob
      read(9,*) ind(pl), z_p(pl) 
    end do
    close(9)

  allocate(Np(Nprob));     Np=0
  allocate(Ump(Nprob));    Ump=0.0
  allocate(Vmp(Nprob));    Vmp=0.0
  allocate(Wmp(Nprob));    Wmp=0.0
  allocate(Rad_mp(Nprob));  Rad_mp=0.0
  allocate(var_1(Nprob));  var_1=0.0
  allocate(var_2(Nprob));  var_2=0.0
  allocate(var_3(Nprob));  var_3=0.0

  allocate(Rad_1(Nprob));  Rad_1=0.0
  allocate(Ncount(Nprob)); Ncount=0
  allocate(Ncount2(Nprob)); Ncount2=0

  call GraPhi(T % n, 3, Uz,.TRUE.)

  count = 0

    if(HOT==YES) then
      allocate(Tmp(Nprob));   Tmp=0.0
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
            if(Rad_2 < Rad_1(i+1) .and. Rad_2 > Rad_1(i).and.zc(c1) < 0.1) then
              R           = (xc(c1)*xc(c1) + yc(c1)*yc(c1))**0.5 + tiny
              Ump(i)   = Ump(i) + U % n(c1) * xc(c1) / R  + V % n(c1) * yc(c1) / R
              Vmp(i)   = Vmp(i) + (-U % n(c1) * yc(c1) / R  + V % n(c1) * xc(c1) / R)
              Wmp(i)   = Wmp(i) + W % n(c1)
              var_1(i) = var_1(i) + zc(c1)
              Rad_mp(i)= Rad_mp(i) + (xc(c1)*xc(c1) + yc(c1)*yc(c1))**0.5
              Tmp(i)   = Tmp(i) + T % n(c2)
              var_2(i) = var_2(i) + T % n(c1) !Uz(c1)
              var_3(i) = var_3(i) + sqrt(U % n(c1)**2 + V % n(c1)**2 + W % n(c1)**2)
              Ncount(i)= Ncount(i) + 1
            end if
          end if
        end if
      end do
    end do
!---- average over all processors
  do pl=1, Nprob
    call IGlSum(Ncount(pl))
    call IGlSum(Ncount2(pl))

    call GloSum(Ump(pl))
    call GloSum(Vmp(pl))
    call GloSum(Wmp(pl))

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
        Vmp(i)    =  Vmp(i)/Ncount(i)
        Ump(i)    =  Ump(i)/Ncount(i)
        Tmp(i)    =  Tmp(i)/Ncount(i)
        var_1(i)  =  var_1(i)/Ncount(i)
        var_2(i)  =  var_2(i)/Ncount(i)
        var_3(i)  =  var_3(i)/Ncount(i)
        Rad_mp(i) =  Rad_mp(i)/Ncount(i)
      end if
    end do
    call wait

    open(3,FILE='Nusselt_Utau.dat')
    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        write(3,'(4E15.7,I6)') Rad_mp(i)/2.0, 0.2/(1.40e-4*(Tmp(i)-20.0)), &
                            (abs((Ump(i)*0.0001/var_1(i))))**0.5,  &
                            Tmp(i), & 
                            Ncount(i)
      end if
    end do
    close(3)

  deallocate(Np)
  deallocate(Ump)
  deallocate(Vmp)
  deallocate(Wmp)

  deallocate(var_1)
  deallocate(var_2)
  deallocate(var_3)
  deallocate(Rad_mp)
  deallocate(Rad_1)
  deallocate(Ncount)
  deallocate(Ncount2)

  if(HOT==YES) then
    deallocate(Tmp)
  end if

  if(this < 2) write(*,*)'Finished with Nusselt file'



  END SUBROUTINE UserCutLines_Nu
