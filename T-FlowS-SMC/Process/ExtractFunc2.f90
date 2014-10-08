!======================================================================!
  SUBROUTINE ExtractFunc2(nameOUT)
!----------------------------------------------------------------------!
! Reads the ".1D" file created by the "Generator" and averages the     !
! results in the planes defined by coordinates in it. Then averages    !
! the values of Un, Vn, Wn, uu, vv, ww, uv, uw and vw and     !
! writes them into file ".1Dr".                                        !
!----------------------------------------------------------------------!
  USE all_mod
  USE allp_mod
  USE les_mod
  USE pro_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
!------------------------------[Calling]-------------------------------!
  INTEGER             :: c, i, j, k, N_x, N_y, count
  CHARACTER           :: nameIN*80
  REAL,ALLOCATABLE    :: x1(:), x2(:), y(:), U_mean(:), V_mean(:), & 
                                 UU_mean(:), VV_mean(:), &
  INTEGER,ALLOCATABLE :: Np(:), Ncount(:)
  CHARACTER, OPTIONAL :: nameOUT*(*)
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $  
!======================================================================!


  open(996, FILE='x.dat')
  if(this < 2) write(6, *) '# Now reading the file: x.dat ' 
  read(996,*) N_x
  allocate(x1(N_x))
  allocate(x2(N_x))
  do i=1, N_x
    read(996,*) x1(i), x2(i)
  end do
  close(996)

!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
    if(this < 2) write(6, *) '# Now reading the file: y.dat ' 
    open(996, FILE='y.dat')

!---- write the number of probes 
    read(996,*) N_y
    allocate(y(N_y))

!---- write the probe coordinates out
    do pl=1,N_y
      read(996,*) y(pl)
    end do
    close(996)

    allocate(Np(N_y));    Np=0 
    allocate(U_mean(N_y));   U_mean=0.0
    allocate(V_mean(N_y));   V_mean=0.0
    allocate(UU_mean(N_y));   UU_mean=0.0
    allocate(VV_mean(N_y));   VV_mean=0.0
    allocate(Ncount(N_y)); Ncount=0
    count = 0

    end if  
!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
  do i = 1, N_x
    do j = 1, N_y-1
      do c=1,NC
        if(xc(c) < x2(i) .and. xc(c) > x1(i)) then
          if(yc(c) > y(j) .and. yc(c) < y(j+1)) then

            U_mean(j)   = U_mean(j) +  U % mean(c)
            V_mean(j)   = V_mean(j) +  V % mean(c)
            UU_mean(j)   = UU_mean(j) + ( UU % mean(c) - U % mean(c) * U % mean(c) )
            VV_mean(j)   = VV_mean(j) + ( VV % mean(c) - V % mean(c) * V % mean(c) )
            Ncount(j) = Ncount(j) + 1

          end if
        end if
      end do 
    end do 
!still inside x loop
!---- average over all processors
    do k=1, N_y
      call kGlSum(Ncount(k))
      call GloSum(U_mean(k))
      call GloSum(V_mean(k))
      call GloSum(UU_mean(k))
      call GloSum(VV_mean(k))

      count =  count + Ncount(k) 

    end do

    !creating files
    nameIN='xx.xxx'
    write(nameIN,'(F6.3)') (x1(i)+x2(i))/2
    write(nameIN(len_trim(nameIN)+1:len_trim(nameIN)+5),'(A)') '.dat'

    open(996, FILE=nameIN)

    do j=1, N_y !writing the results
      if(Ncount(j) /= 0) then ! if there no cells catched there is no point to create file
                            U_mean(j)  = U_mean(j)  / Ncount(j)
                            V_mean(j)  = V_mean(j)  / Ncount(j)
                            UU_mean(j) = UU_mean(j) / Ncount(j)
                            VV_mean(j) = VV_mean(j) / Ncount(j) 
      write(996,'(5E15.7)') y(j),           &
                            U_mean(j),      &
                            V_mean(j),      & 
                            UU_mean(j),     &
                            VV_mean(j)
      end if ! if there no cells catched there is no point to create file
    end do

    end do 
    close(996)
  end do ! end x loop (i=1, N_x)
 
  deallocate(Np)
  deallocate(y)
  deallocate(U_mean)
  deallocate(V_mean)
  deallocate(UU_mean)
  deallocate(VV_mean)
  deallocate(Ncount)


  if(this < 2) write(*,*) 'Finished with UserCutLines_Horiz_LES'

  END SUBROUTINE UserCutLines_Horiz_LES
