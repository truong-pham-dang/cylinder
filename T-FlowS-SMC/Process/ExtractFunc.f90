!======================================================================!
  SUBROUTINE ExtractFunc(nameOUT)
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
!----------------------------------------------------------------------!  
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  CHARACTER, OPTIONAL :: nameOUT*(*)
!------------------------------[Calling]-------------------------------!
  INTEGER              :: c, i, j, N_x, N_y
  INTEGER, ALLOCATABLE :: Ncount(:)
  CHARACTER            :: nameIN*80
  REAL,ALLOCATABLE     :: x1(:), x2(:), y(:), &
                                 U_mean(:), V_mean(:), & 
                                 UU_mean(:), VV_mean(:), UV_mean(:)
  LOGICAL              :: FileExists
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $  
!======================================================================!

  nameIN=trim(nameOUT)
  write(nameIN(len_trim(nameIN)+1:len_trim(nameIN)+4),'(A)') '.2D'
  inquire(FILE=nameIN, EXIST=FileExists)

  if (FileExists) then  ! FileExists will be TRUE if the file
    if(this < 2) write(*, *) ' # Now reading the file: ', nameIN 
    open(996, FILE=nameIN) ! file should consist of cells nodes
  else
    if(this < 2) write(*, *) ' # File ', nameIN, ' wasn`t found'
  end if

  !reading xy data from nameIN.2d file
  !basic structure of file:
  !N_x
  read(996,*) N_x
  allocate(x1(1:N_x)); x1=0.0
  allocate(x2(1:N_x)); x2=0.0
  !x1(1) x2(1) - points of nodes, defining layer of searching
  !..
  !x1(N_x) x2(N_x)
  do i=1, N_x ! natural variable for x_point
    read(996,*) x1(i), x2(i)
  end do
  !N_y
  read(996,*) N_y 
  allocate(y(1:N_y)); y=0.0
  !y(1)
  !..
  !y(N_y) - node coordinates of layer domain
    do j=1,N_y ! natural variable for y_point
      read(996,*) y(j)
    end do
  !y array should in ascending order or sorting function should be use here
  !file *.2D must be finished with empty line!
  close(996)


  allocate(U_mean(1:N_y-1));   U_mean=0.0
  allocate(V_mean(1:N_y-1));   V_mean=0.0
  allocate(UU_mean(1:N_y-1));  UU_mean=0.0
  allocate(VV_mean(1:N_y-1));  VV_mean=0.0
  allocate(UV_mean(1:N_y-1));  UV_mean=0.0

  allocate(Ncount(1:N_y-1));   Ncount = 0.0
!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
  do i = 1, N_x ! for each layer formed by x1(i),x2(i)

    do j = 1, N_y-1 ! define y mesh
      do c=1,NC ! and look in every cell in work domain catched into x1(i),x2(i),y(j),y(j+1) cell irresponding to z
        if(x1(i) < xc(c) .and. xc(c) < x2(i)) then
          if (y(j) < yc(c) .and. yc(c) < y(j+1)) then !gathering statistics
            !1st order means
            U_mean(j)   = U_mean(j) +  U % mean(c)
            V_mean(j)   = V_mean(j) +  V % mean(c)
            !2nd order covarisnces
            UU_mean(j)   = UU_mean(j) + ( UU % mean(c) - U % mean(c) * U % mean(c) )
            VV_mean(j)   = VV_mean(j) + ( VV % mean(c) - V % mean(c) * V % mean(c) )
            UV_mean(j)   = UV_mean(j) + ( UV % mean(c) - U % mean(c) * V % mean(c) )
            !count catched in x1(i),x2(i),y(j),y(j+1) cell
            Ncount(j) = Ncount(j) + 1
          end if
        end if
      end do 
    end do 
!still inside x loop
!---- average over all processors
    do j=1, N_y-1
      call iGlSum(Ncount(j))

      call GloSum(U_mean(j))
      call GloSum(V_mean(j))
      call GloSum(UU_mean(j))
      call GloSum(VV_mean(j))
      call GloSum(UV_mean(j))
    end do

    !creating files
    nameIN='x=x.xxx'
    write(nameIN(3:7),'(F5.3)') (x1(i)+x2(i))/2
    write(nameIN(len_trim(nameIN)+1:len_trim(nameIN)+5),'(A)') '.dat'

    open(996, FILE=nameIN)

    do j=1, N_y-1 !writing the results
      call wait
      if(Ncount(j) > 0) then ! if there no cells catched there is no point to create file
                            U_mean(j)  = U_mean(j)  / Ncount(j)
                            V_mean(j)  = V_mean(j)  / Ncount(j)
                            UU_mean(j) = UU_mean(j) / Ncount(j)
                            VV_mean(j) = VV_mean(j) / Ncount(j) 
                            UV_mean(j) = UV_mean(j) / Ncount(j) 
      write(996,'(6E15.7,I4)') (y(j)+y(j+1))/2,           &
                            U_mean(j),      &
                            V_mean(j),      & 
                            UU_mean(j),     &
                            VV_mean(j),     &
                            UV_mean(j),     &
                            Ncount(j)
      end if ! if there no cells catched there is no point to create file

      U_mean(j)  = 0
      V_mean(j)  = 0
      UU_mean(j) = 0
      VV_mean(j) = 0
      UV_mean(j) = 0
      Ncount(j)  = 0

    end do
  close(996)
  end do ! end x loop (i=1, N_x)
 

  deallocate(x1)
  deallocate(x2)
  deallocate(y)

  deallocate(Ncount)

  deallocate(U_mean)
  deallocate(V_mean)
  deallocate(UU_mean)
  deallocate(VV_mean)
  deallocate(UV_mean)
  


  if(this < 2) write(*,*) 'Finished with ExtractFunc'

  RETURN
  END SUBROUTINE ExtractFunc
