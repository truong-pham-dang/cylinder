!======================================================================!
  SUBROUTINE UserHill_Nu() 
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
  REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), wall_p(:), & 
                                 Tmp(:) 
  INTEGER,ALLOCATABLE :: Np(:), Ncount(:)
  REAL                :: R, dummy, Nu, Umag, H, Utau_p
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $  
!======================================================================!


!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
    if(this < 2) write(6, *) '# Now reading the file: x_dist.1D ' 
    open(9, FILE='x_dist.1D')

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
    allocate(Vmp(Nprob));   Vmp=0.0
    allocate(Wmp(Nprob));   Wmp=0.0
    allocate(wall_p(Nprob));wall_p=0.0

    allocate(Ncount(Nprob)); Ncount=0
    count = 0

    if(HOT==YES) then
      allocate(Tmp(Nprob));   Tmp=0.0
    end if  

!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
    do i = 1, Nprob-1
      do c=-NbC, NC
        if(IsNearWall(c)) then
          if(xc(c) > z_p(i) .and. xc(c) < z_p(i+1)) then
            Ump(i)   = Ump(i) + U % mean(c) 
            Vmp(i)   = Vmp(i) + V % mean(c)
            Wmp(i)   = Wmp(i) + W % mean(c)
            wall_p(i)= WallDs(c)
          end if
        end if
        if(TypeBC(c) == WALLFL) then
          if(xc(c) > z_p(i) .and. xc(c) < z_p(i+1)) then
            if(HOT==YES) then
              Tmp(i)   = Tmp(i) + T % mean(c)
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
      call GloSum(wall_p(pl))

      count =  count + Ncount(pl) 

      if(HOT==YES) then
        call GloSum(Tmp(pl))
      end if
    end do

    do i = 1, Nprob
      Ncount(i) = max(Ncount(i),1) 
    end do
    open(9, FILE='Hill_Nu.dat')
    do i = 1, Nprob-1
      Wmp(i)    =  Wmp(i)/Ncount(i)
      Ump(i)    =  Ump(i)/Ncount(i)
      Vmp(i)    =  Vmp(i)/Ncount(i)
      Tmp(i)    =  Tmp(i)/Ncount(i)
      wall_p(i) =  wall_p(i)/Ncount(i)
      Umag = sqrt( Ump(i)*Ump(i) + Vmp(i)*Vmp(i) + Wmp(i)*Wmp(i))
      Utau_p = sqrt(VISc*Umag/wall_p(i))  
      H    = 1.0
      Nu   = H * 0.1/(0.0000255*(Tmp(i) - 20.0 + tiny))  
      write(9,'(3E15.7,1X,I6)') (z_p(i)+z_p(i+1))/2.0, Utau_p, Nu, &
                                 Ncount(i) 
    end do 
    close(3)

    deallocate(Np)
    deallocate(z_p)
    deallocate(Ump)
    deallocate(Vmp)
    deallocate(Wmp)
    deallocate(wall_p)
    deallocate(Ncount)
    if(HOT==YES) then
      deallocate(Tmp)
    end if

  if(this < 2) write(*,*) 'Finished with UserHill_Nu '

  END SUBROUTINE UserHill_Nu
