!======================================================================!
  SUBROUTINE UserProbe1D_jet_rans(namAut) 
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
  INTEGER             :: Nprob, pl, c, dummy, i, count, k, c1, c2, s
  CHARACTER           :: namCoo*80, namPro*80, answer*80, JetIn*9, namOut*80
  REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:), ind(:),               & 
                                 var_1(:), var_2(:), var_3(:), Rad_mp(:), &
                                 var_4(:), var_5(:)  
  INTEGER,ALLOCATABLE :: Np(:), Ncount(:)
  REAL                :: R, Urad_mean, Utan_mean, R1, R2, Urad, Utan, pi
  CHARACTER, OPTIONAL :: namAut*(*)
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $  
!======================================================================!
    pi = 3.141592
!    call Scale()
    namCoo = name
    namCoo(len_trim(name)+1:len_trim(name)+3) = '.1D'
    if(this < 2) write(6, *) '# Now reading the file:', namCoo
    open(9, FILE=namCoo)

!---- write the number of searching intervals
    read(9,*) Nprob
    allocate(z_p(Nprob*2))
    allocate(ind(Nprob*2))
!---- read the intervals positions
    do pl=1,Nprob
      read(9,*) ind(pl), z_p(pl)
    end do
    close(9)
    call SSORT (z_p, ind, Nprob,0)

    allocate(Np(Nprob));    Np=0 
    allocate(Ump(Nprob));   Ump=0.0
    allocate(Vmp(Nprob));   Vmp=0.0
    allocate(Wmp(Nprob));   Wmp=0.0
    allocate(uup(Nprob));   uup=0.0
    allocate(vvp(Nprob));   vvp=0.0
    allocate(wwp(Nprob));   wwp=0.0
    allocate(uvp(Nprob));   uvp=0.0
    allocate(Rad_mp(Nprob));  Rad_mp=0.0

    allocate(Ncount(Nprob)); Ncount=0
    count = 0

    if(HOT==YES) then
      allocate(Tmp(Nprob));   Tmp=0.0
    end if  

!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!

    do c =  1, NC
      R  = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 + tiny
      PP % n(c)   = (U % n(c) * xc(c) / R  + V % n(c) * yc(c) / R)
    end do
    call GraPhi(PP % n, 1, Ux,.TRUE.)    ! dU/dx
    call GraPhi(PP % n, 2, Uy,.TRUE.)    ! dU/dx
    call GraPhi(W % n, 3, Wz,.TRUE.)    ! dU/dx

    call GradP(P % n,Px,Py,Pz)

   call GraPhi(W % n, 1, Wx,.TRUE.)    ! dU/dx
   call GraPhi(W % n, 2, Wy,.TRUE.)    ! dU/dx
   call CalcShear(U % n, V % n, W % n, Shear)
    do c =  1, NC
      R  = (xc(c)*xc(c) + yc(c)*yc(c))**0.5
      PP % n(c)   = Wx(c) * xc(c) / R  + Wy(c) * yc(c) / R
      PP % n(c)   = (Ux(c) * xc(c) / R + Uy(c) * yc(c) / R)
      PP % n(c)   = -2.0*VISt(c)*PP % n(c) + 2.0/3.0 * Kin % n(c) &
                    -2.0*VISt(c)*Wz(c) + 2.0/3.0 * Kin % n(c)
    end do

  do k = 0, 6
    if(k == 0) then
      R1 = 0.04   
      R2 = 0.0
    else if(k == 1) then
      R1 = 1.0    
      R2 = 0.992  
!      R2 = 0.9986  
    else if(k == 2) then
      R1 = 2.1500 
      R2 = 2.0 
    else if(k == 3) then
      R1 = 3.0684
      R2 = 2.9744
    else if(k == 4) then
      R1 = 4.1433 
      R2 = 3.9098
    else if(k == 5) then
      R1 = 0.5047000E+01 
      R2 = 0.4933200E+01 
    else if(k == 6) then
      R1 = 0.6000000E+01
      R2 = 0.5876600E+01
!    else if(k == 7) then
!      R1 = 1.5281 
!      R2 = 1.4861
!    else if(k == 8) then
!      R1 = 2.5181 
!      R2 = 2.4296
    end if  


    do i = 1, Nprob-1
      do c=1,NC
        Rad_2 = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 + tiny
        if(Rad_2 < R1 .and. Rad_2 > R2) then
          if(zc(c) < z_p(i) .and. zc(c) > z_p(i+1)) then
            R           = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 + tiny
            Urad_mean   = (U % n(c) * xc(c) / R  + V % n(c) * yc(c) / R)
            Utan_mean   = (-U % n(c) * yc(c) / R  + V % n(c) * xc(c) / R) 
            Urad   = (U % n(c) * xc(c) / R  + V % n(c) * yc(c) / R)
            Utan   = (-U % n(c) * yc(c) / R  + V % n(c) * xc(c) / R) 
 
            Ump(i)   = Ump(i) + Urad_mean 
            if(k == 0) then
              Vmp(i)   = Vmp(i) + (U % n(c)**2.0 + V % n(c)** 2.0 + W % n(c)**2.0 + PP % n(c))**0.5 
            else 
              Vmp(i)   = Vmp(i) + (U % n(c)**2.0 + V % n(c)** 2.0 + W % n(c)**2.0)**0.5 
            end if   
            Wmp(i)   = Wmp(i) + W % n(c)
            uuP(i)   = uup(i) + Kin % n(c) 
            vvp(i)   = vvp(i) + Eps % n(c)
            wwp(i)   = wwp(i) + VISt(c)
            uvp(i)   = uvp(i) + & !+ max(0.07,0.22*(sin(pi*0.5*(1.0/max(zc(c),1.0))**1.5)))
                       T % n(c) !max(0.07,0.22*(sin(pi*0.5*(1.0/max(zc(c),1.0))**1.5))) * v_2%n(c)*Kin % n(c)*Tsc(c)


!            uuP(i)   = uup(i) + W % C(c) 
!            vvp(i)   = vvp(i) + W % Do(c)
!            wwp(i)   = wwp(i) - Pz(c) * volume(c)
!            uvp(i)   = uvp(i) + (W % C(c) + W % Do(c) - Pz(c) * volume(c))

!            uuP(i)   = uup(i) + PP % n(c)
!            vvp(i)   = vvp(i) + Wz(c)
!            wwp(i)   = wwp(i) + Shear(c) 
!            uvp(i)   = uvp(i) + (W % C(c) + W % Do(c) - Pz(c) * volume(c))

            if(HOT==YES) then
              Tmp(i)   = Tmp(i) + T % n(c)
            end if
       
            Rad_mp(i) = Rad_mp(i) + (xc(c)*xc(c) + yc(c)*yc(c))**0.5
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
      call GloSum(Rad_mp(pl))

      call GloSum(uup(pl))
      call GloSum(vvp(pl))
      call GloSum(wwp(pl))

      call GloSum(uvp(pl))

      count =  count + Ncount(pl) 

      if(HOT==YES) then
        call GloSum(Tmp(pl))
      end if
    end do

    JetIn  = 'Jet_res_x'
    write(JetIn(9:9),'(I1)') k

    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        Wmp(i)    =  Wmp(i)/Ncount(i)
        Ump(i)    =  Ump(i)/Ncount(i)
        Vmp(i)    =  Vmp(i)/Ncount(i)
        uup(i)    =  uup(i)/Ncount(i)
        vvp(i)    =  vvp(i)/Ncount(i)
        wwp(i)    =  wwp(i)/Ncount(i)
        uvp(i)    =  uvp(i)/Ncount(i)
        Tmp(i)    =  Tmp(i)/Ncount(i)
        Rad_mp(i) =  Rad_mp(i)/Ncount(i)
      end if
    end do

    open(3,FILE=JetIn)
    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        write(3,'(8E15.7)') (z_p(i)+z_p(i+1))/4.0, Ump(i), Vmp(i), Wmp(i) , uup(i), vvp(i), wwp(i), uvp(i) 
      end if
    end do 
    close(3)

    do i = 1, Nprob
      Ncount(i) = 0
      Wmp(i)    = 0.0 
      Ump(i)    = 0.0 
      Vmp(i)    = 0.0 
      uup(i)    = 0.0 
      vvp(i)    = 0.0 
      wwp(i)    = 0.0 
      uvp(i)    = 0.0 
      Tmp(i)    = 0.0 
      Rad_mp(i) =  0.0
    end do
    if(this < 2) write(*,*) 'Finished with cut line  ', k
  end do   !end number of radius

    deallocate(Np)
    deallocate(z_p)
    deallocate(Ump)
    deallocate(Vmp)
    deallocate(Wmp)
    deallocate(uup)
    deallocate(vvp)
    deallocate(wwp)
    deallocate(uvp)
    deallocate(Rad_mp)
    deallocate(Ncount)
    if(HOT==YES) then
      deallocate(Tmp)
    end if


  if(this < 2) write(*,*) 'Finished with UserProbe1D_jet '

  END SUBROUTINE UserProbe1D_jet_rans
