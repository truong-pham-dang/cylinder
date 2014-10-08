!======================================================================!
  SUBROUTINE UserCutLines_annulus() 
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
    INTEGER             :: Nprob, pl, c, i, count, kk
    CHARACTER           :: namCoo*80, namPro*80, answer*80, namRes*80
    CHARACTER           :: namRes_plus*80
    REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:), ind(:),       & 
                                 var_1(:), var_2(:),     &
                                 var_3(:), Wall_p(:), &
                                 Ufric_p(:)
    INTEGER,ALLOCATABLE :: Np(:), Ncount(:)
    REAL                :: R, Urad_mean, Utan_mean, dummy, Lscale, R_max
    REAL                :: b11, b22, b12, b21
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $  
!======================================================================!

    namPro = name

    namRes = name
    namRes(len_trim(name)+1:len_trim(name)+8) = '_res.dat'
    namRes_plus = name
    namRes_plus(len_trim(name)+1:len_trim(name)+13) = '_res_plus.dat'
!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
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

    call SSORT (z_p, ind, Nprob, 0)
    
    call CalcShear(U % n, V % n, W % n, Shear)

    allocate(Np(Nprob));    Np=0 
    allocate(Wall_p(Nprob));Wall_p=0.0
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
    allocate(var_1(Nprob));  var_1=0.0
    allocate(var_2(Nprob));  var_2=0.0
    allocate(Ufric_p(Nprob));  Ufric_p=0.0

    allocate(Ncount(Nprob)); Ncount=0
    count = 0

    if(HOT==YES) then
      allocate(Tmp(Nprob));   Tmp=0.0
      allocate(TTp(Nprob));   TTp=0.0
      allocate(uTp(Nprob));   uTp=0.0
      allocate(vTp(Nprob));   vTp=0.0
      allocate(wTp(Nprob));   wTp=0.0
    end if  
    kk = 0

    Lscale = 0.0 
!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
!    R_max = 0.0

!    do i = 1, Nprob
!      R_max = max(z_p(i),R_max)
!    end do
!    open(12, FILE='inflow.dat') 
!    do c=1, NC
!      if(zc(c) > 0.0.and.zc(c) < 0.16) then
!        write(12,'(13E15.7)') xc(c), yc(c), U%n(c),V%n(c),W%n(c),uu%n(c),vv%n(c),ww%n(c),uv%n(c),uw%n(c),vw%n(c),Eps%n(c),f22%n(c)
!      end if
!    end do
!    close(12)
    do i = 1, Nprob-1
      do c=1, NC
        Lscale = max(Lscale,WallDs(c))
        R = sqrt(xc(c)**2+yc(c)**2)
        if(R > abs(z_p(i+1)) .and. R < abs(z_p(i))) then
        if(zc(c) > 14.5.and.zc(c) < 15.0) then
          Wall_p(i) = Wall_p(i) + R 
            b11 = xc(c)/R
            b12 = yc(c)/R
            b21 = -yc(c)/R
            b22 = xc(c)/R
          if(SIMULA==DNS) then

            Ump(i)   = Ump(i) + U % n(c) * xc(c) / R  + V % n(c) * yc(c) / R
            Vmp(i)   = Vmp(i) + (-U % n(c) * yc(c) / R  + V % n(c) * xc(c) / R)
!           Ump(i)   = Ump(i) + U % n(c) 
!           Vmp(i)   = Vmp(i) + V % n(c) 
            Wmp(i)   = Wmp(i) + W % n(c)
            uup(i)   = uup(i) + b11*b11*uu%n(c) + b11*b12*uv%n(c) &
                              + b12*b11*uv%n(c) + b12*b12*vv%n(c) 
            vvp(i)   = vvp(i) + b21*b21*uu%n(c) + b21*b22*uv%n(c) &
                              + b22*b21*uv%n(c) + b22*b22*vv%n(c)
            wwp(i)   = wwp(i) + ww % n(c) 
            uvp(i)   = uvp(i) + b11*b21*uu%n(c) + b11*b22*uv%n(c) &
                              + b12*b21*uv%n(c) + b12*b22*vv%n(c)
            uwp(i)   = uwp(i) + b11*uw%n(c) + b12*vw%n(c) 
            vwp(i)   = vwp(i) + b21*uw%n(c) + b22*vw%n(c) 
            var_1(i)   = var_1(i) + Kin % n(c) 
            var_2(i)   = var_2(i) + Eps % n(c) 
            if(IsNearWall(c)) then
              Ufric_p(i) = Ufric_p(i) + (VISc * (U % n(c)**2 + V % n(c)**2 + W % n(c)**2)**0.5/WallDs(c))**0.5
            end if
          else if(SIMULA==K_EPS_VV) then
            Ump(i)   = Ump(i) + U % n(c) * xc(c) / R  + V % n(c) * yc(c) / R
            Vmp(i)   = Vmp(i) + (-U % n(c) * yc(c) / R  + V % n(c) * xc(c) / R)
            Wmp(i)   = Wmp(i) + W % n(c)
            uup(i)   = uup(i) + Kin % n(c) 
            vvp(i)   = vvp(i) + Eps % n(c) 
            wwp(i)   = wwp(i) + v_2 % n(c) 
            uvp(i)   = uvp(i) + f22 % n(c) 
            uwp(i)   = uwp(i) + 0.0 
            if(IsNearWall(c)) then
              Ufric_p(i) = Ufric_p(i) + (VISc * (U % n(c)**2 + V % n(c)**2 + W % n(c)**2)**0.5/WallDs(c))**0.5
            end if
          else if(SIMULA==ZETA) then
            Ump(i)   = Ump(i) + U % n(c) * xc(c) / R  + V % n(c) * yc(c) / R
            Vmp(i)   = Vmp(i) + (-U % n(c) * yc(c) / R  + V % n(c) * xc(c) / R)
            Wmp(i)   = Wmp(i) + W % n(c)
            uup(i)   = uup(i) + Kin % n(c) 
            vvp(i)   = vvp(i) + Eps % n(c) 
            wwp(i)   = wwp(i) + v_2 % n(c) 
            uvp(i)   = uvp(i) + VISt(c)*(b11*Wx(c) + b12*Wz(c))
            uwp(i)   = uwp(i) + VISt(c) 
            if(IsNearWall(c)) then
              Ufric_p(i) = Ufric_p(i) + (VISc * (U % n(c)**2 + V % n(c)**2 + W % n(c)**2)**0.5/WallDs(c))**0.5
            end if
          else if(SIMULA==K_EPS) then
            Ump(i)   = Ump(i) + U % n(c) * xc(c) / R  + V % n(c) * yc(c) / R
            Vmp(i)   = Vmp(i) + (-U % n(c) * yc(c) / R  + V % n(c) * xc(c) / R)
            Wmp(i)   = Wmp(i) + W % n(c)
            uup(i)   = uup(i) + Kin % n(c) 
            vvp(i)   = vvp(i) + Eps % n(c) 
            wwp(i)   = wwp(i) + 0.0
            uvp(i)   = uvp(i) + 0.0 
            uwp(i)   = uwp(i) + 0.0
!            if(MODE == WF) then
!              if(IsNearWall(c)) then
!                Ufric_p(i) = Ufric_p(i) + Cmu**0.25*Kin % n(c)**0.5
!              end if
!            else
              if(IsNearWall(c)) then
                Ufric_p(i) = Ufric_p(i) + (VISc * (U % n(c)**2 + V % n(c)**2 + W % n(c)**2)**0.5/WallDs(c))**0.5
              end if
!            end if
          else if(SIMULA==SPA_ALL) then
            Ump(i)   = Ump(i) + U % n(c) * xc(c) / R  + V % n(c) * yc(c) / R
            Vmp(i)   = Vmp(i) + (-U % n(c) * yc(c) / R  + V % n(c) * xc(c) / R)
            Wmp(i)   = Wmp(i) + W % n(c)
            uup(i)   = uup(i) + VISt(c)
            vvp(i)   = vvp(i) + 0.0 
            wwp(i)   = wwp(i) + 0.0
            uvp(i)   = uvp(i) + 0.0 
            uwp(i)   = uwp(i) + 0.0 
            if(IsNearWall(c)) then
              Ufric_p(i) = Ufric_p(i) + (VISc * (U % n(c)**2 + V % n(c)**2 + W % n(c)**2)**0.5/WallDs(c))**0.5
            end if
          else if(SIMULA==LES.or.SIMULA==DES_SPA.or.SIMULA == DNS) then
            Ump(i)   = Ump(i) + U % mean(c) * xc(c) / R  + V % mean(c) * yc(c) / R
            Vmp(i)   = Vmp(i) + (-U % mean(c) * yc(c) / R  + V % mean(c) * xc(c) / R)
            Wmp(i)   = Wmp(i) + W % mean(c)
            uup(i)   = uup(i) + (uu % mean(c)- (U % mean(c) * xc(c) / R  + V % mean(c) * yc(c) / R)**2.0)
            vvp(i)   = vvp(i) + (vv % mean(c)- (-U % mean(c) * yc(c) / R  + V % mean(c) * xc(c) / R)**2.0)
            wwp(i)   = wwp(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
            uvp(i)   = uvp(i) + (uv % mean(c)- (U % mean(c) * xc(c) / R  + V % mean(c) * yc(c) / R)* &
                                (-U % mean(c) * yc(c) / R  + V % mean(c) * xc(c) / R)   )
            uwp(i)   = uwp(i) + (uw % mean(c)- (U % mean(c) * xc(c) / R  + V % mean(c) * yc(c) / R) * W % mean(c))
            vwp(i)   = uwp(i) + (vw % mean(c)- (-U % mean(c) * yc(c) / R  + V % mean(c) * xc(c) / R) * W % mean(c))
            if(IsNearWall(c)) then
              Ufric_p(i) = Ufric_p(i) + (VISc * (U % mean(c)**2 + V % mean(c)**2 + W % mean(c)**2)**0.5/WallDs(c))**0.5
            end if
          end if 
 
          if(HOT==YES) then
            if(SIMULA == LES.or.SIMULA == DES_SPA.or.SIMULA == DNS) then
              Tmp(i)   = Tmp(i) + T % mean(c)
              TTp(i)   = TTp(i) + (TT % mean(c) - T % mean(c) * T % mean(c))
              uTp(i)   = uTp(i) + (uT % mean(c) - (U % mean(c) * xc(c) / R  + V % mean(c) * yc(c) / R) * T % mean(c))
              vTp(i)   = vTp(i) + (vT % mean(c) - (-U % mean(c) * yc(c) / R  + V % mean(c) * xc(c) / R) * T % mean(c))
              wTp(i)   = wTp(i) + (wT % mean(c) - w % mean(c) * T % mean(c))
            else
              Tmp(i)   = Tmp(i) + T % n(c)
            end if
          end if
          Ncount(i) = Ncount(i) + 1
        end if
        end if
      end do 
    end do 

!---- average over all processors
    do pl=1, Nprob-1
      call IGlSum(Ncount(pl))

      call GloSum(Wall_p(pl))

      call GloSum(Ump(pl))
      call GloSum(Vmp(pl))
      call GloSum(Wmp(pl))

      call GloSum(uup(pl))
      call GloSum(vvp(pl))
      call GloSum(wwp(pl))

      call GloSum(uvp(pl))
      call GloSum(uwp(pl))
      call GloSum(vwp(pl))
      call GloSum(var_1(pl))
      call GloSum(var_2(pl))
      call GloSum(Ufric_p(pl))

      count =  count + Ncount(pl) 

      if(HOT==YES) then
        call GloSum(Tmp(pl))
        call GloSum(TTp(pl))
        call GloSum(uTp(pl))
        call GloSum(vTp(pl))
        call GloSum(wTp(pl))
      end if
    end do

    do i = 1, Nprob-1
      if(Ncount(i) /= 0) then
        Wall_p(i) =  Wall_p(i)/Ncount(i)
        Ump(i)    =  Ump(i)/Ncount(i)
        Vmp(i)    =  Vmp(i)/Ncount(i)
        Wmp(i)    =  Wmp(i)/Ncount(i)
        uup(i)    =  uup(i)/Ncount(i)
        vvp(i)    =  vvp(i)/Ncount(i)
        wwp(i)    =  wwp(i)/Ncount(i)
        uvp(i)    =  uvp(i)/Ncount(i)
        uwp(i)    =  uwp(i)/Ncount(i)
        vwp(i)    =  vwp(i)/Ncount(i)
        var_1(i)  =  var_1(i)/Ncount(i)
        var_2(i)  =  var_2(i)/Ncount(i)
        Ufric_p(i)  =  Ufric_p(i)/Ncount(i)
        if(HOT == YES) then
          Tmp(i)    =  Tmp(i)/Ncount(i)
          TTp(i)    =  TTp(i)/Ncount(i)
          uTp(i)    =  uTp(i)/Ncount(i)
          vTp(i)    =  vTp(i)/Ncount(i)
          wTp(i)    =  wTp(i)/Ncount(i)
        end if
      end if
    end do 

    Ufric = 0.0
    do i = 1, Nprob
      Ufric = max(Ufric_p(i),Ufric)       
    end do

10  continue
 
    open(3,FILE=namRes)
    write(3,'(A1,2(A10, F13.5))') '#', 'Utau = ', Ufric, 'Re_tau = ', Ufric*R_max/VISc
    if(SIMULA == LES.or.SIMULA == DES_SPA.or.SIMULA == DNS) then
      write(3,'(A1,2X,A80)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8:uv, 9:uw, 10:vw, 11:Kin' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(11E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i),   &
          uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
          vwp(i), 0.5*(uup(i) + vvp(i) + wwp(i))
        end if
      end do 
    else if(SIMULA == DNS) then
      write(3,'(A1,1X,A50)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:Kin, 6:Eps, 7:zeta, 8:f' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(12E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i),   &
          uup(i), vvp(i), wwp(i), uvp(i), uwp(i), vwp(i), var_1(i), var_2(i)
        end if
      end do 
    else if(SIMULA == ZETA) then
      write(3,'(A1,1X,A50)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:Kin, 6:Eps, 7:zeta, 8:f' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(9E15.7)') (Wall_p(i)-0.5428)/0.4615, abs(Ump(i)), abs(Vmp(i))/2.8, abs(Wmp(i)),   &
          uup(i), vvp(i), wwp(i), uvp(i), uwp(i)
        end if
      end do 
    else if(SIMULA == K_EPS_VV) then
      write(3,'(A1,1X,A50)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:Kin, 6:Eps, 7:v2, 8:f  ' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
!          write(3,'(9E15.7)') Wall_p(i), abs(Ump(i)), abs(Vmp(i)), abs(Wmp(i)),   &
          write(3,'(9E15.7)') 1.0 - Wall_p(i), abs(Ump(i)), abs(Vmp(i)), -abs(Wmp(i)), 20.0,  &
          uup(i), vvp(i), wwp(i), uvp(i) !, uwp(i)
        end if
      end do 
    else if(SIMULA == K_EPS) then
      write(3,'(A1,1X,A50)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:Kin, 6:Eps             ' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(9E15.7)') Wall_p(i), abs(Ump(i)), abs(Vmp(i)), abs(Wmp(i)),   &
          uup(i), vvp(i), wwp(i), uvp(i), uwp(i)
        end if
      end do 
    else if(SIMULA == SPA_ALL) then
      write(3,'(A1,1X,A50)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:VISt                   ' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(9E15.7)') Wall_p(i), abs(Ump(i)), abs(Vmp(i)), abs(Wmp(i)),   &
          uup(i), vvp(i), wwp(i), uvp(i), uwp(i)
        end if
      end do 
    end if
    close(3)

    open(3,FILE=namRes_plus)
    write(3,'(A1,2(A10, F13.5))') '#', 'Utau = ', Ufric, 'Re_tau = ', Ufric*R_max/VISc
    if(SIMULA == LES.or.SIMULA == DES_SPA.or.SIMULA==DNS) then
      write(3,'(A1,2X,A80)') '#', '1:Xrad+, 2:U+, 3:V+, 4:W+, 5:uu+, 6:vv+, 7:ww+, 8:uv+, 9:uw+, 10:vw+, 11:Kin+' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(11E15.7)') Wall_p(i)*Ufric/VISc, abs(Ump(i))/Ufric, abs(Vmp(i))/Ufric, abs(Wmp(i))/Ufric,   &
          uup(i)/Ufric**2.0, vvp(i)/Ufric**2.0, wwp(i)/Ufric**2.0, uvp(i)/Ufric**2.0, uwp(i)/Ufric**2.0,      &
          vwp(i)/Ufric**2.0, (uup(i) + vvp(i) + wwp(i))/Ufric**2.0
        end if
      end do 
    else if(SIMULA == DNS) then
      write(3,'(A1,1X,A80)') '#', '1:Xrad+, 2:U+, 3:V+, 4:W+, 5:Kin+, 6:Eps+, 7:zeta+, 8:f+                     ' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(11E15.7)') (Wall_p(i))*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric,   &
          uup(i)/Ufric**2.0, vvp(i)/Ufric**2.0, wwp(i)/Ufric**2.0, uvp(i)/Ufric**2.0, uwp(i)/Ufric**2.0, &
          vwp(i)/Ufric**2.0, 0.5*(uup(i)/Ufric**2.0+vvp(i)/Ufric**2.0+wwp(i)/Ufric**2.0)
        end if
      end do 
    else if(SIMULA == ZETA) then
      write(3,'(A1,1X,A80)') '#', '1:Xrad+, 2:U+, 3:V+, 4:W+, 5:Kin+, 6:Eps+, 7:zeta+, 8:f+                     ' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(9E15.7)') (Wall_p(i))*Ufric/VISc, abs(Ump(i))/Ufric, abs(Vmp(i))/Ufric, abs(Wmp(i))/Ufric,   &
          uup(i)/Ufric**2.0, vvp(i)*VISc/Ufric**4.0, wwp(i), uvp(i)/Ufric**2.0, uwp(i)/VISc
        end if
      end do 
    else if(SIMULA == K_EPS_VV) then
      write(3,'(A1,1X,A80)') '#', '1:Xrad+, 2:U+, 3:V+, 4:W+, 5:Kin+, 6:Eps+, 7:v2+, 8:f+                       ' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(9E15.7)') (Wall_p(i))*Ufric/VISc, abs(Ump(i))/Ufric, abs(Vmp(i))/Ufric, abs(Wmp(i))/Ufric,   &
          uup(i)/Ufric**2.0, vvp(i)*VISc/Ufric**4.0, wwp(i)/Ufric**2.0, uvp(i)*VISc/Ufric**2.0, uwp(i)/Ufric**2.0
        end if
      end do 
    else if(SIMULA == SPA_ALL) then
      write(3,'(A1,1X,A80)') '#', '1:Xrad+, 2:U+, 3:V+, 4:W+, 5:VISt/VISc                                       ' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(9E15.7)') (Wall_p(i))*Ufric/VISc, abs(Ump(i))/Ufric, abs(Vmp(i))/Ufric, abs(Wmp(i))/Ufric,   &
          uup(i)/Ufric**2.0, vvp(i)*VISc/Ufric**4.0, wwp(i), uvp(i)*VISc/Ufric**2.0, uwp(i)/Ufric**2.0
        end if
      end do 
    end if
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
    deallocate(var_1)
    deallocate(var_2)
    deallocate(Ksgsp)
    deallocate(Ufric_p)
    deallocate(Wall_p)
    if(HOT==YES) then
      deallocate(Tmp)
      deallocate(TTp)
      deallocate(uTp)
      deallocate(vTp)
      deallocate(wTp)
    end if
  END SUBROUTINE UserCutLines_annulus
