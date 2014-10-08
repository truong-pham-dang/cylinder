!======================================================================!
  SUBROUTINE UserCutLines_budgets_cylind_HOT() 
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
  REAL :: Rad_2, Ufric , R_max
!------------------------------[Calling]-------------------------------!
  INTERFACE
    LOGICAL FUNCTION Approx(A,B,tol)
      REAL           :: A,B
      REAL, OPTIONAL :: tol
    END FUNCTION Approx
  END INTERFACE 
!-------------------------------[Locals]-------------------------------!
  INTEGER             :: Nprob, pl, c, i, count, Ncount_wall
  CHARACTER           :: namCoo*80, namPro*80, answer*80
  REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                         uup(:), vvp(:), wwp(:), &
                         uvp(:), uwp(:), vwp(:), &
                         uuup(:), uuvp(:), uuwp(:), &
                         vvup(:), vvvp(:), vvwp(:), &
                         wwup(:), wwvp(:), wwwp(:), &
                         uwup(:), uwvp(:), uwwp(:), &
                         Tmp(:), TTp(:),         &
                         uTp(:), vTp(:), wTp(:), &
                         Ksgsp(:),               & 
                         var_1(:), var_2(:), var_3(:), Rad_1(:), Rad_mp(:), &
                         var_10(:), var_11(:), var_12(:), var_13(:), &
                         Puu(:), Pvv(:), Pww(:), Puw(:), &
                         Diffv_uu(:), Diffv_vv(:), Diffv_ww(:), Diffv_uw(:), &
                         Difft_uu(:), Difft_vv(:), Difft_ww(:), Difft_uw(:), &
                         Diss_uu(:), Diss_vv(:), Diss_ww(:), Diss_uw(:), Diss_sgs(:), &
                         Ruu(:), Rvv(:), Rww(:), Ruw(:), &
                         PDuu(:), PDvv(:), PDww(:), PDuw(:), &
                         Cuu(:), Cvv(:), Cww(:), Cuw(:), &
                         Put_p(:), Pvt_p(:), Pwt_p(:), Ptt(:), &
                         Diffv_ut(:), Diffv_vt(:), Diffv_wt(:), Diffv_tt(:), &
                         Difft_ut(:), Difft_vt(:), Difft_wt(:), Difft_tt(:), &
                         Diss_ut(:), Diss_vt(:), Diss_wt(:), Diss_tt(:), &
                         Rut(:), Rvt(:), Rwt(:), Rtt(:), &
                         PDut(:), PDvt(:), PDwt(:), PDtt(:), &
                         Cut(:), Cvt(:), Cwt(:), Ctt(:), Ufric_p(:)
    
  
  INTEGER,ALLOCATABLE :: Np(:), Ncount(:), int(:)
  REAL                :: R, Urad_mean, Utan_mean, NF, vol1, Tfric, Twall, Twall_p, FFF, dummy
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $  
!======================================================================!


!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
    namCoo = name
    namCoo(len_trim(name)+1:len_trim(name)+3) = '.1D'
    if(this < 2) write(6, *) '# Now reading the file:', namCoo
    open(9, FILE=namCoo)


!---- write the number of probes 
  read(9,*) Nprob
  allocate(z_p(Nprob))
  allocate(int(Nprob))

!---- write the probe coordinates out
  do pl=1,Nprob
    read(9,*) int(pl), z_p(pl) 
  end do
  close(9)

  call SSORT (z_p, int, Nprob, 0)
  call CalcShear(U % n, V % n, W % n, Shear)
  allocate(Ufric_p(Nprob));  Ufric_p=0.0
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
  allocate(Rad_mp(Nprob));  Rad_mp=0.0
  allocate(var_1(Nprob));  var_1=0.0
  allocate(var_10(Nprob));  var_10=0.0
  allocate(var_11(Nprob));  var_11=0.0
  allocate(var_12(Nprob));  var_12=0.0
  allocate(var_13(Nprob));  var_13=0.0

  allocate(uuup(Nprob));   uuup=0.0
  allocate(uuvp(Nprob));   uuvp=0.0
  allocate(uuwp(Nprob));   uuwp=0.0

  allocate(vvup(Nprob));   vvup=0.0
  allocate(vvvp(Nprob));   vvvp=0.0
  allocate(vvwp(Nprob));   vvwp=0.0

  allocate(wwup(Nprob));   wwup=0.0
  allocate(wwvp(Nprob));   wwvp=0.0
  allocate(wwwp(Nprob));   wwwp=0.0

  allocate(uwup(Nprob));   uwup=0.0
  allocate(uwvp(Nprob));   uwvp=0.0
  allocate(uwwp(Nprob));   uwwp=0.0

  allocate(Puu(Nprob));   Puu=0.0
  allocate(Pvv(Nprob));   Pvv=0.0
  allocate(Pww(Nprob));   Pww=0.0
  allocate(Puw(Nprob));   Puw=0.0

  allocate(Diss_uu(Nprob));   Diss_uu=0.0
  allocate(Diss_vv(Nprob));   Diss_vv=0.0
  allocate(Diss_ww(Nprob));   Diss_ww=0.0
  allocate(Diss_uw(Nprob));   Diss_uw=0.0

  allocate(Diss_sgs(Nprob));   Diss_sgs=0.0

  allocate(Diffv_uu(Nprob));   Diffv_uu=0.0
  allocate(Diffv_vv(Nprob));   Diffv_vv=0.0
  allocate(Diffv_ww(Nprob));   Diffv_ww=0.0
  allocate(Diffv_uw(Nprob));   Diffv_uw=0.0

  allocate(Difft_uu(Nprob));   Difft_uu=0.0
  allocate(Difft_vv(Nprob));   Difft_vv=0.0
  allocate(Difft_ww(Nprob));   Difft_ww=0.0
  allocate(Difft_uw(Nprob));   Difft_uw=0.0

  allocate(Ruu(Nprob));   Ruu=0.0
  allocate(Rvv(Nprob));   Rvv=0.0
  allocate(Rww(Nprob));   Rww=0.0
  allocate(Ruw(Nprob));   Ruw=0.0


  allocate(Cuu(Nprob));   Cuu=0.0
  allocate(Cvv(Nprob));   Cvv=0.0
  allocate(Cww(Nprob));   Cww=0.0
  allocate(Cuw(Nprob));   Cuw=0.0

  allocate(PDuu(Nprob));   PDuu=0.0
  allocate(PDvv(Nprob));   PDvv=0.0
  allocate(PDww(Nprob));   PDww=0.0
  allocate(PDuw(Nprob));   PDuw=0.0

  allocate(Rad_1(Nprob));  Rad_1=0.0
  allocate(Ncount(Nprob)); Ncount=0
  allocate(Put_p(Nprob));   Put_p=0.0
  allocate(Pvt_p(Nprob));   Pvt_p=0.0
  allocate(Pwt_p(Nprob));   Pwt_p=0.0
  allocate(Ptt(Nprob));   Ptt=0.0

  allocate(Diss_ut(Nprob));   Diss_ut=0.0
  allocate(Diss_vt(Nprob));   Diss_vt=0.0
  allocate(Diss_wt(Nprob));   Diss_wt=0.0
  allocate(Diss_tt(Nprob));   Diss_tt=0.0

  allocate(Diffv_ut(Nprob));   Diffv_ut=0.0
  allocate(Diffv_vt(Nprob));   Diffv_vt=0.0
  allocate(Diffv_wt(Nprob));   Diffv_wt=0.0
  allocate(Diffv_tt(Nprob));   Diffv_tt=0.0

  allocate(Difft_ut(Nprob));   Difft_ut=0.0
  allocate(Difft_vt(Nprob));   Difft_vt=0.0
  allocate(Difft_wt(Nprob));   Difft_wt=0.0
  allocate(Difft_tt(Nprob));   Difft_tt=0.0

  allocate(Rut(Nprob));   Rut=0.0
  allocate(Rvt(Nprob));   Rvt=0.0
  allocate(Rwt(Nprob));   Rwt=0.0

  allocate(Cut(Nprob));   Cut=0.0
  allocate(Cvt(Nprob));   Cvt=0.0
  allocate(Cwt(Nprob));   Cwt=0.0
  allocate(Ctt(Nprob));   Ctt=0.0

  allocate(PDut(Nprob));   PDut=0.0
  allocate(PDvt(Nprob));   PDvt=0.0
  allocate(PDwt(Nprob));   PDwt=0.0
  
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
    R_max = 0.0

    do i = 1, Nprob
      R_max = max(z_p(i),R_max)
    end do

    do i = 1, Nprob
      Rad_1(i) = R_max - WallDs(c)    
    end do 
    Ncount_wall = 0
    Twall_p     = 0.0
    do c = -Nbc,1
      if(TypeBC(c) == WALLFL) then  
        Twall_p = Twall_p + T % mean(c)
        Ncount_wall = Ncount_wall + 1
      end if
    end do
    call GloSum(Twall_p)
    call IGlSum(Ncount_wall)

    Twall =  Twall_p / Ncount_wall
    do i = 1, Nprob-1
      do c=1,NC
        Rad_2 = (xc(c)**2 + yc(c)**2)**0.5 
        if(Rad_2 > abs(z_p(i+1)) .and. Rad_2 < abs(z_p(i))) then
          R           = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 + tiny
          Urad_mean   = (U % mean(c) * xc(c) / R  + V % mean(c) * yc(c) / R)
          Utan_mean   = (-U % mean(c) * yc(c) / R  + V % mean(c) * xc(c) / R) 

          if(IsNearWall(c)) then
            Ufric_p(i) = Ufric_p(i) + (VISc * (U % mean(c)**2 + V % mean(c)**2 + W % mean(c)**2)**0.5/WallDs(c))**0.5
          end if
 
          Ump(i)   = Ump(i) + U % mean(c)
          Vmp(i)   = Vmp(i) + V % mean(c)
          Wmp(i)   = Wmp(i) + W % mean(c)
          uup(i)   = uup(i) + (uu % mean(c)- Urad_mean * Urad_mean)
          vvp(i)   = vvp(i) + (vv % mean(c)- Utan_mean * Utan_mean)
          wwp(i)   = wwp(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
          uvp(i)   = uvp(i) + (uv % mean(c)- Urad_mean * Utan_mean )
          uwp(i)   = uwp(i) + (uw % mean(c)- Urad_mean * W % mean(c))

          uuup(i)   = uuup(i) + uuu % mean(c)
          uuvp(i)   = uuvp(i) + uuv % mean(c)
          uuwp(i)   = uuwp(i) + uuw % mean(c)

          vvup(i)   = vvup(i) + vvu % mean(c)
          vvvp(i)   = vvvp(i) + vvv % mean(c)
          vvwp(i)   = vvwp(i) + vvw % mean(c)

          wwup(i)   = wwup(i) + wwu % mean(c)
          wwvp(i)   = wwvp(i) + wwv % mean(c)
          wwwp(i)   = wwwp(i) + www % mean(c)

          uwup(i)   = uwup(i) + uwu % mean(c)
          uwvp(i)   = uwvp(i) + uwv % mean(c)
          uwwp(i)   = uwwp(i) + uww % mean(c)

          var_1(i) = var_1(i) + TT % mean(c)


          var_10(i) = var_10(i) + uu % n(c)
          var_11(i) = var_11(i) + vv % n(c)
          var_12(i) = var_12(i) + ww % n(c)
          var_13(i) = var_13(i) + volume(c)

            Puu(i)    = Puu(i) + Puu_mean(c)
            Pvv(i)    = Pvv(i) + Pvv_mean(c)
            Pww(i)    = Pww(i) + Pww_mean(c)
            Puw(i)    = Puw(i) + Puw_mean(c)

            Cuu(i)    = Cuu(i) + C_uu_mean(c)
            Cvv(i)    = Cvv(i) + C_vv_mean(c)
            Cww(i)    = Cww(i) + C_ww_mean(c)
            Cuw(i)    = Cuw(i) + C_uw_mean(c)

            Ruu(i)    = Ruu(i) + PR_uu_mean(c)
            Rvv(i)    = Rvv(i) + PR_vv_mean(c)
            Rww(i)    = Rww(i) + PR_ww_mean(c)
            Ruw(i)    = Ruw(i) + PR_uw_mean(c)

            PDuu(i)    = PDuu(i) + PD_uu_mean(c)
            PDvv(i)    = PDvv(i) + PD_vv_mean(c)
            PDww(i)    = PDww(i) + PD_ww_mean(c)
            PDuw(i)    = PDuw(i) + PD_uw_mean(c)

            Difft_uu(i)  = Difft_uu(i) + Dift_uu_mean(c)
            Difft_vv(i)  = Difft_vv(i) + Dift_vv_mean(c)
            Difft_ww(i)  = Difft_ww(i) + Dift_ww_mean(c)
            Difft_uw(i)  = Difft_uw(i) + Dift_uw_mean(c)

            Diffv_uu(i)  = Diffv_uu(i) + Difv_uu_mean(c)
            Diffv_vv(i)  = Diffv_vv(i) + Difv_vv_mean(c)
            Diffv_ww(i)  = Diffv_ww(i) + Difv_ww_mean(c)
            Diffv_uw(i)  = Diffv_uw(i) + Difv_uw_mean(c)

            Diss_uu(i)  = Diss_uu(i) + Diss_uu_mean(c)
            Diss_vv(i)  = Diss_vv(i) + Diss_vv_mean(c)
            Diss_ww(i)  = Diss_ww(i) + Diss_ww_mean(c)
            Diss_uw(i)  = Diss_uw(i) + Diss_uw_mean(c)

            Diss_sgs(i)  = Diss_sgs(i) + Diss_sgs_mean(c)

            Put_p(i)    = Put_p(i) + Put_mean(c)
            Pvt_p(i)    = Pvt_p(i) + Pvt_mean(c)
            Pwt_p(i)    = Pwt_p(i) + Pwt_mean(c)
            Ptt(i)    = Ptt(i) + Ptt_mean(c)

            Cut(i)    = Cut(i) + C_ut_mean(c)
            Cvt(i)    = Cvt(i) + C_vt_mean(c)
            Cwt(i)    = Cwt(i) + C_wt_mean(c)
            Ctt(i)    = Ctt(i) + C_tt_mean(c)

            Rut(i)    = Rut(i) + PR_ut_mean(c)
            Rvt(i)    = Rvt(i) + PR_vt_mean(c)
            Rwt(i)    = Rwt(i) + PR_wt_mean(c)

            PDut(i)    = PDut(i) + PD_ut_mean(c)
            PDvt(i)    = PDvt(i) + PD_vt_mean(c)
            PDwt(i)    = PDwt(i) + PD_wt_mean(c)

            Difft_ut(i)  = Difft_ut(i) + Dift_ut_mean(c)
            Difft_vt(i)  = Difft_vt(i) + Dift_vt_mean(c)
            Difft_wt(i)  = Difft_wt(i) + Dift_wt_mean(c)
            Difft_tt(i)  = Difft_tt(i) + Dift_tt_mean(c)

            Diffv_ut(i)  = Diffv_ut(i) + Difv_ut_mean(c) + Difv_ut_tot(c)
            Diffv_vt(i)  = Diffv_vt(i) + Difv_vt_mean(c) + Difv_vt_tot(c)
            Diffv_wt(i)  = Diffv_wt(i) + Difv_wt_mean(c) + Difv_wt_tot(c)
            Diffv_tt(i)  = Diffv_tt(i) + Difv_tt_mean(c)

            Diss_ut(i)  = Diss_ut(i) + Diss_ut_mean(c)
            Diss_vt(i)  = Diss_vt(i) + Diss_vt_mean(c)
            Diss_wt(i)  = Diss_wt(i) + Diss_wt_mean(c)
            Diss_tt(i)  = Diss_tt(i) + Diss_tt_mean(c)

          if(HOT==YES) then
            Tmp(i)   = Tmp(i) + T % mean(c)
            TTp(i)   = TTp(i) + (TT % mean(c) - T % mean(c) * T % mean(c))
            uTp(i)   = uTp(i) + (uT % mean(c) - Urad_mean * T % mean(c))
            vTp(i)   = vTp(i) + (vT % mean(c) - Utan_mean * T % mean(c))
            wTp(i)   = wTp(i) + (wT % mean(c) - w % mean(c) * T % mean(c))
          end if
       
          Rad_mp(i) = Rad_mp(i) + WallDs(c)
          Ncount(i) = Ncount(i) + 1
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

    call GloSum(Ufric_p(pl))

    call GloSum(uup(pl))
    call GloSum(vvp(pl))
    call GloSum(wwp(pl))

    call GloSum(uvp(pl))
    call GloSum(uwp(pl))
    call GloSum(vwp(pl))
    call GloSum(var_1(pl))

    call GloSum(uuup(pl))
    call GloSum(uuvp(pl))
    call GloSum(uuwp(pl))

    call GloSum(vvup(pl))
    call GloSum(vvvp(pl))
    call GloSum(vvwp(pl))

    call GloSum(wwup(pl))
    call GloSum(wwvp(pl))
    call GloSum(wwwp(pl))

    call GloSum(uwup(pl))
    call GloSum(uwvp(pl))
    call GloSum(uwwp(pl))

    call GloSum(var_10(pl))
    call GloSum(var_11(pl))
    call GloSum(var_12(pl))
    call GloSum(var_13(pl))

      call GloSum(Puu(pl))
      call GloSum(Pvv(pl))
      call GloSum(Pww(pl))
      call GloSum(Puw(pl))

      call GloSum(PDuu(pl))
      call GloSum(PDvv(pl))
      call GloSum(PDww(pl))
      call GloSum(PDuw(pl))

      call GloSum(Cuu(pl))
      call GloSum(Cvv(pl))
      call GloSum(Cww(pl))
      call GloSum(Cuw(pl))

      call GloSum(Ruu(pl))
      call GloSum(Rvv(pl))
      call GloSum(Rww(pl))
      call GloSum(Ruw(pl))

      call GloSum(Difft_uu(pl))
      call GloSum(Difft_vv(pl))
      call GloSum(Difft_ww(pl))
      call GloSum(Difft_uw(pl))

      call GloSum(Diffv_uu(pl))
      call GloSum(Diffv_vv(pl))
      call GloSum(Diffv_ww(pl))
      call GloSum(Diffv_uw(pl))

      call GloSum(Diss_uu(pl))
      call GloSum(Diss_vv(pl))
      call GloSum(Diss_ww(pl))
      call GloSum(Diss_uw(pl))

      call GloSum(Diss_sgs(pl))

      call GloSum(Put_p(pl))
      call GloSum(Pvt_p(pl))
      call GloSum(Pwt_p(pl))
      call GloSum(Ptt(pl))

      call GloSum(Cut(pl))
      call GloSum(Cvt(pl))
      call GloSum(Cwt(pl))
      call GloSum(Ctt(pl))

      call GloSum(Rut(pl))
      call GloSum(Rvt(pl))
      call GloSum(Rwt(pl))

      call GloSum(PDut(pl))
      call GloSum(PDvt(pl))
      call GloSum(PDwt(pl))

      call GloSum(Difft_ut(pl))
      call GloSum(Difft_vt(pl))
      call GloSum(Difft_wt(pl))
      call GloSum(Difft_tt(pl))

      call GloSum(Diffv_ut(pl))
      call GloSum(Diffv_vt(pl))
      call GloSum(Diffv_wt(pl))
      call GloSum(Diffv_tt(pl))

      call GloSum(Diss_ut(pl))
      call GloSum(Diss_vt(pl))
      call GloSum(Diss_wt(pl))
      call GloSum(Diss_tt(pl))


    count =  count + Ncount(pl) 

    if(HOT==YES) then
      call GloSum(Tmp(pl))
      call GloSum(TTp(pl))
      call GloSum(uTp(pl))
      call GloSum(vTp(pl))
      call GloSum(wTp(pl))
    end if
  end do

    open(3,FILE='pipe_mean.dat')
    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        Ufric_p(i)  =  Ufric_p(i)/Ncount(i)
        Wmp(i)    =  Wmp(i)/Ncount(i)
        Ump(i)    =  Ump(i)/Ncount(i)
        Vmp(i)    =  Vmp(i)/Ncount(i)
        uup(i)    =  uup(i)/Ncount(i)
        vvp(i)    =  vvp(i)/Ncount(i)
        wwp(i)    =  wwp(i)/Ncount(i)
        uvp(i)    =  uvp(i)/Ncount(i)
        uwp(i)    =  uwp(i)/Ncount(i)

        uuup(i)    =  uuup(i)/Ncount(i)
        uuvp(i)    =  uuvp(i)/Ncount(i)
        uuwp(i)    =  uuwp(i)/Ncount(i)

        vvup(i)    =  vvup(i)/Ncount(i)
        vvvp(i)    =  vvvp(i)/Ncount(i)
        vvwp(i)    =  vvwp(i)/Ncount(i)

        wwup(i)    =  wwup(i)/Ncount(i)
        wwvp(i)    =  wwvp(i)/Ncount(i)
        wwwp(i)    =  wwwp(i)/Ncount(i)

        uwup(i)    =  uwup(i)/Ncount(i)
        uwvp(i)    =  uwvp(i)/Ncount(i)
        uwwp(i)    =  uwwp(i)/Ncount(i)

        var_1(i)  =  var_1(i)/Ncount(i)
        var_10(i)  =  var_10(i)/Ncount(i)
        var_11(i)  =  var_11(i)/Ncount(i)
        var_12(i)  =  var_12(i)/Ncount(i)
        var_13(i)  =  var_13(i)/Ncount(i)
      Puu(i)    = Puu(i)/Ncount(i)
      Pvv(i)    = Pvv(i)/Ncount(i)
      Pww(i)    = Pww(i)/Ncount(i)
      Puw(i)    = Puw(i)/Ncount(i)

      PDuu(i)    = PDuu(i)/Ncount(i)
      PDvv(i)    = PDvv(i)/Ncount(i)
      PDww(i)    = PDww(i)/Ncount(i)
      PDuw(i)    = PDuw(i)/Ncount(i)

      Cuu(i)    = Cuu(i)/Ncount(i)
      Cvv(i)    = Cvv(i)/Ncount(i)
      Cww(i)    = Cww(i)/Ncount(i)
      Cuw(i)    = Cuw(i)/Ncount(i)


      Ruu(i)    = Ruu(i)/Ncount(i)
      Rvv(i)    = Rvv(i)/Ncount(i)
      Rww(i)    = Rww(i)/Ncount(i)
      Ruw(i)    = Ruw(i)/Ncount(i)

      Difft_uu(i)    = Difft_uu(i)/Ncount(i)
      Difft_vv(i)    = Difft_vv(i)/Ncount(i)
      Difft_ww(i)    = Difft_ww(i)/Ncount(i)
      Difft_uw(i)    = Difft_uw(i)/Ncount(i)

      Diffv_uu(i)    = Diffv_uu(i)/Ncount(i)
      Diffv_vv(i)    = Diffv_vv(i)/Ncount(i)
      Diffv_ww(i)    = Diffv_ww(i)/Ncount(i)
      Diffv_uw(i)    = Diffv_uw(i)/Ncount(i)

      Diss_sgs(i)    = Diss_sgs(i)/Ncount(i)

      Diss_uu(i)    = Diss_uu(i)/Ncount(i)
      Diss_vv(i)    = Diss_vv(i)/Ncount(i)
      Diss_ww(i)    = Diss_ww(i)/Ncount(i)
      Diss_uw(i)    = Diss_uw(i)/Ncount(i)        

      Rad_mp(i) =  Rad_mp(i)/Ncount(i)
      Put_p(i)    = Put_p(i)/Ncount(i)
      Pvt_p(i)    = Pvt_p(i)/Ncount(i)
      Pwt_p(i)    = Pwt_p(i)/Ncount(i)
      Ptt(i)    = Ptt(i)/Ncount(i)
      Cut(i)    = Cut(i)/Ncount(i)
      Cvt(i)    = Cvt(i)/Ncount(i)
      Cwt(i)    = Cwt(i)/Ncount(i)
      Ctt(i)    = Ctt(i)/Ncount(i)
      PDut(i)    = PDut(i)/Ncount(i)
      PDvt(i)    = PDvt(i)/Ncount(i)
      PDwt(i)    = PDwt(i)/Ncount(i)
      Rut(i)    = Rut(i)/Ncount(i)
      Rvt(i)    = Rvt(i)/Ncount(i)
      Rwt(i)    = Rwt(i)/Ncount(i)
      Difft_ut(i)    = Difft_ut(i)/Ncount(i)
      Difft_vt(i)    = Difft_vt(i)/Ncount(i)
      Difft_wt(i)    = Difft_wt(i)/Ncount(i)
      Difft_tt(i)    = Difft_tt(i)/Ncount(i)
      Diffv_ut(i)    = Diffv_ut(i)/Ncount(i)
      Diffv_vt(i)    = Diffv_vt(i)/Ncount(i)
      Diffv_wt(i)    = Diffv_wt(i)/Ncount(i)
      Diffv_tt(i)    = Diffv_tt(i)/Ncount(i)
      Diss_ut(i)    = Diss_ut(i)/Ncount(i)
      Diss_vt(i)    = Diss_vt(i)/Ncount(i)
      Diss_wt(i)    = Diss_wt(i)/Ncount(i)
      Diss_tt(i)    = Diss_tt(i)/Ncount(i)
      Tmp(i)        = Tmp(i)/Ncount(i)  
      TTp(i)        = TTp(i)/Ncount(i)  
      uTp(i)        = uTp(i)/Ncount(i)  
      vTp(i)        = vTp(i)/Ncount(i)  
      wTp(i)        = wTp(i)/Ncount(i)  

        write(3,'(13E15.7,I8)') Rad_mp(i), Wmp(i) , uup(i), vvp(i), wwp(i), uvp(i), uwp(i),&
                               Tmp(i), TTp(i), uTp(i), vTp(i), wTp(i), var_1(i), Ncount(i)
      end if
    end do 
    close(3)

    Ufric = 0.0
    do i = 1, Nprob
      Ufric = max(Ufric_p(i),Ufric)
    end do


      Tfric = 0.005/Ufric
      if(this < 2)write(*,*) 'mean Utau = ', Ufric, 'mean Ttau = ', Tfric,        &
       'Twall_mean = ', Twall , 'ukupno nasao celija: ', count

    open(3,FILE='pipe_mean_plus.dat')
    write(3,*)'# Utau = ', Ufric, 'Re_tau = ', Ufric/VISc
    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        write(3,'(11E15.7)') (Rad_mp(i))*Ufric/VISc, abs(Wmp(i))/Ufric,  &
         (abs(uup(i))/Ufric**2.0)**0.5, (abs(vvp(i))/Ufric**2.0)**0.5,  &
         (abs(wwp(i))/Ufric**2.0)**0.5, uwp(i)/Ufric**2.0, &
         -(Tmp(i)-Twall)/Tfric, TTp(i)/(Tfric*Tfric), uTp(i)/(Tfric*Ufric), &
          vTp(i)/(Tfric*Ufric), wTp(i)/(Tfric*Ufric) 
      end if
    end do 
    close(3)
    open(4,FILE='Budget_uu.dat')
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    write(4,'(A1,1X,A94)') '#', '1:Xrad, 2:Xrad+, 3:Prod, 4:Conv, 5:PressStrain, 6:ViscDiff, 7:TurbDiff, 8:Diss, 9:PresDiff' 
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    do i = 1, Nprob
      if( Ncount(i) > 2) then
        NF    = VISc/Ufric**4.0
        vol1  = 1.0/var_13(i) 
        write(4,'(9E15.7,1X,I6)') Rad_mp(i), Rad_mp(i)*Ufric/VISc, &
        -Puu(i)*NF, Cuu(i)*NF*vol1, Ruu(i)*NF, Diffv_uu(i)*NF*vol1, -Difft_uu(i)*NF, &
        Diss_uu(i)*NF, -PDuu(i)*NF, &
        Ncount(i)
      end if 
    end do
    close(4)

    open(4,FILE='Budget_vv.dat')
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    write(4,'(A1,1X,A94)') '#', '1:Xrad, 2:Xrad+, 3:Prod, 4:Conv, 5:PressStrain, 6:ViscDiff, 7:TurbDiff, 8:Diss, 9:PresDiff' 
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    do i = 1, Nprob
      if( Ncount(i) > 2) then
        NF    = VISc/Ufric**4.0
        vol1  = 1.0/var_13(i) 
        write(4,'(9E15.7,1X,I6)') Rad_mp(i), Rad_mp(i)*Ufric/VISc, &
        -Pvv(i)*NF, Cvv(i)*NF*vol1, Rvv(i)*NF, Diffv_vv(i)*NF*vol1, -Difft_vv(i)*NF, &
        Diss_vv(i)*NF, -PDvv(i)*NF, &
        Ncount(i)
      end if 
    end do
    close(4)

    open(4,FILE='Budget_ww.dat')
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    write(4,'(A1,1X,A94)') '#', '1:Xrad, 2:Xrad+, 3:Prod, 4:Conv, 5:PressStrain, 6:ViscDiff, 7:TurbDiff, 8:Diss, 9:PresDiff' 
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    do i = 1, Nprob
      if( Ncount(i) > 2) then
        NF    = VISc/Ufric**4.0
        vol1  = 1.0/var_13(i) 
        write(4,'(9E15.7,1X,I6)') Rad_mp(i), Rad_mp(i)*Ufric/VISc, &
        -Pww(i)*NF, Cww(i)*NF*vol1, Rww(i)*NF, Diffv_ww(i)*NF*vol1, -Difft_ww(i)*NF, &
        Diss_ww(i)*NF, -PDww(i)*NF, &
        Ncount(i)
      end if 
    end do
    close(4)

    open(4,FILE='Budget_kin.dat')
    write(4,'(A1,1X,A99)') '#', '=============================================================================================' 
    write(4,'(A1,1X,A99)') '#', '1:Xrad,2:Xrad+,3:Prod,4:Conv,5:PressStrain,6:ViscDiff,7:TurbDiff,8:Diss,9:PresDiff,10:Diss_sgs' 
    write(4,'(A1,1X,A99)') '#', '==============================================================================================' 
    do i = 1, Nprob
      if( Ncount(i) > 2) then
        NF    = VISc/Ufric**4.0
        vol1  = 1.0/var_13(i) 
        write(4,'(10E15.7,1X,I6)') Rad_mp(i), Rad_mp(i)*Ufric/VISc, &
        -0.5*(Puu(i)*NF+Pvv(i)*NF+Pww(i)*NF), 0.5*(Cuu(i)*NF*vol1+Cvv(i)*NF*vol1+Cww(i)*NF*vol1), &
         0.5*(Ruu(i)*NF+Rvv(i)*NF+Rww(i)*NF), 0.5*(Diffv_uu(i)*NF*vol1+Diffv_vv(i)*NF*vol1+Diffv_ww(i)*NF*vol1),&
         -0.5*(Difft_uu(i)*NF+Difft_vv(i)*NF+Difft_ww(i)*NF), &
        0.5*(Diss_uu(i)*NF+Diss_vv(i)*NF+Diss_ww(i)*NF), -0.5*(PDuu(i)*NF+PDvv(i)*NF+PDww(i)*NF), Diss_sgs(i)*NF, &
        Ncount(i)
      end if 
    end do
    close(4)

    open(4,FILE='Budget_uw.dat')
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    write(4,'(A1,1X,A94)') '#', '1:Xrad, 2:Xrad+, 3:Prod, 4:Conv, 5:PressStrain, 6:ViscDiff, 7:TurbDiff, 8:Diss, 9:PresDiff' 
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    do i = 1, Nprob
      if( Ncount(i) > 2) then
        NF    = VISc/Ufric**4.0
        vol1  = 1.0/var_13(i) 
        write(4,'(9E15.7,1X,I6)') Rad_mp(i), Rad_mp(i)*Ufric/VISc, &
        Puw(i)*NF, -Cuw(i)*NF*vol1, -Ruw(i)*NF, -Diffv_uw(i)*NF*vol1, +Difft_uw(i)*NF, &
        -Diss_uw(i)*NF, PDuw(i)*NF, &
        Ncount(i)
      end if 
    end do
    close(4)

    open(4,FILE='Budget_ut.dat')
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    write(4,'(A1,1X,A94)') '#', '1:Xrad, 2:Xrad+, 3:Prod, 4:Conv, 5:PressStrain, 6:ViscDiff, 7:TurbDiff, 8:Diss, 9:PresDiff' 
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    do i = 1, Nprob
!      Ufric = (VISc*abs(Wmp(1))/(Rad_mp(1)))**0.5
      if( Ncount(i) > 2) then
        NF    = VISc/(2.0*0.005**2.0)
        FFF    = VISc/(Ufric**3.0*Tfric)
        vol1  = 1.0/var_13(i)
        write(4,'(9E15.7,1X,I6)') Rad_mp(i), (Rad_mp(i))*Ufric/VISc, &
        -Put_p(i)*FFF, Cut(i)*vol1*FFF, Rut(i)*FFF, Diffv_ut(i)*vol1*FFF, &
        -Difft_ut(i)*FFF, Diss_ut(i)*FFF, -PDut(i)*FFF, &
        Ncount(i)
      end if
    end do
    close(4)

    open(4,FILE='Budget_vt.dat')
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    write(4,'(A1,1X,A94)') '#', '1:Xrad, 2:Xrad+, 3:Prod, 4:Conv, 5:PressStrain, 6:ViscDiff, 7:TurbDiff, 8:Diss, 9:PresDiff' 
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    do i = 1, Nprob
!      Ufric = (VISc*abs(Wmp(1))/(Rad_mp(1)))**0.5
      if( Ncount(i) > 2) then
        NF    = VISc/(2.0*0.005**2.0)
        FFF    = VISc/(Ufric**3.0*Tfric)
        vol1  = 1.0/var_13(i)
        write(4,'(9E15.7,1X,I6)') Rad_mp(i), (Rad_mp(i))*Ufric/VISc, &
        -Pvt_p(i)*FFF, Cvt(i)*vol1*FFF, Rvt(i)*FFF, Diffv_vt(i)*vol1*FFF, &
        -Difft_vt(i)*FFF, Diss_vt(i)*FFF, -PDvt(i)*FFF, &
        Ncount(i)
      end if
    end do
    close(4)

    open(4,FILE='Budget_wt.dat')
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    write(4,'(A1,1X,A94)') '#', '1:Xrad, 2:Xrad+, 3:Prod, 4:Conv, 5:PressStrain, 6:ViscDiff, 7:TurbDiff, 8:Diss, 9:PresDiff' 
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    do i = 1, Nprob
!      Ufric = (VISc*abs(Wmp(1))/(Rad_mp(1)))**0.5
      if( Ncount(i) > 2) then
        NF    = VISc/(2.0*0.005**2.0)
        FFF    = VISc/(Ufric**3.0*Tfric)
        vol1  = 1.0/var_13(i)
        write(4,'(9E15.7,1X,I6)') Rad_mp(i), (Rad_mp(i))*Ufric/VISc, &
        -Pwt_p(i)*FFF, Cwt(i)*vol1*FFF, Rwt(i)*FFF, &
        Diffv_wt(i)*vol1*FFF, -Difft_wt(i)*FFF, Diss_wt(i)*FFF, -PDwt(i)*FFF, &
        Ncount(i)
      end if
    end do
    close(4)


    open(4,FILE='Budget_tt.dat')
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    write(4,'(A1,1X,A74)') '#', '1:Xrad, 2:Xrad+, 3:Prod, 4:Conv, 5:ViscDiff, 6:TurbDiff, 7:Diss' 
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    do i = 1, Nprob
!      Ufric = (VISc*abs(Wmp(1))/(Rad_mp(1)))**0.5
      if( Ncount(i) > 2) then
        NF    = VISc/(2.0*0.005**2.0)
        FFF    = VISc/(Ufric**3.0*Tfric)
        vol1  = 1.0/var_13(i)
        write(4,'(7E15.7,1X,I6)') Rad_mp(i), (Rad_mp(i))*Ufric/VISc, &
        -Ptt(i)*NF, Ctt(i)*NF*vol1, Diffv_tt(i)*NF*vol1, -Difft_tt(i)*NF, Diss_tt(i)*NF, &
        Ncount(i)
      end if
    end do
    close(4)


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
  deallocate(uuup)
  deallocate(uuvp)
  deallocate(uuwp)
  deallocate(vvup)
  deallocate(vvvp)
  deallocate(vvwp)
  deallocate(wwup)
  deallocate(wwvp)
  deallocate(wwwp)
  deallocate(uwup)
  deallocate(uwvp)
  deallocate(uwwp)
  deallocate(var_1)
  deallocate(var_10)
  deallocate(var_11)
  deallocate(var_12)
  deallocate(var_13)
  deallocate(Ksgsp)
    deallocate(Puu)
    deallocate(Pvv)
    deallocate(Pww)
    deallocate(Puw)
    deallocate(Cuu)
    deallocate(Cvv)
    deallocate(Cww)
    deallocate(Cuw)
    deallocate(Ruu)
    deallocate(Rvv)
    deallocate(Rww)
    deallocate(Ruw)
    deallocate(PDuu)
    deallocate(PDvv)
    deallocate(PDww)
    deallocate(PDuw)
    deallocate(Diss_uu)
    deallocate(Diss_vv)
    deallocate(Diss_ww)
    deallocate(Diss_uw)
    deallocate(Diss_sgs)
    deallocate(Difft_uu)
    deallocate(Difft_vv)
    deallocate(Difft_ww)
    deallocate(Difft_uw)
    deallocate(Diffv_uu)
    deallocate(Diffv_vv)
    deallocate(Diffv_ww)
    deallocate(Diffv_uw)
    deallocate(Put_p)
    deallocate(Pvt_p)
    deallocate(Pwt_p)
    deallocate(Ptt)
    deallocate(Cut)
    deallocate(Cvt)
    deallocate(Cwt)
    deallocate(Ctt)
    deallocate(Rut)
    deallocate(Rvt)
    deallocate(Rwt)
    deallocate(Diss_ut)
    deallocate(Diss_vt)
    deallocate(Diss_wt)
    deallocate(Diss_tt)
    deallocate(Difft_ut)
    deallocate(Difft_vt)
    deallocate(Difft_wt)
    deallocate(Difft_tt)
    deallocate(Diffv_ut)
    deallocate(Diffv_vt)
    deallocate(Diffv_wt)
    deallocate(Diffv_tt)
    deallocate(PDut)
    deallocate(PDvt)
    deallocate(PDwt)
 

  if(HOT==YES) then
    deallocate(Tmp)
    deallocate(TTp)
    deallocate(uTp)
    deallocate(vTp)
    deallocate(wTp)
  end if
  END SUBROUTINE UserCutLines_budgets_cylind_HOT
