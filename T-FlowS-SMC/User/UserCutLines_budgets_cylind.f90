!======================================================================!
  SUBROUTINE UserCutLines_budgets_cylind() 
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
  REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), ind(:),& 
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
                         Puu(:), Pvv(:), Pww(:), Puv(:), Puw(:), Pvw(:), &
                         Diffv_uu(:), Diffv_vv(:), Diffv_ww(:), Diffv_uv(:), Diffv_uw(:), Diffv_vw(:),&
                         Difft_uu(:), Difft_vv(:), Difft_ww(:), Difft_uv(:), Difft_uw(:), Difft_vw(:), & 
                         Diss_uu(:), Diss_vv(:), Diss_ww(:), Diss_uv(:), Diss_uw(:), Diss_vw(:), Diss_sgs(:),& 
                         Ruu(:), Rvv(:), Rww(:), Ruv(:), Ruw(:),Rvw(:),&
                         PDuu(:), PDvv(:), PDww(:), PDuv(:), PDuw(:),PDvw(:),&
                         Cuu(:), Cvv(:), Cww(:), Cuv(:), Cuw(:),Cvw(:),&
                         Put_p(:), Pvt_p(:), Pwt_p(:), Ptt(:), &
                         Diffv_ut(:), Diffv_vt(:), Diffv_wt(:), Diffv_tt(:), &
                         Difft_ut(:), Difft_vt(:), Difft_wt(:), Difft_tt(:), &
                         Diss_ut(:), Diss_vt(:), Diss_wt(:), Diss_tt(:), &
                         Rut(:), Rvt(:), Rwt(:), Rtt(:), &
                         PDut(:), PDvt(:), PDwt(:), PDtt(:), &
                         Cut(:), Cvt(:), Cwt(:), Ctt(:), Ufric_p(:)
    
  
  INTEGER,ALLOCATABLE :: Np(:), Ncount(:)
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
  allocate(ind(Nprob))

!---- write the probe coordinates out
  do pl=1,Nprob
    read(9,*) ind(pl), z_p(pl) 
  end do
  close(9)

  call SSORT (z_p, ind, Nprob, 0)
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
  allocate(Puv(Nprob));   Puv=0.0
  allocate(Puw(Nprob));   Puw=0.0
  allocate(Pvw(Nprob));   Pvw=0.0

  allocate(Diss_uu(Nprob));   Diss_uu=0.0
  allocate(Diss_vv(Nprob));   Diss_vv=0.0
  allocate(Diss_ww(Nprob));   Diss_ww=0.0
  allocate(Diss_uv(Nprob));   Diss_uv=0.0
  allocate(Diss_uw(Nprob));   Diss_uw=0.0
  allocate(Diss_vw(Nprob));   Diss_vw=0.0

  allocate(Diss_sgs(Nprob));   Diss_sgs=0.0

  allocate(Diffv_uu(Nprob));   Diffv_uu=0.0
  allocate(Diffv_vv(Nprob));   Diffv_vv=0.0
  allocate(Diffv_ww(Nprob));   Diffv_ww=0.0
  allocate(Diffv_uv(Nprob));   Diffv_uv=0.0
  allocate(Diffv_uw(Nprob));   Diffv_uw=0.0
  allocate(Diffv_vw(Nprob));   Diffv_vw=0.0

  allocate(Difft_uu(Nprob));   Difft_uu=0.0
  allocate(Difft_vv(Nprob));   Difft_vv=0.0
  allocate(Difft_ww(Nprob));   Difft_ww=0.0
  allocate(Difft_uv(Nprob));   Difft_uv=0.0
  allocate(Difft_uw(Nprob));   Difft_uw=0.0
  allocate(Difft_vw(Nprob));   Difft_vw=0.0

  allocate(Ruu(Nprob));   Ruu=0.0
  allocate(Rvv(Nprob));   Rvv=0.0
  allocate(Rww(Nprob));   Rww=0.0
  allocate(Ruv(Nprob));   Ruv=0.0
  allocate(Ruw(Nprob));   Ruw=0.0
  allocate(Rvw(Nprob));   Rvw=0.0


  allocate(Cuu(Nprob));   Cuu=0.0
  allocate(Cvv(Nprob));   Cvv=0.0
  allocate(Cww(Nprob));   Cww=0.0
  allocate(Cuv(Nprob));   Cuv=0.0
  allocate(Cuw(Nprob));   Cuw=0.0
  allocate(Cvw(Nprob));   Cvw=0.0

  allocate(PDuu(Nprob));   PDuu=0.0
  allocate(PDvv(Nprob));   PDvv=0.0
  allocate(PDww(Nprob));   PDww=0.0
  allocate(PDuv(Nprob));   PDuv=0.0
  allocate(PDuw(Nprob));   PDuw=0.0
  allocate(PDvw(Nprob));   PDvw=0.0

  allocate(Rad_1(Nprob));  Rad_1=0.0
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
    R_max = 0.0

    do i = 1, Nprob
      R_max = max(z_p(i),R_max)
    end do

    Ncount_wall = 0
    do c = -Nbc,1
      if(TypeBC(c) == WALLFL) then  
        Ncount_wall = Ncount_wall + 1
      end if
    end do
    call IGlSum(Ncount_wall)
    do i = 1, Nprob-1
      do c=1,NC
        Rad_2 = (xc(c)**2 + yc(c)**2)**0.5
        if(Rad_2 > abs(z_p(i+1)) .and. Rad_2 < abs(z_p(i))) then
 
          R           = (xc(c)*xc(c) + yc(c)*yc(c))**0.5 + tiny
          Urad_mean   = (U % mean(c) * xc(c) / R  + V % mean(c) * yc(c) / R)
          Utan_mean   = (-U % mean(c) * yc(c) / R  + V % mean(c) * xc(c) / R) 

          if(IsNearWall(c)) then
            Ufric_p(i) = Ufric_p(i) + 1.0  !(VISc * (U % mean(c)**2 + V % mean(c)**2 + W % mean(c)**2)**0.5/WallDs(c))**0.5
          end if
 
          Ump(i)   = Ump(i) + U % mean(c)
          Vmp(i)   = Vmp(i) + V % mean(c)
          Wmp(i)   = Wmp(i) + W % mean(c)
          uup(i)   = uup(i) + (uu % mean(c)- Urad_mean * Urad_mean)
          vvp(i)   = vvp(i) + (vv % mean(c)- Utan_mean * Utan_mean)
          wwp(i)   = wwp(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
          uvp(i)   = uvp(i) + (uv % mean(c)- Urad_mean * Utan_mean )
          uwp(i)   = uwp(i) + (uw % mean(c)- Urad_mean * W % mean(c))
          vwp(i)   = vwp(i) + (vw % mean(c)- Utan_mean * W % mean(c))

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



          var_10(i) = var_10(i) + uu % n(c)
          var_11(i) = var_11(i) + vv % n(c)
          var_12(i) = var_12(i) + ww % n(c)
          var_13(i) = var_13(i) + volume(c)

            Puu(i)    = Puu(i) + Puu_mean(c)
            Pvv(i)    = Pvv(i) + Pvv_mean(c)
            Pww(i)    = Pww(i) + Pww_mean(c)
            Puv(i)    = Puv(i) + Puv_mean(c)
            Puw(i)    = Puw(i) + Puw_mean(c)
            Pvw(i)    = Pvw(i) + Pvw_mean(c)

            Cuu(i)    = Cuu(i) + C_uu_mean(c)
            Cvv(i)    = Cvv(i) + C_vv_mean(c)
            Cww(i)    = Cww(i) + C_ww_mean(c)
            Cuw(i)    = Cuw(i) + C_uw_mean(c)

            Ruu(i)    = Ruu(i) + PR_uu_mean(c)
            Rvv(i)    = Rvv(i) + PR_vv_mean(c)
            Rww(i)    = Rww(i) + PR_ww_mean(c)
            Ruv(i)    = Ruv(i) + PR_uv_mean(c)
            Ruw(i)    = Ruw(i) + PR_uw_mean(c)
            Rvw(i)    = Rvw(i) + PR_vw_mean(c)

            PDuu(i)    = PDuu(i) + PD_uu_mean(c)
            PDvv(i)    = PDvv(i) + PD_vv_mean(c)
            PDww(i)    = PDww(i) + PD_ww_mean(c)
            PDuv(i)    = PDuv(i) + PD_uv_mean(c)
            PDuw(i)    = PDuw(i) + PD_uw_mean(c)
            PDvw(i)    = PDvw(i) + PD_vw_mean(c)

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
            Diss_uv(i)  = Diss_uv(i) + Diss_uv_mean(c)
            Diss_uw(i)  = Diss_uw(i) + Diss_uw_mean(c)
            Diss_vw(i)  = Diss_vw(i) + Diss_vw_mean(c)

            Diss_sgs(i)  = Diss_sgs(i) + Diss_sgs_mean(c)

          Rad_mp(i) = Rad_mp(i) + R
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
      call GloSum(Puv(pl))
      call GloSum(Puw(pl))
      call GloSum(Pvw(pl))

      call GloSum(PDuu(pl))
      call GloSum(PDvv(pl))
      call GloSum(PDww(pl))
      call GloSum(PDuv(pl))
      call GloSum(PDuw(pl))
      call GloSum(PDvw(pl))

      call GloSum(Cuu(pl))
      call GloSum(Cvv(pl))
      call GloSum(Cww(pl))
      call GloSum(Cuv(pl))
      call GloSum(Cuw(pl))
      call GloSum(Cvw(pl))

      call GloSum(Ruu(pl))
      call GloSum(Rvv(pl))
      call GloSum(Rww(pl))
      call GloSum(Ruv(pl))
      call GloSum(Ruw(pl))
      call GloSum(Rvw(pl))

      call GloSum(Difft_uu(pl))
      call GloSum(Difft_vv(pl))
      call GloSum(Difft_ww(pl))
      call GloSum(Difft_uv(pl))
      call GloSum(Difft_uw(pl))
      call GloSum(Difft_vw(pl))

      call GloSum(Diffv_uu(pl))
      call GloSum(Diffv_vv(pl))
      call GloSum(Diffv_ww(pl))
      call GloSum(Diffv_uv(pl))
      call GloSum(Diffv_uw(pl))
      call GloSum(Diffv_vw(pl))

      call GloSum(Diss_uu(pl))
      call GloSum(Diss_vv(pl))
      call GloSum(Diss_ww(pl))
      call GloSum(Diss_uv(pl))
      call GloSum(Diss_uw(pl))
      call GloSum(Diss_vw(pl))

      call GloSum(Diss_sgs(pl))

    count =  count + Ncount(pl) 

  end do

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
      Puv(i)    = Puv(i)/Ncount(i)
      Puw(i)    = Puw(i)/Ncount(i)
      Pvw(i)    = Pvw(i)/Ncount(i)

      PDuu(i)    = PDuu(i)/Ncount(i)
      PDvv(i)    = PDvv(i)/Ncount(i)
      PDww(i)    = PDww(i)/Ncount(i)
      PDuv(i)    = PDuv(i)/Ncount(i)
      PDuw(i)    = PDuw(i)/Ncount(i)
      PDvw(i)    = PDvw(i)/Ncount(i)

      Cuu(i)    = Cuu(i)/Ncount(i)
      Cvv(i)    = Cvv(i)/Ncount(i)
      Cww(i)    = Cww(i)/Ncount(i)
      Cuv(i)    = Cuv(i)/Ncount(i)
      Cuw(i)    = Cuw(i)/Ncount(i)
      Cvw(i)    = Cvw(i)/Ncount(i)


      Ruu(i)    = Ruu(i)/Ncount(i)
      Rvv(i)    = Rvv(i)/Ncount(i)
      Rww(i)    = Rww(i)/Ncount(i)
      Ruv(i)    = Ruv(i)/Ncount(i)
      Ruw(i)    = Ruw(i)/Ncount(i)
      Rvw(i)    = Rvw(i)/Ncount(i)

      Difft_uu(i)    = Difft_uu(i)/Ncount(i)
      Difft_vv(i)    = Difft_vv(i)/Ncount(i)
      Difft_ww(i)    = Difft_ww(i)/Ncount(i)
      Difft_uv(i)    = Difft_uv(i)/Ncount(i)
      Difft_uw(i)    = Difft_uw(i)/Ncount(i)
      Difft_vw(i)    = Difft_vw(i)/Ncount(i)

      Diffv_uu(i)    = Diffv_uu(i)/Ncount(i)
      Diffv_vv(i)    = Diffv_vv(i)/Ncount(i)
      Diffv_ww(i)    = Diffv_ww(i)/Ncount(i)
      Diffv_uv(i)    = Diffv_uv(i)/Ncount(i)
      Diffv_uw(i)    = Diffv_uw(i)/Ncount(i)
      Diffv_vw(i)    = Diffv_vw(i)/Ncount(i)

      Diss_sgs(i)    = Diss_sgs(i)/Ncount(i)

      Diss_uu(i)    = Diss_uu(i)/Ncount(i)
      Diss_vv(i)    = Diss_vv(i)/Ncount(i)
      Diss_ww(i)    = Diss_ww(i)/Ncount(i)
      Diss_uv(i)    = Diss_uv(i)/Ncount(i)        
      Diss_uw(i)    = Diss_uw(i)/Ncount(i)        
      Diss_vw(i)    = Diss_vw(i)/Ncount(i)        

      Rad_mp(i) =  Rad_mp(i)/Ncount(i)

      end if
    end do 

    Ufric = 0.0
    do i = 1, Nprob
      Ufric = 1.0  !max(Ufric_p(i),Ufric)
    end do

    NF    = 1.0   !VISc/Ufric**4.0
    open(4,FILE='Budget_uu.dat')
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    write(4,'(A1,1X,A94)') '#', '1:Xrad, 2:Xrad+, 3:Prod, 4:Conv, 5:PressStrain, 6:ViscDiff, 7:TurbDiff, 8:Diss, 9:PresDiff' 
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    do i = 1, Nprob
      if( Ncount(i) > 2) then
        vol1  = 1.0/var_13(i) 
        write(4,'(10E15.7,1X,I6)') Rad_mp(i), Rad_mp(i)*Ufric/VISc, &
        -Puu(i)*NF, -Cuu(i)*NF*vol1, Ruu(i)*NF, Diffv_uu(i)*NF*vol1, -Difft_uu(i)*NF, &
        Diss_uu(i)*NF, -PDuu(i)*NF,&
        (-Puu(i)-Cuu(i)*NF*vol1+Ruu(i)+Diffv_uu(i)*NF*vol1-Difft_uu(i)+Diss_uu(i)-PDuu(i)),&
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
        vol1  = 1.0/var_13(i) 
        write(4,'(10E15.7,1X,I6)') Rad_mp(i), Rad_mp(i)*Ufric/VISc, &
        -Pvv(i)*NF, -Cvv(i)*NF*vol1, Rvv(i)*NF, Diffv_vv(i)*NF*vol1, -Difft_vv(i)*NF, &
        Diss_vv(i)*NF, -PDvv(i)*NF, &
        (-Pvv(i)-Cvv(i)*NF*vol1+Rvv(i)+Diffv_vv(i)*NF*vol1-Difft_vv(i)+Diss_vv(i)-PDvv(i)),&
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
        vol1  = 1.0/var_13(i) 
        write(4,'(10E15.7,1X,I6)') Rad_mp(i), Rad_mp(i)*Ufric/VISc, &
        -Pww(i)*NF, -Cww(i)*NF*vol1, Rww(i)*NF, Diffv_ww(i)*NF*vol1, -Difft_ww(i)*NF, &
        Diss_ww(i)*NF, -PDww(i)*NF, &
        (-Pww(i)-Cww(i)*NF*vol1+Rww(i)+Diffv_ww(i)*NF*vol1-Difft_ww(i)+Diss_ww(i)-PDww(i)),&
        Ncount(i)
      end if 
    end do
    close(4)

    open(4,FILE='Budget_kin.dat')
    write(4,'(A1,1X,A99)') '#','===========================================================================' 
    write(4,'(A1,1X,A99)') '#','1:Xrad,2:Xrad+,3:Prod,4:Conv,5:PressStrain,6:ViscDiff,7:TurbDiff,8:Diss,9:PresDiff,10:Diss_sgs' 
    write(4,'(A1,1X,A99)') '#','=============================================================================' 
    do i = 1, Nprob
      if( Ncount(i) > 2) then
        vol1  = 1.0/var_13(i) 
        write(4,'(10E15.7,1X,I6)') Rad_mp(i), Rad_mp(i)*Ufric/VISc, &
        0.5*(Puu(i)*NF+Pvv(i)*NF+Pww(i)*NF), 0.5*(Cuu(i)*NF*vol1+Cvv(i)*NF*vol1+Cww(i)*NF*vol1), &
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
        vol1  = 1.0/var_13(i) 
        write(4,'(10E15.7,1X,I6)') Rad_mp(i), Rad_mp(i)*Ufric/VISc, &
        -Puw(i)*NF, -Cuw(i)*NF*vol1, Ruw(i)*NF, Diffv_uw(i)*NF*vol1, -Difft_uw(i)*NF, &
        Diss_uw(i)*NF, -PDuw(i)*NF, &
        (-Puw(i)-Cuw(i)*NF*vol1+Ruw(i)+Diffv_uw(i)*NF*vol1-Difft_uw(i)+Diss_uw(i)-PDuw(i)),&
        Ncount(i)
      end if 
    end do
    close(4)

    open(4,FILE='Budget_uv.dat')
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    write(4,'(A1,1X,A94)') '#', '1:Xrad, 2:Xrad+, 3:Prod, 4:Conv, 5:PressStrain, 6:ViscDiff, 7:TurbDiff, 8:Diss, 9:PresDiff' 
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    do i = 1, Nprob
      if( Ncount(i) > 2) then
        vol1  = 1.0/var_13(i) 
        write(4,'(10E15.7,1X,I6)') Rad_mp(i), Rad_mp(i)*Ufric/VISc, &
        -Puv(i)*NF, -Cuv(i)*NF*vol1, Ruv(i)*NF, Diffv_uv(i)*NF*vol1, -Difft_uv(i)*NF, &
        Diss_uv(i)*NF, -PDuv(i)*NF, &
        (-Puv(i)-Cuv(i)*NF*vol1+Ruv(i)+Diffv_uv(i)*NF*vol1-Difft_uv(i)+Diss_uv(i)-PDuv(i)),&
        Ncount(i)
      end if 
    end do
    close(4)

    open(4,FILE='Budget_vw.dat')
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    write(4,'(A1,1X,A94)') '#', '1:Xrad, 2:Xrad+, 3:Prod, 4:Conv, 5:PressStrain, 6:ViscDiff, 7:TurbDiff, 8:Diss, 9:PresDiff' 
    write(4,'(A1,1X,A94)') '#', '==========================================================================================' 
    do i = 1, Nprob
      if( Ncount(i) > 2) then
        vol1  = 1.0/var_13(i) 
        write(4,'(10E15.7,1X,I6)') Rad_mp(i), Rad_mp(i)*Ufric/VISc, &
        -Pvw(i)*NF, -Cvw(i)*NF*vol1, Rvw(i)*NF, Diffv_vw(i)*NF*vol1, -Difft_vw(i)*NF, &
        Diss_vw(i)*NF, -PDvw(i)*NF, &
        (-Pvw(i)-Cvw(i)*NF*vol1+Rvw(i)+Diffv_vw(i)*NF*vol1-Difft_vw(i)+Diss_vw(i)-PDvw(i)),&
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
    deallocate(Puv)
    deallocate(Pvw)
    deallocate(Cuu)
    deallocate(Cvv)
    deallocate(Cww)
    deallocate(Cuw)
    deallocate(Cuv)
    deallocate(Cvw)
    deallocate(Ruu)
    deallocate(Rvv)
    deallocate(Rww)
    deallocate(Ruw)
    deallocate(Ruv)
    deallocate(Rvw)
    deallocate(PDuu)
    deallocate(PDvv)
    deallocate(PDww)
    deallocate(PDuw)
    deallocate(PDuv)
    deallocate(PDvw)
    deallocate(Diss_uu)
    deallocate(Diss_vv)
    deallocate(Diss_ww)
    deallocate(Diss_uw)
    deallocate(Diss_uv)
    deallocate(Diss_vw)
    deallocate(Diss_sgs)
    deallocate(Difft_uu)
    deallocate(Difft_vv)
    deallocate(Difft_ww)
    deallocate(Difft_uw)
    deallocate(Difft_uv)
    deallocate(Difft_vw)
    deallocate(Diffv_uu)
    deallocate(Diffv_vv)
    deallocate(Diffv_ww)
    deallocate(Diffv_uw)
    deallocate(Diffv_uv)
    deallocate(Diffv_vw)

  if(HOT==YES) then
    deallocate(Tmp)
    deallocate(TTp)
    deallocate(uTp)
    deallocate(vTp)
    deallocate(wTp)
  end if
  END SUBROUTINE UserCutLines_budgets_cylind
