!======================================================================!
  SUBROUTINE UserDiffuser() 
!----------------------------------------------------------------------!
! Reads the ".1D" file created by the "Generator" and averages the     !
! results in the planes defined by coordinates in it. Then averages    !
! the values of Umean, Vmean, Wmean, uu, vv, ww, uv, uw and vw and     !
! writes them into file ".1Dr".                                        !
!----------------------------------------------------------------------!
  USE all_mod
  USE allp_mod
  USE pro_mod
  USE rans_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
!------------------------------[Calling]-------------------------------!
  INTERFACE
    LOGICAL FUNCTION Approx(A,B,tol)
      REAL           :: A,B
      REAL, OPTIONAL :: tol
    END FUNCTION Approx
  END INTERFACE 
!-------------------------------[Locals]-------------------------------!
  INTEGER             :: Nprob, pl, c, i, count, k, c1, c2, s, j
  CHARACTER           :: namCoo*80, namPro*80, answer*80, JetIn*31, namOut*80
  REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), xc_p(:), yc_p(:), zc_p(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), tke(:), nu(:)
  INTEGER,ALLOCATABLE :: Np(:), Ncount(:)
  REAL                :: dummy
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exc_p $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $  
!======================================================================!

    call GraPhi(U % n, 1, Ux,.TRUE.)    ! dU/dx
    call GraPhi(V % n, 1, Vx,.TRUE.)    ! dU/dx
    call GraPhi(W % n, 1, Wx,.TRUE.)    ! dU/dx
    call GraPhi(U % n, 2, Uy,.TRUE.)    ! dU/dx
    call GraPhi(V % n, 2, Vy,.TRUE.)    ! dU/dx
    call GraPhi(W % n, 2, Wy,.TRUE.)    ! dU/dx
    call GraPhi(U % n, 3, Uz,.TRUE.)    ! dU/dx
    call GraPhi(V % n, 3, Vz,.TRUE.)    ! dU/dx
    call GraPhi(W % n, 3, Wz,.TRUE.)    ! dU/dx



  do k = 0, 55

!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
    if(k == 0) then
      open(9, FILE='y0.dat')
      j = k + 1
      JetIn  = 'c13.2_z0250_x-2_IUS_ZETA-F.dat'
    else if(k == 1) then 
      open(9, FILE='y0.dat')
      j = k+1  
      JetIn  = 'c13.2_z0500_x-2_IUS_ZETA-F.dat'
    else if(k == 2) then 
      open(9, FILE='y0.dat')
      j = k+1  
      JetIn  = 'c13.2_z0750_x-2_IUS_ZETA-F.dat'
    else if(k == 3) then 
      open(9, FILE='y0.dat')
      j = k+1  
      JetIn  = 'c13.2_z0875_x-2_IUS_ZETA-F.dat'
    else if(k == 4) then 
      open(9, FILE='y0.dat')
      j = k+1  
      JetIn  = 'c13.2_z0250_x0_IUS_ZETA-F.dat'
    else if(k == 5) then 
      open(9, FILE='y0.dat')
      j = k+1  
      JetIn  = 'c13.2_z0500_x0_IUS_ZETA-F.dat'
    else if(k == 6) then 
      open(9, FILE='y0.dat')
      j = k+1  
      JetIn  = 'c13.2_z0750_x0_IUS_ZETA-F.dat'
    else if(k == 7) then 
      open(9, FILE='y0.dat')
      j = k+1  
      JetIn  = 'c13.2_z0875_x0_IUS_ZETA-F.dat'
    else if(k == 8) then 
      open(9, FILE='y2.dat')
      j = k+1 
      JetIn  = 'c13.2_z0250_x2_IUS_ZETA-F.dat' 
    else if(k == 9) then 
      open(9, FILE='y2.dat')
      j = k+1  
      JetIn  = 'c13.2_z0500_x2_IUS_ZETA-F.dat'
    else if(k == 10) then 
      open(9, FILE='y2.dat')
      j = k+1  
      JetIn  = 'c13.2_z0750_x2_IUS_ZETA-F.dat'
    else if(k == 11) then 
      open(9, FILE='y2.dat')
      j = k+1  
      JetIn  = 'c13.2_z0875_x2_IUS_ZETA-F.dat'
    else if(k == 12) then 
      open(9, FILE='y4.dat')
      j = k+1  
      JetIn  = 'c13.2_z0250_x4_IUS_ZETA-F.dat'
    else if(k == 13) then 
      open(9, FILE='y4.dat')
      j = k+1  
      JetIn  = 'c13.2_z0500_x4_IUS_ZETA-F.dat'
    else if(k == 14) then 
      open(9, FILE='y4.dat')
      j = k+1 
      JetIn  = 'c13.2_z0750_x4_IUS_ZETA-F.dat' 
    else if(k == 15) then 
      open(9, FILE='y4.dat')
      j = k+1  
      JetIn  = 'c13.2_z0875_x4_IUS_ZETA-F.dat'
    else if(k == 16) then 
      open(9, FILE='y6.dat')
      j = k+1 
      JetIn  = 'c13.2_z0250_x6_IUS_ZETA-F.dat' 
    else if(k == 17) then 
      open(9, FILE='y6.dat')
      j = k+1  
      JetIn  = 'c13.2_z0500_x6_IUS_ZETA-F.dat'
    else if(k == 18) then 
      open(9, FILE='y6.dat')
      j = k+1  
      JetIn  = 'c13.2_z0750_x6_IUS_ZETA-F.dat'
    else if(k == 19) then 
      open(9, FILE='y6.dat')
      j = k+1  
      JetIn  = 'c13.2_z0875_x6_IUS_ZETA-F.dat'
    else if(k == 20) then 
      open(9, FILE='y8.dat')
      j = k+1  
      JetIn  = 'c13.2_z0250_x8_IUS_ZETA-F.dat'
    else if(k == 21) then 
      open(9, FILE='y8.dat')
      j = k+1  
      JetIn  = 'c13.2_z0500_x8_IUS_ZETA-F.dat'
    else if(k == 22) then 
      open(9, FILE='y8.dat')
      j = k+1  
      JetIn  = 'c13.2_z0750_x8_IUS_ZETA-F.dat'
    else if(k == 23) then 
      open(9, FILE='y8.dat')
      j = k+1  
      JetIn  = 'c13.2_z0875_x8_IUS_ZETA-F.dat'
    else if(k == 24) then 
      open(9, FILE='y10.dat')
      j = k+1  
      JetIn  = 'c13.2_z0250_x10_IUS_ZETA-F.dat'
    else if(k == 25) then 
      open(9, FILE='y10.dat')
      j = k+1  
      JetIn  = 'c13.2_z0500_x10_IUS_ZETA-F.dat'
    else if(k == 26) then 
      open(9, FILE='y10.dat')
      j = k+1  
      JetIn  = 'c13.2_z0750_x10_IUS_ZETA-F.dat'
    else if(k == 27) then 
      open(9, FILE='y10.dat')
      j = k+1  
      JetIn  = 'c13.2_z0875_x10_IUS_ZETA-F.dat'
    else if(k == 28) then 
      open(9, FILE='y12.dat')
      j = k+1  
      JetIn  = 'c13.2_z0250_x12_IUS_ZETA-F.dat'
    else if(k == 29) then 
      open(9, FILE='y12.dat')
      j = k+1  
      JetIn  = 'c13.2_z0500_x12_IUS_ZETA-F.dat'
    else if(k == 30) then 
      open(9, FILE='y12.dat')
      j = k+1  
      JetIn  = 'c13.2_z0750_x12_IUS_ZETA-F.dat'
    else if(k == 31) then 
      open(9, FILE='y12.dat')
      j = k+1  
      JetIn  = 'c13.2_z0875_x12_IUS_ZETA-F.dat'
    else if(k == 32) then 
      open(9, FILE='y14.dat')
      j = k+1  
      JetIn  = 'c13.2_z0250_x14_IUS_ZETA-F.dat'
    else if(k == 33) then 
      open(9, FILE='y14.dat')
      j = k+1  
      JetIn  = 'c13.2_z0500_x14_IUS_ZETA-F.dat'
    else if(k == 34) then 
      open(9, FILE='y14.dat')
      j = k+1  
      JetIn  = 'c13.2_z0750_x14_IUS_ZETA-F.dat'
    else if(k == 35) then 
      open(9, FILE='y14.dat')
      j = k+1  
      JetIn  = 'c13.2_z0875_x14_IUS_ZETA-F.dat'
    else if(k == 36) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0250_x155_IUS_ZETA-F.dat'
    else if(k == 37) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0500_x155_IUS_ZETA-F.dat'
    else if(k == 38) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0750_x155_IUS_ZETA-F.dat'
    else if(k == 39) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0875_x155_IUS_ZETA-F.dat'
    else if(k == 40) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0250_x17_IUS_ZETA-F.dat'
    else if(k == 41) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0500_x17_IUS_ZETA-F.dat'
    else if(k == 42) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0750_x17_IUS_ZETA-F.dat' 
    else if(k == 43) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0875_x17_IUS_ZETA-F.dat'
    else if(k == 44) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0250_x185_IUS_ZETA-F.dat'
    else if(k == 45) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0500_x185_IUS_ZETA-F.dat'
    else if(k == 46) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0750_x185_IUS_ZETA-F.dat'
    else if(k == 47) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0875_x185_IUS_ZETA-F.dat'
    else if(k == 48) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0250_x20_IUS_ZETA-F.dat'
    else if(k == 49) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0500_x20_IUS_ZETA-F.dat'
    else if(k == 50) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0750_x20_IUS_ZETA-F.dat'
    else if(k == 51) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0875_x20_IUS_ZETA-F.dat'
    else if(k == 52) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0250_x215_IUS_ZETA-F.dat'
    else if(k == 53) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0500_x215_IUS_ZETA-F.dat'
    else if(k == 54) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0750_x215_IUS_ZETA-F.dat'
    else if(k == 55) then 
      open(9, FILE='y15.5.dat')
      j = k+1  
      JetIn  = 'c13.2_z0975_x215_IUS_ZETA-F.dat'
    end if

    if(this < 2) write(6, *) '# Now opening file: ', JetIn 
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
    allocate(uup(Nprob));   uup=0.0
    allocate(vvp(Nprob));   vvp=0.0
    allocate(wwp(Nprob));   wwp=0.0
    allocate(uvp(Nprob));   uvp=0.0
    allocate(uwp(Nprob));   uwp=0.0
    allocate(vwp(Nprob));   vwp=0.0
    allocate(tke(Nprob));   tke=0.0
    allocate(nu(Nprob));    nu=0.0
    allocate(xc_p(Nprob));    xc_p=0.0
    allocate(yc_p(Nprob));    yc_p=0.0
    allocate(zc_p(Nprob));    zc_p=0.0

    allocate(Ncount(Nprob)); Ncount=0
    count = 0


!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!

    do i = 1, Nprob
      do c=1,NC
        if( Approx( xc(c),xc(Cm(j)),4.0e-2 ).and.Approx( zc(c),zc(Cm(j)),2.0e-1 )  ) then
          if(yc(c) > z_p(i) .and. yc(c) < z_p(i+1)) then
            xc_p(i)    = xc_p(i) + xc(c) 
            yc_p(i)    = yc_p(i) + yc(c) 
            zc_p(i)    = zc_p(i) + zc(c) 
            Ump(i)   = Ump(i) + U % n(c) 
            Vmp(i)   = Vmp(i) + V % n(c) 
            Wmp(i)   = Wmp(i) + W % n(c)
            uuP(i)   = uup(i) -2.0*VISt(c)*Ux(c) + 2.0/3.0 * Kin % n(c) 
            vvp(i)   = vvp(i) -2.0*VISt(c)*Ux(c) + 2.0/3.0 * Kin % n(c)  
            wwp(i)   = wwp(i) -2.0*VISt(c)*Ux(c) + 2.0/3.0 * Kin % n(c) 
            uvp(i)   = uvp(i) -VISt(c)*(Uy(c) + Vx(c))
            uwp(i)   = uwp(i) -VISt(c)*(Uz(c) + Wx(c))
            vwp(i)   = vwp(i) -VISt(c)*(Vz(c) + Wy(c))
            tke(i)   = tke(i) + Kin % n(c)     
            nu(i)    = nu(i) + VISt(c)
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

      call GloSum(tke(pl))
      call GloSum(nu(pl))
      call GloSum(xc_p(pl))
      call GloSum(yc_p(pl))
      call GloSum(zc_p(pl))
      count =  count + Ncount(pl) 

    end do

!    JetIn  = 'Diff_res_x'
!    write(JetIn(9:10),'(I2)') k

    open(3,FILE=JetIn)

    write(3,'(A)') '# International University of Sarajevo'
    write(3,'(A)') '# Muhamed Hadziabdic'
    write(3,'(A)') '# mhadziabdic@ius.edu.ba'                  
    write(3,'(A)') '# T-FlowS, finite-volume, un-structured' 
    write(3,'(A)') '# 2nd order in space, 2nd order in time' 
    write(3,'(A)') '# grid size:1.25e+6'                    
    write(3,'(A)') '# Model used:  RANS zeta-f model'       
    write(3,'(A)') '# x/H   y/H   z/B   U/Ubulk  V/Ubulk  W/Ubulk  k/Ubulk2  &
                            uu/Ubulk2  vv/Ubulk2  ww/Ubulk2  uv/Ubulk2  uw/Ubulk2  vw/Ubulk2  νt/ν'
    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        Wmp(i)    =  Wmp(i)/Ncount(i)
        Ump(i)    =  Ump(i)/Ncount(i)
        Vmp(i)    =  Vmp(i)/Ncount(i)
        uup(i)    =  uup(i)/Ncount(i)
        vvp(i)    =  vvp(i)/Ncount(i)
        wwp(i)    =  wwp(i)/Ncount(i)
        uvp(i)    =  uvp(i)/Ncount(i)
        uwp(i)    =  uwp(i)/Ncount(i)
        vwp(i)    =  vwp(i)/Ncount(i)
        tke(i)    =  tke(i)/Ncount(i)
        nu(i)     =  nu(i)/Ncount(i)
        xc_p(i)     =  xc_p(i)/Ncount(i)
        yc_p(i)     =  yc_p(i)/Ncount(i)
        zc_p(i)     =  zc_p(i)/Ncount(i)
        write(3,'(14E15.7,I6)') xc_p(i), yc_p(i), zc_p(i)/3.33, Ump(i), Vmp(i), Wmp(i), & 
                               tke(i), uup(i), vvp(i), wwp(i), uvp(i), uwp(i), vwp(i), nu(i)/VISc, Ncount(i) 
      end if
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
    deallocate(tke)
    deallocate(nu)
    deallocate(xc_p)
    deallocate(yc_p)
    deallocate(zc_p)
    deallocate(Ncount)
  end do   !end number of radius


  if(this < 2) write(*,*) 'Finished with UserDiffuser '

  END SUBROUTINE UserDiffuser
