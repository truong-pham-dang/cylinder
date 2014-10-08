!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                    ________  ______                                  !
!                   |        ||      \                                 !
!                   `--.  .--'|  ,-.  \________  ___                   !
!                      |  |___|  |_/  /  _  \  \/  /                   !
!                      |  |___|      /  _____\    /                    !
!                      |  |   |  |\  \  \____/    \                    !
!                      |__|   |__| \__\_____/__/\__\                   !
!                                                                      !
!                   UNSTRUCTURED PRAGMATIC LES SOLVER                  !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
!   Pragmatic means that it can use both FRACTIONAL STEP and SIMPLE    !
!                                                                      !
!                            version: Chatou                           !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
!                                     Bojan NICENO                     !
!                                     Delft University of Technology   !
!                                     Faculty of Applied Sciences      !
!                                     Section Thermofluids             !
!                                     niceno@ws.tn.tudelft.nl          !
!                                                                      !
!======================================================================!
  PROGRAM Processor
!----------------------------------------------------------------------!
!   Unstructured Finite Volume LES/RANS solver.                        !
!   Authors: Bojan NICENO & Muhamed HADZIABDIC                         !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------[Calling]-------------------------------!
  REAL :: CorUVW
!-------------------------------[Locals]-------------------------------!
  INTEGER          :: i, m, n, Ndtt_temp, HOTtemp, SIMULAtemp, c, Nproc,s,c1,c2
  REAL             :: Mres, CPUtim
  REAL             :: start, finish
  CHARACTER        :: namSav*10, converted_INT*80
  LOGICAL          :: restar, multiple
!-----------------------------[Interface]------------------------------!
  INTERFACE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
    SUBROUTINE NewUVW(var, Ui, dUidi, dUidj, dUidk,  &
          Si, Sj, Sk, Di, Dj, Dk, Pi, dUjdi, dUkdi) 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
      USE all_mod
      USE pro_mod
      IMPLICIT NONE
      INTEGER       :: var
      TYPE(Unknown) :: Ui
      REAL          :: dUidi(-NbC:NC), dUidj(-NbC:NC), dUidk(-NbC:NC)
      REAL          :: Si(NS), Sj(NS), Sk(NS) 
      REAL          :: Di(NS), Dj(NS), Dk(NS) 
      REAL          :: Pi(-NbC:NC), dUjdi(-NbC:NC), dUkdi(-NbC:NC) 
    END SUBROUTINE NewUVW
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
    SUBROUTINE CalcSc(var, PHI, dPHIdx, dPHIdy, dPHIdz)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
      USE all_mod
      USE pro_mod
      IMPLICIT NONE
      INTEGER       :: var
      TYPE(Unknown) :: PHI
      REAL          :: dPHIdx(-NbC:NC),dPHIdy(-NbC:NC),dPHIdz(-NbC:NC)
    END SUBROUTINE CalcSc
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
    SUBROUTINE CalcTurb(var, PHI, dPHIdx, dPHIdy, dPHIdz)  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
      USE all_mod
      USE pro_mod
      IMPLICIT NONE
      INTEGER       :: var
      TYPE(Unknown) :: PHI
      REAL          :: dPHIdx(-NbC:NC),dPHIdy(-NbC:NC),dPHIdz(-NbC:NC)
    END SUBROUTINE CalcTurb
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
    SUBROUTINE ProSav(namAut)  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
      USE all_mod
      USE pro_mod
      IMPLICIT NONE
       CHARACTER, OPTIONAL :: namAut*(*)
    END  SUBROUTINE ProSav  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
    SUBROUTINE DatSav(namAut)  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
      USE all_mod
      USE pro_mod
      IMPLICIT NONE
      CHARACTER, OPTIONAL :: namAut*(*)
    END SUBROUTINE DatSav  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
    SUBROUTINE SavRes(namAut)  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
      USE all_mod
      USE pro_mod
      IMPLICIT NONE
      CHARACTER, OPTIONAL :: namAut*(*)
    END SUBROUTINE SavRes  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
!    SUBROUTINE GenLoaXML2_01(sub, NCsub, NamAut)  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
!      USE all_mod
!      USE pro_mod
!      IMPLICIT NONE
!      INTEGER      :: sub, NCsub
!      CHARACTER, OPTIONAL :: namAut*(*)
!    END SUBROUTINE GenLoaXML2_01  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!    SUBROUTINE UserProbe1D()
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!      USE all_mod
!      USE pro_mod
!      IMPLICIT NONE
!    END SUBROUTINE UserProbe1D
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    SUBROUTINE UserProbe2D(namAut)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
      USE all_mod
      USE pro_mod
      IMPLICIT NONE
      CHARACTER, OPTIONAL :: namAut*(*)
    END SUBROUTINE UserProbe2D
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    SUBROUTINE CalcShear(Ui, Vi, Wi, She)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
      USE all_mod
      USE pro_mod
      USE les_mod
      IMPLICIT NONE
      REAL          :: Ui(-NbC:NC), Vi(-NbC:NC), Wi(-NbC:NC)
      REAL          :: She(-NbC:NC)
    END SUBROUTINE CalcShear
  END INTERFACE

!--------------------------------[CVS]---------------------------------!
!  $Id: Processor.f90,v 1.92 2009/01/21 10:10:51 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/Processor.f90,v $  
!======================================================================!

  call cpu_time(start)
!---- Test the precision
  open(90,FORM='UNFORMATTED',FILE='Processor.real');
  write(90) 3.1451592
  close(90)
  open(765, FILE='force.dat',POSITION='APPEND') ! file should consist of cells nodes

!//////////////////////////////////!
!     Start parallel execution     !
!//////////////////////////////////!
  call StaPar  

  call Timex(CPUtim) 

!///////////////////////////////!
!     Open the command file     !
!///////////////////////////////!
  open(7, FILE='T-FlowS.cmn')    

  if(this  < 2) then
    call logo
  endif

  call Wait   

!---- initialize parameters
  call IniPar
  
!---- load the finite volume grid      
  call CnsLoa
  call GeoAloc
  call GeoLoa
  call BufLoa
  call Exchng(volume(-NbC))

  call wait   

  call TopolM

  call Wait   

!  if(JET==yes) call Set_BC   ! IZMJENA
!>>>>>>>>>>>>>>>>>>>!
!                   !
!     TIME LOOP     !
!                   !
!<<<<<<<<<<<<<<<<<<<!
  Ndt  = 0
  Ndtt = 0
  call LoaRes(restar)

  if(restar) then
    call BouLoa(.true.)
  end if
 
!----- Read command file (T-FlowS.cmn) 
1 call ReaCom(restar)

!----- Initialize variables
  if(.not. restar) then
    call BouLoa(.true.)
    call IniVar 
    call Wait
  end if   

!----- Interpolate between diff. meshes
  call LoaIni()

  if(.not. restar) then
    do m=1,Nmat
      if(SHAKE(m) == YES) then
        call UserPerturb2(5.0,0,m)
      end if
    end do  
  end if

!----- Check if there are more materials
  multiple = .FALSE.
  i = StateMat(1)
  do m=1,Nmat
    if(StateMat(m) /= i) multiple = .TRUE.
  end do

  if(this  < 2)             &
    write(*,'(A18,15A11)')  &
    '#        N        ',   &
    ' Mass inb. ',          &
    '    U:     ',          &
    '    V:     ',          &
    '    W:     ',          &
    '    P:     '   

  if(this  < 2)             &
    write(*,'(A18,15A11)')  &
    '  CFL max: ',          &
    '  Pe max:  ',          &
    '  Flux x:  ',          &
    '  Flux y:  ',          &
    '  Flux z:  ',          &
    'iterations ',          &
    '  dp/dx:   ',          &
    '    K:     '

!----- Prepare ...
  call Calc3()
  call FindBad()
  if(SIMULA==LES.and.MODE==SMAG.and..NOT.restar) call NearWallCell()

!----- Prepare the gradient matrix for velocities
  call CalcG(.TRUE.) 

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!        LET THE TIME        !
!     INTEGRATION  BEGIN     !
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

  LinMon0 = ' '
  LinMon1 = '# [1]'
  LinMon2 = '# [2]'
  LinMon3 = '# [3]'
  LineRes = '#'

  if(ALGOR  ==  FRACT) call ForAPF() 

!---- Print the areas of monitoring planes
  if(This < 2) then
    do m=1,Nmat
      write(*,'(A5,I2,A2,1PE12.3)') '# Ax(',m,')=', AreaX(m)
      write(*,'(A5,I2,A2,1PE12.3)') '# Ay(',m,')=', AreaY(m)
      write(*,'(A5,I2,A2,1PE12.3)') '# Az(',m,')=', AreaZ(m)
    end do
  end if

  if(Ndt == 0)  goto 6 

  do n=Ndtt+1, Ndt 


    Time = Time + dt
    write( converted_INT, '(I7)' ) n

    if(Cm(1) > 0 .and. n < 1000000 ) then
      write(LinMon0(1: 6), '(I6)')      n;    
      write(LinMon0(7:18), '(1PE12.3)') Time; 
    else if (Cm(1) > 0 .and. n > 1000000 ) then ! need to print only last 6 digits of current n
      write(LinMon0(1: 6),'(A)') converted_INT(2:7)
      write(LinMon0(7:18), '(1PE12.3)') Time; 
    end if

    if(SIMULA==DES_SPA) then
      call CalcShear(U % n, V % n, W % n, Shear)
      call CalcVort(U % n, V % n, W % n, Vort)
    end if

    if(SIMULA == LES) then
      call CalcShear(U % n, V % n, W % n, Shear)
      if(MODE == DYN) call CalcSGS_Dynamic() 
      if(MODE == WALE) call CalcWALE() 
      call CalcSGS()
    end if  
    If(SIMULA==HYB_ZETA) then  
      call CalcSGS_Dynamic()      
      call CalcSGS_hybrid()
    end if

    call CalcConvect
    if(SIMULA==EBM.or.SIMULA==HJ) call CalcVISt_EBM

    do ini=1,Nini                   !  FRACTION & SIMPLE  

      if(.NOT. multiple) then 
        call GradP(P % n,Px,Py,Pz)
      else 
        call GradP3(P % n,Px,Py,Pz)
      end if
      call GraPhi(U % n, 1, Ux,.TRUE.)    ! dU/dx
      call GraPhi(U % n, 2, Uy,.TRUE.)    ! dU/dy
      call GraPhi(U % n, 3, Uz,.TRUE.)    ! dU/dz
      call GraPhi(V % n, 1, Vx,.TRUE.)    ! dV/dx
      call GraPhi(V % n, 2, Vy,.TRUE.)    ! dV/dy
      call GraPhi(V % n, 3, Vz,.TRUE.)    ! dV/dz
      call GraPhi(W % n, 1, Wx,.TRUE.)    ! dW/dx
      call GraPhi(W % n, 2, Wy,.TRUE.)    ! dW/dy 
      call GraPhi(W % n, 3, Wz,.TRUE.)    ! dW/dz
!---- U velocity component --------------------------------!
      call NewUVW(1, U,                 &   
                     Ux,   Uy,   Uz,    & 
                     Sx,   Sy,   Sz,    &                             
                     Dx,   Dy,   Dz,    &                             
                     Px,   Vx,   Wx)      ! dP/dx, dV/dx, dW/dx


!---- V velocity component --------------------------------!
      call NewUVW(2, V,                 &
                     Vy,   Vx,   Vz,    & 
                     Sy,   Sx,   Sz,    &
                     Dy,   Dx,   Dz,    &
                     Py,   Uy,   Wy)      ! dP/dy, dU/dy, dW/dy

!---- W velocity component --------------------------------!
      call NewUVW(3, W,                 & 
                     Wz,   Wx,   Wy,    & 
                     Sz,   Sx,   Sy,    &                         
                     Dz,   Dx,   Dy,    &                         
                     Pz,   Uz,   Vz)      ! dP/dz, dU/dz, dV/dz
      if(ALGOR == FRACT) then
        call Exchng(Asave)  
        call ModOut()
        call CalcPF()
      endif
      if(ALGOR == SIMPLE) then
        call Exchng(Asave)  
        call ModOut()
        call CalcPS()
      end if

      if(.NOT. multiple) then 
        call GradP(PP % n,Px,Py,Pz)
      else 
        call GradP3(PP % n,Px,Py,Pz)
      end if

      call CalcFlux
      Mres = CorUVW() !  project the velocities

!---- Temperature
      if(HOT==YES) then
        call GraPhi(T % n,1,PHIx,.TRUE.)            ! dT/dx
        call GraPhi(T % n,2,PHIy,.TRUE.)            ! dT/dy
        call GraPhi(T % n,3,PHIz,.TRUE.)            ! dT/dz
!!!!        call GraCorNew(T % n,PHIx,PHIy,PHIz)    ! 2mat
        call CalcSc(5, T, PHIx, PHIy, PHIz)  ! dT/dx, dT/dy, dT/dz
      end if 

!---- Rans models
      if(SIMULA==K_EPS.or.SIMULA == HYB_PITM) then
        call CalBou
        call CalcShear(U % n, V % n, W % n, Shear)
        call GraPhi(Kin % n,1,PHIx,.TRUE.)             ! dK/dx
        call GraPhi(Kin % n,2,PHIy,.TRUE.)             ! dK/dy
        call GraPhi(Kin % n,3,PHIz,.TRUE.)             ! dK/dz
        call CalcTurb(6, Kin, PHIx, PHIy, PHIz) ! dK/dx, dK/dy, dK/dz

        call GraPhi(Eps % n,1,PHIx,.TRUE.)           ! dEps/dx
        call GraPhi(Eps % n,2,PHIy,.TRUE.)           ! dEps/dy
        call GraPhi(Eps % n,3,PHIz,.TRUE.)           ! dEps/dz
        call CalcTurb(7, Eps, PHIx, PHIy, PHIz)

        call CalcVISt_KEps()
      end if 

      if(SIMULA==K_EPS_VV.or.SIMULA==ZETA.or.SIMULA == HYB_ZETA) then
        call CalcShear(U % n, V % n, W % n, Shear)
        call GraPhi(Kin % n,1,PHIx,.TRUE.)             ! dK/dx
        call GraPhi(Kin % n,2,PHIy,.TRUE.)             ! dK/dy
        call GraPhi(Kin % n,3,PHIz,.TRUE.)             ! dK/dz
        call CalcTurb(6, Kin, PHIx, PHIy, PHIz) ! dK/dx, dK/dy, dK/dz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !do c2=1,NC
            if ( this < 2 ) write(*,'(A,E20.10)') 'MaxKin=', maxval(Kin % n) 
        !end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

        call GraPhi(Eps % n,1,PHIx,.TRUE.)             ! dEps/dx
        call GraPhi(Eps % n,2,PHIy,.TRUE.)             ! dEps/dy
        call GraPhi(Eps % n,3,PHIz,.TRUE.)             ! dEps/dz
        call CalcTurb(7, Eps, PHIx, PHIy, PHIz)
         
        call CalBou


        call GraPhi(f22 % n,1,PHIx,.TRUE.)             ! df22/dx
        call GraPhi(f22 % n,2,PHIy,.TRUE.)             ! df22/dy
        call GraPhi(f22 % n,3,PHIz,.TRUE.)             ! df22/dz
        call CalcF22(8, f22, PHIx, PHIy, PHIz) 

        call GraPhi(v_2 % n,1,PHIx,.TRUE.)             ! dv_2/dx
        call GraPhi(v_2 % n,2,PHIy,.TRUE.)             ! dv_2/dy
        call GraPhi(v_2 % n,3,PHIz,.TRUE.)             ! dv_2/dz
        call CalcTurb(9, v_2, PHIx, PHIy, PHIz)  
        call CalcVISt_KepsV2F()
      end if                 

      if(SIMULA==EBM.or.SIMULA==HJ) then

        call CalBou

        if(SIMULA==EBM) call Scale  

        call GraPhi(U % n, 1, Ux,.TRUE.)    ! dU/dx
        call GraPhi(U % n, 2, Uy,.TRUE.)    ! dU/dy
        call GraPhi(U % n, 3, Uz,.TRUE.)    ! dU/dz
 
        call GraPhi(V % n, 1, Vx,.TRUE.)    ! dV/dx
        call GraPhi(V % n, 2, Vy,.TRUE.)    ! dV/dy
        call GraPhi(V % n, 3, Vz,.TRUE.)    ! dV/dz

        call GraPhi(W % n, 1, Wx,.TRUE.)    ! dW/dx
        call GraPhi(W % n, 2, Wy,.TRUE.)    ! dW/dy 
        call GraPhi(W % n, 3, Wz,.TRUE.)    ! dW/dz

        call GraPhi(uu % n,1,PHIx,.TRUE.)             ! dK/dx
        call GraPhi(uu % n,2,PHIy,.TRUE.)             ! dK/dy
        call GraPhi(uu % n,3,PHIz,.TRUE.)             ! dK/dz
        call CalcStresses(6, uu, PHIx, PHIy, PHIz) ! dK/dx, dK/dy, dK/dz

        call GraPhi(vv % n,1,PHIx,.TRUE.)             ! dEps/dx
        call GraPhi(vv % n,2,PHIy,.TRUE.)             ! dEps/dy
        call GraPhi(vv % n,3,PHIz,.TRUE.)             ! dEps/dz
        call CalcStresses(7, vv, PHIx, PHIy, PHIz)
         
        call GraPhi(ww % n,1,PHIx,.TRUE.)             ! df22/dx
        call GraPhi(ww % n,2,PHIy,.TRUE.)             ! df22/dy
        call GraPhi(ww % n,3,PHIz,.TRUE.)             ! df22/dz
        call CalcStresses(8, ww, PHIx, PHIy, PHIz) 


        call GraPhi(uv % n,1,PHIx,.TRUE.)             ! dv_2/dx
        call GraPhi(uv % n,2,PHIy,.TRUE.)             ! dv_2/dy
        call GraPhi(uv % n,3,PHIz,.TRUE.)             ! dv_2/dz
        call CalcStresses(9, uv, PHIx, PHIy, PHIz)  

        call GraPhi(uw % n,1,PHIx,.TRUE.)             ! df22/dx
        call GraPhi(uw % n,2,PHIy,.TRUE.)             ! df22/dy
        call GraPhi(uw % n,3,PHIz,.TRUE.)             ! df22/dz
        call CalcStresses(10, uw, PHIx, PHIy, PHIz) 

        call GraPhi(vw % n,1,PHIx,.TRUE.)             ! df22/dx
        call GraPhi(vw % n,2,PHIy,.TRUE.)             ! df22/dy
        call GraPhi(vw % n,3,PHIz,.TRUE.)             ! df22/dz
        call CalcStresses(11, vw, PHIx, PHIy, PHIz) 

        if(SIMULA==EBM) then
          call GraPhi(f22 % n,1,PHIx,.TRUE.)             ! df22/dx
          call GraPhi(f22 % n,2,PHIy,.TRUE.)             ! df22/dy
          call GraPhi(f22 % n,3,PHIz,.TRUE.)             ! df22/dz
          call CalcF22(12, f22, PHIx, PHIy, PHIz) 
        end if 

        call GraPhi(Eps % n,1,PHIx,.TRUE.)             ! df22/dx
        call GraPhi(Eps % n,2,PHIy,.TRUE.)             ! df22/dy
        call GraPhi(Eps % n,3,PHIz,.TRUE.)             ! df22/dz
        call CalcStresses(13, Eps, PHIx, PHIy, PHIz) 
      end if                 

      if(SIMULA==SPA_ALL.or.SIMULA==DES_SPA) then
        call CalcShear(U % n, V % n, W % n, Shear)
        call CalcVort(U % n, V % n, W % n, Vort)

        call GraPhi(VIS % n,1,PHIx,.TRUE.)             ! dVIS/dx
        call GraPhi(VIS % n,2,PHIy,.TRUE.)             ! dVIS/dy
        call GraPhi(VIS % n,3,PHIz,.TRUE.)             ! dVIS/dz
        call CalcTurb(6, VIS, PHIx, PHIy, PHIz)        ! dVIS/dx, dVIS/dy, dVIS/dz
        call CalcVISt_SPA_ALL(n)
      end if

      
      call CalBou !---- Update the values at boundaries

      write(LineRes(2:4),'(I3)') ini

      if(Cm(1) /= 0) write(*,'(A100)') LineRes !----- End of the current iteration 


      if(ALGOR == SIMPLE) then
        if( res(1) <= SIMTol .and. res(2) <= SIMTol .and. &
           res(3) <= SIMTol .and. res(4) <= SIMTol ) goto 4 
      endif
    end do 

!----- End of the current time step
4   if(Cm(1) /= 0) then
      write(*,'(A138)') LinMon0
      do m=1,Nmat
        if(m .eq.1) write(*,'(A138)') LinMon1
        if(m .eq.2) write(*,'(A138)') LinMon2
        if(m .eq.3) write(*,'(A138)') LinMon3
      end do
    end if

!----- Write the values in monitoring points
    do i=1,Nmon
      if(Cm(i)  > 0) then                               
        if(HOT==NO) then
          write(10+i,'(I9,4E16.6)')                    &
           n, U % n(Cm(i)), V%n(Cm(i)), W%n(Cm(i)), P%n(Cm(i))
        else
          write(10+i,'(I9,5E16.6)')                    &
           n, U % n(Cm(i)), V%n(Cm(i)), W%n(Cm(i)), P%n(Cm(i)), T%n(Cm(i))
        end if
      end if
    end do  
!----- Write the separation point for monitoring
    
  !call TraceSeparationPoint_new(n)
  
   if(PIPE==YES.or.JET==YES) then 
     !call CalcMn_Cylind(Nstat, n)  !  calculate mean values 
     !if(BUDG == YES.and.HOT==YES) call CalcBudgets_cylind(Nbudg, n)
     !if(BUDG == YES.and.HOT==NO)  call CalcBudgets_cylind(Nbudg, n)
   else
     call CalcMn(Nstat, n, restar)  !  calculate mean values 
     call Wait
   end if

!-----------------------------------------------------!  
!                                                     ! 
!   Recalculate the pressure drop                     !
!   to keep the constant mass flux                    !
!                                                     !
!   First Newtons law:                                !
!   ~~~~~~~~~~~~~~~~~~                                !
!   F = m * a                                         !
!                                                     !
!   where:                                            !
!   ~~~~~~                                            !
!   a = dv / dt = dFlux / dt * 1 / (A * rho)          !
!   m = rho * V                                       !
!   F = Pdrop * l * A = Pdrop * V                     !
!                                                     !
!   finally:                                          !
!   ~~~~~~~~                                          !
!   Pdrop * V = rho * V * dFlux / dt * 1 / (A * rho)  !
!                                                     !
!   after cancelling: V and rho, it yields:           !
!                                                     !
!   Pdrop = dFlux/dt/A                                !
!                                                     !
!-----------------------------------------------------!
    do m=1,Nmat
      if(SHAKE(m) == YES) then
        if( n  < SHAKE_PER(m) .and. MOD(n+1,SHAKE_INT(m)) == 0 ) then 
          call UserPerturb2(5.0,n,m)
        endif 
      endif

      if( FLUXoX(m)  /=  0.0 ) then
        PdropX(m) = (FLUXoX(m)-FLUXx(m)) / (dt*AreaX(m)+TINY) 
      end if
      if( FLUXoY(m)  /=  0.0 ) then
        PdropY(m) = (FLUXoY(m)-FLUXy(m)) / (dt*AreaY(m)+TINY) 
      end if
      if( FLUXoZ(m)  /=  0.0 ) then
        PdropZ(m) = (FLUXoZ(m)-FLUXz(m)) / (dt*AreaZ(m)+TINY) 
      end if
    end do

!---- Regular restart files savings, each 1000 time steps            
    if(mod(n,6000) == 0) then
      Ndtt_temp = Ndtt
      Ndtt      = n
      namSav = 'SAVExxxxxx'
      if ( n < 1000000 ) then 
        write(namSav(5:10),'(I6.6)') n
      else
        write(namSav(5:10),'(A)') converted_INT(2:7)
      end if 
        call SavRes(namSav)
        !call SaveYplus(Nstat, n) 
      Ndtt = Ndtt_temp
    end if

!---- Is user forcing it to extract function ?                    
    open(1, file='data_now', err=7, status='old')
      call Wait
      if(This < 2) close(1, status='delete')
      Ndtt_temp = Ndtt
      Ndtt      = n
      namSav = 'SAVExxxxxx'
      if ( n < 1000000 ) then 
        write(namSav(5:10),'(I6.6)') n
      else
        write(namSav(5:10),'(A)') converted_INT(2:7)
      end if
      call ExtractFunc('coarse_far')
      call Wait
      call ExtractFunc('coarse_near')
      call Wait
      !call SaveYplus(Nstat, n) 
7   continue

!---- Is user forcing it to save ?                    
    open(1, file='save_now', err=8, status='old')
      call Wait
      if(This < 2) close(1, status='delete')
      Ndtt_temp = Ndtt
      Ndtt      = n
      namSav = 'SAVExxxxxx'
      if ( n < 1000000 ) then 
        write(namSav(5:10),'(I6.6)') n
      else
        write(namSav(5:10),'(A)') converted_INT(2:7)
      end if
      call SavRes(namSav)
      call SavParView(this,NC,namSav)
      call UserCutLines_channel(yc)
      !call SaveYplus(Nstat, n) 
      Ndtt = Ndtt_temp
8   continue


  end do                    ! n, number of time steps  

  if(This < 2) then
    open(9,FILE='stop')
    close(9)
  end if

5 Ndtt = n - 1                                 

!--------------------------!
!     Save the results     !
!--------------------------!
6 call SavRes                           
  call SavIni
  call ProSav
  call DatSav
  namSav = 'SAVExxxxxx'
  write(namSav(5:10),'(I6.6)') n
  call SavParView(this,NC,namSav)
  !call SaveYplus(Nstat, n) 
  call Wait

!---- Create a cut line with normalise a values
  if(CHANNEL == YES) then
    !call UserCutLines_channel(yc)
  else if(PIPE == YES) then
    call UserCutLines_pipe
    if(BUDG==YES.and.HOT==NO) call UserCutLines_budgets_cylind 
    if(BUDG==YES.and.HOT==YES) call UserCutLines_budgets_cylind_HOT 
!    call UserCutLines_annulus
  else if(JET == YES) then
    call UserCutLines_Nu
    call UserProbe1D_jet_rans
    call UserProbe1D_Nusselt_jet
  end if

 
!  call UserDiffuser 

3 if(this  < 2) write(*,*) '# Exiting !'

!////////////////////////////////!
!     Close the command file     !
!////////////////////////////////!
  close(7)                


!////////////////////////////////////!
!     Close the monitoring files     !
!////////////////////////////////////!
  do n=1,Nmon
    if(Cm(n)  > 0) close(10+n)
  end do

!////////////////////////////////!
!     End parallel execution     !
!////////////////////////////////!
  call endpar            
  call cpu_time(finish)
  if(this < 2) print '("Time = ",f14.3," seconds.")',finish-start

  END PROGRAM Processor  
