!======================================================================!
  SUBROUTINE DatSav(namAut)
!----------------------------------------------------------------------!
! Writes: NAME.dat                                                     !
! ~~~~~~~                                                              !
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE par_mod
  USE les_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER             :: c,  c2,  n, s, Nadd
  INTEGER             :: Nfac(10), NtotFac
  CHARACTER           :: answer*80, namOut*80, namTem*80
  CHARACTER, OPTIONAL :: namAut*(*)
  REAL                :: R, RR
!--------------------------------[CVS]---------------------------------!
!  $Id: DatSav.f90,v 1.30 2009/06/30 12:07:19 IUS\mhadziabdic Exp $
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/DatSav.f90,v $           
!======================================================================!
!   See also: number                                                   !
!----------------------------------------------------------------------!

!---- store the name
  namTem = name     

  if(PRESENT(namAut)) then
    if(this  < 2) write(*,*) namAut
    name = namAut  
  else
    if(this  < 2)  &
      write(*,*) '# Input result file name [skip cancels]:'
    call ReadC(7,inp,tn,ts,te)  
    read(inp(ts(1):te(1)),'(A80)')  name
    answer=name
    call ToUppr(answer) 
    if(answer == 'SKIP') then
      name = namTem  
      return
    end if 
  end if

  call wait 

!<<<<<<<<<<<<<<<<<<<<<<<<<!
!                         !
!     create DAT file     !
!                         !
!<<<<<<<<<<<<<<<<<<<<<<<<<!
  call NamFil(THIS, namOut, '.dat', len_trim('.dat'))
  open(9, FILE=namOut)
  if(this  < 2) write(*,*) '# Now creating the file:', namOut  

    call GraPhi(U % n, 1, Ux,.TRUE.)    ! dU/dx
    call GraPhi(V % n, 1, Vx,.TRUE.)    ! dU/dx
    call GraPhi(W % n, 1, Wx,.TRUE.)    ! dU/dx
    call GraPhi(U % n, 2, Uy,.TRUE.)    ! dU/dx
    call GraPhi(V % n, 2, Vy,.TRUE.)    ! dU/dx
    call GraPhi(W % n, 2, Wy,.TRUE.)    ! dU/dx
    call GraPhi(U % n, 3, Uz,.TRUE.)    ! dU/dx
    call GraPhi(V % n, 3, Vz,.TRUE.)    ! dU/dx
    call GraPhi(W % n, 3, Wz,.TRUE.)    ! dU/dx

!---------------!
!     start     !
!---------------!
  write(9,'(A34)') '(0 "============================")'
  write(9,'(A34)') '(0 "Created by T-Rex - Processor")'
  write(9,'(A34)') '(0 "============================")'
  Nadd = 0
  if(SIMULA == LES .or. SIMULA == DES_SPA) then
    Nadd = Nadd + 9   ! Um, Vm, Wm, uu, vv, ww, uv, uw, vw
  end if
  if(HOT == YES) then
    Nadd = Nadd + 1
    if(SIMULA == LES .or. SIMULA == DES_SPA) then
      Nadd = Nadd + 5 ! Tm, TT, uT, vT, wT
    end if
  endif
  write(9,'(A3, I3, A2)') '(0 ', Nadd, ' )'

!-----------------!
!     results     !
!-----------------+----------------------------------------------!
!  (300 (sub_section_id  zone_id  size  n_time_levels  n_phases  !
!        first_id  last_id)                                      !         
!----------------------------------------------------------------!

!====================!
!     velocities     !
!====================!
  write(9,'(A16)') '(0 "velocities")'
  write(9,'(A16,I9,I9,A2)') '(300 (2 1 3 0 0 ',  1, NC, ')(' 
  do c=1,NC
    write(9,'(3F14.6)') U % n(c), V % n(c), W % n(c) 
  end do  
  write(9,'(A2)') '))'

!==================!
!     pressure     !
!==================!
  call DatSavSc('pressure', 1, P % n)

!===========================!
!     aditional scalars     !
!===========================!
  Nadd = 0 
  R    = 0.0  
  RR   = 0.0  
!----------------------!
!     other scalars    !
!----------------------!
  if(TGV == YES) then
    do c =  1, NC
      R   = R + (U % n(c) - (-sin(xc(c))*cos(yc(c))*exp(-2.0*VISc*Time)))**2.0
      RR  = RR + U % n(c)*U % n(c)
      PP % n(c) = sqrt(R/RR) 
    end do
    Nadd = Nadd + 1
    if(this  < 2) write(*,'(A,I3,A,E18.9)') '# Scalar ', Nadd, ' is L2 error', sqrt(R/RR) 
    call DatSavSc('temperature', 699+Nadd, PP % n)
  end if

  if(HOT == YES) then
    Nadd = Nadd + 1
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is temperature' 
    call DatSavSc('temperature', 699+Nadd, T % n)
  end if
  
  if(SIMULA == K_EPS) then
!-- VIS  
    Nadd = Nadd + 1
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is K' 
    call DatSavSc('Kin', 699+Nadd, Kin % n)
!-- VISt  
    Nadd = Nadd + 1
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is Eps' 
    call DatSavSc('Eps', 699+Nadd, Eps % n)
  end if


  if(SIMULA == K_EPS_VV.or.SIMULA==ZETA.or.SIMULA==HYB_ZETA) then
!-- VIS  
    Nadd = Nadd + 1
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is K' 
    call DatSavSc('Kin', 699+Nadd, Kin % n)
!-- VISt  
    Nadd = Nadd + 1
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is Eps' 
    call DatSavSc('Eps', 699+Nadd, Eps % n)
!-- Vort  
    Nadd = Nadd + 1
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is v_2' 
    call DatSavSc('v_2', 699+Nadd, v_2 % n)
!-- WallDs
    Nadd = Nadd + 1
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is f22' 
    call DatSavSc('f22', 699+Nadd, f22 % n)

  end if ! SIMULA=K_EPS_VV

  if(SIMULA == SPA_ALL) then
!-- VIS  
    Nadd = Nadd + 1
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is VIS' 
    call DatSavSc('VIS', 699+Nadd, VIS % n)
!-- VISt  
    Nadd = Nadd + 1
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is VISt' 
    call DatSavSc('VISt', 699+Nadd, VISt)
!-- Vort  
    Nadd = Nadd + 1
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is Vort' 
    call DatSavSc('Vort', 699+Nadd, Vort)
!-- WallDs
    Nadd = Nadd + 1
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is WallDs' 
    call DatSavSc('WallDs', 699+Nadd, WallDs)
  end if ! SIMULA=SPA_ALL

  if(SIMULA == LES .or. SIMULA == DES_SPA) then
!-- Umean
    Nadd = Nadd + 1
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is Umean' 
    call DatSavSc('Umean', 699+Nadd, U % mean)
!-- Vmean
    Nadd = Nadd + 1
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is Vmean' 
    call DatSavSc('Vmean', 699+Nadd, V % mean)
!-- Wmean
    Nadd = Nadd + 1
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is Wmean' 
    call DatSavSc('Wmean', 699+Nadd, W % mean)
!-- uu
    Nadd = Nadd + 1
    PP % n = uu % mean - U % mean * U % mean
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is uu' 
    call DatSavSc('uu', 699+Nadd, PP % n)
!-- vv
    Nadd = Nadd + 1
    PP % n = vv % mean - V % mean * V % mean
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is vv' 
    call DatSavSc('vv', 699+Nadd, PP % n) 
!-- ww
    Nadd = Nadd + 1
    PP % n = ww % mean - W % mean * W % mean
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is ww' 
    call DatSavSc('ww', 699+Nadd, PP % n)
!-- uv
    Nadd = Nadd + 1
    PP % n = uv % mean - U % mean * V % mean
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is uv' 
    call DatSavSc('uv', 699+Nadd, PP % n)
!-- uw
    Nadd = Nadd + 1
    PP % n = uw % mean - U % mean * W % mean
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is uw' 
    call DatSavSc('uw', 699+Nadd, PP % n) 
!-- vw
    Nadd = Nadd + 1
    PP % n = vw % mean - V % mean * W % mean
    if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is ww' 
    call DatSavSc('vw', 699+Nadd, PP % n)

    if(HOT == YES) then
!---- Tmean
      Nadd = Nadd + 1
      if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is Tmean' 
      call DatSavSc('Tmean', 699+Nadd, T % mean)
!---- TT    
      Nadd = Nadd + 1
      PP % n = TT % mean - T % mean * T % mean
      if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is TT' 
      call DatSavSc('TT', 699+Nadd, PP % n)
!---- uT    
      Nadd = Nadd + 1
      PP % n = uT % mean - U % mean * T % mean
      if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is uT' 
      call DatSavSc('uT', 699+Nadd, PP % n)
!---- vT    
      Nadd = Nadd + 1
      PP % n = vT % mean - V % mean * T % mean
      if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is vT' 
      call DatSavSc('vT', 699+Nadd, PP % n)
!---- wT    
      Nadd = Nadd + 1
      PP % n = wT % mean - W % mean * T % mean
      if(this  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is wT' 
      call DatSavSc('wT', 699+Nadd, PP % n)
    end if
  end if ! SIMULA == LES .or. SIMULA == DES_SPA

!---------------------------------------!
!     this is important for connect     !
!---------------------------------------!
  NtotFac = 0
  do n=1,10   ! Browse through boundary condition types
    Nfac(n) = 0
    do s=1,NS   ! Count the faces with boundary condition "n"
      c2 = SideC(2,s)
      if(c2 < 0) then
        if(BCmark(c2) == n) Nfac(n)=Nfac(n)+1
      end if
    end do    ! sides

!---- prepare for next boundary
    if(this  < 2) write(*,*) 'Number of faces:', Nfac(n), NtotFac+1, NtotFac+Nfac(n)
    NtotFac = NtotFac+Nfac(n)
  end do   ! n -> boundary condition types

  write(9,'(I8)') NC
  do n=1,10
    write(9,'(I8)') Nfac(n)
  end do

  close(9)

!---- restore the name
  name = namTem  

  END SUBROUTINE DatSav
