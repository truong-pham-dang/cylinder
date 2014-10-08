!======================================================================!
  SUBROUTINE ReadFLUENT(FLUfile)
!----------------------------------------------------------------------!
! Purpose: Reads NAME.gmv file to initialize variables                 !
! ~~~~~~~~                                                             !
!-----------------------------[Modules]--------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!---------------------------[Parameters]-------------------------------!
  CHARACTER :: FLUfile*20
!----------------------------[Calling]---------------------------------!
  REAL      :: Dist
!-----------------------------[Locals]---------------------------------!
  CHARACTER :: answer*20 !, Line*300, token*130
  INTEGER   :: i, j, c, s, c1, c2, dum1, dum2, error
  INTEGER   :: Ncell, Nnode, CellNodes, Node(10)
  REAL      :: Mres, T1, T2, L1, L2        !,dummy
  REAL,ALLOCATABLE :: Xnode(:), Ynode(:), Znode(:)
  REAL,ALLOCATABLE :: Xcent(:), Ycent(:), Zcent(:),     &
                      Ucent(:), Vcent(:), Wcent(:),     &
                      Pcent(:), Tcent(:), Kcent(:),     &
                      Ecent(:), Pkcent(:), VVcent(:),   &
                      F22cent(:), ZETAcent(:), VIScent(:)
!------------------------------[CVS]-----------------------------------!
!  $Id: ReadGMV.f90,v 1.3 2003/01/08 16:22:11 mirza Exp $
!  $Source: /home/mirza/.CVSROOT_mirza/T-Rex/Process/ReadGMV.f90,v $
!======================================================================!



!----------------------!
!-- Coordinates Loop --!
!----------------------!
  if(this < 2) write(*,*) '# NOW READING FILE: ', FLUfile
  open (5, file=FLUfile)

  do i=1,6                             !header
    read(5,*) inp                      !header
  end do                               !header

!---- read the line with usefull information
  call ReadC(5,inp,tn,ts,te)
  read(inp(ts(1):te(1)),*) Nnode
  read(inp(ts(2):te(2)),*) Ncell
  read(5,*) inp                        !ENDOFSECTION
!---- allocation
  allocate(Xnode(Nnode)) ; Xnode = 0.0
  allocate(Ynode(Nnode)) ; Ynode = 0.0
  allocate(Znode(Nnode)) ; Znode = 0.0

  allocate(Xcent(Ncell)) ; Xcent = 0.0
  allocate(Ycent(Ncell)) ; Ycent = 0.0
  allocate(Zcent(Ncell)) ; Zcent = 0.0

!---- read the nodal coordinates
  read(5,*) inp                        !NODAL COORDINATES
  do i=1,Nnode
    call ReadC(5,inp,tn,ts,te)
    read(inp(ts(2):te(2)),*) Xnode(i)
    read(inp(ts(3):te(3)),*) Ynode(i)
    read(inp(ts(4):te(4)),*) Znode(i)
  end do
  read(5,*) inp                        !ENDOFSECTION

!---- read nodes of each cell
  read(5,*) inp                        !ELEMENTS/CELLS
  do i=1,Ncell
    read(5,*) dum1, dum2, CellNodes, (Node(j), j=1,CellNodes)
    do j = 1, CellNodes
      Xcent(i) = Xcent(i) + Xnode(Node(j)) / CellNodes
      Ycent(i) = Ycent(i) + Ynode(Node(j)) / CellNodes
      Zcent(i) = Zcent(i) + Znode(Node(j)) / CellNodes
    end do
  end do
  read(5,*) inp                        !ENDOFSECTION

  close(5)

!--------------------!
!-- Variables Loop --!
!--------------------!
  FLUfile(len_trim(FLUfile)-3 : len_trim(FLUfile)) = ".dat"
  if(this < 2) write(*,*) '# NOW READING FILE: ', FLUfile
  open (5, file=FLUfile)

  do i=1,4                             !header
    read(5,*) inp                      !header
  end do                               !header

!-----Velocities
  allocate(Ucent(Ncell)); Ucent = 0.0
  allocate(Vcent(Ncell)); Vcent = 0.0
  allocate(Wcent(Ncell)); Wcent = 0.0
  read(5,*) inp                        !velocities
  read(5,*) inp                        !header begin
  do i=1,Ncell
    call ReadC(5,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) Ucent(i)
    read(inp(ts(2):te(2)),*) Vcent(i)
    read(inp(ts(3):te(3)),*) Wcent(i)
  end do
  read(5,*) inp                        !header end

!-----Other variables
2 CONTINUE
  read(UNIT=5, IOSTAT=error, FMT=*) inp

!    read(5,*)                       !header
!    read(5,*) inp                      !header
!    read(5,*)                       !header
  if(error == 0) then
    backspace(5)
    call ReadC(5,inp,tn,ts,te)
!   if(inp(ts(tn):te(tn)) /= '")') then
!   read(inp(ts(tn):te(tn)),*) answer
!   call ToUppr(answer)

!--Pressure
    if(inp == '(0 "pressure")') then
      allocate(Pcent(Ncell)); Pcent = 0.0
      read(5,*) inp                    !header begin
      do i=1,Ncell
        call ReadC(5,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) Pcent(i)
      end do
      read(5,*) inp                    !header end

!--Temperature
    else if(inp == '(0 "temperature")') then
      allocate(Tcent(Ncell)); Tcent = 0.0
      read(5,*) inp                    !header begin
      do i=1,Ncell
        call ReadC(5,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) Tcent(i)
      end do
      read(5,*) inp                    !header end

!--Kin
    else if(inp == '(0 "Kin")') then
      allocate(Kcent(Ncell)); Kcent = 0.0
      read(5,*) inp                    !header begin
      do i=1,Ncell
        call ReadC(5,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) Kcent(i)
      end do
      read(5,*) inp                    !header end

!--Eps
    else if(inp == '(0 "Eps")') then
      allocate(Ecent(Ncell)); Ecent = 0.0
      read(5,*) inp                    !header begin
      do i=1,Ncell
        call ReadC(5,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) Ecent(i)
      end do
      read(5,*) inp                    !header end

!--Pk
    else if(inp == '(0 "Pk")') then
      allocate(Pkcent(Ncell)); Pkcent = 0.0
      read(5,*) inp                    !header begin
      do i=1,Ncell
        call ReadC(5,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) Pkcent(i)
      end do
      read(5,*) inp                    !header end

!--v_2
    else if(inp == '(0 "v_2")') then
      allocate(VVcent(Ncell)); VVcent = 0.0
      read(5,*) inp                    !header begin
      do i=1,Ncell
        call ReadC(5,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) VVcent(i)
      end do
      read(5,*) inp                    !header end

!--zeta
    else if(inp == '(0 "zeta")') then
      allocate(ZETAcent(Ncell)); ZETAcent = 0.0
      read(5,*) inp                    !header begin
      do i=1,Ncell
        call ReadC(5,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) ZETAcent(i)
      end do
      read(5,*) inp                    !header end

!--f22
    else if(inp == '(0 "f22")') then
      allocate(F22cent(Ncell)); F22cent = 0.0
      read(5,*) inp                    !header begin
      do i=1,Ncell
        call ReadC(5,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) F22cent(i)
      end do
      read(5,*) inp                    !header end

!--Vis
    else if(inp == '(0 "VIS")') then
      allocate(VIScent(Ncell)); VIScent = 0.0
      read(5,*) inp                    !header begin
      do i=1,Ncell
        call ReadC(5,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) VIScent(i)
      end do
      read(5,*) inp                    !header end

    end if !end if answer==K,E,V2...
!   end if   !end if inp=='")'

    GOTO 2   !go back to read new variables

  else       !else if err/=0

    write(*,*) '# Finished reading file: ', FLUfile

    GOTO 4   !go to close variables loop

  end if     !end if err==0

4 CONTINUE
  close(5)



!-----Initialization
  do c=1, NC
    Mres = HUGE
    do i = 1, Ncell
      if(Dist(Xcent(i),Ycent(i),Zcent(i),xc(c),yc(c),zc(c)) < Mres) then
        Mres = Dist(Xcent(i),Ycent(i),Zcent(i),xc(c),yc(c),zc(c))
        j = i  !closest cell
      end if
    end do

    U % mean(c) = 0.0
    V % mean(c) = 0.0
    W % mean(c) = 0.0
    U % n(c)    = Ucent(j)
    U % o(c)    = U % n(c)
    U % oo(c)   = U % n(c)
    V % n(c)    = Vcent(j)
    V % o(c)    = V % n(c)
    V % oo(c)   = V % n(c)
    W % n(c)    = Wcent(j)
    W % o(c)    = W % n(c)
    W % oo(c)   = W % n(c)

    if(ALLOCATED(Pcent)) then
      P % n(c) = Pcent(j)
    else
      if(c==2) write(*,*) '@ReadGMV: P not found in ', FLUfile,  &
                          ' Initialization with constant (0.0)!'
      P % n(c) = 0.0
    end if

    if(HOT==YES) then
      if(ALLOCATED(Tcent)) then
        T % n(c) = Tcent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: T not found in ', FLUfile, &
                            ' Initialization with constant (20.0)!'
        T % n(c) = 20.0  !or 20.0
      end if
      T % o(c)  = T % n(c)
      T % oo(c) = T % n(c)
    end if

    if(SIMULA==K_EPS.or.SIMULA==HYB_PITM) then
      if(ALLOCATED(Kcent)) then
        Kin % n(c) = Kcent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: Kin not found in ', FLUfile, &
                            ' Initialization with constant (1.0E-3)!'
        Kin % n(c) = 1.0E-3
      end if
      Kin % o(c)  = Kin % n(c)
      Kin % oo(c) = Kin % n(c)

      if(ALLOCATED(Ecent)) then
        Eps % n(c) = Ecent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: Eps not found in ', FLUfile, &
                            ' Initialization with constant (1.0E-4)!'
        Eps % n(c)  = 1.0E-4
      end if
      Eps % o(c)  = Eps % n(c)
      Eps % oo(c) = Eps % n(c)

      if(ALLOCATED(Pkcent)) then
        Pk(c) = Pkcent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: Pk not found in ', FLUfile,  &
                            ' Initialization with constant (1.0E-4)!'
        Pk(c)  = 1.0E-4
      end if

      VISt(c) = Cmu * Kin%n(c)*Kin%n(c) / Eps%n(c)
    end if

    if(SIMULA==K_EPS_VV.or.SIMULA==ZETA.or.SIMULA==HYB_ZETA) then
      if(ALLOCATED(Kcent)) then
        Kin % n(c) = Kcent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: Kin not found in ', FLUfile, &
                            ' Initialization with constant (1.0E-3)!'
        Kin % n(c) = 1.0E-3
      end if
      Kin % o(c)  = Kin % n(c)
      Kin % oo(c) = Kin % n(c)

      if(ALLOCATED(Ecent)) then
        Eps % n(c) = Ecent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: Eps not found in ', FLUfile, &
                            ' Initialization with constant (1.0E-4)!'
        Eps % n(c) = 1.0E-4
      end if
      Eps % o(c)  = Eps % n(c)
      Eps % oo(c) = Eps % n(c)

      if(ALLOCATED(Pkcent)) then
        Pk(c) = Pkcent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: Pk not found in ', FLUfile,  &
                            ' Initialization with constant (1.0E-4)!'
        Pk(c)  = 1.0E-4
      end if

      if(ALLOCATED(VVcent)) then
        v_2 % n(c) = VVcent(j) 
      else
        if(c==2) write(*,*) '@ReadGMV: zeta not found in ', FLUfile,  &
                            ' Initialization with constant (1.0E-2)!'
        v_2 % n(c) = 1.0E-2
      end if
      v_2 % o(c)  = v_2 % n(c)
      v_2 % oo(c) = v_2 % n(c)

      if(ALLOCATED(F22cent)) then
        f22 % n(c) = F22cent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: f22 not found in ', FLUfile, &
                            ' Initialization with constant (1.0E-2)!'
        f22 % n(c) = 1.0E-2
      end if
      f22 % o(c)  = f22 % n(c)
      f22 % oo(c) = f22 % n(c)

!      Uf(c)  = 0.047
!      Ynd(c) = 30.0
      T1 = Kin%n(c) / (Eps%n(c) + TINY)                                          !ZETA!--New scales
      T2 = 1.0 * Ct * ( VISc / (Eps%n(c) + TINY) )**0.5
      L1 = Cl * Kin%n(c)**1.5 / (Eps%n(c)+TINY)                                  !ZETA!--New scales
      L2 = Cl * Cni * ( VISc**3.0 / (Eps%n(c)+TINY) )**0.25
      Tsc(c) = max(T1, T2)
      Lsc(c) = max(L1, L2)
      if(SIMULA==ZETA) then
        VISt(c) = CmuD * v_2%n(c) * Kin%n(c) * Tsc(c)
      else if(SIMULA==HYB_ZETA) then
        VISt(c) = CmuD*v_2%n(c)*Kin % n(c)*Tsc(c)
        VISt_eff(c) = max(VISt(c),VISt_sgs(c))
      else if(SIMULA==K_EPS_VV) then
        VISt(c) = CmuD * v_2%n(c) * Tsc(c)
      end if
    end if

    if(SIMULA == SPA_ALL .or. SIMULA==DES_SPA) then
      if(ALLOCATED(VIScent)) then
        VIS % n(c)  = VIScent(j)
        VIS % o(c)  = VIScent(j)
        VIS % oo(c) = VIScent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: viscosity not found in ', FLUfile
        if(c==2) write(*,*) '#         Viscosity initialized with constant!!'
        VIS % n(c)  = 1.0
        VIS % o(c)  = 1.0
        VIS % oo(c) = 1.0
      end if
    end if

!-----Status printout
    if(mod(c,5000) == 0) write(*,*) 100.0 * c / (1.0 * NC), '% complete...'

  end do   !end do=1,NC
!-----Fluxes, convective and diffusive terms
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    Flux(s) = f(s) * DENc(material(c1)) * (U%n(c1)*Sx(s)+V%n(c1)*Sy(s)+W%n(c1)*Sz(s))     &
            + (1.0-f(s)) * DENc(material(c2)) * (U%n(c2)*Sx(s)+V%n(c2)*Sy(s)+W%n(c2)*Sz(s))
  end do
  
  call IniCoDoXo(U)
  call IniCoDoXo(V)
  call IniCoDoXo(W)
  if(SIMULA==K_EPS.or.SIMULA==HYB_PITM) then
    call IniCoDoXo(Kin)
    call IniCoDoXo(Eps)
  end if
  if(SIMULA==K_EPS_VV.or.SIMULA==ZETA.or.SIMULA==HYB_ZETA) then
    call IniCoDoXo(Kin)
    call IniCoDoXo(Eps)
    call IniCoDoXo(V_2)
    call IniCoDoXo(F22)
    f22%C = 0.0
    f22%Co = 0.0
  end if
!----------------------------------------------------------------------!

!-----Deallocation
  deallocate(Xnode)
  deallocate(Ynode)
  deallocate(Znode)
  deallocate(Xcent)
  deallocate(Ycent)
  deallocate(Zcent)
  if(ALLOCATED(Ucent)) deallocate(Ucent)
  if(ALLOCATED(Vcent)) deallocate(Vcent)
  if(ALLOCATED(Wcent)) deallocate(Wcent)
  if(ALLOCATED(Pcent)) deallocate(Pcent)
  if(ALLOCATED(Tcent)) deallocate(Tcent)
  if(ALLOCATED(Kcent)) deallocate(Kcent)
  if(ALLOCATED(Ecent)) deallocate(Ecent)
  if(ALLOCATED(Pkcent)) deallocate(Pkcent)
  if(ALLOCATED(VVcent)) deallocate(VVcent)
  if(ALLOCATED(F22cent)) deallocate(F22cent)
  if(ALLOCATED(VIScent)) deallocate(VIScent)
  call wait
6 CONTINUE

  RETURN
  END SUBROUTINE ReadFLUENT
!======================================================================!
