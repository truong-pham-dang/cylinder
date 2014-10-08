!======================================================================!
  SUBROUTINE ReadGMV(GMVfile)
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
  CHARACTER :: GMVfile*20
!----------------------------[Calling]---------------------------------!
  REAL      :: Dist
!-----------------------------[Locals]---------------------------------!
  CHARACTER :: answer*20
  INTEGER   :: i, j, c, s, c1,c2
  INTEGER   :: CellNumb, Nnodes, CellNodes, Node(10)
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


!-----Function---------------------------------------------------------!
!  do c=1,NC
!    U % mean(c) = 0.0
!    V % mean(c) = 0.0
!    W % mean(c) = 0.0
!    U % n(c)    = sqrt( abs((0.0345/zc(c))*(0.0345/zc(c))-yc(c)*yc(c)) )
!    U % o(c)    = U % n(c)
!    U % oo(c)   = U % n(c)
!    V % n(c)    = sqrt( abs((0.0345/zc(c))*(0.0345/zc(c))-xc(c)*xc(c)) )
!    V % o(c)    = V % n(c)
!    V % oo(c)   = V % n(c)
!    W % n(c)    = -0.0345/sqrt(xc(c)*xc(c)+yc(c)*yc(c))
!    W % o(c)    = W % n(c)
!    W % oo(c)   = W % n(c)
!    Kin % n(c)  = 1.0E-3
!    Kin % o(c)  = 1.0E-3
!    Kin % oo(c) = 1.0E-3
!    Eps % n(c)  = 1.0E-4
!    Eps % o(c)  = 1.0E-4
!    Eps % oo(c) = 1.0E-4
!  end do
!  GOTO 6
!----------------------------------------------------------------------!


!-----Input file-------------------------------------------------------!
  if(this < 2) write(*,*) '# NOW READING FILE: ', GMVfile
  open (5, file=GMVfile)

!-----Header
  call ReadC(5,inp,tn,ts,te)

!-----Nodes
  call ReadC(5,inp,tn,ts,te)
  read(inp(ts(2):te(2)),*) Nnodes

  allocate(Xnode(Nnodes)) ; Xnode = 0.0
  allocate(Ynode(Nnodes)) ; Ynode = 0.0
  allocate(Znode(Nnodes)) ; Znode = 0.0

  do i = 1, Nnodes
    call ReadC(5,inp,tn,ts,te)
    read(inp,*) Xnode(i)
  end do
  do i = 1, Nnodes
    call ReadC(5,inp,tn,ts,te)
    read(inp,*) Ynode(i)
  end do
  do i = 1, Nnodes
    call ReadC(5,inp,tn,ts,te)
    read(inp,*) Znode(i)
  end do

!-----Cells
  call ReadC(5,inp,tn,ts,te)
  read(inp(ts(2):te(2)),*) CellNumb

  allocate(Xcent(CellNumb)) ; Xcent = 0.0
  allocate(Ycent(CellNumb)) ; Ycent = 0.0
  allocate(Zcent(CellNumb)) ; Zcent = 0.0

  do i = 1, CellNumb
    call ReadC(5,inp,tn,ts,te)
    read(inp(ts(2):te(2)),*) CellNodes

    call ReadC(5,inp,tn,ts,te)
    do j = 1, CellNodes
      read(inp(ts(j):te(j)),*) Node(j)
    end do

    do j = 1, CellNodes
      Xcent(i) = Xcent(i) + Xnode(Node(j)) / CellNodes
      Ycent(i) = Ycent(i) + Ynode(Node(j)) / CellNodes
      Zcent(i) = Zcent(i) + Znode(Node(j)) / CellNodes
    end do
  end do

!-----Variables Loop
2 CONTINUE
  call ReadC(5,inp,tn,ts,te)
  read(inp(ts(1):te(1)),*) answer
  call ToUppr(answer)
  if(answer == 'ENDVARS') GOTO 4

!-----Velocities
  if(answer == 'VELOCITY') then
    allocate(Ucent(CellNumb)); Ucent = 0.0
    allocate(Vcent(CellNumb)); Vcent = 0.0
    allocate(Wcent(CellNumb)); Wcent = 0.0

    do i = 1, CellNumb
      call ReadC(5,inp,tn,ts,te)
      read(inp,*) Ucent(i)
    end do
    do i = 1, CellNumb
      call ReadC(5,inp,tn,ts,te)
      read(inp,*) Vcent(i)
    end do
    do i = 1, CellNumb
      call ReadC(5,inp,tn,ts,te)
      read(inp,*) Wcent(i)
    end do
  end if

!-----Pressure
  if(answer == 'P') then
    allocate(Pcent(CellNumb)); Pcent = 0.0
    do i = 1, CellNumb
      call ReadC(5,inp,tn,ts,te)
      read(inp,*) Pcent(i)
    end do
  end if

!-----Temperature
  if(answer == 'T') then
    allocate(Tcent(CellNumb)); Tcent = 0.0
    do i = 1, CellNumb
      call ReadC(5,inp,tn,ts,te)
      read(inp,*) Tcent(i)
    end do
  end if

!-----Kin
  if(answer == 'KIN') then
    allocate(Kcent(CellNumb)); Kcent = 0.0
    do i = 1, CellNumb
      call ReadC(5,inp,tn,ts,te)
      read(inp,*) Kcent(i)
    end do
  end if

!-----Eps
  if(answer == 'EPS') then
    allocate(Ecent(CellNumb)); Ecent = 0.0
    do i = 1, CellNumb
      call ReadC(5,inp,tn,ts,te)
      read(inp,*) Ecent(i)
    end do
  end if

!-----Pk
  if(answer == 'PK') then
    allocate(Pkcent(CellNumb)); Pkcent = 0.0
    do i = 1, CellNumb
      call ReadC(5,inp,tn,ts,te)
      read(inp,*) Pkcent(i)
    end do
  end if

!-----v_2
  if(answer == 'VV') then
    allocate(VVcent(CellNumb)); VVcent = 0.0
    do i = 1, CellNumb
      call ReadC(5,inp,tn,ts,te)
      read(inp,*) VVcent(i)
    end do
  end if

!-----zeta
  if(answer == 'ZETA') then
    allocate(ZETAcent(CellNumb)); ZETAcent = 0.0
    do i = 1, CellNumb
      call ReadC(5,inp,tn,ts,te)
      read(inp,*) ZETAcent(i)
    end do
  end if

!-----f22
  if(answer == 'F22') then
    allocate(F22cent(CellNumb)); F22cent = 0.0
    do i = 1, CellNumb
      call ReadC(5,inp,tn,ts,te)
      read(inp,*) F22cent(i)
    end do
  end if

!-----Vis
  if(answer == 'VIS') then
    allocate(VIScent(CellNumb)); VIScent = 0.0
    do i = 1, CellNumb
      call ReadC(5,inp,tn,ts,te)
      read(inp,*) VIScent(i)
    end do
  end if

  GOTO 2   !go to close variables loop
4 CONTINUE !variables loop exit

  close(5)


!-----Initialization
  if(this < 2) write(*,*) '# Initialization of flow field from file:'

  do c=1,NC
    Mres = HUGE
    if(material(c)==1) then
    do i = 1, CellNumb
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
      if(c==2) write(*,*) '@ReadGMV: P not found in ', GMVfile,  &
                          ' Initialization with constant (0.0)!'
      P % n(c) = 0.0
    end if

    if(HOT==YES) then
      if(ALLOCATED(Tcent)) then
        T % n(c) = Tcent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: T not found in ', GMVfile, &
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
        if(c==2) write(*,*) '@ReadGMV: Kin not found in ', GMVfile, &
                            ' Initialization with constant (1.0E-3)!'
        Kin % n(c) = 1.0E-3
      end if
      Kin % o(c)  = Kin % n(c)
      Kin % oo(c) = Kin % n(c)

      if(ALLOCATED(Ecent)) then
        Eps % n(c) = Ecent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: Eps not found in ', GMVfile, &
                            ' Initialization with constant (1.0E-4)!'
        Eps % n(c)  = 1.0E-4
      end if
      Eps % o(c)  = Eps % n(c)
      Eps % oo(c) = Eps % n(c)

      if(ALLOCATED(Pkcent)) then
        Pk(c) = Pkcent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: Pk not found in ', GMVfile,  &
                            ' Initialization with constant (1.0E-4)!'
        Pk(c)  = 1.0E-4
      end if

      VISt(c) = Cmu * Kin%n(c)*Kin%n(c) / Eps%n(c)
    end if

    if(SIMULA==K_EPS_VV.or.SIMULA==ZETA.or.SIMULA==HYB_ZETA) then
      if(ALLOCATED(Kcent)) then
        Kin % n(c) = Kcent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: Kin not found in ', GMVfile, &
                            ' Initialization with constant (1.0E-3)!'
        Kin % n(c) = 1.0E-3
      end if
      Kin % o(c)  = Kin % n(c)
      Kin % oo(c) = Kin % n(c)

      if(ALLOCATED(Ecent)) then
        Eps % n(c) = Ecent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: Eps not found in ', GMVfile, &
                            ' Initialization with constant (1.0E-4)!'
        Eps % n(c) = 1.0E-4
      end if
      Eps % o(c)  = Eps % n(c)
      Eps % oo(c) = Eps % n(c)

      if(ALLOCATED(Pkcent)) then
        Pk(c) = Pkcent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: Pk not found in ', GMVfile,  &
                            ' Initialization with constant (1.0E-3)!'
        Pk(c)  = 1.0E-3
      end if

      if(ALLOCATED(VVcent)) then
        if(ZETA==YES) then                          !Zeta!--inverted reading
          v_2 % n(c) = VVcent(j) / (Kcent(j)+TINY)  !Zeta!--inverted reading
        else
          v_2 % n(c) = VVcent(j)
        end if
      else if(ALLOCATED(ZETAcent)) then
        if(ZETA==YES) then
          v_2 % n(c) = ZETAcent(j)
        else                                        !Zeta!--inverted reading
          v_2 % n(c) = ZETAcent(j) * Kcent(j)       !Zeta!--inverted reading
        end if
      else
        if(c==2) write(*,*) '@ReadGMV: zeta not found in ', GMVfile,  &
                            ' Initialization with constant (1.0E-3)!'
        v_2 % n(c) = 1.0E-3  !1.0E-2
      end if
      v_2 % o(c)  = v_2 % n(c)
      v_2 % oo(c) = v_2 % n(c)

      if(ALLOCATED(F22cent)) then
        f22 % n(c) = F22cent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: f22 not found in ', GMVfile, &
                            ' Initialization with constant (1.0E-3)!'
        f22 % n(c) = 1.0E-4  !1.0E-3
      end if
      f22 % o(c)  = f22 % n(c)
      f22 % oo(c) = f22 % n(c)

      Uf(c)  = 0.047
      Ynd(c) = 30.0
      T1 = Kin%n(c) / (Eps%n(c) + TINY)                                          !ZETA!--New scales
      T2 = 1.0 * Ct * ( VISc / (Eps%n(c) + TINY) )**0.5
      L1 = Cl * Kin%n(c)**1.5 / (Eps%n(c)+TINY)                                  !ZETA!--New scales
      L2 = Cl * Cni * ( VISc**3.0 / (Eps%n(c)+TINY) )**0.25
      Tsc(c) = max(T1, T2)
      Lsc(c) = max(L1, L2)
!      if(ZETA==YES) then
!        VISt(c) = CmuD * v_2%n(c) * Kin%n(c) * Tsc(c)
!      else
!        VISt(c) = CmuD * v_2%n(c) * Tsc(c)
!      end if
       VISt(c) = Cmu * Kin%n(c)*Kin%n(c) / (Eps%n(c) + TINY)
    end if

    if(SIMULA == SPA_ALL .or. SIMULA==DES_SPA) then
      if(ALLOCATED(VIScent)) then
        VIS % n(c)  = VIScent(j)
        VIS % o(c)  = VIScent(j)
        VIS % oo(c) = VIScent(j)
      else
        if(c==2) write(*,*) '@ReadGMV: viscosity not found in ', GMVfile
        if(c==2) write(*,*) '#         Viscosity initialized with constant!!'
        VIS % n(c)  = 1.0
        VIS % o(c)  = 1.0
        VIS % oo(c) = 1.0
      end if
    end if

!-----Status printout
    if(mod(c,5000) == 0) write(*,*) 100.0 * c / (1.0 * NC), '% complete...'
    else 
      U % mean(c) = 0.0
      V % mean(c) = 0.0
      W % mean(c) = 0.0
      U % n(c)    = 0.0     
      U % o(c)    = 0.0     
      U % oo(c)   = 0.0     
      V % n(c)    = 0.0      
      V % o(c)    = 0.0     
      V % oo(c)   = 0.0      
      W % n(c)    = 0.0     
      W % o(c)    = 0.0      
      W % oo(c)   = 0.0     

      P % n(c) = 0.0

      T % o(c)  = 20.0
      T % oo(c) = 20.0
        
    end if
  end do   !end do=1,NC

!-----Fluxes, convective and diffusive terms
  do s=1,NS
  c1=SideC(1,s)
  c2=SideC(2,s)
    Flux(s) = f(s) * DENc(material(c1)) * (U%n(c1)*Sx(s)+V%n(c1)*Sy(s)+W%n(c1)*Sz(s))     &
            + (1.0-f(s)) * DENc(material(c2)) * (U%n(c2)*Sx(s)+V%n(c2)*Sy(s)+W%n(c2)*Sz(s))
  end do
!  call Exchng(Flux)    !moze li se razmjenjivati flux?!?!
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

6 CONTINUE

  RETURN
  END SUBROUTINE ReadGMV
!======================================================================!
