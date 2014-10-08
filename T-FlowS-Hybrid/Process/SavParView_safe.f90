!======================================================================!
  SUBROUTINE SavParView(sub, NCsub, namAut)
!----------------------------------------------------------------------!
! Reads: NAME.gmv and generates NAME.vti Paraview XML output file      !
! ~~~~~~~                                                              ! 
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE par_mod
  USE allp_mod
  USE les_mod
  USE pro_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER ::  NCsub, sub
  CHARACTER :: storename*80, namTem*80, namXML*80
  CHARACTER :: namAut*(*)
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: c,  c1,  c2,  n, s, contauai
  CHARACTER :: namOut*80, Line*300, stringadummy*100, nameIn*80
  REAL,ALLOCATABLE :: x(:), y(:), z(:)    ! self evident
  INTEGER,ALLOCATABLE :: connessione(:,:) ! connection
  INTEGER :: celleconnessione
  INTEGER :: NNsub, NmaterBC, NNsub_new, NCsub_new
  INTEGER :: off_set_connection
  INTEGER :: i
  INTEGER :: numprocessi
  REAL,ALLOCATABLE :: vettore (:)         !local vector for postprox
  REAL,ALLOCATABLE :: vettore2 (:)        !local vector for postprox
  REAL,ALLOCATABLE :: vettore3 (:)        !local vector for postprox
  REAL    :: Nx, Ny, Nz
  REAL    :: Cs, R
  REAL    :: Stot, lf, UtauL, Uff
  REAL    :: Utot, Unor, Utan, Apow, Bpow, nu, dely, yPlus
  REAL    :: frictionv(NC)


!======================================================================!

!  write(*,*) 'writing paraview XML data file: ', namAut

!  if (this < 2) then
!      write(*,*) 'writing paraview XML data file: ', namAut
!  end if


  namTem = name
  storename = namAut
!<<<<<<<<<<<<<<<<<<<<<<<<<!
!                         !
!     reads GMV file      !
!                         !
!<<<<<<<<<<<<<<<<<<<<<<<<<!
  call NamFil(sub, namOut, '.gmv', len_trim('.gmv'))
  open(9, FILE=namOut)
  if (this <2) then
    write(*,*) 'Now reading the file: ', namOut
  end if

!---------------!
!     start     !
!---------------!
  call ReadC(9,inp,tn,ts,te)  !read 'gmvinput ascii' line

!---------------!
!     nodes     !
!---------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp(ts(1):te(1)),*) stringadummy
  read(inp(ts(2):te(2)),*) NNsub
  allocate(x(NNsub)); x=0.
  allocate(y(NNsub)); y=0.
  allocate(z(NNsub)); z=0.
  
  do n=1,NNsub
    call ReadC(9,inp,tn,ts,te)                           
    read(inp(ts(1):te(1)),*) x(n)
  end do
  do n=1,NNsub
    call ReadC(9,inp,tn,ts,te)                           
    read(inp(ts(1):te(1)),*) y(n)
  end do
  do n=1,NNsub
    call ReadC(9,inp,tn,ts,te)                           
    read(inp(ts(1):te(1)),*) z(n)
  end do

!----------------------!
!     cell section     !
!----------------------!

  celleconnessione = 0
  off_set_connection = 0

  call ReadC(9,inp,tn,ts,te)                           
  read(inp(ts(1):te(1)),*) stringadummy
  read(inp(ts(2):te(2)),*) NCsub_new

  if (NCsub_new.ne.NCsub) then 
     write(*,*) 'number of cells read and processed is different, exiting!'
     stop
  end if
  
  allocate(connessione(NCsub,9)); connessione=0

  do n=1,NCsub
    call ReadC(9,inp,tn,ts,te)                           
    read(inp(ts(1):te(1)),*) stringadummy
    read(inp(ts(2):te(2)),*) off_set_connection
    if (n==1) then 
       connessione(n,1) = off_set_connection
    else
       connessione(n,1) = connessione (n-1,1)+ off_set_connection
    end if
    
    celleconnessione = celleconnessione + off_set_connection

    call ReadC(9,inp,tn,ts,te)
    do c=1,off_set_connection
      read(inp(ts(c):te(c)),*) connessione(n,c+1)
    end do
  end do

  close(9)

  name = namAut
  call NamFil(sub, namXML, '.vtu', len_trim('.vtu'))

  open(9, FILE=namXML)
  if (this <2) then
  write(6, *) 'Now writing the file:', namXML
  end if


  write(9,'(A21)') '<?xml version="1.0"?>'
  write(9,'(A73)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(9,'(A20)') '  <UnstructuredGrid>'

  write(9,*) '    <Piece NumberOfPoints="', NNsub,'" NumberOfCells="', NCsub, '" >' 


  write(9,*) '       <CellData Scalars="scalars" vectors="velocity">'

!------------------------------------------------------------------------------------------------------------
! S C A L A R S  S C A L A R S  S C A L A R S   S C A L A R S   S C A L A R S   S C A L A R S   S C A L A R S
!
! NOTICE: please remember to put thermal field - related variables in their own section
!------------------------------------------------------------------------------------------------------------

  !--scalar: pressure
  write(9,*) '        <DataArray type="Float32" Name="pressure" format="ascii">'
                             do c=1,NCsub
                                write(9,*) P % n(c)
                             end do  
  write(9,*) '        </DataArray>'
  
!--scalar: k and eps     
if ((SIMULA == K_EPS).or.(SIMULA == K_EPS_VV).or.(SIMULA == ZETA).or.(SIMULA==HYB_ZETA)) then

  write(9,*) '        <DataArray type="Float32" Name="k" format="ascii">'
                             do c=1,NCsub
                                write(9,*) Kin % n(c)
                             end do  
  write(9,*) '        </DataArray>'
  
  write(9,*) '        <DataArray type="Float32" Name="eps" format="ascii">'
                             do c=1,NCsub
                                write(9,*) Eps % n(c)
                             end do  
  write(9,'(A20)') '        </DataArray>'

end if

!--scalar: zeta or v2 and eventually f   
if ((SIMULA == K_EPS_VV).or.(SIMULA == ZETA).or.(SIMULA==HYB_ZETA)) then
  
  if (SIMULA == K_EPS_VV) then
    write(9,*) '        <DataArray type="Float32" Name="v2" format="ascii">'
  else 
    write(9,*) '        <DataArray type="Float32" Name="zeta" format="ascii">'
  end if
                             do c=1,NCsub
                                write(9,*) v_2 % n(c)
                             end do  
  write(9,'(A20)') '        </DataArray>'
  
!  if (ZETA_EFF == YES ) then
!       write(9,*) '        <DataArray type="Float32" Name="f" format="ascii">'
!                             do c=1,NCsub
!                                write(9,*) f22 % n(c)
!                             end do
!       write(9,'(A20)') '        </DataArray>'
!  end if !ZETA_EFF

end if 

!--scalar: alpha field
if (SIMULA == HYB_ZETA) then
  write(9,*) '        <DataArray type="Float32" Name="alpha" format="ascii">'
                             do c=1,NCsub
                                write(9,*) VISt(c)
                             end do
  write(9,*) '        </DataArray>'
end if


!--scalar: viscosities for HYB_ZETA
if (SIMULA == HYB_ZETA) then
  write(9,'(A65)') '        <DataArray type="Float32" Name="visT" format="ascii">'
                             do c=1,NCsub
                                write(9,*) VISt(c)/VISC
                             end do
  write(9,'(A20)') '        </DataArray>'
!  write(9,'(A65)') '        <DataArray type="Float32" Name="visSGS" format="ascii">'
!                             do c=1,NCsub
!                                write(9,*) VISt_sgs(c)/VISC
!                             end do
!  write(9,'(A20)') '        </DataArray>'
!  write(9,'(A65)') '        <DataArray type="Float32" Name="visEFF" format="ascii">'
!                             do c=1,NCsub
!                                write(9,*) VISt_eff(c)/VISC
!                             end do
!  write(9,'(A20)') '        </DataArray>'
!  write(9,'(A65)') '<DataArray type="Float32" Name="visT_mean" format="ascii">'
!                             do c=1,NCsub
!                                write(9,*) VISt_mean(c)/VISC
!                             end do
!  write(9,'(A20)') '        </DataArray>'
!  write(9,'(A75)') '<DataArray type="Float32" Name="visT_SGS_mean" format="ascii">'
!                             do c=1,NCsub
!                                write(9,*) VISt_sgs_mean(c)/VISC
!                             end do
!  write(9,'(A20)') '        </DataArray>'
end if  
if (SIMULA==HYB_ZETA) then
  write(9,'(A65)') '        <DataArray type="Float32" Name="K_mean" format="ascii">'
                             do c=1,NCsub
                                write(9,*) Kin % mean(c)
                             end do
  write(9,'(A20)') '        </DataArray>'
  write(9,'(A65)') '        <DataArray type="Float32" Name="Eps_mean" format="ascii">'
                             do c=1,NCsub
                                write(9,*) Eps % mean(c)
                             end do
  write(9,'(A20)') '        </DataArray>'
  write(9,'(A65)') '<DataArray type="Float32" Name="Zeta_mean" format="ascii">'
                             do c=1,NCsub
                                write(9,*) v_2 % mean(c)
                             end do
  write(9,'(A20)') '        </DataArray>'
end if




!--scalar: wall distance  
!if(MATTONE==YES) then
  write(9,'(A65)') '<DataArray type="Float32" Name="wall distance" format="ascii">'
                             do c=1,NCsub
                                write(9,*) WallDs(c)
                             end do  
  write(9,'(A20)') '        </DataArray>'
!end if

!--scalar: mean pressure
if(SIMULA == LES) then
  write(9,'(A65)') '        <DataArray type="Float32" Name="Pmean" format="ascii">'
                             do c=1,NCsub
                                write(9,*) P % mean(c)
                             end do  
  write(9,'(A20)') '        </DataArray>'
end if

!--scalar: HYB_ZETA kLES
if (SIMULA == HYB_ZETA) then
  write(9,'(A65)') '        <DataArray type="Float32" Name="kLES" format="ascii">'
                             do c=1,NCsub
                                write(9,*) 0.5 *((uu % mean(c) - U % mean(c) * U % mean(c)) + &
                                           (vv % mean(c) - V % mean(c) * V % mean(c)) +                  &
                                           (ww % mean(c) - W % mean(c) * W % mean(c)))
                             end do
  write(9,'(A20)') '        </DataArray>'
end if
!if (SIMULA == HYB_ZETA) then
!  write(9,'(A65)') '<DataArray type="Float32" Name="LRANS" format="ascii">'
!                             do c=1,NCsub
!                                write(9,*) Lrans(c)
!                             end do
!  write(9,'(A20)') '        </DataArray>'
!  write(9,'(A65)') '        <DataArray type="Float32" Name="LLES" format="ascii">'
!                             do c=1,NCsub
!                                write(9,*) LLES(c)
!                             end do
!  write(9,'(A20)') '        </DataArray>'
!end if

!--scalar: reynolds stresses  
!                             per domenico: qualora mi volessi far notare che e' un tensore, pensa a quanto puo'
!                             essere mortalmente noioso scrivere sta routine e fatti due conti se ti conviene
!                             scassare i cabasisi :oP

if ( (SIMULA == YES).and.(DNS == YES) ) then
  write(9,'(A65)') '        <DataArray type="Float32" Name="uu" format="ascii">'
                             do c=1,NCsub
                                write(9,*) uu % mean(c) - U % mean(c) * U % mean(c)
                             end do  
  write(9,'(A20)') '        </DataArray>'
  write(9,'(A65)') '        <DataArray type="Float32" Name="vv" format="ascii">'
                             do c=1,NCsub
                                write(9,*) vv % mean(c) - V % mean(c) * V % mean(c)
                             end do  
  write(9,'(A20)') '        </DataArray>'
  write(9,'(A65)') '        <DataArray type="Float32" Name="ww" format="ascii">'
                             do c=1,NCsub
                                write(9,*) ww % mean(c) - W % mean(c) * W % mean(c)
                             end do  
  write(9,'(A20)') '        </DataArray>'
  write(9,'(A65)') '        <DataArray type="Float32" Name="uv" format="ascii">'
                             do c=1,NCsub
                                write(9,*) uv % mean(c) - U % mean(c) * V % mean(c)
                             end do  
  write(9,'(A20)') '        </DataArray>'
  write(9,'(A65)') '        <DataArray type="Float32" Name="uw" format="ascii">'
                             do c=1,NCsub
                                write(9,*) uw % mean(c) - U % mean(c) * W % mean(c)
                             end do  
  write(9,'(A20)') '        </DataArray>'
  write(9,'(A65)') '        <DataArray type="Float32" Name="vw" format="ascii">'
                             do c=1,NCsub
                                write(9,*) vw % mean(c) - V % mean(c) * W % mean(c)
                             end do  
  write(9,'(A20)') '        </DataArray>'
end if !ReStresses

!--scalar: pressure laplacian 
!if (PLOT_LAPP == YES) then
!  allocate (vettore(-NbC:NC)); vettore = 0.0
!  call CalcLapP(vettore)
!  write(9,'(A65)') '        <DataArray type="Float32" Name="LapP" format="ascii">'
!                             do c=1,NCsub
!                                write(9,*) vettore(c)
!                             end do  
!  write(9,'(A20)') '        </DataArray>'
!  deallocate (vettore)
!end if
!--scalar: Q 
!if (PLOT_Q == YES) then
!  allocate (vettore(-NbC:NC)); vettore = 0.0
!  call CalcQ(vettore)
!  write(9,'(A65)') '        <DataArray type="Float32" Name="Q" format="ascii">'
!                             do c=1,NCsub
!                                write(9,*) vettore(c)
!                             end do  
!  write(9,'(A20)') '        </DataArray>'
!  deallocate (vettore)
!end if
!--scalar: Cdyn - notice: there should be an if to check the model, but don't want to bother with this
!if(PLOT_CDYN==YES) then
!  write(9,'(A65)') '        <DataArray type="Float32" Name="Cdyn" format="ascii">'
!                             do c=1,NCsub
!                                write(9,*) Cdyn(c)
!                             end do  
!  write(9,'(A20)') '        </DataArray>'
!end if
!if(PLOT_CELLV==YES) then
!  write(9,'(A65)') '<DataArray type="Float32" Name="cell volume" format="ascii">'
!                             do c=1,NCsub
!                                write(9,*) volume(c)
!                             end do  
!  write(9,'(A20)') '        </DataArray>'
!end if

if(SIMULA==HYB_ZETA) then
      allocate (vettore(-NbC:NC)); vettore = 0.0
      call GraPhi(U % n, 1, Ux,.TRUE.)    ! dU/dx
      call GraPhi(U % n, 2, Uy,.TRUE.)    ! dU/dy
      call GraPhi(U % n, 3, Uz,.TRUE.)    ! dU/dz
      call GraPhi(V % n, 1, Vx,.TRUE.)    ! dV/dx
      call GraPhi(V % n, 2, Vy,.TRUE.)    ! dV/dy
      call GraPhi(V % n, 3, Vz,.TRUE.)    ! dV/dz
      call GraPhi(W % n, 1, Wx,.TRUE.)    ! dW/dx
      call GraPhi(W % n, 2, Wy,.TRUE.)    ! dW/dy
      call GraPhi(W % n, 3, Wz,.TRUE.)    ! dW/dz

!  if(SIMULA == HYB_ZETA) then
    do c=1,NC
      vettore(c)=0.
      vettore(c)=vettore(c)+(Ux(c))**2
      vettore(c)=vettore(c)+(Vy(c))**2
      vettore(c)=vettore(c)+(Wz(c))**2
      vettore(c)=vettore(c)+(0.25*(Uy(c)+Vx(c)))
      vettore(c)=vettore(c)+(0.25*(Uy(c)+Vx(c)))
      vettore(c)=vettore(c)+(0.25*(Uz(c)+Wx(c)))
      vettore(c)=vettore(c)+(0.25*(Uz(c)+Wx(c)))
      vettore(c)=vettore(c)+(0.25*(Vz(c)+Wy(c)))
      vettore(c)=vettore(c)+(0.25*(Vz(c)+Wy(c)))
      vettore(c)=max(0.0,vettore(c))
      vettore(c)=vettore(c)**.5
      vettore(c)=vettore(c)**3  
      vettore(c)=((eps %n(c)*eps%n(c))**.5)/(vettore(c)+TINY)
      vettore(c)=max(0.0,vettore(c))
      vettore(c)=(vettore(c))**.5
    end do
    write(9,'(A65)') '<DataArray type="Float32" Name="LRANS-P" format="ascii">'
                             do c=1,NCsub
                                write(9,*) vettore(c)
                             end do
    write(9,'(A20)') '        </DataArray>'
    deallocate (vettore)
end if

call wait

1290 continue
!------------------------------------------------------------------------------------------------------------
! V E C T O R S  V E C T O R S  V E C T O R S   V E C T O R S   V E C T O R S   V E C T O R S   V E C T O R S
!------------------------------------------------------------------------------------------------------------

!--vector: velocity
  write(9,'(A160)') '<DataArray type="Float32" Name="velocity" NumberOfComponents="3" format="ascii">'
                             do c=1,NCsub
                                write(9,*) U % n(c), V % n(c), W % n(c) 
                             end do  
  write(9,'(A20)') '        </DataArray>'

if (SIMULA == LES) then
  write(9,'(A99)') '<DataArray type="Float32" Name="mean velocity" NumberOfComponents="3" format="ascii">'
                             do c=1,NCsub
                                write(9,*) U % mean(c), V % mean(c), W % mean(c) 
                             end do  
  write(9,'(A20)') '        </DataArray>'
end if   

!------------------------------------------------------------------------------------------------------------
! T H E R M A L  F I E L D  R E L A T E D   V A R I A B L E S (HOT == YES)
!------------------------------------------------------------------------------------------------------------

if(HOT == YES) then

! scalar: temperature
  write(9,'(A65)') '        <DataArray type="Float32" Name="T" format="ascii">'
                             do c=1,NCsub
                                write(9,*) T % n(c)
                             end do  
  write(9,'(A20)') '        </DataArray>'

! scalar: mean temperature



  if (SIMULA == LES) then
    write(9,*) '        <DataArray type="Float32" Name="Tmean" format="ascii">'
                             do c=1,NCsub
                                write(9,*) T % mean(c)
                             end do  
    write(9,'(A20)') '        </DataArray>'
  end if

  if (SIMULA==LES) then
    write(9,'(A65)') '        <DataArray type="Float32" Name="TT" format="ascii">'
                             do c=1,NCsub
                                write(9,*) TT % mean(c) - T % mean(c) * T % mean(c)
                             end do  
    write(9,'(A20)') '        </DataArray>'
    write(9,'(A65)') '        <DataArray type="Float32" Name="uT" format="ascii">'
                             do c=1,NCsub
                                write(9,*) uT % mean(c) - U % mean(c) * T % mean(c)
                             end do  
    write(9,'(A20)') '        </DataArray>'
    write(9,'(A65)') '        <DataArray type="Float32" Name="vT" format="ascii">'
                             do c=1,NCsub
                                write(9,*) vT % mean(c) - V % mean(c) * T % mean(c)
                             end do  
    write(9,'(A20)') '        </DataArray>'
    write(9,'(A65)') '        <DataArray type="Float32" Name="wT" format="ascii">'
                             do c=1,NCsub
                                write(9,*) wT % mean(c) - W % mean(c) * T % mean(c)
                             end do  
    write(9,'(A20)') '        </DataArray>'
  end if


end if !HOT == YES




  write(9,'(A17)') '      </CellData>'
  
  write(9,'(A14)') '      <Points>'
  write(9,*) '         <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
                             do c=1,NNsub
                                 write(9,*) x(c), y(c), z(c)
                             end do
  write(9,'(A20)') '        </DataArray>'
  write(9,'(A15)') '      </Points>'
  write(9,'(A13)') '      <Cells>'
  write(9,'(A71)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
  do n=1,NCsub
                                   !hexa
                                   if (connessione(n,9) /= 0) then
                                     do i=2,9
                                            write(9,*) connessione(n,i)-1
                                     end do
                                   else
                                       !prism
                                       if ((connessione(n,8) == 0).and.(connessione(n,7) /= 0)) then
                                         write(9,*) connessione(n,2)-1, connessione(n,4)-1, & 
                                         connessione(n,3)-1, connessione(n,5)-1, connessione(n,7)-1, connessione(n,6)-1
                                       else
                                          if ((connessione(n,7) == 0).and.(connessione(n,6) /= 0)) then
                                              write(9,*) connessione(n,5)-1, connessione(n,4)-1,&
                                              connessione(n,3)-1, connessione(n,6)-1, connessione(n,2)-1
                                          else
                                              write(9,*) connessione(n,5)-1, connessione(n,4)-1,&
                                              connessione(n,3)-1, connessione(n,2)-1
                                          end if
                                       end if
                                   end if

                              end do 
 write(9,'(A20)') '        </DataArray>'
  write(9,'(A62)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
                              do n=1,NCsub
                                 write(9,*) connessione(n,1)
                              end do
  write(9,'(A20)') '        </DataArray>'
  write(9,'(A60)') '        <DataArray type="UInt8" Name="types" format="ascii">'
                              do n=1,NCsub
                                 if ((connessione(n,6)==0)) then !thetra
                                      write(9,*) '10'
                                 else 
                                 if ((connessione(n,7)==0)) then
                                      write(9,*) '14'          !pyr
                                        else
                                        if ((connessione(n,8)==0)) then
                                            write(9,*) '13'      !prism
                                        else
                                            write(9,*) '12'      !hexa
                                        end if
                                     end if
                                 end if
                              end do
  write(9,'(A20)') '        </DataArray>'
  write(9,'(A14)') '      </Cells>'
  write(9,'(A14)') '</Piece>'
  write(9,'(A19)') '</UnstructuredGrid>'
  write(9,'(A14)') '</VTKFile>'

! write down the master file for parallel jobs
! actually it writes it down also for sequential process but it's not
! necessary for paraview
  numprocessi = this
  call GloMax(numprocessi)
  if (numprocessi==0) then
     numprocessi=1
  end if

  if (this < 2) then
    name = storename
    call NamFil(0, name, '.pvtu', len_trim('.pvtu'))
    open(12, FILE=name)
!     write(*,*) 'writing parallel master file: ', name
     write(12,'(A21)') '<?xml version="1.0"?>'
     write(12,'(A74)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
     write(12,*) '<PUnstructuredGrid GhostLevel="0">'
     write(12,*) '       <PCellData Scalars="scalars" vectors="velocity">'
     write(12,*) '        <PDataArray type="Float32" Name="pressure"/>'
     !--scalar: k and eps
     if ((SIMULA == K_EPS).or.(SIMULA == K_EPS_VV).or.(SIMULA == ZETA).or.(SIMULA==HYB_ZETA)) then
       write(12,*) '        <PDataArray type="Float32" Name="k"/>'
       write(12,*) '        <PDataArray type="Float32" Name="eps"/>'
     end if
     !..scalar: zeta (or v2) and eventually f
     if ((SIMULA == K_EPS_VV).or.(SIMULA == ZETA).or.(SIMULA==HYB_ZETA)) then
        if (SIMULA == K_EPS_VV) then
           write(12,*) '        <PDataArray type="Float32" Name="v2"/>' 
        else
           write(12,*) '        <PDataArray type="Float32" Name="zeta"/>'
        end if
!       if (ZETA_EFF == YES ) then
!         write(12,*) '        <PDataArray type="Float32" Name="f"/>'
!       end if
     end if
     !--scalar: alpha field
     if (SIMULA == HYB_ZETA) then
        write(12,*)      '        <PDataArray type="Float32" Name="alpha"/>'
     end if
     !--scalar: viscosities for HYB_ZETA
     if (SIMULA == HYB_ZETA) then
        write(12,'(A65)') '        <PDataArray type="Float32" Name="visT"/>'
        write(12,'(A65)') '        <PDataArray type="Float32" Name="visSGS"/>'
        write(12,'(A65)') '        <PDataArray type="Float32" Name="visEFF"/>'
        write(12,'(A65)') '        <PDataArray type="Float32" Name="visT_mean"/>'
        write(12,'(A65)') '        <PDataArray type="Float32" Name="visT_SGS_mean"/>'
     end if
     !-- scalar: kin, eps and zeta mean for HYB_ZETA simulation
     if (SIMULA==HYB_ZETA) then
         write(12,'(A65)') '        <PDataArray type="Float32" Name="K_mean"/>'
         write(12,'(A65)') '        <PDataArray type="Float32" Name="Eps_mean"/>'
         write(12,'(A65)') '        <PDataArray type="Float32" Name="Zeta_mean"/>'
     end if
     !--scalar: wall distance
!     if(MATTONE==YES) then
       write(12,'(A65)') '        <PDataArray type="Float32" Name="wall distance"/>'
!     end if
     !--scalar: mean pressure
     if(SIMULA==LES) then
       write(12,'(A65)') '        <PDataArray type="Float32" Name="Pmean"/>'
     end if
     !--scalar: HYB_ZETA kLES
     if (SIMULA == HYB_ZETA) then
       write(12,'(A65)') '        <PDataArray type="Float32" Name="kLES"/>'
     end if
     if (SIMULA == HYB_ZETA) then
         write(12,'(A65)') '<PDataArray type="Float32" Name="LRANS"/>'
         write(12,'(A65)') '<PDataArray type="Float32" Name="LLES"/>'
     end if
     !--scalar: Reynolds Stresses
     if ( SIMULA == LES ) then
       write(12,'(A65)') '        <PDataArray type="Float32" Name="uu"/>'
       write(12,'(A65)') '        <PDataArray type="Float32" Name="vv"/>'
       write(12,'(A65)') '        <PDataArray type="Float32" Name="ww"/>'
       write(12,'(A65)') '        <PDataArray type="Float32" Name="uv"/>'
       write(12,'(A65)') '        <PDataArray type="Float32" Name="uw"/>'
       write(12,'(A65)') '        <PDataArray type="Float32" Name="vw"/>'
     end if !ReStresses 
     !--scalar: pressure laplacian
     if(SIMULA==HYB_ZETA) then
       write(12,'(A65)') '<PDataArray type="Float32" Name="LRANS-P"/>'
     end if
     !------------------------------------------------------------------------------------------------------------
     ! NOTICE: IF YOU WANT TO ADD ANOTHER SCALAR DO IT BEFORE THIS COMMENT OR AFTER VORTICITY MAGNITUDE AS
     !         VORTICITY PLOT HAS A TRICKY TREATMENT YOU DON'T WANT TO MESS WITH!
     !------------------------------------------------------------------------------------------------------------

     !if user has decided to plot all the vorticity components go to vectors
!     if ((PLOT_VORT_X==YES).and.(PLOT_VORT_Y==YES).and.(PLOT_VORT_Z==YES)) goto 2290
     !--scalar: x-vorticity
!     if (PLOT_VORT_X == YES) then
!       write(12,'(A65)') '        <PDataArray type="Float32" Name="x-vorticity"/>'
!     end if
     !--scalar: y-vorticity
!     if (PLOT_VORT_Y == YES) then
!       write(12,'(A65)') '        <PDataArray type="Float32" Name="y-vorticity"/>'
!     end if
     !--scalar: z-vorticity
!     if (PLOT_VORT_Z == YES) then
!       write(12,'(A65)') '        <PDataArray type="Float32" Name="z-vorticity"/>'
!     end if
     !--scalar: vorticity magnitude
2290 continue 
!     if (PLOT_VORT_MAGN == YES) then
!       write(12,*) '        <PDataArray type="Float32" Name="vorticity magnitude"/>' 
!     end if

!------------------------------------------------------------------------------------------------------------
! V E C T O R S  V E C T O R S  V E C T O R S   V E C T O R S   V E C T O R S   V E C T O R S   V E C T O R S
!------------------------------------------------------------------------------------------------------------
     write(12,*) '        <PDataArray type="Float32" Name="velocity" NumberOfComponents="3"/>'

     if (SIMULA==LES) then
       write(12,'(A80)') '<PDataArray type="Float32" Name="mean velocity" NumberOfComponents="3"/>'
     end if

!     if ((PLOT_VORT_X==YES).and.(PLOT_VORT_Y==YES).and.(PLOT_VORT_Z==YES))then
!       write(12,*) '        <PDataArray type="Float32" Name="vorticity" NumberofComponents="3"/>'
!     end if
!------------------------------------------------------------------------------------------------------------
! T H E R M A L  F I E L D  R E L A T E D   V A R I A B L E S (HOT == YES)
!------------------------------------------------------------------------------------------------------------

     if(HOT == YES) then
       ! scalar: temperature
         write(12,'(A65)') '        <PDataArray type="Float32" Name="T"/>'
       ! scalar: mean temperature
       if (SIMULA==LES) then
         write(12,*) '        <PDataArray type="Float32" Name="Tmean"/>'
       end if
       if (SIMULA==LES) then
         write(12,'(A65)') '        <PDataArray type="Float32" Name="TT"/>'
         write(12,'(A65)') '        <PDataArray type="Float32" Name="uT"/>'
         write(12,'(A65)') '        <PDataArray type="Float32" Name="vT"/>'
         write(12,'(A65)') '        <PDataArray type="Float32" Name="wT"/>'
       end if
     end if !HOT == YES
     write(12,*) '       </PCellData>'
     write(12,*) '       <PPoints>'
     write(12,*) '         <PDataArray type="Float32" NumberOfComponents="3"/>'
     write(12,*) '       </PPoints>'
     ! qui ci devo mettere la stampa dei vari files di sottoprocesso  <Piece Source="chan1_180.vtu"/>
        name = namAut
        do i=1,numprocessi
           call NamFil(i, nameIn, '.vtu"/>', len_trim('.vtu"/>'))
!           write(12,*) '<Piece Source="',nameIn, '"/>'
           write(12,'(A)',advance="no") '<Piece Source="'
           write(12,'(A)') nameIn
        end do
     
     write(12,*) '</PUnstructuredGrid>'
     write(12,*) '</VTKFile>'
  end if !(this <2)

  name = namTem

  call wait

  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(connessione)

  END SUBROUTINE SavParView
