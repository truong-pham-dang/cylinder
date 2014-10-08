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

!=============================================================================================================================
!! HEAD of .VTU
!=============================================================================================================================
  write(9,'(A21)') '<?xml version="1.0"?>'
  write(9,'(A73)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(9,'(A20)') '  <UnstructuredGrid>'
  write(9,      *) '    <Piece NumberOfPoints="', NNsub,'" NumberOfCells="', NCsub, '">' 
  write(9,'(A16)') '      <CellData>'
!=============================================================================================================================
!! Some solely values 
!=============================================================================================================================

!--scalar: alpha field
  write(9,'(A65)') '        <DataArray type="Float32" Name="Kin_mean" format="ascii">'
                             do c=1,NCsub
                                write(9,*) Kin % mean(c)
                             end do
  write(9,'(A20)') '        </DataArray>'

  write(9,'(A65)') '        <DataArray type="Float32" Name="Eps_mean" format="ascii">'
                             do c=1,NCsub
                                write(9,*) Eps % mean(c)
                             end do
  write(9,'(A20)') '        </DataArray>'

  !write(9,'(A62)') '        <DataArray type="Float32" Name="alpha" format="ascii">'
  !                           do c=1,NCsub
  !                              write(9,*) max(1.0, (Kin % n(c)**1.5) / (Eps % n(c) * 0.8 *(volume(c)**(1/3)) )  )
  !                           end do
  !write(9,'(A20)') '        </DataArray>'

  write(9,'(A62)') '        <DataArray type="Float32" Name="LRANS" format="ascii">'
                             do c=1,NCsub
                                write(9,*) (Kin % mean(c)**1.5) / Eps % mean(c)
                             end do
  write(9,'(A20)') '        </DataArray>'

  write(9,'(A61)') '        <DataArray type="Float32" Name="LSGS" format="ascii">'
                             do c=1,NCsub
                                write(9,*) 0.8*volume(c)**(0.33333)
                             end do
  write(9,'(A20)') '        </DataArray>'

  write(9,'(A59)') '        <DataArray type="Float32" Name="Ls" format="ascii">'
                             do c=1,NCsub
                                write(9,*) Ls % mean(c)
                             end do
  write(9,'(A20)') '        </DataArray>'


  !write(9,'(A61)') '        <DataArray type="Float32" Name="VISt" format="ascii">'
  !                           do c=1,NCsub
  !                              write(9,*) VISt(c)
  !                           end do
  !write(9,'(A20)') '        </DataArray>'

  write(9,'(A66)') '        <DataArray type="Float32" Name="VISt_mean" format="ascii">'
                             do c=1,NCsub
                                write(9,*) VISt_mean(c)
                             end do  
  write(9,'(A20)') '        </DataArray>'

  !write(9,'(A65)') '        <DataArray type="Float32" Name="VISt_sgs" format="ascii">'
  !                           do c=1,NCsub
  !                              write(9,*) VISt_sgs(c)
  !                           end do
  !write(9,'(A20)') '        </DataArray>'

  !write(9,'(A70)') '        <DataArray type="Float32" Name="Ceps1RNGratio" format="ascii">'
  !                           do c=1,NCsub
  !                              write(9,'(ES22.15)') Ce1RNGratio(c)
  !                           end do
  !write(9,'(A20)') '        </DataArray>'

  !write(9,'(A65)') '        <DataArray type="Float32" Name="VISt_eff" format="ascii">'
  !                           do c=1,NCsub
  !                              write(9,*) VISt_eff(c)
  !                           end do
  !write(9,'(A20)') '        </DataArray>'
  
  !write(9,'(A70)') '        <DataArray type="Float32" Name="VISt_eff_mean" format="ascii">'
  !                           do c=1,NCsub
  !                              write(9,*) VISt_eff_mean(c)
  !                           end do
  !write(9,'(A20)') '        </DataArray>'



!=============================================================================================================================
!-- mean uu,vv,ww,uv,uw,vw stresses-------------------------------------------------------------------------------------------

  write(9,'(A110)') '        <DataArray type="Float32" Name="Velocity2ndOrderCovarianceMean" NumberOfComponents="6" format="ascii">'
                              do c=1,NCsub
                     write(9,*) (UU % mean(c) - U % mean(c) * U % mean(c) ), &!+ (2/3)* Kin % mean(c) - 2*VISt_mean(c)*( Ux(c) ) ), &
                                (VV % mean(c) - V % mean(c) * V % mean(c) ), &!+ (2/3)* Kin % mean(c) - 2*VISt_mean(c)*( Vx(c) )  ), &
                                (WW % mean(c) - W % mean(c) * W % mean(c) ), &!+ (2/3)* Kin % mean(c) - 2*VISt_mean(c)*( Wx(c) ) ), &
                                (UV % mean(c) - U % mean(c) * V % mean(c) ), &!- VISt_mean(c)*( Uy(c) + Vx(c) ) ), &
                                (UW % mean(c) - U % mean(c) * W % mean(c) ), &!- VISt_mean(c)*( Uz(c) + Wx(c) ) ), &
                                (VW % mean(c) - V % mean(c) * W % mean(c) )   !- VISt_mean(c)*( Vz(c) + Wy(c) ) )
                              end do  

  write(9,'(A20)') '        </DataArray>'
call wait

1290 continue
!------------------------------------------------------------------------------------------------------------
! V E C T O R S  V E C T O R S  V E C T O R S   V E C T O R S   V E C T O R S   V E C T O R S   V E C T O R S
!------------------------------------------------------------------------------------------------------------


!--vector: velocity
  write(9,'(A)')'        <DataArray type="Float32" Name="VelocityInstant" NumberOfComponents="3" format="ascii">'
                             do c=1,NCsub
                                write(9,*) U % n(c), V % n(c), W % n(c) 
                             end do  
  write(9,'(A)') '        </DataArray>'

  write(9,'(A110)') '        <DataArray type="Float32" Name="VelocityMeanAccumulated" NumberOfComponents="3" format="ascii">'
                             do c=1,NCsub
                                write(9,*) U % mean(c), V % mean(c), W % mean(c) 
                             end do  
  write(9,'(A20)') '        </DataArray>'


  write(9,'(A17)') '      </CellData>'
  
  write(9,'(A14)') '      <Points>'
  write(9,'(A72)') '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
                             do c=1,NNsub
                                 write(9,*) x(c), y(c), z(c)
                             end do
  write(9,'(A20)') '        </DataArray>'
  write(9,'(A15)') '      </Points>'
  write(9,'(A13)') '      <Cells>'
  write(9,'(A67)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
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
  write(9,'(A12)') '    </Piece>'
  write(9,'(A21)') '  </UnstructuredGrid>'
  write(9,'(A10)') '</VTKFile>'

! write down the master file for parallel jobs
! actually it writes it down also for sequential process but it's not
! necessary for paraview

  numprocessi = Npro

  call GloMax(numprocessi)
  if (numprocessi==0) then
     numprocessi=1
  end if

  close(9)


  if (this < 2) then
!!_-------------------------------------------------------_!!    
    name = storename
    call NamFil(0, name, '.pvtu', len_trim('.pvtu'))
    open(12, FILE=name)
!     write(*,*) 'writing parallel master file: ', name
     write(12,'(A21)')   '<?xml version="1.0"?>'
     write(12,'(A74)')   '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
     write(12,'(A36)')   '  <PUnstructuredGrid GhostLevel="0">'

     write(12,'(A15)')   '    <PCellData>'

  !--scalar: k and eps
     write(12,'(A52)') '        <PDataArray type="Float32" Name="Kin_mean"/>'
     write(12,'(A52)') '        <PDataArray type="Float32" Name="Eps_mean"/>'
     write(12,'(A49)') '        <PDataArray type="Float32" Name="LRANS"/>'
     write(12,'(A48)') '        <PDataArray type="Float32" Name="LSGS"/>'
     write(12,'(A46)') '        <PDataArray type="Float32" Name="Ls"/>'
     write(12,'(A53)') '        <PDataArray type="Float32" Name="VISt_mean"/>'
     write(12,'(A97)') '        <PDataArray type="Float32" Name="Velocity2ndOrderCovarianceMean" NumberOfComponents="6"/>'
     write(12,'(A82)') '        <PDataArray type="Float32" Name="VelocityInstant" NumberOfComponents="3"/>'
     write(12,'(A90)') '        <PDataArray type="Float32" Name="VelocityMeanAccumulated" NumberOfComponents="3"/>'
     
     
     write(12,'(A16)') '    </PCellData>'
     write(12,'(A13)') '    <PPoints>'
     write(12,'(A57)') '      <PDataArray type="Float32" NumberOfComponents="3"/>' !x,y,z
     write(12,'(A14)') '    </PPoints>'
     ! qui ci devo mettere la stampa dei vari files di sottoprocesso  <Piece Source="chan1_180.vtu"/>
        name = namAut
        do i=1,numprocessi
           call NamFil(i, nameIn, '.vtu"/>', len_trim('.vtu"/>'))
           write(12,'(A)',advance="no") '<Piece Source="'
           write(12,'(A)') nameIn
        end do
     
     write(12,'(A22)') '  </PUnstructuredGrid>'
     write(12,'(A10)') '</VTKFile>'
 !!_-------------------------------------------------------_!!
  end if !(this <2)

  name = namTem

  call wait

  close(12)

  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(connessione)

  END SUBROUTINE SavParView
