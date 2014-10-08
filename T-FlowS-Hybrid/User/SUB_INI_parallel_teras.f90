!======================================================================*
      PROGRAM SUB_INI
!----------------------------------------------------------------------*
!  This program reads files *.xyz, *.U__, *.V__, *.W__, *.P__ and *.T__
!  and creates the similar files for each subdomain. In these files 
!  the results obtained on one mesh are stored. 
!  The LoaIni.f90 reads these files and interpolate the previous solution
!  into new mesh.
!----------------------------------------------------------------------*
      IMPLICIT NONE
!======================================================================*

  REAL,ALLOCATABLE :: xc(:),yc(:),zc(:)
  REAL,ALLOCATABLE :: Sx(:),Sy(:),Sz(:)
  REAL,ALLOCATABLE :: volume(:)            ! cell's volume
  REAL,ALLOCATABLE :: delta(:)             ! delta (max(dx,dy,dz))
  REAL,ALLOCATABLE :: Dx(:),Dy(:),Dz(:)
  REAL,ALLOCATABLE :: xsp(:),ysp(:),zsp(:) ! face coordinates
  REAL,ALLOCATABLE :: WallDs(:), f(:)
  REAL             :: Xmax, Xmin, Ymin, Ymax, Zmin, Zmax

  INTEGER   :: NC, NS, ND, NN                    ! num. of nodes and cells
  INTEGER   :: NbC, Ncopy, NSsh, Nmat

  INTEGER,ALLOCATABLE :: material(:)     ! material markers
  INTEGER,ALLOCATABLE :: SideC(:,:)      !  c0, c1, c2

  INTEGER,ALLOCATABLE :: TypeBC(:)       ! type of boundary condition

  INTEGER,ALLOCATABLE :: CopyC(:)        !  might be shorter
  INTEGER,ALLOCATABLE :: CopyS(:,:)      !  similar to SideC

  INTEGER          :: i, l1, n, IND
  INTEGER          :: j, k,  c, nearest, var, Nvar, c1, c2, s
  INTEGER          :: NCold
  CHARACTER*80 nameIn
  CHARACTER*80 namOut, namSav
  REAL,ALLOCATABLE :: Xold(:),Yold(:),Zold(:)
  REAL,ALLOCATABLE :: Uold(:),Vold(:),Wold(:),Told(:)
  REAL,ALLOCATABLE :: UCold(:),VCold(:),WCold(:),TCold(:)
  REAL,ALLOCATABLE :: UCoold(:),VCoold(:),WCoold(:),TCoold(:)
  REAL,ALLOCATABLE :: Uoold(:),Voold(:),Woold(:),Toold(:)
  REAL,ALLOCATABLE :: UDoold(:),VDoold(:),WDoold(:),TDoold(:)
  REAL,ALLOCATABLE :: UXold(:),VXold(:),WXold(:),TXold(:)
  REAL,ALLOCATABLE :: UXoold(:),VXoold(:),WXoold(:),TXoold(:)
  REAL,ALLOCATABLE :: Pold(:)
  REAL,ALLOCATABLE :: PPold(:)
  REAL,ALLOCATABLE :: Pxold(:),Pyold(:),Pzold(:)
!---- Variables for ReadC:
  CHARACTER  :: namU*38, namV*38, namW*38, naOut*38, naIn*80, namFin*80
  CHARACTER  :: answer*80, NameOut*80, namP*38, namT*38, name*80
!----------------------------------------------------------------------!
! The answer name is case dependent
!----------------------------------------------------------------------!

  write(*,*)' Enter the name of data files: '
  read(*,*) name 
  write(*,*)' Enter the name of case files: '
  read(*,*) namSav 
  write(*,*)' Enter the number of subdomains: '
  read(*,*) ND 
  write(*,*)' Enter 0 if case is not HOT and 1 if case is HOT: '
  read(*,*) IND


  answer = name
  naOut = answer
  naOut(len_trim(answer)+1:len_trim(answer)+4)='.xyz'

  answer = name
  namU = answer
  namU(len_trim(answer)+1:len_trim(answer)+4)='.U__'

  answer = name
  namV = answer
  namV(len_trim(answer)+1:len_trim(answer)+4)='.V__'

  answer = name
  namW = answer
  namW(len_trim(answer)+1:len_trim(answer)+4)='.W__'

  answer = name
  namP = answer
  namP(len_trim(answer)+1:len_trim(answer)+4)='.P__'

  answer = name
  namT = answer
  namT(len_trim(answer)+1:len_trim(answer)+4)='.T__'

  write(*,*)'Files to open are: ', naOut, namU, namV, namW, namP, namT
  open(5, FILE=naOut)
  read(5,*) NCold

  allocate (Xold(NCold)); Xold = 0.0
  allocate (Yold(NCold)); Yold = 0.0
  allocate (Zold(NCold)); Zold = 0.0
  allocate (Uold(NCold)); Uold = 0.0
  allocate (Vold(NCold)); Vold = 0.0
  allocate (Wold(NCold)); Wold = 0.0
  allocate (Told(NCold)); Told = 0.0
  allocate (Uoold(NCold)); Uoold = 0.0
  allocate (Voold(NCold)); Voold = 0.0
  allocate (Woold(NCold)); Woold = 0.0
  allocate (Toold(NCold)); Toold = 0.0
  allocate (UDoold(NCold)); UDoold = 0.0
  allocate (VDoold(NCold)); VDoold = 0.0
  allocate (WDoold(NCold)); WDoold = 0.0
  allocate (TDoold(NCold)); TDoold = 0.0
  allocate (UCold(NCold)); UCold = 0.0
  allocate (VCold(NCold)); VCold = 0.0
  allocate (WCold(NCold)); WCold = 0.0
  allocate (TCold(NCold)); TCold = 0.0
  allocate (UCoold(NCold)); UCoold = 0.0
  allocate (VCoold(NCold)); VCoold = 0.0
  allocate (WCoold(NCold)); WCoold = 0.0
  allocate (TCoold(NCold)); TCoold = 0.0
  allocate (UXold(NCold)); UXold = 0.0
  allocate (VXold(NCold)); VXold = 0.0
  allocate (WXold(NCold)); WXold = 0.0
  allocate (TXold(NCold)); TXold = 0.0
  allocate (UXoold(NCold)); UXoold = 0.0
  allocate (VXoold(NCold)); VXoold = 0.0
  allocate (WXoold(NCold)); WXoold = 0.0
  allocate (TXoold(NCold)); TXoold = 0.0
  allocate (Pold(NCold)); Pold = 0.0
  allocate (PPold(NCold)); PPold = 0.0
  allocate (Pxold(NCold)); Pxold = 0.0
  allocate (Pyold(NCold)); Pyold = 0.0
  allocate (Pzold(NCold)); Pzold = 0.0

  j = NCold
  do k = 1, j
    if(mod(k,500000) == 0) write(*,*) (100.*k/(1.*j)), '% complete...'
    read(5,*) Xold(k), Yold(k), Zold(k)
  end do
  close(5)
  write(*,*) ' Finished with reading ulaz.xyz file'
  open(5, FILE=namU)
  do k = 1, j
    if(mod(k,500000) == 0) write(*,*) (100.*k/(1.*j)), '% complete...'
    read(5,*) Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k)
  end do
  close(5)
  write(*,*) ' Finished with reading ulaz.U__ file'
  open(5, FILE=namV)
  do k = 1, j
    if(mod(k,500000) == 0) write(*,*) (100.*k/(1.*j)), '% complete...'
    read(5,*) Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k)
  end do
  close(5)
  write(*,*) ' Finished with reading ulaz.V__ file'
  open(5, FILE=namW)
  do k = 1, j
    if(mod(k,500000) == 0) write(*,*) (100.*k/(1.*j)), '% complete...'
      read(5,*) Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k)
    end do
  close(5)
  write(*,*) ' Finished with reading ulaz.W__ file'
  open(5, FILE=namP)
  do k = 1, j
    if(mod(k,500000) == 0) write(*,*) (100.*k/(1.*j)), '% complete...'
      read(5,*) Pold(k), PPold(k), Pxold(k), Pyold(k), Pzold(k)
    end do
  close(5)
  write(*,*) ' Finished with reading ulaz.P__ file'

  if(IND == 1) then
    open(5, FILE=namT)
    do k = 1, j
      if(mod(k,500000) == 0) write(*,*) (100.*k/(1.*j)), '% complete...'
        read(5,*) Told(k), Toold(k), TCold(k), TCoold(k), TDoold(k), TXold(k), TXoold(k) 
      end do
    close(5)
    write(*,*) ' Finished with reading ulaz.T__ file'
    write(*,*) 'LoaInI: finished with reading the files'
  end if

  NameOut = namSav
  if(ND > 1) then
    NameOut(len_trim(namSav)+1:len_trim(namSav)+5)="-0000"
    l1=len_trim(NameOut)

!----------------------------------------------------------------------!
! 48 is the number of subdomain of the new mesh,it is case dependent
!----------------------------------------------------------------------!
    do j = 1, ND 
      if(j  <  10) then
        write(NameOut(l1  :l1),'(I1)') j
      else if(j  < 100) then
        write(NameOut(l1-1:l1),'(I2)') j
      else
        write(NameOut(l1-2:l1),'(I3)') j
      end if
      write(*,*) NameOut

      name = NameOut

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
!     Read the binary file with the     *
!       connections between cells       *
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
      nameIn=name
      nameIn(len_trim(name)+1:len_trim(name)+4)='.cns'
      open(9, FILE=nameIn,FORM='UNFORMATTED')
      write(*,*) '# Now reading the binary .cns file:', nameIn

!///// number of cells, boundary cells and sides
      read(9) NC
      read(9) NbC
      read(9) NS
      read(9) NSsh
      read(9) Nmat

!///// cell materials
      allocate (material(-NbC:NC))
      read(9) (material(c), c=1,NC)
      read(9) (material(c), c=-1,-NBC,-1)

      close(9)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
!     Read the binary file with     *
!       geometrical quantities      *
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
      nameIn = name
      nameIn(len_trim(name)+1:len_trim(name)+4)='.geo'
      open(9, FILE=nameIn, FORM='UNFORMATTED')
      write(*,*) '# Now reading the binary .geo file:', nameIn

      allocate (xc(-NbC:NC))
      allocate (yc(-NbC:NC))
      allocate (zc(-NbC:NC))
      allocate (volume(-NbC:NC))
      allocate (delta(-NbC:NC))
      allocate (WallDs(NS))
      allocate (Sx(NS))
      allocate (Sy(NS))
      allocate (Sz(NS))
      allocate (Dx(NS))
      allocate (Dy(NS))
      allocate (Dz(NS))
      allocate (f(NS))
      allocate (xsp(NS))
      allocate (ysp(NS))
      allocate (zsp(NS))

      read(9) (xc(c), c=1,NC)
      read(9) (yc(c), c=1,NC)
      read(9) (zc(c), c=1,NC)
      read(9) (xc(c), c=-1,-NBC,-1)
      read(9) (yc(c), c=-1,-NBC,-1)
      read(9) (zc(c), c=-1,-NBC,-1)
      close(9)
 
      Xmax = -1.0e+6
      Xmin = 1.0e+6
      Ymax = -1.0e+6
      Ymin = 1.0e+6
      Zmax = -1.0e+6
      Zmin = 1.0e+6

      do c = 1, NC
        Xmax = max(xc(c),Xmax)
        Xmin = min(xc(c),Xmin)
        Ymax = max(yc(c),Ymax)
        Ymin = min(yc(c),Ymin)
        Zmax = max(zc(c),Zmax)
        Zmin = min(zc(c),Zmin)
      end do
    
!
! max will be increased for 5% in order to cover all cells in doman
!

      Xmax = Xmax * 1.01 
      Ymax = Ymax * 1.01 
      Zmax = Zmax * 1.01 
    
      namFin=name
      namFin(len_trim(name)+1:len_trim(name)+4)='.ini'
      write(*,*) namFin
      open(9,file = namFin)
      n = 0
      k = 0
      do k = 1, NCold
        if(Xold(k) <= Xmax.and.Xold(k) >= Xmin) then
          if(Yold(k) <= Ymax.and.Yold(k) >= Ymin) then
            if(Zold(k) <= Zmax.and.Zold(k) >= Zmin) then
              n = n + 1
            end if
          end if
        end if
      end do
      write(9,*) n
      NN = 0
      if(IND == 1) then
        do k = 1, NCold
          if(Xold(k) <= Xmax.and.Xold(k) >= Xmin) then
            if(Yold(k) <= Ymax.and.Yold(k) >= Ymin) then
              if(Zold(k) <= Zmax.and.Zold(k) >= Zmin) then
                NN = NN + 1
                write(9,*) Xold(k), Yold(k), Zold(k), &
                                    Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                                    Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                                    Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                                    Told(k), Toold(k), TCold(k), TCoold(k), TDoold(k), TXold(k), TXoold(k), &
                                    Pold(k), PPold(k), Pxold(k), Pyold(k), Pzold(k)
              end if
            end if
          end if
        end do
      else if(IND == 0) then
        do k = 1, NCold
          if(Xold(k) <= Xmax.and.Xold(k) >= Xmin) then
            if(Yold(k) <= Ymax.and.Yold(k) >= Ymin) then
              if(Zold(k) <= Zmax.and.Zold(k) >= Zmin) then
                NN = NN + 1
                write(9,*) Xold(k), Yold(k), Zold(k), &
                                    Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                                    Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                                    Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                                    Pold(k), PPold(k), Pxold(k), Pyold(k), Pzold(k)
              end if
            end if
          end if
        end do
      end if
      write(*,*) 'The number of cells in subdomain is: ', NC, 'Found cells: ', NN
      close(9)
      deallocate(material)
      deallocate (xc)
      deallocate (yc)
      deallocate (zc)
      deallocate (volume)
      deallocate (delta)
      deallocate (WallDs)
      deallocate (Sx)
      deallocate (Sy)
      deallocate (Sz)
      deallocate (Dx)
      deallocate (Dy)
      deallocate (Dz)
      deallocate (f)
      deallocate (xsp)
      deallocate (ysp)
      deallocate (zsp)
    end do
  else
    write(*,*) NameOut

    name = NameOut

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
!     Read the binary file with the     *
!       connections between cells       *
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
    nameIn=name
    nameIn(len_trim(name)+1:len_trim(name)+4)='.cns'
    open(9, FILE=nameIn,FORM='UNFORMATTED')
    write(*,*) '# Now reading the binary .cns file:', nameIn

!///// number of cells, boundary cells and sides
    read(9) NC
    read(9) NbC
    read(9) NS
    read(9) NSsh
    read(9) Nmat

!///// cell materials
    allocate (material(-NbC:NC))
    read(9) (material(c), c=1,NC)
    read(9) (material(c), c=-1,-NBC,-1)

    close(9)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
!     Read the binary file with     *
!       geometrical quantities      *
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
    nameIn = name
    nameIn(len_trim(name)+1:len_trim(name)+4)='.geo'
    open(9, FILE=nameIn, FORM='UNFORMATTED')
    write(*,*) '# Now reading the binary .geo file:', nameIn

    allocate (xc(-NbC:NC))
    allocate (yc(-NbC:NC))
    allocate (zc(-NbC:NC))
    allocate (volume(-NbC:NC))
    allocate (delta(-NbC:NC))
    allocate (WallDs(NS))
    allocate (Sx(NS))
    allocate (Sy(NS))
    allocate (Sz(NS))
    allocate (Dx(NS))
    allocate (Dy(NS))
    allocate (Dz(NS))
    allocate (f(NS))
    allocate (xsp(NS))
    allocate (ysp(NS))
    allocate (zsp(NS))

    read(9) (xc(c), c=1,NC)
    read(9) (yc(c), c=1,NC)
    read(9) (zc(c), c=1,NC)
    read(9) (xc(c), c=-1,-NBC,-1)
    read(9) (yc(c), c=-1,-NBC,-1)
    read(9) (zc(c), c=-1,-NBC,-1)
    close(9)
 
    Xmax = -1.0e+6
    Xmin = 1.0e+6
    Ymax = -1.0e+6
    Ymin = 1.0e+6
    Zmax = -1.0e+6
    Zmin = 1.0e+6

    do c = 1, NC
      Xmax = max(xc(c),Xmax)
      Xmin = min(xc(c),Xmin)
      Ymax = max(yc(c),Ymax)
      Ymin = min(yc(c),Ymin)
      Zmax = max(zc(c),Zmax)
      Zmin = min(zc(c),Zmin)
    end do

!
! max will be increased for 5% in order to cover all cells in doman
!

    Xmax = Xmax * 1.01 
    Ymax = Ymax * 1.01 
    Zmax = Zmax * 1.01 
    Xmin = Xmin * 0.99 
    Ymin = Ymin * 0.99 
    Zmin = Zmin * 0.99 
    
    namFin=name
    namFin(len_trim(name)+1:len_trim(name)+4)='.ini'
    write(*,*) namFin
    open(9,file = namFin)
    n = 0
    k = 0
    do k = 1, NCold
      if(Xold(k) <= Xmax.and.Xold(k) >= Xmin) then
        if(Yold(k) <= Ymax.and.Yold(k) >= Ymin) then
          if(Zold(k) <= Zmax.and.Zold(k) >= Zmin) then
            n = n + 1
          end if
        end if
      end if
    end do
    write(9,*) n
    NN = 0
    if(IND == 1) then
      do k = 1, NCold
        if(Xold(k) <= Xmax.and.Xold(k) >= Xmin) then
          if(Yold(k) <= Ymax.and.Yold(k) >= Ymin) then
            if(Zold(k) <= Zmax.and.Zold(k) >= Zmin) then
              NN = NN + 1
              write(9,*) Xold(k), Yold(k), Zold(k), &
                                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                                  Told(k), Toold(k), TCold(k), TCoold(k), TDoold(k), TXold(k), TXoold(k), &
                                  Pold(k), PPold(k), Pxold(k), Pyold(k), Pzold(k)
            end if
          end if
        end if
      end do
    else if(IND == 0) then
      do k = 1, NCold
        if(Xold(k) <= Xmax.and.Xold(k) >= Xmin) then
          if(Yold(k) <= Ymax.and.Yold(k) >= Ymin) then
            if(Zold(k) <= Zmax.and.Zold(k) >= Zmin) then
              NN = NN + 1
              write(9,*) Xold(k), Yold(k), Zold(k), &
                                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                                  Pold(k), PPold(k), Pxold(k), Pyold(k), Pzold(k)
            end if
          end if
        end if
      end do
    end if
    write(*,*) 'The number of cells in subdomain is: ', NC, 'Found cells: ', NN, Xmax, Ymax, Zmax
    close(9)
    deallocate(material)
    deallocate (xc)
    deallocate (yc)
    deallocate (zc)
    deallocate (volume)
    deallocate (delta)
    deallocate (WallDs)
    deallocate (Sx)
    deallocate (Sy)
    deallocate (Sz)
    deallocate (Dx)
    deallocate (Dy)
    deallocate (Dz)
    deallocate (f)
    deallocate (xsp)
    deallocate (ysp)
    deallocate (zsp)
  end if
  deallocate(Xold)
  deallocate(Yold)
  deallocate(Zold)
  deallocate(Uold)
  deallocate(Vold)
  deallocate(Wold)
  deallocate(Told)
  deallocate(Uoold)
  deallocate(Voold)
  deallocate(Woold)
  deallocate(Toold)
  deallocate(UDoold)
  deallocate(VDoold)
  deallocate(WDoold)
  deallocate(TDoold)
  deallocate(UCold)
  deallocate(VCold)
  deallocate(WCold)
  deallocate(TCold)
  deallocate(UCoold)
  deallocate(VCoold)
  deallocate(WCoold)
  deallocate(TCoold)
  deallocate(UXold)
  deallocate(VXold)
  deallocate(WXold)
  deallocate(TXold)
  deallocate(UXoold)
  deallocate(VXoold)
  deallocate(WXoold)
  deallocate(TXoold)
  deallocate(Pold)
  deallocate(PPold)
  deallocate(Pxold)
  deallocate(Pyold)
  deallocate(Pzold)

  END PROGRAM

