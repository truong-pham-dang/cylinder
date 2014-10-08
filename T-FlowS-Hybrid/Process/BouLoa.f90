!======================================================================!
  SUBROUTINE BouLoa(in_out)
!----------------------------------------------------------------------!
! Reads: NAME.b                                                        !
! ~~~~~~                                                               !
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE rans_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  LOGICAL       :: in_out
!-------------------------------[Locals]-------------------------------!
  INTEGER       :: c, n, dum1, NB, NP, Ninit, m 
  CHARACTER*80  :: namBou, namPro(128), dir
  INTEGER       :: typBou(128)
  REAL          :: xyz(1024)
  REAL          :: wi
  LOGICAL       :: here
!--------------------------------[CVS]---------------------------------!
!  $Id: BouLoa.f90,v 1.38 2009/06/30 11:42:00 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/BouLoa.f90,v $  
!======================================================================!

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Read the file with boundary conditions     !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  namBou = name
  namBou(len_trim(name)+1:len_trim(name)+2) = '.b'
  open(9, FILE=namBou)
  if(this < 2) write(*,*) '# Now reading the file:', namBou

!---------------------!
! Phisical properties !
!---------------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) Nmat
  do n=1,Nmat
    call ReadC(9,inp,tn,ts,te)
    call ToUppr(  inp(ts(2):te(2))  )
    if( inp(ts(2):te(2))  ==  'FLUID') then 
      StateMat(n)=FLUID
    else if( inp(ts(2):te(2))  ==  'SOLID') then 
      StateMat(n)=SOLID
    else 
      if(this < 2) write(*,*) 'BouLoa: Unknown material state'
      stop  
    end if
    read(inp(ts(3):te(3)),*) VISc
    read(inp(ts(4):te(4)),*) DENc(n)
    if(HOT==YES) read(inp(ts(5):te(5)),*) CONc(n)
    if(HOT==YES) read(inp(ts(6):te(6)),*) CAPc(n)
  end do
!---------------------!
! Boundary conditions !
!---------------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) NB
  do n=1,NB
    call ReadC(9,inp,tn,ts,te)
    call ToUppr(  inp(ts(2):te(2))  )
    call ToUppr(  inp(ts(3):te(3))  )
    read(inp(ts(1):te(1)),*) dum1
    if( inp(ts(2):te(2)) == 'INFLOW') then 
      typBou(n)=INFLOW
    else if( inp(ts(2):te(2)) == 'WALL') then 
      typBou(n)=WALL
    else if( inp(ts(2):te(2)) == 'OUTFLOW') then 
      typBou(n)=OUTFLOW
    else if( inp(ts(2):te(2)) == 'SYMMETRY') then 
      typBou(n)=SYMMETRY
    else if( inp(ts(2):te(2)) == 'WALLFLUX') then 
      typBou(n)=WALLFL
    else if( inp(ts(2):te(2)) == 'CONVECTIVE') then 
      typBou(n)=CONVECT
    else if( inp(ts(2):te(2)) == 'PRESSURE') then 
      typBou(n)=PRESSURE
    else
      if(this < 2) write(*,*) 'BouLoa: Unknown boundary condition type: ', inp(ts(2):te(2))
      stop  
    end if
    if( inp(ts(3):te(3))  ==  'FILE') then
      read(inp(ts(4):te(4)),'(A80)') namPro(n)
    else
      read(inp(ts(3):te(3)),*) U % bound(n)
      read(inp(ts(4):te(4)),*) V % bound(n)
      read(inp(ts(5):te(5)),*) W % bound(n)
      if(typBou(n)==PRESSURE) then
        read(inp(ts(6):te(6)),*) P % bound(n)
        if(HOT==YES) then 
          read(inp(ts(7):te(7)),*) T % bound(n)
          if(SIMULA==EBM.or.SIMULA==HJ) then
            read(inp(ts(8):te(8)),*)   uu % bound(n)
            read(inp(ts(9):te(9)),*)   vv % bound(n)
            read(inp(ts(10):te(10)),*) ww % bound(n)
            read(inp(ts(11):te(11)),*) uv % bound(n)
            read(inp(ts(12):te(12)),*) uw % bound(n)
            read(inp(ts(13):te(13)),*) vw % bound(n)
            read(inp(ts(14):te(14)),*) Eps% bound(n)
            if(SIMULA==EBM) read(inp(ts(15):te(15)),*) f22 % bound(n)
          end if
          if(SIMULA==K_EPS) then
            read(inp(ts(8):te(8)),*) Kin % bound(n)
            read(inp(ts(9):te(9)),*) Eps % bound(n)
          end if
          if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
            read(inp(ts(8):te(8)),*) Kin % bound(n)
            read(inp(ts(9):te(9)),*) Eps % bound(n)
            read(inp(ts(10):te(10)),*) v_2 % bound(n)
            read(inp(ts(11):te(11)),*) f22 % bound(n)
          end if
          if(SIMULA == SPA_ALL) then
            read(inp(ts(8):te(8)),*) VIS % bound(n)
          end if
          if(SIMULA == DES_SPA) then
            read(inp(ts(8):te(8)),*) VIS % bound(n)
          end if
        else
          if(SIMULA==EBM.or.SIMULA==HJ) then
            read(inp(ts(7):te(7)),*)   uu % bound(n)
            read(inp(ts(8):te(8)),*)   vv % bound(n)
            read(inp(ts(9):te(9)),*) ww % bound(n)
            read(inp(ts(10):te(10)),*) uv % bound(n)
            read(inp(ts(11):te(11)),*) uw % bound(n)
            read(inp(ts(12):te(12)),*) vw % bound(n)
            read(inp(ts(13):te(13)),*) Eps% bound(n)
            if(SIMULA==EBM) read(inp(ts(14):te(14)),*) f22 % bound(n)
          end if
          if(SIMULA==K_EPS) then
            read(inp(ts(7):te(7)),*) Kin % bound(n)
            read(inp(ts(8):te(8)),*) Eps % bound(n)
          end if
          if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
            read(inp(ts(7):te(7)),*) Kin % bound(n)
            read(inp(ts(8):te(8)),*) Eps % bound(n)
            read(inp(ts(9):te(9)),*) v_2  % bound(n)
            read(inp(ts(10):te(10)),*) f22 % bound(n)
          end if
          if(SIMULA == SPA_ALL) then
            read(inp(ts(7):te(7)),*) VIS % bound(n)
          end if
          if(SIMULA == DES_SPA) then
            read(inp(ts(7):te(7)),*) VIS % bound(n)
          end if
        end if
        namPro(n)=''
      else  
        if(HOT==YES) then 
          read(inp(ts(6):te(6)),*) T % bound(n)
          if(SIMULA==EBM.or.SIMULA==HJ) then
            read(inp(ts(7):te(7)),*) uu % bound(n)
            read(inp(ts(8):te(8)),*) vv % bound(n)
            read(inp(ts(9):te(9)),*) ww % bound(n)
            read(inp(ts(10):te(10)),*) uv % bound(n)
            read(inp(ts(11):te(11)),*) uw % bound(n)
            read(inp(ts(12):te(12)),*) vw % bound(n)
            read(inp(ts(13):te(13)),*) Eps% bound(n)
            if(SIMULA==EBM) read(inp(ts(14):te(14)),*) f22 % bound(n)
          end if
          if(SIMULA==K_EPS) then
            read(inp(ts(7):te(7)),*) Kin % bound(n)
            read(inp(ts(8):te(8)),*) Eps % bound(n)
          end if
          if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
            read(inp(ts(7):te(7)),*) Kin % bound(n)
            read(inp(ts(8):te(8)),*) Eps % bound(n)
            read(inp(ts(9):te(9)),*) v_2 % bound(n)
            read(inp(ts(10):te(10)),*) f22 % bound(n)
          end if
          if(SIMULA == SPA_ALL) then
            read(inp(ts(7):te(7)),*) VIS % bound(n)
          end if
          if(SIMULA == DES_SPA) then
            read(inp(ts(7):te(7)),*) VIS % bound(n)
          end if
        else
          if(SIMULA==EBM.or.SIMULA==HJ) then
            read(inp(ts(6):te(6)),*) uu % bound(n)
            read(inp(ts(7):te(7)),*) vv % bound(n)
            read(inp(ts(8):te(8)),*) ww % bound(n)
            read(inp(ts(9):te(9)),*) uv % bound(n)
            read(inp(ts(10):te(10)),*) uw % bound(n)
            read(inp(ts(11):te(11)),*) vw % bound(n)
            read(inp(ts(12):te(12)),*) Eps% bound(n)
            if(SIMULA==EBM) read(inp(ts(13):te(13)),*) f22 % bound(n)
          end if
          if(SIMULA==K_EPS) then
            read(inp(ts(6):te(6)),*) Kin % bound(n)
            read(inp(ts(7):te(7)),*) Eps % bound(n)
          end if
          if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
            read(inp(ts(6):te(6)),*) Kin % bound(n)
            read(inp(ts(7):te(7)),*) Eps % bound(n)
            read(inp(ts(8):te(8)),*) v_2  % bound(n)
            read(inp(ts(9):te(9)),*) f22 % bound(n)
          end if
          if(SIMULA == SPA_ALL) then
            read(inp(ts(6):te(6)),*) VIS % bound(n)
          end if
          if(SIMULA == DES_SPA) then
            read(inp(ts(6):te(6)),*) VIS % bound(n)
          end if
        end if
        namPro(n)=''
      end if
    end if
  end do      

!--------------------!
! Initial conditions !
!--------------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) Ninit
  if(Ninit > Nmat) then
    if(this < 2) write(*,*) 'Warning: there are more initial conditions then materials'
  end if

  do n=1,Ninit
    call ReadC(9,inp,tn,ts,te)
    call ToUppr(inp(ts(2):te(2)))

!-----Initial conditions given in GMV file
    if(inp(ts(2):te(2)) == 'FILE') then
      read(inp(ts(3):te(3)),'(A80)') namIni(n)
      write(*,*) '@BouLoa: material ', n, '; init. cond. given by file: ', namIni(n)
    else
      namIni(n) = ''

!-----Initial conditions given by constant
    read(inp(ts(2):te(2)),*) U % init(n)
    read(inp(ts(3):te(3)),*) V % init(n)
    read(inp(ts(4):te(4)),*) W % init(n)
 
    if(HOT==YES) then
      read(inp(ts(5):te(5)),*) T % init(n)
      if(SIMULA==EBM.or.SIMULA==HJ) then
        read(inp(ts(6):te(6)),*) uu % init(n)
        read(inp(ts(7):te(7)),*) vv % init(n)
        read(inp(ts(8):te(8)),*) ww % init(n)
        read(inp(ts(9):te(9)),*) uv % init(n)
        read(inp(ts(10):te(10)),*) uw % init(n)
        read(inp(ts(11):te(11)),*) vw % init(n)
        read(inp(ts(12):te(12)),*) Eps% init(n)
        if(SIMULA==EBM) read(inp(ts(13):te(13)),*) f22 % init(n)
      end if
      if(SIMULA==K_EPS) then
        read(inp(ts(6):te(6)),*) Kin % init(n)
        read(inp(ts(7):te(7)),*) Eps % init(n)
      end if
      if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
        read(inp(ts(6):te(6)),*) Kin % init(n)
        read(inp(ts(7):te(7)),*) Eps % init(n)
        read(inp(ts(8):te(8)),*) v_2  % init(n)
        read(inp(ts(9):te(9)),*) f22 % init(n)
      end if
      if(SIMULA == SPA_ALL) then
        read(inp(ts(6):te(6)),*) VIS % init(n)
      end if
      if(SIMULA == DES_SPA) then
        read(inp(ts(6):te(6)),*) VIS % init(n)
      end if
    else ! HOT /= YES
      if(SIMULA==EBM.or.SIMULA==HJ) then
        read(inp(ts(5):te(5)),*) uu % init(n)
        read(inp(ts(6):te(6)),*) vv % init(n)
        read(inp(ts(7):te(7)),*) ww % init(n)
        read(inp(ts(8):te(8)),*) uv % init(n)
        read(inp(ts(9):te(9)),*) uw % init(n)
        read(inp(ts(10):te(10)),*) vw % init(n)
        read(inp(ts(11):te(11)),*) Eps% init(n)
        if(SIMULA==EBM) read(inp(ts(12):te(12)),*) f22 % init(n)
      end if
      if(SIMULA==K_EPS) then
        read(inp(ts(5):te(5)),*) Kin % init(n)
        read(inp(ts(6):te(6)),*) Eps % init(n)
      end if
      if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
        read(inp(ts(5):te(5)),*) Kin % init(n)
        read(inp(ts(6):te(6)),*) Eps % init(n)
        read(inp(ts(7):te(7)),*) v_2  % init(n)
        read(inp(ts(8):te(8)),*) f22 % init(n)
      end if
      if(SIMULA == SPA_ALL) then
        read(inp(ts(5):te(5)),*) VIS % init(n)
      end if
      if(SIMULA == DES_SPA) then
        read(inp(ts(5):te(5)),*) VIS % init(n)
      end if
    end if
    end if
  end do  

  close(9)

!----------------------------------!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!----------------------------------!
  do n=1,NB

!---- Boundary condition is given by a single constant
    if(namPro(n) == '') then 
      do c=-1,-NbC,-1
	if(bcmark(c) == n) then
	  TypeBC(c) = typBou(n)

!===== if in_out is set to true, set boundary values,
!===== otherwise, just the TypeBC remains set.

          if(in_out) then
	    U % n(c) = U % bound(n) 
	    V % n(c) = V % bound(n)
	    W % n(c) = W % bound(n)
            P % n(c) = P % bound(n) 
            if(HOT == YES) then
              if(TypeBC(c).eq.WALLFL) then
                T % q(c) =  T % bound(n)
              else
                T % n(c) =  T % bound(n)
              endif
            end if  ! for HOT==YES
            if(SIMULA==EBM.or.SIMULA==HJ) then
              uu % n(c) = uu % bound(n)
              vv % n(c) = vv % bound(n)
              ww % n(c) = ww % bound(n)
              uv % n(c) = uv % bound(n)
              uw % n(c) = uw % bound(n)
              vw % n(c) = vw % bound(n)
              Eps % n(c) = Eps % bound(n)
              if(SIMULA==EBM) f22 % n(c)   = f22 % bound(n)
            end if
            if(SIMULA==K_EPS) then
              Kin % n(c) = Kin % bound(n)
              Eps % n(c) = Eps % bound(n)
              Uf(c)        = 0.047
              Ynd(c)       = 30.0
            end if
            if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
              Kin % n(c)   = Kin % bound(n)
              Eps % n(c)   = Eps % bound(n)
              f22 % n(c)   = f22 % bound(n)
              v_2 % n(c)   = v_2 % bound(n)
            end if
            if(SIMULA == SPA_ALL) then
              VIS % n(c)   = VIS % bound(n)
            end if
            if(SIMULA == DES_SPA) then
              VIS % n(c)   = VIS % bound(n)
            end if
          end if
	end if 
      end do
!---- Boundary condition is prescribed in a file 
    else
      open(9, FILE=namPro(n))
      if(this < 2) write(*,*) '# Now reading the file:', namPro(n)
      call ReadC(9,inp,tn,ts,te)
      read(inp(ts(1):te(1)),*) NP                  ! number of points
      if(NP  > 1000) then
	if(this < 2) write(*,*) 'BouLoa: Too much points in a profile !'
	stop  
      end if
      call ReadC(9,inp,tn,ts,te)
      read(inp(ts(1):te(1)),*) dir  ! direction
      do m=1,NP
	call ReadC(9,inp,tn,ts,te)
	read(inp(ts(1):te(1)),*) xyz(m)
	read(inp(ts(2):te(2)),*) U % pro(m)
	read(inp(ts(3):te(3)),*) V % pro(m)
	read(inp(ts(4):te(4)),*) W % pro(m)
	if(HOT==YES) then
          read(inp(ts(5):te(5)),*) T % pro(m)
          if(SIMULA==K_EPS) then
            read(inp(ts(6):te(6)),*) Kin % pro(m)
            read(inp(ts(7):te(7)),*) Eps % pro(m)
          end if
          if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
            read(inp(ts(6):te(6)),*) Kin % pro(m)
            read(inp(ts(7):te(7)),*) Eps % pro(m)
            read(inp(ts(8):te(8)),*) v_2 % pro(m)
            read(inp(ts(9):te(9)),*) f22 % pro(m)
          end if
          if(SIMULA == SPA_ALL) then
            read(inp(ts(6):te(6)),*) VIS % pro(m)
          end if
          if(SIMULA == DES_SPA) then
            read(inp(ts(6):te(6)),*) VIS % pro(m)
          end if
        else
          if(SIMULA==K_EPS) then
            read(inp(ts(5):te(5)),*) Kin % pro(m)
            read(inp(ts(6):te(6)),*) Eps % pro(m)
          end if
          if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
            read(inp(ts(5):te(5)),*) Kin % pro(m)
            read(inp(ts(6):te(6)),*) Eps % pro(m)
            read(inp(ts(7):te(7)),*) v_2 % pro(m)
            read(inp(ts(8):te(8)),*) f22 % pro(m)
          end if
          if(SIMULA == SPA_ALL) then
            read(inp(ts(5):te(5)),*) VIS % pro(m)
          end if
          if(SIMULA == DES_SPA) then
            read(inp(ts(5):te(5)),*) VIS % pro(m)
          end if
        end if
      end do

      do c=-1,-NbC,-1
	if(bcmark(c) == n) then
	  TypeBC(c) = typBou(n)
          
!===== if in_out is set to true, set boundary values,
!===== otherwise, just the TypeBC remains set.

          if(in_out) then
	    do m=1,NP-1
              here = .FALSE. 
!----- compute the weight factors
              if( (dir == 'X' .or. dir == 'x') .and.                  &
                   xc(c) >= xyz(m) .and. xc(c) <= xyz(m+1) ) then
                wi = ( xyz(m+1)-xc(c) ) / ( xyz(m+1) - xyz(m) )
                here = .TRUE.
              else if( (dir == 'Y' .or. dir == 'y') .and.             &
                   yc(c) >= xyz(m) .and. yc(c) <= xyz(m+1) ) then
                wi = ( xyz(m+1)-yc(c) ) / ( xyz(m+1) - xyz(m) )
                here = .TRUE.
              else if( (dir == 'Z' .or. dir == 'z') .and.             &
                   zc(c) >= xyz(m) .and. zc(c) <= xyz(m+1) ) then
                wi = ( xyz(m+1)-zc(c) ) / ( xyz(m+1) - xyz(m) )
                here = .TRUE.
              else if( (dir == 'RX' .or. dir == 'rx') .and.           &
                   sqrt(yc(c)*yc(c)+zc(c)*zc(c)) >= xyz(m) .and.      &
                   sqrt(yc(c)*yc(c)+zc(c)*zc(c)) <= xyz(m+1) ) then
                wi = ( xyz(m+1) - sqrt(yc(c)*yc(c)+zc(c)*zc(c)) )     &
                   / ( xyz(m+1) - xyz(m) )
                here = .TRUE.
              else if( (dir == 'RY' .or. dir == 'ry') .and.           &
                   sqrt(xc(c)*xc(c)+zc(c)*zc(c)) >= xyz(m) .and.      &
                   sqrt(xc(c)*xc(c)+zc(c)*zc(c)) <= xyz(m+1) ) then
                wi = ( xyz(m+1) - sqrt(xc(c)*xc(c)+zc(c)*zc(c)) )     &
                   / ( xyz(m+1) - xyz(m) )
                here = .TRUE.
              else if( (dir == 'RZ' .or. dir == 'rz') .and.           &
                   sqrt(xc(c)*xc(c)+yc(c)*yc(c)) <= xyz(m) .and.      &
                   sqrt(xc(c)*xc(c)+yc(c)*yc(c)) >= xyz(m+1) ) then
                wi = ( xyz(m+1) - sqrt(xc(c)*xc(c)+yc(c)*yc(c)) )     &
                   / ( xyz(m+1) - xyz(m) )
                here = .TRUE.
              end if

!----- interpolate the profiles     
              if(here) then
                U % n(c) = wi*U % pro(m) + (1.-wi)*U % pro(m+1)
                V % n(c) = wi*V % pro(m) + (1.-wi)*V % pro(m+1)
                W % n(c) = wi*W % pro(m) + (1.-wi)*W % pro(m+1)
                if(HOT==YES) &
                  T % n(c) = wi*T % pro(m) + (1.-wi)*T % pro(m+1)
                if(SIMULA==K_EPS) then
                  Kin % n(c) = wi*Kin % pro(m) + (1.-wi)*Kin % pro(m+1)
                  Eps % n(c) = wi*Eps % pro(m) + (1.-wi)*Eps % pro(m+1)
                end if
                if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
                  Kin % n(c) = wi*Kin % pro(m) + (1.-wi)*Kin % pro(m+1)
                  Eps % n(c) = wi*Eps % pro(m) + (1.-wi)*Eps % pro(m+1)
                  f22 % n(c) = wi*f22 % pro(m) + (1.-wi)*f22 % pro(m+1)
                  v_2 % n(c) = wi*v_2 % pro(m)  + (1.-wi)*v_2% pro(m+1)
                end if
                if(SIMULA == SPA_ALL) then
                  VIS % n(c) = wi*VIS % pro(m) + (1.-wi)*VIS % pro(m+1)
                end if
                if(SIMULA == DES_SPA) then
                  VIS % n(c) = wi*VIS % pro(m) + (1.-wi)*VIS % pro(m+1)
                end if
              end if
	    end do
          end if  ! if(in_out)
	end if 
      end do
      close(9)
    end if
  end do 

!---------------------------------!
! Finally handle the buffer cells !
!---------------------------------!
  do c=-1,-NbC,-1
    if(bcmark(c) == BUFFER) TypeBC(c)=BUFFER 
  end do

  RETURN 

  END SUBROUTINE BouLoa
