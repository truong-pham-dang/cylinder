!======================================================================!
  SUBROUTINE LoaIni()
!----------------------------------------------------------------------!
! This version of LoaIni is optimised for very large meshes
! Program SUB_INI needs to be used to create files needed by this
! subroutine 
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER          :: j, k,  c, nearest, var, Nvar, c1, c2, s 
  INTEGER          :: NCold  
  TYPE(Unknown)    :: PHI
  REAL,ALLOCATABLE :: Xold(:),Yold(:),Zold(:)
  REAL,ALLOCATABLE :: Uold(:),Vold(:),Wold(:),Told(:), Kold(:), Eold(:), v_2old(:), f22old(:)
  REAL,ALLOCATABLE :: UCold(:),VCold(:),WCold(:),TCold(:), KCold(:), ECold(:), v_2Cold(:), f22Cold(:)
  REAL,ALLOCATABLE :: UCoold(:),VCoold(:),WCoold(:),TCoold(:), KCoold(:), ECoold(:), v_2Coold(:), f22Coold(:)
  REAL,ALLOCATABLE :: Uoold(:),Voold(:),Woold(:),Toold(:), Koold(:), Eoold(:), v_2oold(:), f22oold(:)
  REAL,ALLOCATABLE :: UDoold(:),VDoold(:),WDoold(:),TDoold(:), KDoold(:), EDoold(:), v_2Doold(:), f22Doold(:)
  REAL,ALLOCATABLE :: UXold(:),VXold(:),WXold(:),TXold(:), KXold(:), EXold(:), v_2Xold(:), f22Xold(:)
  REAL,ALLOCATABLE :: UXoold(:),VXoold(:),WXoold(:),TXoold(:), KXoold(:), EXoold(:), v_2Xoold(:), f22Xoold(:)
  REAL,ALLOCATABLE :: Pold(:)
  REAL,ALLOCATABLE :: PPold(:)
  REAL,ALLOCATABLE :: Pxold(:),Pyold(:),Pzold(:)

  REAL             :: DISTnew, DISTold
  REAL             :: Dist, Us, Ws, Vs, R, Rnew, Rold, DistR

!---- Variables for ReadC:
  CHARACTER  :: namCoo*80, ext*4, answer*80, answer_hot*80
  CHARACTER*80 :: nameIn
!--------------------------------[CVS]---------------------------------!
  character*80 rcs1,rcs2
  data rcs1/                                                        &
  '$Id: LoaIni.f90,v 1.4 2008/12/10 14:44:15 IUS\mhadziabdic Exp $'/
  data rcs2/                                                        &
  '$Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/LoaIni.f90,v $'/
!======================================================================!  

  call ReadC(7,inp,tn,ts,te)
!  read(inp, '(A80)')  namCoo
  read(inp, *)  namCoo
  answer=namCoo
  call ToUppr(answer)
  if(answer == 'SKIP') then
    return
  end if  
!---- save the name
  answer = name
  name = namCoo 

  if(tn==2) read(inp(ts(2):te(2)),'(A8)') answer
  call ToUppr(answer)
  if(answer == 'HOT') then
    HOTini = YES
  end if

  call NamFil(THIS, nameIn, '.ini', len_trim('.ini'))

  if(this < 2) write(*,*)'now reading file:', nameIn 

  open(5, FILE=nameIn) 
  read(5,*) NCold

  allocate (Xold(NCold)); Xold = 0.0
  allocate (Yold(NCold)); Yold = 0.0
  allocate (Zold(NCold)); Zold = 0.0
  allocate (Uold(NCold)); Uold = 0.0
  allocate (Vold(NCold)); Vold = 0.0
  allocate (Wold(NCold)); Wold = 0.0
  allocate (Uoold(NCold)); Uoold = 0.0
  allocate (Voold(NCold)); Voold = 0.0
  allocate (Woold(NCold)); Woold = 0.0
  allocate (UDoold(NCold)); UDoold = 0.0
  allocate (VDoold(NCold)); VDoold = 0.0
  allocate (WDoold(NCold)); WDoold = 0.0
  allocate (UCold(NCold)); UCold = 0.0
  allocate (VCold(NCold)); VCold = 0.0
  allocate (WCold(NCold)); WCold = 0.0
  allocate (UCoold(NCold)); UCoold = 0.0
  allocate (VCoold(NCold)); VCoold = 0.0
  allocate (WCoold(NCold)); WCoold = 0.0
  allocate (UXold(NCold)); UXold = 0.0
  allocate (VXold(NCold)); VXold = 0.0
  allocate (WXold(NCold)); WXold = 0.0
  allocate (UXoold(NCold)); UXoold = 0.0
  allocate (VXoold(NCold)); VXoold = 0.0
  allocate (WXoold(NCold)); WXoold = 0.0
  allocate (Pold(NCold)); Pold = 0.0
  allocate (PPold(NCold)); PPold = 0.0
  allocate (Pxold(NCold)); Pxold = 0.0
  allocate (Pyold(NCold)); Pyold = 0.0
  allocate (Pzold(NCold)); Pzold = 0.0
  if(HOTini == YES) then
    allocate (Told(NCold)); Told = 0.0
    allocate (Toold(NCold)); Toold = 0.0
    allocate (TDoold(NCold)); TDoold = 0.0
    allocate (TCold(NCold)); TCold = 0.0
    allocate (TCoold(NCold)); TCoold = 0.0
    allocate (TXold(NCold)); TXold = 0.0
    allocate (TXoold(NCold)); TXoold = 0.0
  end if
  allocate (Kold(NCold)); Kold = 0.0
  allocate (Koold(NCold)); Koold = 0.0
  allocate (KDoold(NCold)); KDoold = 0.0
  allocate (KCold(NCold)); KCold = 0.0
  allocate (KCoold(NCold)); KCoold = 0.0
  allocate (KXold(NCold)); KXold = 0.0
  allocate (KXoold(NCold)); KXoold = 0.0

  allocate (Eold(NCold)); Eold = 0.0
  allocate (Eoold(NCold)); Eoold = 0.0
  allocate (EDoold(NCold)); EDoold = 0.0
  allocate (ECold(NCold)); ECold = 0.0
  allocate (ECoold(NCold)); ECoold = 0.0
  allocate (EXold(NCold)); EXold = 0.0
  allocate (EXoold(NCold)); EXoold = 0.0

  allocate (v_2old(NCold)); v_2old = 0.0
  allocate (v_2oold(NCold)); v_2oold = 0.0
  allocate (v_2Doold(NCold)); v_2Doold = 0.0
  allocate (v_2Cold(NCold)); v_2Cold = 0.0
  allocate (v_2Coold(NCold)); v_2Coold = 0.0
  allocate (v_2Xold(NCold)); v_2Xold = 0.0
  allocate (v_2Xoold(NCold)); v_2Xoold = 0.0

  allocate (F22old(NCold)); F22old = 0.0
  allocate (F22oold(NCold)); F22oold = 0.0
  allocate (F22Doold(NCold)); F22Doold = 0.0
  allocate (F22Xold(NCold)); F22Xold = 0.0
  allocate (F22Xoold(NCold)); F22Xoold = 0.0

  j = NCold
  do k = 1, j
    if(this < 2) then
      if(mod(k,20000) == 0) write(*,*) (100.*k/(1.*j)), '% complete...'  
    end if
    if(SIMULA == LES) then 
      if(HOTini==YES) then
        read(5,*) Xold(k), Yold(k), Zold(k), &
                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                  Told(k), Toold(k), TCold(k), TCoold(k), TDoold(k), TXold(k), TXoold(k), &
                  Pold(k), PPold(k), Pxold(k), Pyold(k), Pzold(k)   
      else
        read(5,*) Xold(k), Yold(k), Zold(k), &
                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                  Pold(k), PPold(k), Pxold(k), Pyold(k), Pzold(k)   
      end if
    end if 
    if(SIMULA == ZETA.or.SIMULA == K_EPS_VV) then 
      if(HOTini==YES) then
        read(5,*) Xold(k), Yold(k), Zold(k), &
                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                  Told(k), Toold(k), TCold(k), TCoold(k), TDoold(k), TXold(k), TXoold(k), &
                  Kold(k), Koold(k), KCold(k), KCoold(k), KDoold(k), KXold(k), KXoold(k), &
                  Eold(k), Eoold(k), ECold(k), ECoold(k), EDoold(k), EXold(k), EXoold(k), &
                  v_2old(k), v_2oold(k), v_2Cold(k), v_2Coold(k), v_2Doold(k), v_2Xold(k), v_2Xoold(k), &
                  F22old(k), F22oold(k), F22Doold(k), F22Xold(k), F22Xoold(k), &
                  Pold(k), PPold(k), Pxold(k), Pyold(k), Pzold(k)   
      else
        read(5,*) Xold(k), Yold(k), Zold(k), &
                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                  Kold(k), Koold(k), KCold(k), KCoold(k), KDoold(k), KXold(k), KXoold(k), &
                  Eold(k), Eoold(k), ECold(k), ECoold(k), EDoold(k), EXold(k), EXoold(k), &
                  v_2old(k), v_2oold(k), v_2Cold(k), v_2Coold(k), v_2Doold(k), v_2Xold(k), v_2Xoold(k), &
                  F22old(k), F22oold(k), F22Doold(k), F22Xold(k), F22Xoold(k), &
                  Pold(k), PPold(k), Pxold(k), Pyold(k), Pzold(k)   
      end if
    end if 
    if(SIMULA == K_EPS) then 
      if(HOTini==YES) then
        read(5,*) Xold(k), Yold(k), Zold(k), &
                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                  Told(k), Toold(k), TCold(k), TCoold(k), TDoold(k), TXold(k), TXoold(k), &
                  Kold(k), Koold(k), KCold(k), KCoold(k), KDoold(k), KXold(k), KXoold(k), &
                  Eold(k), Eoold(k), ECold(k), ECoold(k), EDoold(k), EXold(k), EXoold(k), &
                  Pold(k), PPold(k), Pxold(k), Pyold(k), Pzold(k)   
      else
        read(5,*) Xold(k), Yold(k), Zold(k), &
                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                  Kold(k), Koold(k), KCold(k), KCoold(k), KDoold(k), KXold(k), KXoold(k), &
                  Eold(k), Eoold(k), ECold(k), ECoold(k), EDoold(k), EXold(k), EXoold(k), &
                  Pold(k), PPold(k), Pxold(k), Pyold(k), Pzold(k)   
      end if
    end if 
  end do
  close(5)
  if(this < 2) write(*,*) 'LoaInI: finished with reading the files'

  nearest = 0
  near = 0
  DISTold = HUGE
  if(JET == YES.or.PIPE == YES) then
    do c = 1, NC
      if(this < 2) then
        if(mod(c,20000) == 0) write(*,*) (100.*c/(1.*NC)), '% complete...'  
      end if
      DISTold = HUGE
      do k = 1, j
        Rold = sqrt(Xold(k)*Xold(k) + Yold(k)*Yold(k))
        Rnew = sqrt(xc(c)*xc(c) + yc(c)*yc(c))
        DISTnew = DistR(Rold,Zold(k),Rnew,zc(c))
        if(DISTnew <= DISTold) then
          nearest =  k
          DISTold = DISTnew
          if(DISTold<1.0e-8) go to 100
        end if 
      end do
100   continue   
      near(c) = nearest 
    end do
  else 
    do c = 1, NC
      if(this < 2) then
        if(mod(c,20000) == 0) write(*,*) (100.*c/(1.*NC)), '% complete...'  
      end if
      DISTold = HUGE
      do k = 1, j
        DISTnew = Dist(Xold(k),Yold(k),Zold(k),xc(c),yc(c),zc(c))
        if(DISTnew <= DISTold) then
          nearest =  k
          DISTold = DISTnew
        end if 
      end do
      near(c) = nearest 
    end do
  end if 


  if(this < 2) write(*,*) 'LoaInI: finished with sorting the cells'

  do c = 1, NC
    U % n(c)  = Uold(near(c)) 
    U % o(c)  = Uoold(near(c)) 
    U % C(c)  = UCold(near(c)) 
    U % Co(c) = UCoold(near(c)) 
    U % Do(c) = UDoold(near(c)) 
    U % X(c)  = UXold(near(c)) 
    U % Xo(c) = UXoold(near(c)) 
    V % n(c)  = Vold(near(c)) 
    V % o(c)  = Voold(near(c)) 
    V % C(c)  = VCold(near(c)) 
    V % Co(c) = VCoold(near(c)) 
    V % Do(c) = VDoold(near(c)) 
    V % X(c)  = VXold(near(c)) 
    V % Xo(c) = VXoold(near(c)) 
    W % n(c)  = Wold(near(c)) 
    W % o(c)  = Woold(near(c)) 
    W % C(c)  = WCold(near(c)) 
    W % Co(c) = WCoold(near(c)) 
    W % Do(c) = WDoold(near(c)) 
    W % X(c)  = WXold(near(c)) 
    W % Xo(c) = WXoold(near(c)) 
    P % n(c)  = Pold(near(c)) 
    PP % n(c) = PPold(near(c)) 
    Px(c)     = Pxold(near(c)) 
    Py(c)     = Pyold(near(c)) 
    Pz(c)     = Pzold(near(c)) 
  end do
  if(HOTini==YES) then
    do c = 1, NC  
      T % n(c)  = Told(near(c)) 
      T % o(c)  = Toold(near(c)) 
    end do
  end if
  if(SIMULA==ZETA.or.SIMULA==K_EPS_VV) then
    do c = 1, NC  
      Kin % n(c)  = Kold(near(c)) 
      Kin % o(c)  = Koold(near(c)) 
      Kin % C(c)  = KCold(near(c)) 
      Kin % Co(c) = KCoold(near(c)) 
      Kin % Do(c) = KDoold(near(c)) 
      Kin % X(c)  = KXold(near(c)) 
      Kin % Xo(c) = KXoold(near(c)) 

      Eps % n(c)  = Eold(near(c)) 
      Eps % o(c)  = Eoold(near(c)) 
      Eps % C(c)  = ECold(near(c)) 
      Eps % Co(c) = ECoold(near(c)) 
      Eps % Do(c) = EDoold(near(c)) 
      Eps % X(c)  = EXold(near(c)) 
      Eps % Xo(c) = EXoold(near(c)) 

      v_2 % n(c)  = v_2old(near(c)) 
      v_2 % o(c)  = v_2oold(near(c)) 
      v_2 % C(c)  = v_2Cold(near(c)) 
      v_2 % Co(c) = v_2Coold(near(c)) 
      v_2 % Do(c) = v_2Doold(near(c)) 
      v_2 % X(c)  = v_2Xold(near(c)) 
      v_2 % Xo(c) = v_2Xoold(near(c)) 

      F22 % n(c)  = F22old(near(c)) 
      F22 % o(c)  = F22oold(near(c)) 
      F22 % Do(c) = F22Doold(near(c)) 
      F22 % X(c)  = F22Xold(near(c)) 
      F22 % Xo(c) = F22Xoold(near(c)) 
    end do
  end if
  if(SIMULA==K_EPS) then
    do c = 1, NC  
      Kin % n(c)  = Kold(near(c)) 
      Kin % o(c)  = Koold(near(c)) 
      Kin % C(c)  = KCold(near(c)) 
      Kin % Co(c) = KCoold(near(c)) 
      Kin % Do(c) = KDoold(near(c)) 
      Kin % X(c)  = KXold(near(c)) 
      Kin % Xo(c) = KXoold(near(c)) 

      Eps % n(c)  = Eold(near(c)) 
      Eps % o(c)  = Eoold(near(c)) 
      Eps % C(c)  = ECold(near(c)) 
      Eps % Co(c) = ECoold(near(c)) 
      Eps % Do(c) = EDoold(near(c)) 
      Eps % X(c)  = EXold(near(c)) 
      Eps % Xo(c) = EXoold(near(c)) 
    end do
  end if
  do s=1, NS
    c1=SideC(1,s)
    c2=SideC(2,s)

!---- interpoliraj gustocu i brzine
    Us = f(s) * U % n(c1) + (1.0-f(s)) * U % n(c2)
    Vs = f(s) * V % n(c1) + (1.0-f(s)) * V % n(c2)
    Ws = f(s) * W % n(c1) + (1.0-f(s)) * W % n(c2)
    Flux(s) = ( Us*Sx(s) + Vs*Sy(s) + Ws*Sz(s) )
  end do 

  name = answer

  deallocate(Xold)
  deallocate(Yold)
  deallocate(Zold)
  deallocate(Uold)
  deallocate(Vold)
  deallocate(Wold)
  deallocate(Kold)
  deallocate(Eold)
  deallocate(v_2old)
  deallocate(F22old)
  deallocate(Uoold)
  deallocate(Voold)
  deallocate(Woold)
  deallocate(Koold)
  deallocate(Eoold)
  deallocate(v_2oold)
  deallocate(F22oold)
  deallocate(UDoold)
  deallocate(VDoold)
  deallocate(WDoold)
  deallocate(KDoold)
  deallocate(EDoold)
  deallocate(v_2Doold)
  deallocate(F22Doold)
  deallocate(UCold)
  deallocate(VCold)
  deallocate(WCold)
  deallocate(KCold)
  deallocate(ECold)
  deallocate(v_2Cold)
  deallocate(UCoold)
  deallocate(VCoold)
  deallocate(WCoold)
  deallocate(KCoold)
  deallocate(ECoold)
  deallocate(v_2Coold)
  deallocate(UXold)
  deallocate(VXold)
  deallocate(WXold)
  deallocate(KXold)
  deallocate(EXold)
  deallocate(v_2Xold)
  deallocate(F22Xold)
  deallocate(UXoold)
  deallocate(VXoold)
  deallocate(WXoold)
  deallocate(KXoold)
  deallocate(EXoold)
  deallocate(v_2Xoold)
  deallocate(F22Xoold)
  deallocate(Pold)
  deallocate(PPold)
  deallocate(Pxold)
  deallocate(Pyold)
  deallocate(Pzold)
  if(HOTini==YES) then
    deallocate(Told)
    deallocate(Toold)
    deallocate(TDoold)
    deallocate(TCold)
    deallocate(TCoold)
    deallocate(TXold)
    deallocate(TXoold)
  end if


  write(*,*) 'Finished with LoaIni  Processor: ', this

  END SUBROUTINE LoaIni

