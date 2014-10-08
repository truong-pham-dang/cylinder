!======================================================================!
  SUBROUTINE ProSav(namAut)
!----------------------------------------------------------------------!
! Writes: NAME.r.gmv                                                   !
! ~~~~~~~                                                              !
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER             :: c, i 
  CHARACTER           :: namOut*80, answer*80, namTem*80
  CHARACTER, OPTIONAL :: namAut*(*)
!--------------------------------[CVS]---------------------------------!
!  $Id: ProSav.f90,v 1.32 2008/12/10 14:49:32 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/ProSav.f90,v $  
!======================================================================!

!---- store the name
  namTem = name     

  if(PRESENT(namAut)) then
    write(*,*) namAut
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

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!                                 !
!     Create GMV results file     !
!                                 !
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
    call NamFil(THIS, namOut, '.r.gmv', len_trim('.r.gmv'))
    open(9, FILE=namOut)
    write(*,*) '# Now creating the file:', namOut

!----------------!
!    Geometry    !
!----------------!
!  if(this < 2) then
!  answer = namTem
!  answer(len_trim(namTem)+1 : len_trim(namTem)+4) = ".gmv"
!  open(19, FILE=answer)
!
!--rewriting NAME.gmv
!  do i=1,10*NC
!    call ReadC(19,inp,tn,ts,te)
!    read(inp(ts(1):te(1)),'(A80)') answer
!    call ToUppr(answer)
!    if(answer == 'ENDGMV') then
!      GOTO 5
!    else
!      if(answer == 'GMVINPUT') then
!        write(9,'(A14)') 'gmvinput ascii'
!      else
!        write(9,*) inp(1 : len_trim(inp))
!      end if
!    end if
!  end do

!--no end of NAME.gmv found
!  if(THIS < 2) write(*,*) '@ProSav: not found end of file', answer
!  if(THIS < 2) write(*,*) '         Canceled making ', namOut
!  RETURN

!--closing NAME.gmv
!5 CONTINUE
!  close(19)
!  end if

!----- velocity
     write(9, *) 'velocity 0'  !  0 is for cell data
     write(9,'(1F16.6)') (U % n(c),c=1,NC)
     write(9,'(1F16.6)') (V % n(c),c=1,NC)
     write(9,'(1F16.6)') (W % n(c),c=1,NC)

!----- other variables
     write(9, *) 'variable' 

!----- pressure
     write(9, *) 'P 0'
     write(9,'(1F16.6)') (P % n(c),c=1,NC)
!----- pressure corrections
     write(9, *) 'PP 0'
     write(9,'(1F16.6)') (PP % n(c),c=1,NC)

!----- temperature
     if(HOT==YES) then
       write(9, *) 'T 0'
       write(9,'(1F16.6)') (T % n(c),c=1,NC)
     end if

!=============!
!    K_EPS    !
!=============!
     if(SIMULA == K_EPS.or.SIMULA == HYB_PITM) then 
       write(9, *) 'k 0'
       write(9,'(1F16.6)') (Kin % n(c),c=1,NC)
       write(9, *) 'eps 0'
       write(9,'(1F16.6)') (Eps % n(c),c=1,NC)
       write(9, *) 'VISt 0'
       write(9,'(1F16.6)') (VISt(c),c=1,NC)
       write(9, *) 'Pk 0'
       write(9,'(1F16.6)') (Pk(c),c=1,NC)
       write(9, *) 'Yplus 0'
       write(9,'(1F16.6)') (Ynd(c),c=1,NC)
       write(9, *) 'Uf 0'
       write(9,'(1F16.6)') (Uf(c),c=1,NC)
       write(9, *) 'uv 0'
       write(9,'(1F16.6)') (VISt(c)*Uz(c),c=1,NC)
     end if 

!================!
!    K_EPS_VV    !
!================!

     if(SIMULA == K_EPS_VV.or.SIMULA==ZETA.or.SIMULA==HYB_ZETA) then 
       write(9, *) 'k 0'
       write(9,'(1F16.6)') (Kin % n(c),c=1,NC)
       write(9, *) 'eps 0'
       write(9,'(1F16.6)') (Eps % n(c),c=1,NC)
       write(9, *) 'vv 0'
       write(9,'(1F16.6)') (v_2 % n(c),c=1,NC)
       write(9, *) 'f22 0'
       write(9,'(1F16.6)') (f22 % n(c),c=1,NC)
       write(9, *) 'Pk 0'
       write(9,'(1F16.6)') ( Pk(c),c=1,NC)
       write(9, *) 'VISt 0'
       write(9,'(1F16.6)') (VISt(c),c=1,NC)
       write(9, *) 'uw 0'
       write(9,'(1F16.6)') (-VISt(c)*Uz(c),c=1,NC)  
!       write(9, *) 'Yplus 0'
!       write(9,'(1F16.6)') (Ynd(c),c=1,NC)  
!       write(9, *) 'Utan 0'
!       write(9,'(1F16.6)') ((-U % n(c) * yc(c) / sqrt(xc(c)**2 + yc(c)**2) &
!           + V % n(c) * xc(c) / sqrt(xc(c)**2 + yc(c)**2)),c=1,NC)  
!       write(9, *) 'Urad 0'
!       write(9,'(1F16.6)') ((U % n(c) * xc(c) / sqrt(xc(c)**2 + yc(c)**2) &
!           + V % n(c) * yc(c) / sqrt(xc(c)**2 + yc(c)**2)),c=1,NC)  
     end if 

!===============!
!    SPA_ALL    !
!===============!
     if(SIMULA == SPA_ALL) then
       write(9, *) 'VIS 0'
       write(9,'(1F16.6)') (VIS % n(c),c=1,NC)
       write(9, *) 'Vort 0'
       write(9,'(1F16.6)') (Vort(c),c=1,NC)
       write(9, *) 'VISt 0'
       write(9,'(1F16.6)') (VISt(c),c=1,NC)
       write(9, *) 'WallDs 0'
       write(9,'(1F16.6)') (WallDs(c),c=1,NC)
       write(9, *) 'delta 0'
       write(9,'(1F16.6)') (delta(c),c=1,NC)
     end if

!===============!
!    DES_SPA    !
!===============!
     if(SIMULA == DES_SPA) then

!----- mean velocities 
       write(9, *) 'Umean 0'
       write(9,'(1F16.6)') (U % mean(c),c=1,NC)
       write(9, *) 'Vmean 0'
       write(9,'(1F16.6)') (V % mean(c),c=1,NC)
       write(9, *) 'Wmean 0'
       write(9,'(1F16.6)') (W % mean(c),c=1,NC)
       if(HOT == YES) then
         write(9, *) 'Tmean 0'
         write(9,'(1F16.6)') (T % mean(c),c=1,NC)
       end if

!----- velocity fluctuations
       write(9, *) 'uu 0'  !  0 is for cell data
       write(9,'(1F16.6)')(UU % mean(c)-U % mean(c)*U % mean(c),c=1,NC)
       write(9, *) 'vv 0'  !  0 is for cell data
       write(9,'(1F16.6)')(VV % mean(c)-V % mean(c)*V % mean(c),c=1,NC)
       write(9, *) 'ww 0'  !  0 is for cell data
       write(9,'(1F16.6)')(WW % mean(c)-W % mean(c)*W % mean(c),c=1,NC)
       write(9, *) 'uv 0'  !  0 is for cell data
       write(9,'(1F16.6)')(UV % mean(c)-U % mean(c)*V % mean(c),c=1,NC)
       write(9, *) 'uw 0'  !  0 is for cell data
       write(9,'(1F16.6)')(UW % mean(c)-U % mean(c)*W % mean(c),c=1,NC)
       write(9, *) 'vw 0'  !  0 is for cell data
       write(9,'(1F16.6)')(VW % mean(c)-V % mean(c)*W % mean(c),c=1,NC)

       write(9, *) 'VIS 0'
       write(9,'(1F16.6)') (VIS % n(c),c=1,NC)
       write(9, *) 'Vort 0'
       write(9,'(1F16.6)') (Vort(c),c=1,NC)
       write(9, *) 'VISt 0'
       write(9,'(1F16.6)') (VISt(c),c=1,NC)
       write(9, *) 'WallDs 0'
       write(9,'(1F16.6)') (WallDs(c),c=1,NC)
       write(9, *) 'delta 0'
       write(9,'(1F16.6)') (delta(c),c=1,NC)
     end if

!===========!
!    LES    !
!===========!
     if(SIMULA == LES) then 
!----- mean velocities 
       write(9, *) 'Umean 0'
       write(9,'(1F16.6)') (U % mean(c),c=1,NC)
       write(9, *) 'Vmean 0'
       write(9,'(1F16.6)') (V % mean(c),c=1,NC)
       write(9, *) 'Wmean 0'
       write(9,'(1F16.6)') (W % mean(c),c=1,NC)
       if(HOT == YES) then
         write(9, *) 'Tmean 0'
         write(9,'(1F16.6)') (T % mean(c),c=1,NC)
       end if
!----- velocity fluctuations
       write(9, *) 'uu 0'  !  0 is for cell data
       write(9,'(1F16.6)')(UU % mean(c)-U % mean(c)*U % mean(c),c=1,NC)
       write(9, *) 'vv 0'  !  0 is for cell data
       write(9,'(1F16.6)')(VV % mean(c)-V % mean(c)*V % mean(c),c=1,NC)
       write(9, *) 'ww 0'  !  0 is for cell data
       write(9,'(1F16.6)')(WW % mean(c)-W % mean(c)*W % mean(c),c=1,NC)
       write(9, *) 'uv 0'  !  0 is for cell data
       write(9,'(1F16.6)')(UV % mean(c)-U % mean(c)*V % mean(c),c=1,NC)
       write(9, *) 'uw 0'  !  0 is for cell data
       write(9,'(1F16.6)')(UW % mean(c)-U % mean(c)*W % mean(c),c=1,NC)
       write(9, *) 'vw 0'  !  0 is for cell data
       write(9,'(1F16.6)')(VW % mean(c)-V % mean(c)*W % mean(c),c=1,NC)

!----- turbulent viscosity
       write(9, *) 'muT 0'
       write(9,'(1F16.6)') (VISt(c),c=1,NC)

!----- turbulent viscosity
       write(9, *) 'ShearMean 0'
       write(9,'(1F16.6)') (ShearMean(c),c=1,NC)

!----- wall distance            
       write(9, *) 'wall 0'
       write(9,'(1F16.6)') (WallDs(c),c=1,NC)
     end if

     write(9, *) 'endvars' 

     write(9, *) 'endgmv'

!----- this is important for parallel version
     write(9,'(I8)') NC    

     close(9)

!---- restore the name
   name = namTem

   END SUBROUTINE ProSav
