!======================================================================!
  SUBROUTINE SavIni()
!----------------------------------------------------------------------!
!   
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER   ::  c, Nvar, var, c1, c2, s, SC, nn
  TYPE(Unknown) :: PHI
  CHARACTER :: namOut*80, answer*80, ext*4
!--------------------------------[CVS]---------------------------------!
  character*80 rcs1,rcs2
  data rcs1/                                                        &
  '$Id: SavIni.f90,v 1.4 2008/12/10 14:57:50 IUS\mhadziabdic Exp $'/
  data rcs2/                                                        &
  '$Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/SavIni.f90,v $'/
!======================================================================!

!---- store the name
   if(this  < 2)                                                     &
   write(*,*) '# Now saving initial files [skip cancels]:'  
   call ReadC(7,inp,tn,ts,te)
  
   read(inp(ts(1):te(1)), '(A80)')  namOut
   answer=namOut
   call ToUppr(answer)
  
   if(answer == 'SKIP') return
  
!---- save the name
   answer = name
   name = namOut
   nn = 0
! For axicyllindrical coord. only 10 degree segment is used. 
  if(JET == YES.or.PIPE == YES) then
    ext = '.xyz'
    call NamFil(THIS, namOut, ext, len_trim(ext))
    open(9,FILE=namOut)
    do c= 1, NC
      if(yc(c)/xc(c) < 0.1763.and. yc(c)/xc(c) > 0.0) then
        nn = nn + 1
      end if
    end do    ! through centers 
    write(9,'(I10)') nn
    do c= 1, NC
      if(yc(c)/xc(c) < 0.1763.and. yc(c)/xc(c) > 0.0) then
        nn = nn + 1
        write(9,'(3E25.8)') xc(c),yc(c),zc(c)
      end if
    end do    ! through centers 
    close(9)

    ext = '.U__'
    call NamFil(THIS, namOut, ext, len_trim(ext))
    open(9,FILE=namOut)
    do c= 1, NC
      if(yc(c)/xc(c) < 0.1763.and. yc(c)/xc(c) > 0.0) then
        write(9,'(7E18.8)') U % n(c), U % o(c), U % C(c), U % Co(c),  &
                            U % Do(c), U % X(c), U % Xo(c)
      end if
    end do    ! through centers 
    close(9)

    ext = '.V__'
    call NamFil(THIS, namOut, ext, len_trim(ext))
    open(9,FILE=namOut)
    do c= 1, NC
      if(yc(c)/xc(c) < 0.1763.and. yc(c)/xc(c) > 0.0) then
        write(9,'(7E18.8)') V % n(c), V % o(c), V % C(c), V % Co(c),  &
                            V % Do(c), V % X(c), V % Xo(c)
      end if
    end do    ! through centers 
    close(9)

    ext = '.W__'
    call NamFil(THIS, namOut, ext, len_trim(ext))
    open(9,FILE=namOut)
    do c= 1, NC
      if(yc(c)/xc(c) < 0.1763.and. yc(c)/xc(c) > 0.0) then
        write(9,'(7E18.8)') W % n(c), W % o(c), W % C(c), W % Co(c),  &
                            W % Do(c), W % X(c), W % Xo(c)
      end if
    end do    ! through centers 
    close(9)

    ext = '.P__'
    call NamFil(THIS, namOut, ext, len_trim(ext))
    open(9,FILE=namOut)
    do c= 1, NC
      if(yc(c)/xc(c) < 0.1763.and. yc(c)/xc(c) > 0.0) then
        write(9,'(5E18.8)') P % n(c), PP % n(c), Px(c), Py(c),  &
                            Pz(c)
      end if
    end do    ! through centers 
    close(9)
 
    if(HOT == YES) then 
      ext = '.T__'
      call NamFil(THIS, namOut, ext, len_trim(ext))
      open(9,FILE=namOut)
      do c= 1, NC
        if(yc(c)/xc(c) < 0.1763.and. yc(c)/xc(c) > 0.0) then
          write(9,'(7E18.8)') T % n(c), T % o(c), T % C(c), T % Co(c),  &
                            T % Do(c), T % X(c), T % Xo(c)
        end if
      end do    ! through centers 
      close(9)
    end if

    if(SIMULA == ZETA.or.SIMULA == K_EPS_VV) then
      ext = '.Kin'
      call NamFil(THIS, namOut, ext, len_trim(ext))
      open(9,FILE=namOut)
      do c= 1, NC
        if(yc(c)/xc(c) < 0.1763.and. yc(c)/xc(c) > 0.0) then
          write(9,'(7E18.8)') Kin % n(c), Kin % o(c), Kin % C(c), Kin % Co(c),  &
                              Kin % Do(c), Kin % X(c), Kin % Xo(c)
        end if
      end do    ! through centers 
      close(9)

      ext = '.Eps'
      call NamFil(THIS, namOut, ext, len_trim(ext))
      open(9,FILE=namOut)
      do c= 1, NC
        if(yc(c)/xc(c) < 0.1763.and. yc(c)/xc(c) > 0.0) then
          write(9,'(7E18.8)') Eps % n(c), Eps % o(c), Eps % C(c), Eps % Co(c),  &
                              Eps % Do(c), Eps % X(c), Eps % Xo(c)
        end if
      end do    ! through centers 
      close(9)

      ext = '.v_2'
      call NamFil(THIS, namOut, ext, len_trim(ext))
      open(9,FILE=namOut)
      do c= 1, NC
        if(yc(c)/xc(c) < 0.1763.and. yc(c)/xc(c) > 0.0) then
          write(9,'(7E18.8)') v_2 % n(c), v_2 % o(c), v_2 % C(c), v_2 % Co(c),  &
                              v_2 % Do(c), v_2 % X(c), v_2 % Xo(c)
        end if
      end do    ! through centers 
      close(9)

      ext = '.f22'
      call NamFil(THIS, namOut, ext, len_trim(ext))
      open(9,FILE=namOut)
      do c= 1, NC
        if(yc(c)/xc(c) < 0.1763.and. yc(c)/xc(c) > 0.0) then
          write(9,'(7E18.8)') f22 % n(c), f22 % o(c),  &
                              f22 % Do(c), f22 % X(c), f22 % Xo(c)
        end if
      end do    ! through centers 
      close(9)
    end if
    if(SIMULA == K_EPS) then
      ext = '.Kin'
      call NamFil(THIS, namOut, ext, len_trim(ext))
      open(9,FILE=namOut)
      do c= 1, NC
        if(yc(c)/xc(c) < 0.1763.and. yc(c)/xc(c) > 0.0) then
          write(9,'(7E18.8)') Kin % n(c), Kin % o(c), Kin % C(c), Kin % Co(c),  &
                              Kin % Do(c), Kin % X(c), Kin % Xo(c)
        end if
      end do    ! through centers 
      close(9)

      ext = '.Eps'
      call NamFil(THIS, namOut, ext, len_trim(ext))
      open(9,FILE=namOut)
      do c= 1, NC
        if(yc(c)/xc(c) < 0.1763.and. yc(c)/xc(c) > 0.0) then
          write(9,'(7E18.8)') Eps % n(c), Eps % o(c), Eps % C(c), Eps % Co(c),  &
                              Eps % Do(c), Eps % X(c), Eps % Xo(c)
        end if
      end do
      close(9)
    end if 
  else
    ext = '.xyz'
    call NamFil(THIS, namOut, ext, len_trim(ext))
    open(9,FILE=namOut)
    do c= 1, NC
      nn = nn + 1
    end do    ! through centers 
    write(9,'(I10)') nn
    do c= 1, NC
      write(9,'(3E25.8)') xc(c),yc(c),zc(c)
    end do    ! through centers 
    close(9)

    ext = '.U__'
    call NamFil(THIS, namOut, ext, len_trim(ext))
    open(9,FILE=namOut)
    do c= 1, NC
      write(9,'(7E18.8)') U % n(c), U % o(c), U % C(c), U % Co(c),  &
                          U % Do(c), U % X(c), U % Xo(c)
    end do    ! through centers 
    close(9)

    ext = '.V__'
    call NamFil(THIS, namOut, ext, len_trim(ext))
    open(9,FILE=namOut)
    do c= 1, NC
      write(9,'(7E18.8)') V % n(c), V % o(c), V % C(c), V % Co(c),  &
                          V % Do(c), V % X(c), V % Xo(c)
    end do    ! through centers 
    close(9)

    ext = '.W__'
    call NamFil(THIS, namOut, ext, len_trim(ext))
    open(9,FILE=namOut)
    do c= 1, NC
      write(9,'(7E18.8)') W % n(c), W % o(c), W % C(c), W % Co(c),  &
                          W % Do(c), W % X(c), W % Xo(c)
    end do    ! through centers 
    close(9)

    ext = '.P__'
    call NamFil(THIS, namOut, ext, len_trim(ext))
    open(9,FILE=namOut)
    do c= 1, NC
      write(9,'(5E18.8)') P % n(c), PP % n(c), Px(c), Py(c),  &
                          Pz(c)
    end do    ! through centers 
    close(9)
 
    if(HOT == YES) then 
      ext = '.T__'
      call NamFil(THIS, namOut, ext, len_trim(ext))
      open(9,FILE=namOut)
      do c= 1, NC
        write(9,'(7E18.8)') T % n(c), T % o(c), T % C(c), T % Co(c),  &
                            T % Do(c), T % X(c), T % Xo(c)
      end do    ! through centers 
      close(9)
    end if 
 
    if(SIMULA == ZETA.or.SIMULA==K_EPS_VV) then
      ext = '.Kin'
      call NamFil(THIS, namOut, ext, len_trim(ext))
      open(9,FILE=namOut)
      do c= 1, NC
        write(9,'(7E18.8)') Kin % n(c), Kin % o(c), Kin % C(c), Kin % Co(c),  &
                            Kin % Do(c), Kin % X(c), Kin % Xo(c)
      end do    ! through centers 
      close(9)

      ext = '.Eps'
      call NamFil(THIS, namOut, ext, len_trim(ext))
      open(9,FILE=namOut)
      do c= 1, NC
        write(9,'(7E18.8)') Eps % n(c), Eps % o(c), Eps % C(c), Eps % Co(c),  &
                            Eps % Do(c), Eps % X(c), Eps % Xo(c)
      end do    ! through centers 
      close(9)

      ext = '.v_2'
      call NamFil(THIS, namOut, ext, len_trim(ext))
      open(9,FILE=namOut)
      do c= 1, NC
        write(9,'(7E18.8)') v_2 % n(c), v_2 % o(c), v_2 % C(c), v_2 % Co(c),  &
                            v_2 % Do(c), v_2 % X(c), v_2 % Xo(c)
      end do    ! through centers 
      close(9)

      ext = '.f22'
      call NamFil(THIS, namOut, ext, len_trim(ext))
      open(9,FILE=namOut)
      do c= 1, NC
        write(9,'(7E18.8)') f22 % n(c), f22 % o(c),  &
                            f22 % Do(c), f22 % X(c), f22 % Xo(c)
      end do    ! through centers 
      close(9)
    end if
    if(SIMULA == K_EPS) then
      ext = '.Kin'
      call NamFil(THIS, namOut, ext, len_trim(ext))
      open(9,FILE=namOut)
      do c= 1, NC
        write(9,'(7E18.8)') Kin % n(c), Kin % o(c), Kin % C(c), Kin % Co(c),  &
                            Kin % Do(c), Kin % X(c), Kin % Xo(c)
      end do    ! through centers 
      close(9)

      ext = '.Eps'
      call NamFil(THIS, namOut, ext, len_trim(ext))
      open(9,FILE=namOut)
      do c= 1, NC
        write(9,'(7E18.8)') Eps % n(c), Eps % o(c), Eps % C(c), Eps % Co(c),  &
                            Eps % Do(c), Eps % X(c), Eps % Xo(c)
      end do    ! through centers 
      close(9)
    end if
  end if

  name = answer 

  END SUBROUTINE SavIni
