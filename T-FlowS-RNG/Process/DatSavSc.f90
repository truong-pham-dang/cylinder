!======================================================================!
  SUBROUTINE DatSavSc(namSc, idSc, PHI)
!----------------------------------------------------------------------!
! Writes: NAME.dat                                                     !
! ~~~~~~~                                                              !
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER   :: idSc
  CHARACTER :: namSc*(*)
  REAL      :: PHI(-NbC:NC)
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: N, s, c, c1, c2, Nfac(10), NtotFac
!--------------------------------[CVS]---------------------------------!
!  $Id: DatSavSc.f90,v 1.4 2008/12/10 14:34:23 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/DatSavSc.f90,v $             
!======================================================================!
!----------------!
!     inside     !
!----------------!
  write(9,'(A4,A,A2)') '(0 "', namSc, '")'
  write(9,'(A6,I3,A9,I8,2X,I9,A2)') '(300 (', idSc, ' 1 1 0 0 ',  1, NC, ')(' 
  do c=1,NC
    write(9,'(F14.6)') PHI(c)
  end do  
  write(9,'(A2)') '))'

!-------------------------!
!     on the boundary     !
!-------------------------!
  NtotFac = 0
  do n=1,10   ! Browse through boundary condition types
    Nfac(n) = 0
    do s=1,NS   ! Count the faces with boundary condition "n"
      c2 = SideC(2,s)
      if(c2 < 0) then
        if(BCmark(c2) == n) Nfac(n)=Nfac(n)+1
      end if
    end do    ! sides

    if(Nfac(n) > 0) then
      write(9,'(A4,A,A17,I3,A3)') '(0 "', namSc, ' on the boundary ', n, ' ")'
      write(9,'(A6,I3,I4,A7,I7,I7,A4)') '(300 (', idSc, 100+n, ' 1 0 0 ', &
                                        NtotFac+1, NtotFac+Nfac(n), ')('
      do s=1,NS !@@@ +NSsh
        c1 = SideC(1,s)
        c2 = SideC(2,s)
        if(c2 < 0) then
          if(BCmark(c2) == n) then
            if(TypeBC(c2) == SYMMETRY) then
              write(9,'(F14.6)') PHI(c1)
            else
              write(9,'(F14.6)') PHI(c2)
            end if
          end if
        end if
      end do
      write(9,'(A2)') '))'
    end if

!---- prepare for next boundary
    NtotFac = NtotFac+Nfac(n)

  end do   ! n -> boundary condition types

  END SUBROUTINE DatSavSc
