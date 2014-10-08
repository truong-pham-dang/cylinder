!======================================================================!
  SUBROUTINE Inertia(sub) 
!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod 
  USE div_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: sub                           ! subdomain 
!-------------------------------[Locals]-------------------------------!
  INTEGER :: i, NCloc
  REAL    :: xm, ym, zm
  REAL    :: Im(3,3), d(3), v(3,3), d_max(3)    
!--------------------------------[CVS]---------------------------------!
!  $Id: Inertia.f90,v 1.2 2002/10/30 16:29:18 niceno Exp $    
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Divide/Inertia.f90,v $     
!======================================================================!

  xm=0.0
  ym=0.0
  zm=0.0
  NCloc=0
  do i=1,NC
    if(proces(i)==sub) then
      xm = xm + xc(i)
      ym = ym + yc(i)
      zm = zm + zc(i)
      NCloc=NCloc+1 
    end if
  end do 
  xm=xm/NCloc
  ym=ym/NCloc
  zm=zm/NCloc

  write(*,*) 'Center of mass for subdomain ', sub, ' is: ', xm, ym, zm

  Im = 0.
  do i=1,NC
    if(proces(i)==sub) then
      Im(1,1) = Im(1,1) + (yc(i)-ym)**2 + (zc(i)-zm)**2
      Im(2,2) = Im(2,2) + (xc(i)-xm)**2 + (zc(i)-zm)**2
      Im(3,3) = Im(3,3) + (xc(i)-xm)**2 + (yc(i)-ym)**2

      Im(1,2) = Im(1,2) - (xc(i)-xm)*(yc(i)-ym) 
      Im(1,3) = Im(1,3) - (xc(i)-xm)*(zc(i)-zm) 
      Im(2,3) = Im(2,3) - (yc(i)-ym)*(zc(i)-zm) 
    end if
  end do 
  Im(2,1) = Im(1,2)
  Im(3,1) = Im(1,3)
  Im(3,2) = Im(2,3)

  call Jacobi(Im, 3, 3, d, v, i)

  write(*,*) 'd=',d(1), d(2), d(3)

  write(*,*) 'v=', (v(1,i), i=1,3)
  write(*,*) '  ', (v(2,i), i=1,3)
  write(*,*) '  ', (v(3,i), i=1,3)

  if(min(d(1),d(2),d(3)) == d(1)) then
    d_max(1)=v(1,1)
    d_max(2)=v(2,1)
    d_max(3)=v(3,1)
  else if(min(d(1),d(2),d(3)) == d(2)) then
    d_max(1)=v(1,2)
    d_max(2)=v(2,2)
    d_max(3)=v(3,2)
  else if(min(d(1),d(2),d(3)) == d(3)) then
    d_max(1)=v(1,3)
    d_max(2)=v(2,3)
    d_max(3)=v(3,3)
  end if 

  write(*,*) 'Sorting the cells'
  do i=1,NC
    iin(i) = i
    criter(i) = xc(i)*d_max(1) + yc(i)*d_max(2) + zc(i)*d_max(3)
  end do
  call RISort(criter(1),iin(1),NC,2)

  return

  END SUBROUTINE Inertia
