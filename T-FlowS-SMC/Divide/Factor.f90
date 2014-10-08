!======================================================================!
  SUBROUTINE Factor(i, number, factrs)
!----------------------------------------------------------------------!
!   Splits the number into it's primal factors. ( 30 = 5 x 3 x 2 )     !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: number, factrs(128)
!-------------------------------[Locals]-------------------------------!
  INTEGER :: i, j, primes(15)
!--------------------------------[CVS]---------------------------------!
!  $Id: Factor.f90,v 1.6 2002/10/30 16:29:18 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Divide/Factor.f90,v $   
!======================================================================!
  DATA primes /47,43,41,37,31,29,23,19,17,13,11, 7, 5, 3, 2/   
!----------------------------------------------------------------------!
!   For example if factor is called with number = 30 it will split it  !
!   into 5,3 and 2, and make the following output:                     !
!   i=3, factrs(1)=5, factrs(2)=3, factrs(3)=2                         !
!   See also: Divide                                                   !
!----------------------------------------------------------------------!

  i=0
  do j=1,15

1   if( mod(number,primes(j)) == 0 ) then
      number = number/primes(j)
      i = i + 1
      factrs(i) = primes(j)
      goto 1
    endif 

  end do

  END SUBROUTINE Factor
