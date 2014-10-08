PROGRAM WallDemo

IMPLICIT NONE

INTEGER i
REAL    u_pow, u_log, D, z
REAL    A_pow, B_pow,  &
        A_log, B_log,  &
        A_dam

A_pow = 8.3
B_pow = 1.0/7.0

A_log = 1.0/0.41
B_log = 5.2     

A_dam = 25.0

do i=0,1000
  z = 1.0 * i

  D = 1.0 - exp( -z / A_dam )

  if(z <= 11.81) then
    u_pow = z
    u_log = z
  else
    u_pow = A_pow * z**B_pow
    u_log = A_log * log(z) + 5.5
  end if  

  write(*,'(4F12.6)') z, u_pow, u_log, D

end do

END PROGRAM WallDemo
