
subroutine tcross(n)
! code time at which the first caustic appears

use affine
use numservicef90, only: vintegrate
implicit none

integer,       intent(in) :: n
real(kind=wp)             :: tX, uhat_minus, uhat_cross, uhat, t, du_hat 
real(kind=wp), dimension(:), allocatable :: dt_duhat_vect, t_vect
integer                   :: i

allocate(dt_duhat_vect(n), t_vect(n))

uhat_cross = -sigma(i_equator)*0.5_wp
uhat_minus = th_minus(i_equator) + u_zero(i_equator) - t_zero(i_equator)
du_hat     = (uhat_cross - uhat_minus)/n

open (unit = 89, file = 'estimateX.dat', status = 'unknown')
t = t_minus
do i = 1, n-1
   uhat             = uhat_minus + i*du_hat
   t                = t + uprime_equ(uhat)*du_hat
   dt_duhat_vect(i) = uprime_equ(uhat)
   write(89,99) t,  uhat
end do
tX = t
close (unit = 89)

call vintegrate(dt_duhat_vect, t_vect, uhat_minus, uhat_cross, 0)  

tX = t_vect(n) + t_minus ! t_minus is the integration constant

write(10, *) 'estimated crossing time          =', tX 
write(*,  *) 'estimated crossing time          =', tX 
write(*,  *) 'tX(2nd order)/tX(1st order) - 1  =', tX/t - 1.0
write(*,  *)

deallocate(dt_duhat_vect, t_vect)
return
99 format(2E13.5)
end subroutine tcross
