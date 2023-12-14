module null_grid

   use null_params

   implicit none

   double precision, save :: dd, dx, dt, qmin, qmax
   integer,          save :: it = 0

   double complex,   parameter :: ii = (0., 1.)
   double precision, parameter :: pi = 3.1415926535897932385d0

   double precision, dimension (:),   pointer, save :: x,  xh, rb, rbh
   double precision, dimension (:,:), pointer, save :: qs, ps, pp
   double precision, dimension (:,:), pointer, save :: qsh, psh, pph
   double complex,   dimension(:,:),  pointer, save :: z, zb
   double precision, save :: mloss(0:3) = (/0., 0., 0., 0./)
contains

  subroutine null_grid_allocate

    allocate (x(nx), xh(nx), rb(nx), rbh(nx))
    allocate (qs(nn,nn), ps(nn,nn), pp(nn,nn))
    allocate (qsh(nn,nn), psh(nn,nn), pph(nn,nn))
    allocate (z(nn,nn),zb(nn,nn))

end subroutine null_grid_allocate

subroutine null_setup_grid

   use null_params
   implicit none

   integer :: i

   if (qsize < 1.0d0) then
      write (*,*) 'qsize MUST BE >= 1.0 --- qsize = ', qsize
      stop
   end if

   dd = (2.0 * qsize) / dble(nn - 5)
   qmin = -qsize - 2. * dd
   qmax =  qsize + 2. * dd

   do i = 1, nn
      qs(i,:) = -qsize + (i-3) * dd
      ps(:,i) = -qsize + (i-3) * dd
   end do

   do i = 1, nn
      qsh(i,:) = -qsize + (i-3) * dd
      psh(:,i) = -qsize + (i-3) * dd
   end do

   pp = 1. + qs ** 2 + ps ** 2
   pph = 1. + qsh ** 2 + psh ** 2

   z  = qs + ii * ps
   zb = qs - ii * ps

   dx = 0.5 / dble(nx - 1 - NIntPts)

   do i = 1, nx
      x(i)  = 0.5 + (i - 1 - NIntPts ) * dx
      xh(i) = 0.5 + (i - 0.5 - NIntPts ) * dx
   end do

   do i = 1, nx - 1
      rb(i)  = rwt * x(i) / (1. - x(i))
   end do
   rb(nx) = rb(nx-1)

   do i = 1, nx
      rbh(i)  = rwt * xh(i) / (1. - xh(i))
   end do

   open (unit = 21, file = 'rgrid.dat', status = 'unknown')
   write (unit = 21, fmt = '(a)') '# I    x(I)    rb(I)'
   do i = 1, nx-1
      write (unit = 21, fmt = '(i4,2(1x,e16.5))') i, x(i), rb(i)
   end do
   close (unit = 21)

   ! Note: the CFL condition here is the requirement
   ! dt <= 2 dr AND dt <= r dd, where r=rwt*x/(1-x), and
   ! we take x=1/2.

   !dt = cfl * rwt * min(8. * dx, dd)

end subroutine null_setup_grid

end module null_grid
