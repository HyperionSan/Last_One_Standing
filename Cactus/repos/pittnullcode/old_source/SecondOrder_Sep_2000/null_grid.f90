module null_grid

   use null_params

   implicit none

   double precision, save :: dd, dx, dt
   integer,          save :: it = 0

   double complex,   parameter :: ii = (0., 1.)
   double precision, parameter :: pi = 3.1415926535897932385

   double precision, dimension (:),   pointer, save :: x,  xh, rb, rbh
   double precision, dimension (:,:), pointer, save :: qs, ps, pp
   double complex,   dimension(:,:),  pointer, save :: z, zb

contains

  subroutine null_grid_allocate

    allocate (x(nx), xh(nx), rb(nx), rbh(nx), qs(nn,nn), ps(nn,nn), pp(nn,nn))
    allocate (z(nn,nn),zb(nn,nn))

end subroutine null_grid_allocate

subroutine null_setup_grid

   use null_params

   implicit none

   integer :: i

   dd = 2.0 / dble(nn - 5)

   do i = 1, nn
      qs(i,:) = -1.0 + (i-3) * dd
      ps(:,i) = -1.0 + (i-3) * dd
   end do

   pp = 1. + qs ** 2 + ps ** 2

   z  = qs + ii * ps
   zb = qs - ii * ps

   dx = 0.5 / dble(nx - 1)

   do i = 1, nx
      x(i)  = 0.5 + (i - 1  ) * dx
      xh(i) = 0.5 + (i - 0.5) * dx
   end do

   do i = 1, nx - 1
      rb(i)  = rwt * x(i) / (1. - x(i))
   end do
   rb(nx) = rb(nx-1)

   do i = 1, nx
      rbh(i)  = rwt * xh(i) / (1. - xh(i))
   end do

   dt = cfl * min(dx, dd)


end subroutine null_setup_grid
end module null_grid
