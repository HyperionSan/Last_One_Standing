module null_eth

   implicit none

   double precision, dimension (:,:), pointer, save, private :: pp
   double complex,   dimension (:,:), pointer, save, private :: p_x, p_y
   double complex,   dimension (:,:), pointer, save, private :: p_xy, p_xx, p_yy

contains

subroutine null_eth_allocate

   use null_grid, only : nn

   allocate (pp(nn,nn), p_x(nn,nn), p_y(nn,nn))
   allocate (p_xy(nn,nn), p_xx(nn,nn), p_yy(nn,nn))

end subroutine null_eth_allocate

subroutine null_d1 (out, in, spin, e1)

   use null_grid, only : ii, nn, dd, qs, ps

   implicit none

   double complex, dimension (:,:), intent (out) :: out
   double complex, dimension (:,:), intent (in)  :: in
   integer,                         intent (in)  :: spin, e1

   double precision,    dimension (:,:), pointer     :: x, y

   x => qs
   y => ps

   pp = 0.5 * (1. + x * x + y * y)

   p_x(1,:) = (-3. * in(1,:) + 4. * in(2,:) - in(3,:)) / (2. * dd)
   p_x(2:nn-1,:) = (in(3:nn,:) - in(1:nn-2,:)) / (2. * dd)
   p_x(nn,:) = (3. * in(nn,:) - 4. * in(nn-1,:) + in(nn-2,:)) / (2. * dd)

   p_y(:,1) = (-3. * in(:,1) + 4. * in(:,2) - in(:,3)) / (2. * dd)
   p_y(:,2:nn-1) = (in(:,3:nn) - in(:,1:nn-2)) / (2. * dd)
   p_y(:,nn) = (3. * in(:,nn) - 4. * in(:,nn-1) + in(:,nn-2)) / (2. * dd)

   out = pp * (p_x + ii * e1 * p_y) + spin * (e1 * x + ii * y) * in

end subroutine null_d1

subroutine null_d2 (out, in, spin, e1, e2)

   use null_grid, only : ii, nn, dd, qs, ps

   implicit none

   double complex, dimension (:,:), intent (out) :: out
   double complex, dimension (:,:), intent (in)  :: in
   integer,                         intent (in)  :: spin, e1, e2

   double precision :: h1, h2, c1, c2, c3, c4, c5, c6
   double precision,    dimension (:,:), pointer     :: x, y

   x => qs
   y => ps

   h1 = 1. / (2. * dd)
   h2 = 1. / (dd * dd)

   c1 = 1 + e1 * e2 + spin * (e1 + e2)
   c2 = 2 * spin + e1 + e2
   c3 = 2 * spin * e1 * e2 + e1 + e2
   c4 = e1 * e2
   c5 = e1 + e2
   c6 = e2 - e1

   pp = 0.5 * (1. + x * x + y * y)

   p_xx(1,:)  = (0., 0.)
   p_xx(nn,:) = (0., 0.)
   p_xx(:,1)  = (0., 0.)
   p_xx(:,nn) = (0., 0.)

   p_xy(1,:)  = (0., 0.)
   p_xy(nn,:) = (0., 0.)
   p_xy(:,1)  = (0., 0.)
   p_xy(:,nn) = (0., 0.)

   p_yy(1,:)  = (0., 0.)
   p_yy(nn,:) = (0., 0.)
   p_yy(:,1)  = (0., 0.)
   p_yy(:,nn) = (0., 0.)

   p_x(1,:) = (-3. * in(1,:) + 4. * in(2,:) - in(3,:)) * h1
   p_x(2:nn-1,:) = (in(3:nn,:) - in(1:nn-2,:)) * h1
   p_x(nn,:) = (3. * in(nn,:) - 4. * in(nn-1,:) + in(nn-2,:)) * h1

   p_y(:,1) = (-3. * in(:,1) + 4. * in(:,2) - in(:,3)) * h1
   p_y(:,2:nn-1) = (in(:,3:nn) - in(:,1:nn-2)) * h1
   p_y(:,nn) = (3. * in(:,nn) - 4. * in(:,nn-1) + in(:,nn-2)) * h1

   p_xx(2:nn-1,:) = (in(3:nn,:)-2.*in(2:nn-1,:)+in(1:nn-2,:)) * h2

   p_yy(:,2:nn-1) = (in(:,3:nn)-2.*in(:,2:nn-1)+in(:,1:nn-2)) * h2

   p_xy(2:nn-1,2:nn-1) = (in(3:nn,3:nn) - in(3:nn,1:nn-2) &
        - in(1:nn-2,3:nn) + in(1:nn-2,1:nn-2)) / (4. * dd * dd)

   p_xy(1,2:nn-1) = (-3. * in(1,3:nn) + 4. * in(2,3:nn) - in(3,3:nn) &
        + 3. * in(1,1:nn-2) - 4. * in(2,1:nn-2) + in(3,1:nn-2)) * h1 * h1

   out = pp * pp * (p_xx - c4 * p_yy + ii * c5 * p_xy) &
        + pp * ((c1 * x + ii * c2 * y) * p_x &
        + (c2 * x + ii * c1 * y) * ii * c4 * p_y) &
        + spin * in * 0.5 * (c6 + c3 * (x * x - c4 * y * y) &
        + 2 * c1 * ii * x * y)

end subroutine null_d2

end module null_eth
