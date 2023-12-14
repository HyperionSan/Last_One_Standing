module null_eth

contains

subroutine null_d1 (out, in, spin, e1)

   use null_grid, only : nn, ii, dd, x=>qs, y=>ps, pp

   implicit none

   double complex,   dimension (:,:), intent (out) :: out
   double complex,   dimension (:,:), intent (in)  :: in
   integer,                           intent (in)  :: spin, e1

   double complex, dimension (nn,nn) :: p_x, p_y

   p_x(1,:) = (-3. * in(1,:) + 4. * in(2,:) - in(3,:)) / (2. * dd)
   p_x(2:nn-1,:) = (in(3:nn,:) - in(1:nn-2,:)) / (2. * dd)
   p_x(nn,:) = (3. * in(nn,:) - 4. * in(nn-1,:) + in(nn-2,:)) / (2. * dd)

   p_y(:,1) = (-3. * in(:,1) + 4. * in(:,2) - in(:,3)) / (2. * dd)
   p_y(:,2:nn-1) = (in(:,3:nn) - in(:,1:nn-2)) / (2. * dd)
   p_y(:,nn) = (3. * in(:,nn) - 4. * in(:,nn-1) + in(:,nn-2)) / (2. * dd)

   out = 0.5 * pp * (p_x + ii * e1 * p_y) + spin * (e1 * x + ii * y) * in

end subroutine null_d1

subroutine null_d2 (out, in, spin, e1, e2)

   use null_grid, only : nn, ii, dd, x=>qs, y=>ps, pp

   implicit none

   double complex,   dimension (:,:), intent (out) :: out
   double complex,   dimension (:,:), intent (in)  :: in
   integer,                           intent (in)  :: spin, e1, e2

   double complex, dimension (nn,nn) :: p_x, p_y, p_xx, p_xy, p_yy
   double precision :: h1, h2, c1, c2, c3, c4, c5, c6

   h1 = 1. / (2. * dd)
   h2 = 1. / (dd * dd)

   c1 = 1 + e1 * e2 + spin * (e1 + e2)
   c2 = 2 * spin + e1 + e2
   c3 = 2 * spin * e1 * e2 + e1 + e2
   c4 = e1 * e2
   c5 = e1 + e2
   c6 = e2 - e1

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

   out = 0.25 * pp * pp * (p_xx - c4 * p_yy + ii * c5 * p_xy) &
        + 0.5 * pp * ((c1 * x + ii * c2 * y) * p_x &
        + (c2 * x + ii * c1 * y) * ii * c4 * p_y) &
        + spin * in * 0.5 * (c6 + c3 * (x * x - c4 * y * y) &
        + 2 * c1 * ii * x * y)

end subroutine null_d2

end module null_eth
