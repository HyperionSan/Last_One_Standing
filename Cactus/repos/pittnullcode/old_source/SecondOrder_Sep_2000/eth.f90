module eth

contains

subroutine d1 (nn, out, in, spin, e1)

   implicit none

   double complex, dimension (:,:), intent (out) :: out
   double complex, dimension (:,:), intent (in)  :: in
   integer,                         intent (in)  :: nn, spin, e1

   double precision, dimension (nn,nn) :: x, y, pp
   double complex,   dimension (nn,nn) :: p_x, p_y

   integer :: i
   double precision :: dd
   double complex, parameter :: ii = (0., 1.)

   dd = 2.0 / dble(nn - 5)

   do i = 1, nn
      x(i,:) = -1.0 + (i-3) * dd
      y(:,i) = -1.0 + (i-3) * dd
   end do

   pp = 0.5 * (1. + x * x + y * y)

   p_x(1,:) = (-3. * in(1,:) + 4. * in(2,:) - in(3,:)) / (2. * dd)
   p_x(2:nn-1,:) = (in(3:nn,:) - in(1:nn-2,:)) / (2. * dd)
   p_x(nn,:) = (3. * in(nn,:) - 4. * in(nn-1,:) + in(nn-2,:)) / (2. * dd)

   p_y(:,1) = (-3. * in(:,1) + 4. * in(:,2) - in(:,3)) / (2. * dd)
   p_y(:,2:nn-1) = (in(:,3:nn) - in(:,1:nn-2)) / (2. * dd)
   p_y(:,nn) = (3. * in(:,nn) - 4. * in(:,nn-1) + in(:,nn-2)) / (2. * dd)

   out = pp * (p_x + ii * e1 * p_y) + spin * (e1 * x + ii * y) * in

end subroutine d1

subroutine d2 (nn, out, in, spin, e1, e2)

   implicit none

   double complex, dimension (:,:), intent (out) :: out
   double complex, dimension (:,:), intent (in)  :: in
   integer,                         intent (in)  :: nn, spin, e1, e2

   double precision, dimension (nn,nn) :: x, y, pp
   double complex,   dimension (nn,nn) :: p_x, p_y
   double complex,   dimension (nn,nn) :: p_xy, p_xx, p_yy

   double precision :: h1, h2, c1, c2, c3, c4, c5, c6

   integer :: i
   double precision :: dd
   double complex, parameter :: ii = (0., 1.)

   dd = 2.0 / dble(nn - 5)

   h1 = 1. / (2. * dd)
   h2 = 1. / (dd * dd)

   c1 = 1 + e1 * e2 + spin * (e1 + e2)
   c2 = 2 * spin + e1 + e2
   c3 = 2 * spin * e1 * e2 + e1 + e2
   c4 = e1 * e2
   c5 = e1 + e2
   c6 = e2 - e1

   do i = 1, nn
      x(i,:) = -1.0 + (i-3) * dd
      y(:,i) = -1.0 + (i-3) * dd
   end do

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

end subroutine d2

subroutine cnsint (nn, jnn, jns, spin)

   implicit none

   integer,                         intent (in)    :: nn, spin
   double complex, dimension (:,:), intent (inout) :: jnn
   double complex, dimension (:,:), intent (in)    :: jns

   integer :: k, l, kk, ll
   double precision :: xn, yn, xs, ys
   double precision :: xk, xd, x1, x2, x3, x4
   double precision :: yl, yd, y1, y2, y3, y4
   double precision :: dd, factor
   double complex :: zeta, zetab, jnsinterp
   double complex, parameter :: ii = (0., 1.)

   dd = 2.0 / dble(nn - 5)
   factor = 1. / (36. * (dd * dd) * (dd * dd) * (dd * dd))

   k = 1
   do l = 1, nn
      xn = -1. + (k-3) * dd
      yn = -1. + (l-3) * dd

      xs =  xn / (xn * xn + yn * yn)
      ys = -yn / (xn * xn + yn * yn)

      kk = int((xs + 1.) / dd) + 3
      xk = -1. + (kk-3) * dd
      xd = xs - xk
      x1 = - xd * (xd - 2 * dd) * (xd - dd)
      x2 = 3. * (xd - 2 * dd) * (xd - dd) * (xd + dd)
      x3 = - 3. * xd * (xd -2 * dd) * (xd + dd)
      x4 = xd * (xd - dd) * (xd + dd)

      ll = int((ys + 1.) / dd) + 3
      yl = -1. + (ll-3) * dd
      yd = ys - yl
      y1 = - yd * (yd - 2 * dd) * (yd - dd)
      y2 = 3. * (yd - 2 * dd) * (yd - dd) * (yd + dd)
      y3 = - 3. * yd * (yd - 2 * dd) * (yd + dd)
      y4 = yd * (yd - dd) * (yd + dd)

      jnsinterp = factor * (y1 * ( x1 * jns(kk-1,ll-1) + x2 * jns(kk&
           & ,ll-1) + x3 * jns(kk+1,ll-1) + x4 * jns(kk+2,ll-1)) +&
           & y2 * ( x1 * jns(kk-1,ll)   + x2 * jns(kk,ll) + x3 *&
           & jns(kk+1,ll)   + x4 * jns(kk+2,ll)) + y3 * ( x1 *&
           & jns(kk-1,ll+1) + x2 * jns(kk,ll+1) + x3 * jns(kk+1,ll&
           & +1) + x4 * jns(kk+2,ll+1)) + y4 * ( x1 * jns(kk-1,ll+2)&
           & + x2 * jns(kk,ll+2) + x3 * jns(kk+1,ll+2) + x4 * jns(kk&
           & +2,ll+2)))

      zeta  = xs + ii * ys
      zetab = xs - ii * ys
      if (spin .ne. 0) then
         jnn(k,l) = jnsinterp * ((-zetab / zeta) ** spin)
      else
         jnn(k,l) = jnsinterp
      end if
   end do

   k = nn
   do l = 1, nn
      xn = -1. + (k-3) * dd
      yn = -1. + (l-3) * dd

      xs =  xn / (xn * xn + yn * yn)
      ys = -yn / (xn * xn + yn * yn)

      kk = int((xs + 1.) / dd) + 3
      xk = -1. + (kk-3) * dd
      xd = xs - xk
      x1 = - xd * (xd - 2 * dd) * (xd - dd)
      x2 = 3. * (xd - 2 * dd) * (xd - dd) * (xd + dd)
      x3 = - 3. * xd * (xd -2 * dd) * (xd + dd)
      x4 = xd * (xd - dd) * (xd + dd)

      ll = int((ys + 1.) / dd) + 3
      yl = -1. + (ll-3) * dd
      yd = ys - yl
      y1 = - yd * (yd - 2 * dd) * (yd - dd)
      y2 = 3. * (yd - 2 * dd) * (yd - dd) * (yd + dd)
      y3 = - 3. * yd * (yd - 2 * dd) * (yd + dd)
      y4 = yd * (yd - dd) * (yd + dd)

      jnsinterp = factor * (y1 * ( x1 * jns(kk-1,ll-1) + x2 * jns(kk&
           & ,ll-1) + x3 * jns(kk+1,ll-1) + x4 * jns(kk+2,ll-1)) +&
           & y2 * ( x1 * jns(kk-1,ll)   + x2 * jns(kk,ll) + x3 *&
           & jns(kk+1,ll)   + x4 * jns(kk+2,ll)) + y3 * ( x1 *&
           & jns(kk-1,ll+1) + x2 * jns(kk,ll+1) + x3 * jns(kk+1,ll&
           & +1) + x4 * jns(kk+2,ll+1)) + y4 * ( x1 * jns(kk-1,ll+2)&
           & + x2 * jns(kk,ll+2) + x3 * jns(kk+1,ll+2) + x4 * jns(kk&
           & +2,ll+2)))

      zeta  = xs + ii * ys
      zetab = xs - ii * ys
      if (spin .ne. 0) then
         jnn(k,l) = jnsinterp * ((-zetab / zeta) ** spin)
      else
         jnn(k,l) = jnsinterp
      end if
   end do

   l = 1
   do k = 1, nn
      xn = -1. + (k-3) * dd
      yn = -1. + (l-3) * dd

      xs =  xn / (xn * xn + yn * yn)
      ys = -yn / (xn * xn + yn * yn)

      kk = int((xs + 1.) / dd) + 3
      xk = -1. + (kk-3) * dd
      xd = xs - xk
      x1 = - xd * (xd - 2 * dd) * (xd - dd)
      x2 = 3. * (xd - 2 * dd) * (xd - dd) * (xd + dd)
      x3 = - 3. * xd * (xd -2 * dd) * (xd + dd)
      x4 = xd * (xd - dd) * (xd + dd)

      ll = int((ys + 1.) / dd) + 3
      yl = -1. + (ll-3) * dd
      yd = ys - yl
      y1 = - yd * (yd - 2 * dd) * (yd - dd)
      y2 = 3. * (yd - 2 * dd) * (yd - dd) * (yd + dd)
      y3 = - 3. * yd * (yd - 2 * dd) * (yd + dd)
      y4 = yd * (yd - dd) * (yd + dd)

      jnsinterp = factor * (y1 * ( x1 * jns(kk-1,ll-1) + x2 * jns(kk&
           & ,ll-1) + x3 * jns(kk+1,ll-1) + x4 * jns(kk+2,ll-1)) +&
           & y2 * ( x1 * jns(kk-1,ll)   + x2 * jns(kk,ll) + x3 *&
           & jns(kk+1,ll)   + x4 * jns(kk+2,ll)) + y3 * ( x1 *&
           & jns(kk-1,ll+1) + x2 * jns(kk,ll+1) + x3 * jns(kk+1,ll&
           & +1) + x4 * jns(kk+2,ll+1)) + y4 * ( x1 * jns(kk-1,ll+2)&
           & + x2 * jns(kk,ll+2) + x3 * jns(kk+1,ll+2) + x4 * jns(kk&
           & +2,ll+2)))

      zeta  = xs + ii * ys
      zetab = xs - ii * ys
      if (spin .ne. 0) then
         jnn(k,l) = jnsinterp * ((-zetab / zeta) ** spin)
      else
         jnn(k,l) = jnsinterp
      end if
   end do

   l = nn
   do k = 1, nn
      xn = -1. + (k-3) * dd
      yn = -1. + (l-3) * dd

      xs =  xn / (xn * xn + yn * yn)
      ys = -yn / (xn * xn + yn * yn)

      kk = int((xs + 1.) / dd) + 3
      xk = -1. + (kk-3) * dd
      xd = xs - xk
      x1 = - xd * (xd - 2 * dd) * (xd - dd)
      x2 = 3. * (xd - 2 * dd) * (xd - dd) * (xd + dd)
      x3 = - 3. * xd * (xd -2 * dd) * (xd + dd)
      x4 = xd * (xd - dd) * (xd + dd)

      ll = int((ys + 1.) / dd) + 3
      yl = -1. + (ll-3) * dd
      yd = ys - yl
      y1 = - yd * (yd - 2 * dd) * (yd - dd)
      y2 = 3. * (yd - 2 * dd) * (yd - dd) * (yd + dd)
      y3 = - 3. * yd * (yd - 2 * dd) * (yd + dd)
      y4 = yd * (yd - dd) * (yd + dd)

      jnsinterp = factor * (y1 * ( x1 * jns(kk-1,ll-1) + x2 * jns(kk&
           & ,ll-1) + x3 * jns(kk+1,ll-1) + x4 * jns(kk+2,ll-1)) +&
           & y2 * ( x1 * jns(kk-1,ll)   + x2 * jns(kk,ll) + x3 *&
           & jns(kk+1,ll)   + x4 * jns(kk+2,ll)) + y3 * ( x1 *&
           & jns(kk-1,ll+1) + x2 * jns(kk,ll+1) + x3 * jns(kk+1,ll&
           & +1) + x4 * jns(kk+2,ll+1)) + y4 * ( x1 * jns(kk-1,ll+2)&
           & + x2 * jns(kk,ll+2) + x3 * jns(kk+1,ll+2) + x4 * jns(kk&
           & +2,ll+2)))

      zeta  = xs + ii * ys
      zetab = xs - ii * ys
      if (spin .ne. 0) then
         jnn(k,l) = jnsinterp * ((-zetab / zeta) ** spin)
      else
         jnn(k,l) = jnsinterp
      end if
   end do

end subroutine cnsint

end module eth
