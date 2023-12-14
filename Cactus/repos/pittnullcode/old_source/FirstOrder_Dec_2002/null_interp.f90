module null_interp

contains

   subroutine null_rnsint (jnn, jns)

   use null_grid, only : nn, dd, qs, ps, qmin

   implicit none

   double precision, dimension (:,:), intent (inout) :: jnn
   double precision, dimension (:,:), intent (in)    :: jns

   integer k, l, kk, ll
   double precision  xn, yn, xs, ys
   double precision  xk, xd, x1, x2, x3, x4
   double precision  yl, yd, y1, y2, y3, y4
   double precision  factor
   double precision  jnsinterp

   factor = 1. / (36. * (dd * dd) * (dd * dd) * (dd * dd))

   do k = 1, nn, nn-1
      do l = 1, nn
         xn = qs(k,l)
         yn = ps(k,l)
   
         xs =  xn / (xn * xn + yn * yn)
         ys = -yn / (xn * xn + yn * yn)
   
!        kk = int((xs - qmin) / dd) + 1
!        ll = int((ys - qmin) / dd) + 1
         kk = min(max(int((xs - qmin) / dd) + 1, 3), nn-3)
         ll = min(max(int((ys - qmin) / dd) + 1, 3), nn-3)
   
         xk = qs(kk,ll)
         yl = ps(kk,ll)
   
         xd = xs - xk
         x1 = - xd * (xd - 2 * dd) * (xd - dd)
         x2 = 3. * (xd - 2 * dd) * (xd - dd) * (xd + dd)
         x3 = - 3. * xd * (xd -2 * dd) * (xd + dd)
         x4 = xd * (xd - dd) * (xd + dd)
   
         yd = ys - yl
         y1 = - yd * (yd - 2 * dd) * (yd - dd)
         y2 = 3. * (yd - 2 * dd) * (yd - dd) * (yd + dd)
         y3 = - 3. * yd * (yd - 2 * dd) * (yd + dd)
         y4 = yd * (yd - dd) * (yd + dd)
   
         jnsinterp = factor * (y1 * ( x1 * jns(kk-1,ll-1) + x2 * jns(kk  ,ll-1)  &
                                    + x3 * jns(kk+1,ll-1) + x4 * jns(kk+2,ll-1)) &
                             + y2 * ( x1 * jns(kk-1,ll)   + x2 * jns(kk  ,ll)    &
                                    + x3 * jns(kk+1,ll)   + x4 * jns(kk+2,ll))   &
                             + y3 * ( x1 * jns(kk-1,ll+1) + x2 * jns(kk  ,ll+1)  &
                                    + x3 * jns(kk+1,ll+1) + x4 * jns(kk+2,ll+1)) &
                             + y4 * ( x1 * jns(kk-1,ll+2) + x2 * jns(kk  ,ll+2)  &
                                    + x3 * jns(kk+1,ll+2) + x4 * jns(kk+2,ll+2)))
   
         jnn(k,l) = jnsinterp
      end do
   end do
   
   do l = 1, nn, nn-1
      do k = 1, nn
         xn = qs(k,l)
         yn = ps(k,l)
   
         xs =  xn / (xn * xn + yn * yn)
         ys = -yn / (xn * xn + yn * yn)
   
!        kk = int((xs - qmin) / dd) + 1
!        ll = int((ys - qmin) / dd) + 1
         kk = min(max(int((xs - qmin) / dd) + 1, 3), nn-3)
         ll = min(max(int((ys - qmin) / dd) + 1, 3), nn-3)
   
         xk = qs(kk,ll)
         yl = ps(kk,ll)
   
         xd = xs - xk
         x1 = - xd * (xd - 2 * dd) * (xd - dd)
         x2 = 3. * (xd - 2 * dd) * (xd - dd) * (xd + dd)
         x3 = - 3. * xd * (xd -2 * dd) * (xd + dd)
         x4 = xd * (xd - dd) * (xd + dd)
   
         yd = ys - yl
         y1 = - yd * (yd - 2 * dd) * (yd - dd)
         y2 = 3. * (yd - 2 * dd) * (yd - dd) * (yd + dd)
         y3 = - 3. * yd * (yd - 2 * dd) * (yd + dd)
         y4 = yd * (yd - dd) * (yd + dd)
   
         jnsinterp = factor * (y1 * ( x1 * jns(kk-1,ll-1) + x2 * jns(kk  ,ll-1)  &
                                    + x3 * jns(kk+1,ll-1) + x4 * jns(kk+2,ll-1)) &
                             + y2 * ( x1 * jns(kk-1,ll)   + x2 * jns(kk  ,ll)    &
                                    + x3 * jns(kk+1,ll)   + x4 * jns(kk+2,ll))   &
                             + y3 * ( x1 * jns(kk-1,ll+1) + x2 * jns(kk  ,ll+1)  &
                                    + x3 * jns(kk+1,ll+1) + x4 * jns(kk+2,ll+1)) &
                             + y4 * ( x1 * jns(kk-1,ll+2) + x2 * jns(kk  ,ll+2)  &
                                    + x3 * jns(kk+1,ll+2) + x4 * jns(kk+2,ll+2)))
   
         jnn(k,l) = jnsinterp
      end do
   end do

end subroutine null_rnsint

subroutine null_cnsint (jnn, jns, spin)

   use null_grid, only : ii, nn, dd, qs, ps, qmin

   implicit none

   double complex, dimension (:,:), intent (inout) :: jnn
   double complex, dimension (:,:), intent (in)    :: jns
   integer,                         intent (in)    :: spin

   integer k, l, kk, ll
   double precision  xn, yn, xs, ys
   double precision  xk, xd, x1, x2, x3, x4
   double precision  yl, yd, y1, y2, y3, y4
   double precision  factor
   double complex zeta, zetab, jnsinterp

   factor = 1. / (36. * (dd * dd) * (dd * dd) * (dd * dd))

   do k = 1, nn, nn-1
      do l = 1, nn
         xn = qs(k,l)
         yn = ps(k,l)
   
         xs =  xn / (xn * xn + yn * yn)
         ys = -yn / (xn * xn + yn * yn)
   
!        kk = int((xs - qmin) / dd) + 1
!        ll = int((ys - qmin) / dd) + 1
         kk = min(max(int((xs - qmin) / dd) + 1, 3), nn-3)
         ll = min(max(int((ys - qmin) / dd) + 1, 3), nn-3)
   
         xk = qs(kk,ll)
         yl = ps(kk,ll)
   
         xd = xs - xk
         x1 = - xd * (xd - 2 * dd) * (xd - dd)
         x2 = 3. * (xd - 2 * dd) * (xd - dd) * (xd + dd)
         x3 = - 3. * xd * (xd -2 * dd) * (xd + dd)
         x4 = xd * (xd - dd) * (xd + dd)
   
         yd = ys - yl
         y1 = - yd * (yd - 2 * dd) * (yd - dd)
         y2 = 3. * (yd - 2 * dd) * (yd - dd) * (yd + dd)
         y3 = - 3. * yd * (yd - 2 * dd) * (yd + dd)
         y4 = yd * (yd - dd) * (yd + dd)
   
         jnsinterp = factor * (y1 * ( x1 * jns(kk-1,ll-1) + x2 * jns(kk  ,ll-1)  &
                                    + x3 * jns(kk+1,ll-1) + x4 * jns(kk+2,ll-1)) &
                             + y2 * ( x1 * jns(kk-1,ll)   + x2 * jns(kk  ,ll)    &
                                    + x3 * jns(kk+1,ll)   + x4 * jns(kk+2,ll))   &
                             + y3 * ( x1 * jns(kk-1,ll+1) + x2 * jns(kk  ,ll+1)  &
                                    + x3 * jns(kk+1,ll+1) + x4 * jns(kk+2,ll+1)) &
                             + y4 * ( x1 * jns(kk-1,ll+2) + x2 * jns(kk  ,ll+2)  &
                                    + x3 * jns(kk+1,ll+2) + x4 * jns(kk+2,ll+2)))
   
         zeta  = xs + ii * ys
         zetab = xs - ii * ys
         if (spin .ne. 0) then
            jnn(k,l) = jnsinterp * (-zetab / zeta) ** spin
         else
            jnn(k,l) = jnsinterp
         end if
      end do
   end do

   do l = 1, nn, nn-1
      do k = 1, nn
         xn = qs(k,l)
         yn = ps(k,l)
   
         xs =  xn / (xn * xn + yn * yn)
         ys = -yn / (xn * xn + yn * yn)
   
!        kk = int((xs - qmin) / dd) + 1
!        ll = int((ys - qmin) / dd) + 1
         kk = min(max(int((xs - qmin) / dd) + 1, 3), nn-3)
         ll = min(max(int((ys - qmin) / dd) + 1, 3), nn-3)
   
         xk = qs(kk,ll)
         yl = ps(kk,ll)
   
         xd = xs - xk
         x1 = - xd * (xd - 2 * dd) * (xd - dd)
         x2 = 3. * (xd - 2 * dd) * (xd - dd) * (xd + dd)
         x3 = - 3. * xd * (xd -2 * dd) * (xd + dd)
         x4 = xd * (xd - dd) * (xd + dd)
   
         yd = ys - yl
         y1 = - yd * (yd - 2 * dd) * (yd - dd)
         y2 = 3. * (yd - 2 * dd) * (yd - dd) * (yd + dd)
         y3 = - 3. * yd * (yd - 2 * dd) * (yd + dd)
         y4 = yd * (yd - dd) * (yd + dd)
   
         jnsinterp = factor * (y1 * ( x1 * jns(kk-1,ll-1) + x2 * jns(kk  ,ll-1)  &
                                    + x3 * jns(kk+1,ll-1) + x4 * jns(kk+2,ll-1)) &
                             + y2 * ( x1 * jns(kk-1,ll)   + x2 * jns(kk  ,ll)    &
                                    + x3 * jns(kk+1,ll)   + x4 * jns(kk+2,ll))   &
                             + y3 * ( x1 * jns(kk-1,ll+1) + x2 * jns(kk  ,ll+1)  &
                                    + x3 * jns(kk+1,ll+1) + x4 * jns(kk+2,ll+1)) &
                             + y4 * ( x1 * jns(kk-1,ll+2) + x2 * jns(kk  ,ll+2)  &
                                    + x3 * jns(kk+1,ll+2) + x4 * jns(kk+2,ll+2)))
   
         zeta  = xs + ii * ys
         zetab = xs - ii * ys
         if (spin .ne. 0) then
            jnn(k,l) = jnsinterp * (-zetab / zeta) ** spin
         else
            jnn(k,l) = jnsinterp
         end if
      end do
   end do

end subroutine null_cnsint

subroutine null_radintp (xp, fin, fout)

   use null_grid

   implicit none

   double precision, dimension (:,:),   intent (in)  :: xp
   double precision, dimension (:,:,:), intent (in)  :: fin
   double precision, dimension (:,:),   intent (out) :: fout

   integer i, j, k
   double precision  factor, xl, xk, xd, xw1, xw2, xw3, xw4

   factor = 1. / (6. * dd * dd * dd)
   xl = x(1)

   do i = 1, nn
      do j = 1, nn

         k = int((xp(i,j) - xl) / dx) + 1
         xk = xl + (k - 1) * dx
         xd = xp(i,j) - xk
         xw1 = - xd * ( - xd + dx) * ( - xd + 2 * dx)
         xw2 = 3 * (xd + dx) * ( - xd + dx) * ( - xd + 2 * dx)
         xw3 = 3 * (xd + dx) * xd * ( - xd + 2 * dx)
         xw4 = - (xd + dx) * xd * ( - xd + dx)

         fout(i,j) = factor * (xw1 * fin(i,j,k-1) + xw2 * fin(i,j,k)&
              & + xw3 * fin(i,j,k+1) + xw4 * fin(i,j,k+2))
      end do
   end do

end subroutine null_radintp

end module null_interp
