module null_interp2

contains

subroutine rinterp(f, fout, xs, ys)

   use null_grid, only : nn, dd, qs, ps, qmin
   implicit none

   double precision, dimension (nn,nn), intent (in) :: f
   double precision fout
   double precision xs, ys

   double precision x1, x2, x3, x4, y1, y2,y3,y4, xk, yl
   double precision xd, yd, factor
   integer kk, ll

   factor = 1. / (36. * (dd * dd) * (dd * dd) * (dd * dd))

   kk = min(max(int((xs - qmin) / dd) + 1, 3), nn-3)
   ll = min(max(int((ys - qmin) / dd) + 1, 3), nn-3)

   xk = qs(kk,ll)
   yl = ps(kk,ll)

   xd = xs - xk
   x1 = - xd * (xd - 2 * dd) * (xd - dd)
   x2 = 3 * (xd - 2 * dd) * (xd - dd) * (xd + dd)
   x3 = - 3 * xd * (xd -2 * dd) * (xd + dd)
   x4 = xd * (xd - dd) * (xd + dd)
      
   yd = ys - yl
   y1 = - yd * (yd - 2 * dd) * (yd - dd)
   y2 = 3 * (yd - 2 * dd) * (yd - dd) * (yd + dd)
   y3 = - 3 * yd * (yd - 2 * dd) * (yd + dd)
   y4 = yd * (yd - dd) * (yd + dd)
      
   fout = factor * (y1 * ( x1 * f(kk-1,ll-1) + x2 * f(kk,ll-1)         &   
                         + x3 * f(kk+1,ll-1) + x4 * f(kk+2,ll-1))      &  
                  + y2 * ( x1 * f(kk-1,ll)   + x2 * f(kk,ll)           &  
                         + x3 * f(kk+1,ll)   + x4 * f(kk+2,ll))        & 
                  + y3 * ( x1 * f(kk-1,ll+1) + x2 * f(kk,ll+1)         &
                         + x3 * f(kk+1,ll+1) + x4 * f(kk+2,ll+1))      &
                  + y4 * ( x1 * f(kk-1,ll+2) + x2 * f(kk,ll+2)         &
                         + x3 * f(kk+1,ll+2) + x4 * f(kk+2,ll+2)))

end subroutine rinterp

subroutine cinterp(f, fout, xs, ys)

   use null_grid, only : nn, dd, qs, ps, qmin
   implicit none

   double complex, dimension (nn,nn), intent (in) :: f
   double complex fout
   double precision xs, ys

   double precision x1, x2, x3, x4, y1, y2,y3,y4, xk, yl
   double precision xd, yd, factor
   integer kk, ll

   kk = min(max(int((xs - qmin) / dd) + 1, 3), nn-3)
   ll = min(max(int((ys - qmin) / dd) + 1, 3), nn-3)

   xk = qs(kk,ll)
   yl = ps(kk,ll)

   factor = 1. / (36. * (dd * dd) * (dd * dd) * (dd * dd))

   xd = xs - xk
   x1 = - xd * (xd - 2 * dd) * (xd - dd)
   x2 = 3 * (xd - 2 * dd) * (xd - dd) * (xd + dd)
   x3 = - 3 * xd * (xd -2 * dd) * (xd + dd)
   x4 = xd * (xd - dd) * (xd + dd)

   yd = ys - yl
   y1 = - yd * (yd - 2 * dd) * (yd - dd)
   y2 = 3 * (yd - 2 * dd) * (yd - dd) * (yd + dd)
   y3 = - 3 * yd * (yd - 2 * dd) * (yd + dd)
   y4 = yd * (yd - dd) * (yd + dd)
   
   fout = factor * (y1 * ( x1 * f(kk-1,ll-1) + x2 * f(kk,ll-1)        &     
                         + x3 * f(kk+1,ll-1) + x4 * f(kk+2,ll-1))     &    
                  + y2 * ( x1 * f(kk-1,ll)   + x2 * f(kk,ll)          &    
                         + x3 * f(kk+1,ll)   + x4 * f(kk+2,ll))       &   
                  + y3 * ( x1 * f(kk-1,ll+1) + x2 * f(kk,ll+1)        &   
                         + x3 * f(kk+1,ll+1) + x4 * f(kk+2,ll+1))     &   
                  + y4 * ( x1 * f(kk-1,ll+2) + x2 * f(kk,ll+2)        &   
                         + x3 * f(kk+1,ll+2) + x4 * f(kk+2,ll+2)))

end subroutine cinterp

subroutine null_rnsint2 (fn, fs)

   use null_grid, only : nn, qs, ps
   implicit none
   double precision, dimension (nn,nn), intent (inout) :: fn, fs

   integer i, j
   double precision xn, yn, xs, ys

   do j = 1, nn
      do i = 1, nn

         if (i==1 .or. i==nn .or. j==1 .or. j==nn) then

            xn = qs(i,j)
            yn = ps(i,j)

            xs =  xn / (xn * xn + yn * yn)
            ys = -yn / (xn * xn + yn * yn)

            call rinterp(fs, fn(i,j), xs, ys)
            call rinterp(fn, fs(i,j), xs, ys)

         end if

      end do
   end do

end subroutine null_rnsint2

subroutine null_cnsint2 (fn, fs, spin)

   use null_grid, only : nn, ii, qs, ps
   implicit none
   double complex, dimension (nn,nn), intent (inout) :: fn, fs
   integer, intent (in) :: spin

   integer i, j
   double precision xn, yn, xs, ys
   double complex zs, zsb

   do j = 1, nn
      do i = 1, nn

         if (i==1 .or. i==nn .or. j==1 .or. j==nn) then

            xn = qs(i,j)
            yn = ps(i,j)

            xs =  xn / (xn * xn + yn * yn)
            ys = -yn / (xn * xn + yn * yn)

            zs  = xs + ii * ys
            zsb = xs - ii * ys

            call cinterp(fs, fn(i,j), xs, ys)
            call cinterp(fn, fs(i,j), xs, ys)

            if (spin .ne. 0) then
               fn(i,j) = fn(i,j) * (-zsb / zs) ** spin
               fs(i,j) = fs(i,j) * (-zsb / zs) ** spin
            end if

         end if

      end do
   end do

end subroutine null_cnsint2

subroutine null_ccircint (fn, fs, spin)

   use null_grid, only : nn, ii, qs, ps, pp
   implicit none
   double complex, dimension (nn,nn), intent (inout) :: fn, fs
   integer, intent (in) :: spin

   integer i, j
   double precision xn, yn, xs, ys
   double complex zs, zsb
   double complex, dimension (nn,nn) :: tn, ts

   tn = fn
   ts = fs
   do j = 1, nn
      do i = 1, nn

         if (pp(i,j) .gt. 2.0d0) then

            xn = qs(i,j)
            yn = ps(i,j)

            xs =  xn / (xn * xn + yn * yn)
            ys = -yn / (xn * xn + yn * yn)

            zs  = xs + ii * ys
            zsb = xs - ii * ys

            call cinterp(ts, fn(i,j), xs, ys)
            call cinterp(tn, fs(i,j), xs, ys)

            if (spin .ne. 0) then
               fn(i,j) = fn(i,j) * (-zsb / zs) ** spin
               fs(i,j) = fs(i,j) * (-zsb / zs) ** spin
            end if

         end if

      end do
   end do

end subroutine null_ccircint

subroutine cJtoU(Jn, Un)

   use null_grid, only : nn
   implicit none

   double complex, dimension (:,:), intent (in) :: Jn
   double complex, dimension (:,:), intent (out) :: Un

   integer i, j

   do j = 2, nn
      do i = 2, nn
         Un(i,j) = 0.25 * (Jn(i-1,j-1) + Jn(i,j-1) + Jn(i-1,j) + Jn(i,j))
      end do
   end do

end subroutine cJtoU

end module null_interp2
