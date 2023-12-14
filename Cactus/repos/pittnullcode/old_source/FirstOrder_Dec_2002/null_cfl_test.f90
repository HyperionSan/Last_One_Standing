! Revised CFL. To remove all possible complications, let's
! think about it this way. Being at the point p = (n+1,j,k,i)
! we want to know if dt was too large. So we 'find the null rays'
! connecting p to those in the previous level that were
! involved in the calculation of the fields at p. The ones we care
! about are (n,j,k,i+1) and (n,j+A,k+B,i-1) where A,B=-1,1
! Certainly we can't evaluate everything cause we do not have the
! metric complete at n+1. We'll proceed for now doing a 'first order
! extrapolation, ie setting metric at n+1 = metric at n, and then
! we'll take a smaller value for dt to be safe....
! also note that du is 'negative' in the calculation, we'll absorb
! the sign in the metric quantities to deal with a positive dt

module null_cfl_test

contains

  subroutine null_cfl0(jos, bos, uos, wos, dt_old, maskb)

    use null_grid
    use null_eth

    implicit none

    double complex,   dimension (nn,nn,nx), intent(in)  :: jos, uos
    double precision, dimension (nn,nn,nx), intent(in)  :: bos, wos
    double precision,                       intent(inout) :: dt_old
    integer, dimension(nn,nn), intent(in) :: maskb
    double precision, dimension (9)     :: dtry1, dtry2, dis, dq, dp
    integer,          dimension (nn,nn) :: mask

    double precision :: du, dfor, dforU, dfor2, dste2
    double precision :: guu, gup, guq, gur, gpp, gqq, gpq, a
    double complex   :: JJ, Jb, U, Ub

    double precision :: r, P0, dr1, dr, KK, B, W, q, p, dmin1, dmin2, &
         dutime, dutime1, dutime2
    integer i,j,k


    dmin2 = 100.

    !this condition is to have a 'floor' for how small
    !we'll be willing to go. I'm setting it to 1e-10 but
    ! if we get to even 1e-5 we're in trouble if we want
    ! to evolve for times of order 1....

    if(dt_old.lt.1e-10) then
       print*, 'WARNING!!!! TOO SMALL DT!!!!'
       STOP
    end if

    !now set dutime to something big, it will be redefined
    dutime = 100. * dt_old

    !loop over all angular points involved in the algorithm

    dq = (/ dd, dd,  dd, 0.0d0, 0.0d0,  0.0d0, -dd, -dd, -dd/)
    dp = (/ dd, 0.0d0, -dd, dd, 0.0d0, -dd,  dd,  0.0d0, -dd/)

    !set all these equal to dt, ie the previous one obtained in this routine.
    !however we need to divide by cfl as when we leave this routine the main
    !code multiplies the dt obtained here by a 'safe factor' cfl. Since we'll
    !require never to go 'larger' than the dt previously used, if we don't 
    !'recover' the old dt obtained here by cfl, the net effect will be to
    ! keep multimpliying by "cfl factors" everytime we leave and we'll have
    ! dt converge to zero.

    dforU = dt_old/cfl
    dfor2 = dt_old/cfl
    dfor  = dt_old/cfl
    dmin1 = dt_old/cfl

    !evaluate the null 'distances' over all points
    !but masking the first 4 and the last one. just
    !in case, soon we'll use the 'real mask' from the
    !code to go over valid points only.

    do k = 2, nx - 1

       r = rb(k)
       dr1 = dx * rwt * 1. / (1.-x(k) )**2
       dr  = dx * rwt * 1. / (1.-x(k) )**2

       do j = 2, nn-1
          do i = 2, nn-1

             !use the mask here to evaluate only valid points
             if ( k .gt. maskb(i,j) ) then
               JJ = jos(i,j,k)
                Jb = conjg(JJ)
                B = bos(i,j,k)
                U = 0.5 * ( uos(i,j,k) + uos(i,j,k-1) )
                Ub = conjg(U)
                KK = sqrt( 1. + JJ*Jb )
                W = wos(i,j,k)

                q = qs(i,j)
                p = ps(i,j)
                P0 = 1. / (1. + q**2 + p**2)

                !define metric components
                guu = -exp(2.*B) * (1. + r * W )+r**2 &
                     * ( KK * U * Ub + dble(JJ*Ub**2) )
                gur = exp(2.*B) !note sign changed!
                guq = r**2 * 2.* P0 * dble( JJ*Ub + KK * U ) !sign changed
                gup = r**2 * 2.* P0 * dimag((JJ*Ub+KK*U) ) !sign changed
                gqq = (r*P0)**2 * 4. * ( KK + dble( JJ) )
                gpp = (r*P0)**2 * 4. * ( KK - dble( JJ) )
                gpq = (r*P0)**2 * 4. * dimag(JJ)

                !evaluate 'distance' of point in the radial direction. ignoring
                ! U for the moment
                dfor2 = 2. * dr1 * gur  / abs( exp(2.*B) * (1. + r * W) )

                !set to 10^-8 to avoid taking insanely small steps1
                if (dfor2 .lt. 1.e-8) then
                   dfor2 = 100.
                end if

                !calculate 'distance' i with the full guu

                if (abs(guu) .ge. 1.e-8) then
                   !radial driection only 
                   dfor = 2. * dr1 * gur  / (-guu)

                   if(dfor .le. 0.) then
                      dfor = 100.
                   end if

                   !with all directions.
	           dr = -dr 

                   dis(:) =  ( gur * dr + guq *dq(:) + gup*dp(:) )**2 &
                        - guu * ( 2.* gpq * dp(:) * dq(:) + gpp * dp(:)**2 &
                        + gqq * dq(:)**2 )

                   !check that we take square root of a positive quantity
                   where(dis .ge. 0.0d0)
                      dtry1(:) = (- ( gur * dr + guq * dq(:) + gup * dp(:) ) &
                           - dsqrt(dis(:) )) / guu

                      dtry2(:) = (- ( gur * dr + guq * dq(:) + gup * dp(:) ) &
                           + sqrt(dis(:) ))/ guu
                   elsewhere
                      dtry1(:) = 100.
                      dtry2(:) = 100.
                   end where

                   where(dtry1.lt.-1.e-8)
                      dtry1 = -dtry1
                   elsewhere
                      dtry1 = 100.
                   end where

                   where(dtry2.lt.-1.e-8)
                      dtry2 = -dtry2
                   elsewhere
                      dtry2 = 100.
                   end where
                   dmin1 = min( minval(dtry1), minval(dtry2))

                else  !(i.e. guu = 0.0)
	           dr = -dr 
                   where(abs(gur * dr + guq *dq(:) + gup*dp(:)).gt.1.e-8)
                      dtry1(:) = (2.* gpq * dp(:) * dq(:) + gpp * dp(:)**2  &
                           + gqq * dq(:)**2 ) / &
                           ( gur * dr + guq *dq(:) + gup*dp(:) ) * 0.5
                   elsewhere
                      dtry1(:) = 1.e-8
                   end where

		   dmin1 = max(1.e-8,minval(dtry1))

                end if

                dutime = min(dmin1, dfor, dfor2, dutime)
             endif
          end do
       end do
!       write(*,*) dutime, dmin1, dfor, dfor2
    end do

    if(dutime.lt.dt_old) dt_old =  dutime

    !FINAL WARNING IN CASE IS TOO SMALL

    if(dt_old.le.1.e-8) print*, 'WAY TOO SMALL DT FROM CFL ROUTINE!'


  end subroutine null_cfl0



  !! Note this 3 level subroutine is not being used at the moment
  !! we'll keep it here just in case we ever feel the urge to use it
  subroutine null_cfl(jns, jos, bns, bos, uns, uos, wns, wos, x_wt, dutime_old)

    use null_grid
    use null_eth

    implicit none

    double complex,   dimension (nn,nn,nx), intent(in)  :: jos, uos, jns, uns
    double precision, dimension (nn,nn,nx), intent(in)  :: bos, bns, wos, wns
    double precision, dimension (nn,nn),    intent(in)  :: x_wt
    double precision                        :: dutime, dutime_old

    double precision, dimension (9)     :: dtry1, dtry2, dis, dq, dp
    integer,          dimension (nn,nn) :: mask

    double precision :: du, dfor, dforU, dfor2, dste2
    double precision :: guu, gup, guq, gur, gpp, gqq, gpq, a
    double complex   :: JJ, Jb, U, Ub
    double precision :: r, P0, dr1, dr, KK, B, W, q, p, dmin1, dmin2,&
         dutime1, dutime2
    integer          :: i,j,k


    print*, 'WARNING YOU SHOULD NOT BE USING THIS ONE'

    dutime1 = dt
    dutime2 = dt

    dutime = dutime2
    dmin2  = 100.
    a      = 100.
    dmin1  = dutime2

    dfor  = dutime1
    dforU = dutime1
    dfor2 = dutime1

    dutime = 100. * dt


    dq = (/ dd, dd,     dd, 0.0d0, 0.0d0,   0.0d0, -dd, -dd,    -dd/)
    dp = (/ dd, 0.0d0, -dd, dd,    0.0d0,  -dd,     dd,  0.0d0, -dd/)

    do k = 4, nx - 1
       where ( x(k) .gt. x_wt )
          mask = 1
       elsewhere
          mask = 0
       end where

       r   = rb(k)
       dr1 = dx * rwt *  1. / (1.-x(k) )**2
       dr  = dx * rwt *  1. / (1.-x(k) )**2

       do j = 2, nn-1
          do i = 2, nn-1

             !here same as previous routine, except that we extrapolate
             !to obtain guess of metric at level n+3/2  (assuming we have at n+1 and n)

             JJ = 1.5 * jns(i,j,k) - 0.5 * jos(i,j,k)
             Jb = conjg(JJ)
             B = 1.5 * bns(i,j,k) - 0.5 * bos(i,j,k)
             U = 0.75 * (uns(i,j,k) + uns(i,j,k-1)) &
                  - 0.25 * (uos(i,j,k) + uos(i,j,k-1))

             Ub = conjg(U)
             KK = sqrt( 1. + JJ*Jb )
             W = 1.5 * wns(i,j,k) - 0.5 * wos(i,j,k)
             q = qs(i,j)
             p = ps(i,j)
             P0 = 1. / (1. + q**2 + p**2)

             !define metric components
             guu = -exp(2.*B) * (1. + r * W )+r**2 &
                  * ( KK * U * Ub + dble(JJ*Ub**2) )
             gur = exp(2.*B) !note sign changed!
             guq = r**2 * 2.* P0 * dble( JJ*Ub + KK * U ) !sign changed
             gup = r**2 * 2.* P0 * dimag((JJ*Ub+KK*U) ) !sign changed
             gqq = (r*P0)**2 * 4. * ( KK + dble( JJ) )
             gpp = (r*P0)**2 * 4. * ( KK - dble( JJ) )
             gpq = (r*P0)**2 * 4. * dimag(JJ)

             !evaluate 'distance' of point in the radial direction. ignoring
             ! U for the moment
             dfor2 = 2. * dr1 * gur  / abs( exp(2.*B) * (1. + r * W) )

             !set to 10^-8 to avoid taking insanely small steps1
             if (dfor2 .lt. 1.e-8) then
                dfor2 = 100.
             end if

             !calculate 'distance' i with the full guu

             if (guu .lt. 0) then
                !radial driection only 
                dfor = 2. * dr1 * gur  / (-guu)

                if(dfor .le. 0.) then
                   dfor = 100.
                end if

                !with all directions.
                dr = -dr

                dis(:) =  ( gur * dr + guq *dq(:) + gup*dp(:) )**2 &
                     - guu * ( 2.* gpq * dp(:) * dq(:) + gpp * dp(:)**2 &
                     + gqq * dq(:)**2 )

                !check that we take square root of a positive quantity
                where(dis .ge. 0.0d0)
                   dtry1(:) = (- ( gur * dr + guq * dq(:) + gup * dp(:) ) &
                        - dsqrt(dis(:) )) / guu

                   dtry2(:) = (- ( gur * dr + guq * dq(:) + gup * dp(:) ) &
                        + sqrt(dis(:) ))/ guu
                elsewhere
                   dtry1(:) = 100.
                   dtry2(:) = 100.
                end where

                where(dtry1.lt.-1.e-8)
                   dtry1 = -dtry1
                elsewhere
                   dtry1 = 100.
                end where

                where(dtry2.lt.-1.e-8)
                   dtry2 = -dtry2
                elsewhere
                   dtry2 = 100.
                end where
                dmin1 = min( minval(dtry1), minval(dtry2))
             else
                dmin1 = dt
                dfor  = dt
             end if

             dutime = min(dmin1, dfor, dfor2, dutime)

          end do
       end do
    end do

    if(dutime.lt.dutime_old) dutime_old = dutime
  end subroutine null_cfl

end module null_cfl_test
