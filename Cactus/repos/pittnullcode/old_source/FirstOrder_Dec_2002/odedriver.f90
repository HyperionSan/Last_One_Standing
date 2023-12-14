subroutine solve_focuseq(time, it, nt, dt, debug)
!===============================================
!purpose: call rksuite solver of the focussing equation
!         call embedding routine
use prec,           only: wp
use affine,         only: udot, sigma, sigma2, b,                        &
                          t_zero, t_end, th_minus, u_zero,               &
                          i_equator, p_omega, omega, init
use gridtranslator, only: theta, n_x, n_u, uhatofu, omvec,               &
                          stereogrid_to_thetavect
use Model,          only: eps, t_
use rksuite_90,  wp_rks => wp

implicit none

interface

 subroutine plotne(parvec, uhmat, omM)
   ! purpose:
   use affine
   real(kind=wp), intent(in), dimension(:)   :: parvec
   real(kind=wp), intent(in), dimension(:)   :: uhmat, omM
 end subroutine plotne

end interface

!arguments
integer,       intent(in) :: it, nt, debug
real(kind=wp), intent(in) :: time, dt

!local variables defining the numerical setup
integer                                  :: totf, flag, rksuite77method, lenwrk
real(kind=wp)                            :: tol
real(kind=wp), dimension(:), allocatable, save :: y_start, thres, ymax, &
                                            y_maxvals, uderiv, work
! local loop variables
real(kind=wp) :: t_want, t_got
integer       :: i

! local variables involved in checking for caustics
real(kind=wp)                                       :: d_th_max

! we do a 1d system - use rksuite's appropriate type
type(rk_comm_real_1d), save :: comm

! other local quantities
integer, save             :: times_called
logical                   :: fulloutput = .false., rksuitef90 = .true.
!< === executable statemets === >!

times_called = times_called + 1
if (debug > 0) then
   write(*,*) 'in solve_focuseq: it =', it, 'dt=',dt, 'time=',time
endif

n_u    = nt !! was: + 2

! get variables from module Model
b      = sqrt(1.0_wp + eps)
t_end  = t_ + nt * dt

if (times_called == 1) then
   write(*,*)
   if (it .ne. 0) stop '"it" not zero @ initialization step'
   open (unit = 10, file = "spheroid.log")

   ! construct the grid, get n_x
   call  stereogrid_to_thetavect(n_x, debug)

   ! allocate arrays
   allocate(uhatofu(n_x), uderiv(n_x), ymax(n_x), omvec(n_x), &
        &   y_start(n_x), thres(n_x))
   
   ! initialize sigma, p_omega, etc.
   call init

   !< prepare the ODE solution >!
   tol   = 1.0e-10_wp
   thres = 1.0e-14_wp

   !    We solve the  Equation d\hat u / du = \Lambda (\hat u) 
   !    for \hat u(u) with the integration constant
   !    that \hat u = ystart on the initial slice.
   ! NOTE: the code uses different conventions than the paper:
   !       the next comment translates these - with th the  \hat t of the paper

   y_start = u_zero + th_minus - t_zero ! th = th_minus - t_zero

   if (rksuitef90) then
      call setup(comm,t_,y_start,t_end,tol,thres,'M','R',.false.,0.0_wp,.true.)
   else
      rksuite77method = 1  ! medium accuracy
      lenwrk = 32 * n_x
      allocate(work(lenwrk))
      call setup77(n_x, t_, y_start, t_end, tol, thres, rksuite77method, &
           &        'U', .false., 0.0d0, work, lenwrk, .true.)
   endif
   uhatofu = y_start

   !< write configuration >!
   write(10, *)
   write(10, *) 'ODE tolerance  [tol]            = ', tol
   write(10, *) 'ODE threshold  [min/max(thres)] = ', minval(thres), &
        &                                             maxval(thres) 
   write(10, *) 'ODE start time [t_]             = ', t_ 
   write(10, *) 'ODE end time   [t_end]          = ', t_end 
   write(10, *) 'ODE output dt  [dt]             = ', dt
   write(10, *)

   write( *, *) 'RKSUITE initialized in solve_focuseq'
   write( *, *)
deallocate(y_start, thres)
endif

if (it > 0) then

!< solve the ODEs from it = 1, n_u >!

   !t_want = t_end - (n_u - it)*dt
   !t_want = time + dt
   t_want = time
   write(10, *) '*** call range_integrate @ time lvl',it,'new time',t_want

   if (rksuitef90) then
      call range_integrate(comm, udot, t_want, t_got, uhatofu, uderiv, flag)
   else
      call ut(udot, t_want, t_got, uhatofu, uderiv, ymax, work, flag)
   endif
   d_th_max = maxval(uhatofu + 0.5_wp*sigma) 
   write(10, *) '*** CALLED range_integrate @ time lvl', it, &
        & 'new time =', t_want, 'flag = ', flag
   if ( d_th_max < 0.0) then
      write(10,98) '       max(u^ + sigma/2) =', d_th_max
   else
      write(10,99) '       max(u^ + sigma/2) =', d_th_max, '  CAUSTIC!'
      write(*,  *) 'we went past a caustic at time step', it
   endif
endif

! construct an array for the omega-values
do i = 1, n_x
   omvec(i) = omega(sigma2(i), p_omega(i), uhatofu(i))
end do

if(fulloutput) then
   call plotflat(0.0_wp, 1.0_wp, theta(1), theta(n_x), n_x, n_u)
   call plotne(theta, uhatofu, omvec)
endif


if (it == n_u) then
   if (rksuitef90) then
      
      write(*, *) 'cleaning up RKSUITE'
      allocate(y_maxvals(n_x))
      ! collect data, clean up
      call statistics(comm,y_maxvals=y_maxvals,total_f_calls=totf)
      !write (*,*) ' y_maxvals=', y_maxvals
      deallocate(y_maxvals)
      
      write (10,*)
      write (10,*) ' The cost of the integration in calls to RHS was', totf
      write (10,*)
      
      call collect_garbage(comm)
   else      
      deallocate(work, ymax)
   end if

endif

return
98 format(A,1E13.5)
99 format(A,1E13.5,A)
end subroutine solve_focuseq


!module flatsurface
!contains

subroutine getIP(x,v)
! purpose

use affine

real(kind=wp), intent(in)                :: x
real(kind=wp), dimension(2), intent(out) :: v

v(1) =     sin(x) ! r0
v(2) = b * cos(x) ! z0
return
end subroutine getIP


subroutine tangent(th,tang)
! purpose

use affine

real(kind=wp), intent(in)                :: th
real(kind=wp), dimension(2), intent(out) :: tang

tang(1) =      cos(th) ! r0'
tang(2) = -b * sin(th) ! z0'
return
end subroutine tangent


subroutine normal(th, nvec)
! purpose

use prec

real(kind=wp), intent(in)   :: th
real(kind=wp), dimension(2) :: nvec, tang

double precision :: norm

call tangent(th, tang)
norm = 1.0_wp/sqrt(tang(1)**2 + tang(2)**2)

nvec(1) =  tang(2)*norm
nvec(2) = -tang(1)*norm

return
end subroutine normal


subroutine flatpoint(th,lambda,point)
! purpose

use prec

real(kind=wp), intent(in)                :: th, lambda
real(kind=wp), dimension(2), intent(out) :: point
real(kind=wp), dimension(2)              :: x, nvec

call getIP(th,x)
call normal(th,nvec)

point = x + lambda * nvec

return
end subroutine flatpoint

!end module flatsurface

subroutine plotflat(t_start, t_end, x_start, x_end, x_points, t_points)
! purpose:

use prec

implicit none

interface
subroutine flatpoint(th,lambda,point)
! purpose
use prec

real(kind=wp), intent(in)                :: th, lambda
real(kind=wp), dimension(2), intent(out) :: point
end subroutine flatpoint
end interface

real(kind=wp), intent(in) :: t_start, t_end, x_start, x_end
integer,       intent(in) :: x_points, t_points
integer :: flatfile1 = 11, flatfile2 = 12
integer :: ix, it
real(kind=wp)               :: t, th, dt, dth
real(kind=wp), dimension(2) :: point


dt  = (t_end - t_start)/dble(t_points-1) 
dth = (x_end - x_start)/dble(x_points-1)
  point = 0.0_wp
  open(flatfile1, file='flatfile1.mesh')
  open(flatfile2, file='flatfile2.mesh')
  write(flatfile1,*) '{CMESH'
  write(flatfile2,*) '{CMESH'
  write(flatfile1,*) x_points, t_points
  write(flatfile2,*) x_points, t_points

  t = t_start
  do it = 1, t_points
    th = 0.0_wp
    do ix = 1, x_points
           call flatpoint(th,t,point)
           write(flatfile1,*)  point(1), point(2), t, 1.0, 0.0, 0.0, 1.0
           write(flatfile2,*) -point(1), point(2), t, 0.0, 0.0, 1.0, 1.0
           th = th + dth
    end do
    t = t + dt
  end do

write(flatfile1,*) '}'
write(flatfile2,*) '}'
close(flatfile1)
close(flatfile2)

return
end subroutine plotflat





subroutine plotne(parvec, uhmat, omM)
  ! purpose:

  use affine
  use gridtranslator, only: n_x, n_u

  implicit none

  interface
     subroutine flatpoint(th,lambda,point)
       ! purpose
       use prec

       real(kind=wp),                intent(in)  :: th, lambda
       real(kind=wp), dimension(2),  intent(out) :: point
     end subroutine flatpoint
  end interface

  real(kind=wp), intent(in), dimension(:)   :: parvec
  real(kind=wp), intent(in), dimension(:)   :: uhmat, omM

  integer :: tfile1 = 31, tfile2 = 32, pfile1 = 41, smovf = 51, omfile = 40
  integer :: it, ix, iphi, nphi = 20, icoll
  integer, save :: times_called
  real(kind=wp), dimension(2) :: point
  real(kind=wp) :: t, th, uh, u0, om, maxom, phi, dphi

  times_called = times_called + 1
  it = times_called

  if(b < 1.0_wp) then
     icoll = 2
  else
     icoll = 1
  endif

  dphi = 4.0d0 * asin(1.0d0) / dble(nphi-1) ! 4 pi / dble(nphi-1)

  open (tfile1, file='Lomega.cmesh')
  open (tfile2, file='Romega.cmesh')
  open (pfile1, file='omega_movie.cmesh')
  open (smovf,  file='smov.sva')
  open (omfile, file='omega.dat')

  write(tfile1,*) '{CMESH'
  write(tfile2,*) '{CMESH'
  write(tfile1,*) n_x, n_u
  write(tfile2,*) n_x, n_u

  maxom = maxval(omM)

  t = 0.0_wp

  write(pfile1,*) '{CMESH'
  write(pfile1,*) nphi, n_x
  write(smovf,*) t

  do ix = 1, n_x
     th = parvec(ix)
     uh = uhmat(ix)
     u0 = uzero_t(th)
     uh = uh - u0
     call flatpoint(th,uh,point)
     om = omM(ix)/maxom
     write(omfile,*) it, ix, om
     if (point(icoll).ge.0.0_wp) then
        write(tfile1,*)  point(1), point(2), t, om,  0.5, 0.0, 1.0
        write(tfile2,*) -point(1), point(2), t, 0.0, 0.5, om,  1.0
     else
        write(tfile1,*)  point(1), point(2), t, 0.0, om, om, 1.0
        write(tfile2,*) -point(1), point(2), t, 0.0, om, om, 1.0
     endif

     phi = 0.0d0
     do iphi = 1, nphi
        write(pfile1,*) point(1)*cos(phi),point(1)*sin(phi),point(2), &
             om,  0.5, 0.0, 1.0 

        write(smovf,*) point(1)*cos(phi),point(1)*sin(phi),point(2)
        phi = phi + dphi
     end do
  end do
  write(pfile1,*) '}'

  write(tfile1,*) '}'
  write(tfile2,*) '}'
  close(tfile1)
  close(tfile2)
  close(pfile1)
  close(smovf)
  close(omfile)

  return
end subroutine plotne
