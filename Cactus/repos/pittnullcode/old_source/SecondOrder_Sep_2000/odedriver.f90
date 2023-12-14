subroutine solve_focuseq(utime, it, nt, dt, debug)
!=========================================
!purpose: call rksuite solver of the focussing equation
!         call embedding routine

use prec
use affine, theta_affine => theta
use gridtranslator
use Model
use numservicef90
use rksuite_90, wp_rks => wp
implicit none

interface

 subroutine plotne(parvec, uhmat, omM)
   ! purpose:
   use affine
   real(kind=wp), intent(in), dimension(:)   :: parvec
   real(kind=wp), intent(in), dimension(:)   :: uhmat, omM
 end subroutine plotne

end interface

!local variables defining the numerical setup
integer                                 :: totf
real(kind=wp)                           :: t_start, tol=1.0e-14_wp
real(kind=wp), dimension(:), allocatable   :: y_start, thres, &
                                           y_got, y_maxvals
! local loop variables
real(kind=wp) :: t_want, t_got
integer       :: i

! local variables involved in checking for caustics
real(kind=wp)                                       :: d_th_max

! we do a 1d system
type(rk_comm_real_1d), save :: comm

! other local quantities
integer,       intent(in) :: it, nt, debug
real(kind=wp), intent(in) :: utime, dt

integer                   :: i, ierr
logical                   :: fulloutput = .false.
!! === executable statemets === !!

write(*,*)  'odedriver called with it = ', it

n_u    = nt !! was: + 2
t_inc  = dt

! get variables from module Model
b      = sqrt(1.0_wp + eps) ! b^2 - 1 = eps
t_end  = t_ + nt * dt

if (.NOT. allocated(uhatofu)) then
   if (it .NE. 0) STOP '"it" not zero @ initialization step'
   open (unit = 10, file = "spheroid.log")

   ! construct the grid, get n_x
   call  stereogrid_to_thetavect(n_x, debug)

   ! allocate arrays
   allocate(uhatofu(n_x), omvec(n_x), y_start(n_x), thres(n_x))
   
   ! initialize sigma, p_omega, etc.
   call init(theta)

   !< prepare the ODE solution >!
   thres   = 1.0e-24_wp
   t_start = t_minus
   y_start = u_zero + th_minus - t_zero ! th = th_minus - t_zero

   call setup(comm,t_start,y_start,t_end,tol,thres,method='M')

   uhatofu = y_start

   !< write configuration >!
   write(10, *)
   write(10, *) 'ODE tolerance   [tol]     = ', tol
   write(10, *) 'ODE threshold   [thres]   = ', thres
   write(10, *) 'ODE start time  [t_start] = ', t_start 
   write(10, *) 'ODE end time    [t_end]   = ', t_end 
   write(10, *) 'ODE output dt   [t_inc]   = ', t_inc
   write(10, *)

deallocate(y_start, thres)
endif
 
if (it > 0) then

!< solve the ODEs from it = 1, n_u >!

   t_want = t_end - (n_u - it)*t_inc
   write(10, *) '*** call range_integrate @ time lvl', it,'new time', &
        &                                                 t_want, utime
   call range_integrate(comm, udot, t_want, t_got, y_got=uhatofu)
   d_th_max = maxval(uhatofu + 0.5_wp*sigma) 
 
   if ( d_th_max < 0.0) then
      write(10,98) '       max(u^ + sigma/2) =', d_th_max
   else
      write(10,99) '       max(u^ + sigma/2) =', d_th_max, '  CAUSTIC!'
      write(*,  *) 'we went past a caustic at time step', it
   endif
   write(87, *) t_want, uhatofu(i_equator)
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
   allocate(y_maxvals(n_x))
   ! collect data, clean up
   call statistics(comm,y_maxvals=y_maxvals,total_f_calls=totf)
   !write (*,*) ' y_maxvals=', y_maxvals
   deallocate(y_maxvals)

   write (10,*)
   write (10,*) ' The cost of the integration in calls to RHS was', totf
   write (10,*)

   call collect_garbage(comm)
end if

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
use gridtranslator, ONLY: n_x, n_u

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
integer :: ix, it, iphi, nphi = 20, icoll
real(kind=wp), dimension(2) :: point
real(kind=wp) :: t, th, dt, uh, u0, om, maxom, phi, dphi

if(b < 1.0_wp) then
  icoll = 2
 else
  icoll = 1
endif

dt   = t_end/dble(n_u-1) 
dphi = 6.28/dble(nphi-1)

  open (tfile1, file='Lomega.cmesh')
  open (tfile2, file='Romega.cmesh')
  open (pfile1, file='omega_movie.cmesh')
  open (smovf,  file='smov.sva')
  open (omfile, file='omega.dat')

  write(tfile1,*) '{CMESH'
  write(tfile2,*) '{CMESH'
  write(tfile1,*) n_x, n_u
  write(tfile2,*) n_x, n_u

  maxom= maxval(omM)

  t = 0.0_wp
! the time loop started here previuosly
!  do it = 1, n_u

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
