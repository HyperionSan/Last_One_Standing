! WARNING  - When I switched the code to an it_start aware form I made some
! assumptions as to what the original author thought it_start should be.
! For example: I changed
!  if( mod(it -1 , ....
! to  if(mod(it - it_start   however the original author may have meant
!  if(mod(it - it_start -1    
! watch out for this
! in each case I kept the original statement as a comment
! yosef
module axihorizon

use prec
implicit none

integer, private :: index = 0
integer :: mx  ! number of grid-points - will usuall be mx = n_x
integer :: i   ! the omnipresent counter
! global arrays
real(kind=wp), dimension (:), allocatable, save, public :: &
     r_new,           r,           rx,        &
     rdot_new,        rdot,        rdotx,     &
     gamma_new,       gamma,       gammax,    &
     gammadot_new,    gammadot,               &
     Cgamma_new,      Cgamma,                 &
     Cgammadot_new,   Cgammadot,              &
     r2gamma_new,     r2gamma,                &
     r2Cgamma_new,    r2Cgamma,               &
     lngammadot_new,  lngammadot,             &
     gugg_new,        gugg,                   &
     rlambda_new,     rlambda,                &
     Ulambda_new,     Ulambda,     Ulambdax,  &
     r2Ulambda_new,   r2Ulambda,              &
     Delta_new,       Delta,                  &                   
     gammalambda_new, gammalambda,            &
     exp_axi,         uhat_old
contains

subroutine allocate_axivars(mx)

implicit none
integer, intent(in) :: mx

allocate(r_new(mx),           r(mx),           rx(mx),        &
         rdot_new(mx),        rdot(mx),        rdotx(mx),     &
         gamma_new(mx),       gamma(mx),       gammax(mx),    &
         gammadot_new(mx),    gammadot(mx),                   &
         Cgamma_new(mx),      Cgamma(mx),                     &
         Cgammadot_new(mx),   Cgammadot(mx),                  &
         r2gamma_new(mx),     r2gamma(mx),                    &
         r2Cgamma_new(mx),    r2Cgamma(mx),                   &
         lngammadot_new(mx),  lngammadot(mx),                 &
         gugg_new(mx),        gugg(mx),                       &
         rlambda_new(mx),     rlambda(mx),                    &
         Ulambda_new(mx),     Ulambda(mx),     Ulambdax(mx),  &
         r2Ulambda_new(mx),   r2Ulambda(mx),                  &
         Delta_new(mx),       Delta(mx),                      &
         gammalambda_new(mx), gammalambda(mx),                &
         exp_axi(mx),         uhat_old(mx)  )

end subroutine allocate_axivars


subroutine parity_check (index, time, vect, parity, filnam) 
!   
use prec,            only: wp
use numservicef90,   only: vdiff
use null_params,     only: iot2
use axisymmetricio,  only: axi_outr
implicit none

! input variables:
integer,       intent(in) :: index, parity
real(kind=wp), intent(in) :: time
character*(*), intent(in) :: filnam  
real(kind=wp), dimension(:), intent(in) :: vect

! local variables
integer                                  :: ieq, i
character*(128)                          :: filnam_local
real(kind=wp), dimension(:), allocatable :: halfvect
!<= executable statements =>!
write(*,*) 'parity check for  ', filnam   
filnam_local = filnam // '_paritycheck'

ieq = (size(vect)-1) / 2 + 1
allocate(halfvect(ieq))

do i = 1, ieq
   halfvect(i) = vect(i) - parity * vect(size(vect) + 1 - i) 
enddo

write(*,*) 'max violation = ', maxval(abs(halfvect)), 'at time =', time
write(*,*) 'max relative violation =', &
     &                         maxval(abs(halfvect))/maxval(abs(vect))
write(*,*)

call axi_outr(ieq, index, trim(filnam_local), halfvect/maxval(abs(vect)))

deallocate(halfvect)
end subroutine parity_check


subroutine step_ext_curv_axi(it, dt, time) 
! this was previously called step_ext_curv_axi_beta

use prec,            only: wp
use numservicef90,   only: vdiff, vdif2, pvdiff, pvdif2
use affine,          only: sigma
use gridtranslator,  only: spheroid_RRdotCgamma, cost, uhatofu
use null_params,     only: iot2, it_start
use axisymmetricio,  only: axi_outr 

implicit none

! input variables:
integer,       intent(in) :: it
real(kind=wp), intent(in) :: dt, time

! output variables: no direct output

! local scalars
real(kind=wp) :: xlo = -1.0_wp, xhi = 1.0_wp, time_new !  -1 <= x <= 1
integer       :: debug = 0

! local arrays
real(kind=wp), dimension (:), allocatable ::                                &
   dotUlambda_RHS, dotRlambda_RHS, dotgammalambda_RHS0, dotgammalambda_RHS, &
   rxxS, gammaxxS

real(kind=wp), dimension (:), allocatable  ::                          &
     rS, rxS, rdotS, rdotxS, gammaS, gammaxS, gammadotS,               &
     r2gammaS, lngammadotS, lngammadotxS, rlambdaS,                    &
     UlambdaS, UlambdaxS, r2UlambdaS, gammalambdaS, gammaxugammaS,     &
     guggS, tmparr, sin2theta

!< === executable statements === >!

time_new = time + dt

if (mx > 0) then
allocate(dotUlambda_RHS(mx), dotRlambda_RHS(mx), dotgammalambda_RHS0(mx), &
         dotgammalambda_RHS(mx), rxxS(mx), gammaxxS(mx))

allocate(rS(mx), rxS(mx), rdotS(mx), rdotxS(mx), gammaS(mx), gammaxS(mx), & 
      &  gammadotS(mx), r2gammaS(mx), rlambdaS(mx),                       &
      &  UlambdaS(mx), UlambdaxS(mx), r2UlambdaS(mx),                     &
      &  lngammadotS(mx), lngammadotxS(mx), gammalambdaS(mx),             &
      &  gammaxugammaS(mx), guggS(mx), tmparr(mx), sin2theta(mx))
else
   write(*,*) 'it seems that mx has not been set!'
   stop 'exiting - check code to set mx correctly'
endif

sin2theta     = 1.0_wp - cost**2
sin2theta(1)  = 0.0_wp
sin2theta(mx) = 0.0_wp

call spheroid_RRdotCgamma(r_new, rdot_new, Cgamma_new, r2Cgamma_new, &
                         lngammadot_new, Cgammadot_new, gugg_new)

gamma_new    = sin2theta * Cgamma_new
r2gamma_new  = sin2theta * r2Cgamma_new
gammadot_new = sin2theta * Cgammadot_new

! compute values at the halfstep
call midval(r,          r_new,          rS) 
call midval(rdot,       rdot_new,       rdotS)
call midval(gamma,      gamma_new,      gammaS)
call midval(gammadot,   gammadot_new,   gammadotS)
call midval(r2gamma,    r2gamma_new,    r2gammaS)
call midval(lngammadot, lngammadot_new, lngammadotS)
call midval(gugg,       gugg_new,       guggS)

call pvdiff(1, gammaS,      gammaxS,      xlo, xhi, debug)
call pvdif2(1, gammaS,      gammaxxS,     xlo, xhi, debug)
call pvdiff(1, rdotS,       rdotxS,       xlo, xhi, debug)
call pvdiff(1, rS,          rxS,          xlo, xhi, debug)
call pvdif2(1, rS,          rxxS,         xlo, xhi, debug)
call pvdiff(1, lngammadotS, lngammadotxS, xlo, xhi, debug)

gammaxugammaS = lngammadotxS + guggS*gammaxS ! gamma,xu / gamma

! do the hierarchy of RK2-steps

! for Ulambda
dotUlambda_RHS = -gammaxugammaS*rS**2 + 2.0_wp*rdotS*rxS &
     &  + rS*(-2.0_wp*lngammadotS*rxS - 2.0_wp*rdotxS) ! dot(r^2*Ulambda)

r2Ulambda_new   = r2Ulambda   + dotUlambda_RHS * dt
Ulambda_new     = r2Ulambda_new/r_new**2

! for rlambda
call midval(Ulambda, Ulambda_new, UlambdaS)
call vdiff(UlambdaS, UlambdaxS, xlo, xhi, debug)

tmparr = 1.0_wp + 0.5_wp*gammaxxS
dotRlambda_RHS = gammaxS*rxS/rS + rS*R2gammaS*rxxS  & ! dot(Delta) 
     &  + tmparr - R2gammaS*rxS**2 - 0.5_wp*gammaxS*UlambdaS &
     &  + gammaS*(0.25*UlambdaS**2 - 0.5_wp*UlambdaxS)

Delta_new       = Delta       + dotRlambda_RHS     * dt
rlambda_new     = (Delta_new - time_new)/(2.0_wp*r_new)
exp_axi         = rlambda_new/r_new

! for gammalambda
call midval(rlambda, rlambda_new, rlambdaS)

dotgammalambda_RHS = R2gammaS * (rS*(-(lngammadotS*rlambdaS)        &
     & - R2gammaS*rxS*UlambdaS) - &
     &   (R2gammaS*rS**2*(UlambdaS**2 - 2*UlambdaxS))/4. )

gammalambda_new = gammalambda + dotgammalambda_RHS * dt

if (mod(it - it_start,iot2) == 0) then

   gammaS = gamma_new
   rS = r_new
   rdotS = rdot_new
   lngammadotS = lngammadot_new
   
   call pvdiff(1, gammaS,      gammaxS,      xlo, xhi, debug)
   call pvdif2(1, gammaS,      gammaxxS,     xlo, xhi, debug)
   call pvdiff(1, rdotS,       rdotxS,       xlo, xhi, debug)
   call pvdiff(1, rS,          rxS,          xlo, xhi, debug)
   call pvdif2(1, rS,          rxxS,         xlo, xhi, debug)
   call pvdiff(1, lngammadotS, lngammadotxS, xlo, xhi, debug)
   !   call parity_check (index, time, dotUlambda_RHS, -1, 'dotUlambda_RHS') 
   !   call parity_check (index, time, gammaxugammaS, -1, 'gammaxugammaS')
   !   call parity_check (index, time, R2gammaS, 1, 'R2gammaS')
   call parity_check (index, time, rxS, -1, 'rxS')
   call parity_check (index, time, rxxS, 1, 'rxxS')
   !   call parity_check (index, time, guggS, 1, 'guggS')
   !   call parity_check (index, time, rdotS, 1, 'rdotS')
   !   call parity_check (index, time, sin2theta, 1, 'sin2theta')
   !   call parity_check (index, time, rS, 1, 'rS')
   !   call parity_check (index, time, gammaS, 1, 'gammaS')

   !   call parity_check (index, time, rxS * gammaxS, 1, 'prod')
   
   call parity_check (index, time, dotRlambda_RHS,    1, 'dotRlambda_RHS')
   call parity_check (index, time, UlambdaxS,    1, 'Ulambdax')
   
   call parity_check (index, time, gammalambda_new, 1, 'gammalambda')
   call parity_check (index, time, Ulambda_new,    -1, 'Ulambda')
   call parity_check (index, time, rlambda_new,     1, 'rlambda')
   call parity_check (index, time, uhat_old,     1, 'uhat')
   call parity_check (index, time, sigma,        1, 'sigma')

   call axi_outr(mx, index, 'uhat',     uhat_old)
   call axi_outr(mx, index, 'sigma',    sigma)
   call axi_outr(mx, index, 'dotRlambda_RHS',     dotRlambda_RHS)
   !   call axi_outr(mx, index, 'gammaxugamma', gammaxugammaS)
   !   call axi_outr(mx, index, 'R2gammaS',  R2gammaS)   
   call axi_outr(mx, index, 'rS', rS)
   call axi_outr(mx, index, 'rxS', rxS)
   call axi_outr(mx, index, 'rxxS', rxxS)
   call axi_outr(mx, index, 'gammaxS', gammaxS)
   call axi_outr(mx, index, 'gammaxxS',  gammaxxS)
   !   call axi_outr(mx, index, 'guggS',  guggS)
   call axi_outr(mx, index, 'rdotS',  rdotS)
   index = index + 1
endif

uhat_old = uhatofu

deallocate(dotUlambda_RHS, dotRlambda_RHS, dotgammalambda_RHS0,        &
           dotgammalambda_RHS, rxxS, gammaxxS)
deallocate(rS, rxS, rdotS, rdotxS, gammaS, gammaxS,                    &
           gammadotS, r2gammaS, lngammadotS, rlambdaS,                 &
           UlambdaS, UlambdaxS, r2UlambdaS, gammalambdaS,              &
           lngammadotxS, gammaxugammaS, guggS, tmparr, sin2theta)
return
end subroutine step_ext_curv_axi


subroutine interpolate_axi2stereo (it, time, mass) 
! purpose: interpolate the follwoing fields to the stereographic (p,q)-grid:
!          J, r, omega, rlambda  - probably more (everything needed in HtoR!
use prec,           only: wp
use null_grid,      only: q=>qs, p=>ps, pp
use horizon,        only: rhonew, rhodotnew, Jnew, rholnew, omeganew, Jlnew
use gridtranslator, only: uhatofu, cost, flipper
use numservicef90,  only: vdiff
implicit none

integer, intent(in)       :: it
real(kind=wp), intent(in) :: time, mass

real(kind=wp) :: xlo = -1.0, xhi = 1.0 !  -1 <= x <= 1
integer       :: debug = 0
!     VARIABLES for INTERPOLATION:
REAL(KIND=WP) :: WK(0:1)
INTEGER :: NCD, LWK, IENDC, IER
LOGICAL :: UNIFRM,PER

real(kind=wp),  dimension(:), allocatable :: SIGMA, yp
REAL(KIND=WP) :: HVAL ! evaluation function for interpolation 
real(kind=wp) :: x_eval, tmp
real(kind=wp), dimension(:), allocatable :: x, thetaJ, thetaJx, &
                  omsym, omsymx, rholsym, rholsymx, ro_new, rmlambda, &
                  thetaJl, thetaJlx, tmparr
integer          :: i, i1 ! loop counters
integer          :: nn, middle_nn
complex(kind=wp) :: phase
real(kind=wp)    :: rad, re, im, sin_phi, cos_phi
!< executable statements >!

if (debug > 0) print *, 'entered interpolate_axi2stereo at time step', it

allocate(x(mx), yp(mx), sigma(mx), thetaJ(mx), thetaJx(mx),             &
         omsym(mx), omsymx(mx), rholsym(mx), rholsymx(mx), ro_new(mx),  &
         rmlambda(mx), thetaJl(mx), thetaJlx(mx), tmparr(mx))

nn = size(rhonew, 1)

! interpolate for:
!   rhonew, rhodotnew, jnew,
!   omeganew, rholnew, jlnew

x = cost

NCD    = 1        ! Number of continuous derivatives at the knots.
PER    = .false.  ! periodic is false
UNIFRM = .true.   ! use uniform tension
IENDC  = 0
LWK    = 0

! r
YP     = 0.0d0
SIGMA  = 0.0d0
call TSPSI(mx, x, r_new, NCD, IENDC, PER, UNIFRM, LWK, &
         & WK, YP, SIGMA, IER)
if (IER .ne. 0) then 
   write(*,*) 'TSPSI exited with ifail ', IER, ' for r.'
endif
call vdiff(r_new, rx, xlo, xhi, debug) ! make sure the rx is taken @ right
                                        ! time step
! r_dot
YP     = 0.0d0
SIGMA  = 0.0d0
call TSPSI(mx, x, r_new, NCD, IENDC, PER, UNIFRM, LWK, &
         & WK, YP, SIGMA, IER)
if (IER .ne. 0) then
   write(*,*) 'TSPSI exited with ifail ', IER, ' for r_dot.'
endif
call vdiff(r_new, rx, xlo, xhi, debug) ! make sure the rx is taken @ right
                                        ! time step


! interpolate J
YP     = 0.0d0
SIGMA  = 0.0d0
thetaJ =  0.5_wp*(1.0_wp/Cgamma_new - Cgamma_new)*flipper
call TSPSI(mx, x, thetaJ, NCD, IENDC, PER, UNIFRM, LWK, &
         & WK, YP, SIGMA, IER)
if (IER .ne. 0) then 
   write(*,*) 'TSPSI exited with ifail ', IER, ' for J.'
endif
call vdiff(thetaJ, thetaJx, xlo, xhi, debug)


! interpolate omeganew
YP     = 0.0d0
SIGMA  = 0.0d0
omsym = -0.5_wp*sqrt((1.0_wp - x**2))*Ulambda_new ! =: omega
call TSPSI(mx, x, omsym, NCD, IENDC, PER, UNIFRM, LWK, &
         & WK, YP, SIGMA, IER)
if (IER .ne. 0) then 
   write(*,*) 'TSPSI exited with ifail ', IER, ' for omega.'
endif
call vdiff(omsym, omsymx, xlo, xhi, debug)

! interpolate rholnew
YP     = 0.0d0
SIGMA  = 0.0d0
ro_new = r_new/(2.0_wp*mass) ! Rm = 2*mass = Schwarzschild horizon radius
rmlambda = -time/(4.0_wp*mass)
rholsym = (rlambda_new - ro_new*rmlambda)/(2.0_wp*mass)  ! =: rho,lambda
call TSPSI(mx, x, rholsym, NCD, IENDC, PER, UNIFRM, LWK, &
         & WK, YP, SIGMA, IER)
if (IER .ne. 0) then 
   write(*,*) 'TSPSI exited with ifail ', IER, ' for rholnew.'
endif
call vdiff(rholsym, rholsymx, xlo, xhi, debug)

! interpolate Jl
YP     = 0.0d0
SIGMA  = 0.0d0
!thetaJl = 0.5_wp*(-Cgammanewl/Cgammanew**2 - Cgammanewl)
!        = - 0.5 * gammanewl*(1/Cgammanew**2 + 1)/sin2theta
thetaJl(2:mx-1) =  (-0.5_wp*gammalambda_new(2:mx-1) * &
     &    (1.0_wp/Cgamma_new(2:mx-1)**2 + 1.0_wp)   / &
     &    (1.0_wp - cost(2:mx-1)**2))*flipper
thetaJl(1)  = thetaJl(2)     *flipper  ! \__ FIXME
thetaJl(mx) = thetaJl(mx - 1)*flipper  ! /
call TSPSI(mx, x, thetaJl, NCD, IENDC, PER, UNIFRM, LWK, &
         & WK, YP, SIGMA, IER)
if (IER .ne. 0) then 
   write(*,*) 'TSPSI exited with ifail ', IER, ' for Jl.'
endif
call vdiff(thetaJl, thetaJlx, xlo, xhi, debug)


SIGMA = 0.0

if (modulo(nn,2) == 1) then
   middle_nn = (nn + 1)/2 
else
   middle_nn = -42 ! never reached!
endif

Do i1 = 1, nn
   Do i = 1, nn


      if ((i == middle_nn).and.(i1 == middle_nn)) then
         phase = 1.0_wp ! for spinweight 2
         rad   = 1.0_wp
         ! angle phi   = 0.0_wp
      else 
         phase = cmplx(q(i,i1), p(i,i1))/cmplx(q(i,i1), -p(i,i1))
         ! note: N/S phase difference comes from the different semantics of
         !       p/q w.r.t. theta ??
         re = q(i,i1) ! dble (phase)
         im = p(i,i1) ! aimag(phase)
         rad  = sqrt(q(i,i1)**2 + p(i,i1)**2) ! abs  (phase)
         sin_phi = im/rad  ! y/rad
         cos_phi = re/rad  ! x/rad
         if (i == middle_nn) then
            sin_phi = im/rad  ! y/rad
            cos_phi = re/rad  ! x/rad
         endif
         if (i1 == middle_nn) then
            sin_phi = im/rad  ! y/rad
            cos_phi = re/rad  ! x/rad
         endif
      endif

      ! North
      x_eval = (2.0_wp - pp(i,i1))/pp(i,i1)
      rhonew(i,i1,1) = HVAL(x_eval, mx, x, r_new, rx, SIGMA, IER)/(2.*mass)
      rhodotnew(i,i1,1) = 0.0  ! FIXME
      STOP ' rhodotnew currently filled with dummy values'
      tmp = HVAL(x_eval, mx, x, thetaJ, thetaJx, SIGMA, IER)
      Jnew(i,i1,1) = cmplx(tmp, 0, kind(phase)) * phase
      tmp = HVAL(x_eval, mx, x, omsym, omsymx, SIGMA, IER) 
      omeganew(i,i1,1) = cmplx(tmp, 0, kind(phase)) &
                       * cmplx(cos_phi, -sin_phi, kind(phase)) 
      rholnew(i,i1,1) = HVAL(x_eval, mx, x, rholsym, rholsymx, SIGMA, IER)
      tmp = HVAL(x_eval, mx, x, thetaJl, thetaJlx, SIGMA, IER)
      Jlnew(i,i1,1) = cmplx(tmp, 0, kind(phase)) * phase
      ! South
      x_eval = (pp(i,i1) - 2.0_wp)/pp(i,i1)
      rhonew(i,i1,2) = HVAL(x_eval, mx, x, r_new, rx, SIGMA, IER)/(2.*mass)
      rhodotnew(i,i1,1) = 0.0  ! FIXME
      tmp = HVAL(x_eval, mx, x, thetaJ, thetaJx, SIGMA, IER)
      Jnew(i,i1,2) = cmplx(tmp, 0, kind(phase)) * phase
      tmp = HVAL(x_eval, mx, x, omsym, omsymx, SIGMA, IER) 
      omeganew(i,i1,2) = cmplx(tmp, 0, kind(phase)) &
                       * cmplx(cos_phi, -sin_phi, kind(phase)) 
      rholnew(i,i1,2) = HVAL(x_eval, mx, x, rholsym, rholsymx, SIGMA, IER)
      tmp = HVAL(x_eval, mx, x, thetaJl, thetaJlx, SIGMA, IER)
      Jlnew(i,i1,2) = cmplx(tmp, 0, kind(phase)) * phase 
   End Do
End Do


deallocate(x, yp, sigma, thetaJ, thetaJx, omsym, omsymx, &
         & rholsym, rholsymx, ro_new, rmlambda, thetaJl, thetaJlx, tmparr)

return
end subroutine interpolate_axi2stereo



subroutine axishifter 
!-----------------------------------------------------------------------
! shift down time levels as needed
!-----------------------------------------------------------------------

   implicit none

   r           = r_new
   rdot        = rdot_new
   gamma       = gamma_new
   gammadot    = gammadot_new
   Cgamma      = Cgamma_new
   Cgammadot   = Cgammadot_new
   r2gamma     = r2gamma_new
   r2Cgamma    = r2Cgamma_new
   lngammadot  = lngammadot_new
   rlambda     = rlambda_new
   Ulambda     = Ulambda_new
   r2Ulambda   = r2Ulambda_new
   Delta       = Delta_new
   gammalambda = gammalambda_new

end subroutine axishifter

subroutine axishifter_up
!-----------------------------------------------------------------------
! shift up time levels, this is needed tp produce correct
! output at the initial step
!-----------------------------------------------------------------------

   implicit none

   r_new           = r
   rdot_new        = rdot
   gamma_new       = gamma
   gammadot_new    = gammadot
   Cgamma_new      = Cgamma
   Cgammadot_new   = Cgammadot
   r2gamma_new     = r2gamma
   r2Cgamma_new    = r2Cgamma
   lngammadot_new  = lngammadot
   rlambda_new     = rlambda
   Ulambda_new     = Ulambda
   r2Ulambda_new   = r2Ulambda
   Delta_new       = Delta
   gammalambda_new = gammalambda

end subroutine axishifter_up


subroutine midval(x,y,mid)
use prec,   only: wp

implicit none

real(kind=wp), dimension(:), intent(in)    :: x, y
real(kind=wp), dimension(:), intent(inout) :: mid

mid = 0.5_wp*(x + y) ! room for improvement here, use Jonathan's trick!

end subroutine midval

subroutine axio (it, time, nt, mass)

   use prec,           only: wp
   use null_params,    only: iot0, iot2, it_start
   use axisymmetricio, only: axi_outr 
   implicit none

   integer, intent (in) :: it, nt ! # time step, #total time steps
   double precision, intent (in) :: time, mass
   double precision :: avexp ! average of exp_out over the slice
   double precision, dimension(mx) :: tmparr

   !if (it == 1) then
   if (it == it_start) then
      Open ( unit = 37, file = "ax_expmin.dat",   status = "unknown" )
      Open ( unit = 38, file = "ax_rmin.dat",     status = "unknown" )
      Open ( unit = 39, file = "ax_av_exp.dat",   status = "unknown" )
   end if

   !if (mod(it - 1,iot2) == 0) then
   if (mod(it - it_start,iot2) == 0) then
      write(*,*)'write grid functions in axio @ time step', it, 'plot #', index

      call embed_intp(it, index, nt/iot2, 0)  ! debug level = 0
!      call embed_constx(it, index, nt/iot2, 0)  ! debug level = 0

      call axi_outr(mx, index, 'gamma',   gamma_new)
      tmparr = r_new/(2.*mass) - 1.0
      call axi_outr(mx, index, 'R',       tmparr)
      call axi_outr(mx, index, 'Ulambda', Ulambda_new)
      call axi_outr(mx, index, 'Rl',      rlambda_new)
      call axi_outr(mx, index, 'gammal',  gammalambda_new)
      call axi_outr(mx, index, 'exp_out', exp_axi)
   end if

   !if (mod(it - 1, iot0) == 0) then
   if (mod(it - it_start, iot0) == 0) then
     
     Write (unit = 37, fmt = '(3E17.8)') time, minval(exp_axi),          &
                                               maxval(exp_axi)
     Call Flush (37)
     
     Write (unit = 38, fmt = '(3E17.8)') time, minval(r_new/(2.*mass) ), &
                                               maxval(r_new/(2.*mass) )
     Call Flush (38)

     avexp = (sum(exp_axi * r_new**2)/size(exp_axi))/(-time/(2.0_wp))
     Write (unit = 39, fmt = '(2E17.8)') time, avexp
     Call Flush (39) 

   end if

end subroutine axio

subroutine embed(it, index, timesteps, debug)
!===========================================
! purpose: do embedding
use affine
use gridtranslator, only: uhatofu, omvec, cost, parvec => theta
use numservicef90,  only: vdiff, vintegrate, add2cmesh_2Dto3D, addslice2cmesh
use axisymmetricio, only: axi_outr 

implicit none

integer,       intent(in)                  :: it, index, timesteps, debug

real(kind=wp), allocatable, dimension(:)   :: sigs, xvec, yvec
real(kind=wp), allocatable, dimension(:)   :: allt, Rvec, Gvec, Bvec, Ovec
real(kind=wp), allocatable, dimension(:)   :: intvec, dX
real(kind=wp)   :: uh, om, temp
integer         :: i, embfile=97, nphi = 8
real            :: red, green, blue, opacity, ehtest, caustic

!! === executable statemets === !!

mx = size(cost)
! allocate arrays
allocate(intvec(mx), dX(mx))
allocate(sigs(mx), xvec(mx), yvec(mx) )
allocate(Ovec(mx), allt(mx), Rvec(mx), Gvec(mx), Bvec(mx))

! initialize arrays
xvec   = 0.0d0
yvec   = 0.0d0
intvec = 0.0d0

do i = 1, mx
   allt = dble(it)/dble(timesteps)
   uh = uhatofu(i)
   om = omvec(i)
   xvec(i) = - om*sint(i)*(uh + 0.5d0*sigma(i)) ! (-) for x > 0 @ EH 
end do

! differentiate
call vdiff(xvec, dX, theta(1), theta(mx), 0)   ! 0 -> no debugging

! compute integrand
do i = 1, mx
   uh   = uhatofu(i)
   om   = omvec  (i)
   temp = om**2*(uh - 0.5d0*sigma(i))**2 - dX(i)**2
   sigs(i)   = sign(1.0d0,temp)
   intvec(i) = sqrt(abs(temp))

   if(debug == 1) then
      write(96, 9995) theta(i), om
      write(95, 9995) theta(i), dX(i)
      write(94, 9995) theta(i), -(uh + 0.5_wp*sigma(i))
   endif
end do

! integrate
call vintegrate(intvec,yvec, theta(1), theta(mx), 0)
yvec = yvec - yvec((mx-1)/2) ! shift coord origin

! write 3D embedding data
do i = 1, mx 
   ehtest  = 0.35_wp * (1.0_wp - sign(1.0_wp, xvec(i)))
   caustic = 0.47_wp * (1.0_wp - sign(1.0_wp, xvec(i)))

   red   = 0.94_wp - ehtest
   green = 0.70_wp - caustic
   blue  = 0.15_wp + ehtest

   opacity = 1.0_wp
   if( xvec(i) < 0) opacity = 0.0_wp

   Rvec(i) = red
   Gvec(i) = green
   Bvec(i) = blue
   Ovec(i) = opacity
end do

call add2cmesh_2Dto3D(xvec, yvec, Rvec, Gvec, Bvec, Ovec,        &
        &             mx, nphi, embfile)

call addslice2cmesh(xvec, yvec, Rvec, Gvec, Bvec, Ovec,          &
        &           mx, index, timesteps, 'pairofpants.cmesh')

deallocate(xvec, yvec, allt, sigs, Rvec, Gvec, Bvec, Ovec)
deallocate(intvec, dX)

write(*, *) 'finished embedding'
write(10,*) 'finished embedding'

return
9995 format(2E13.5)
end subroutine embed


subroutine embed_constx(it, index, timesteps, debug)
!===================================================
! purpose: do embedding and more output

use affine
use gridtranslator, only: cost
use numservicef90,  only: vdiff, vintegrate, add2cmesh_2Dto3D, addslice2cmesh
use axisymmetricio, only: axi_outr 

implicit none

integer,       intent(in)                  :: it, index, timesteps, debug

real(kind=wp), allocatable, dimension(:)   :: sigs, xvec, yvec
real(kind=wp), allocatable, dimension(:)   :: Rvec, Gvec, Bvec, Ovec         
real(kind=wp), allocatable, dimension(:)   :: intvec, dX
integer         :: i, ierr, embfile=97, nphi = 16
real            :: red, green, blue, opacity, ehtest, caustic

!! === executable statemets === !!

! allocate arrays
mx = size(r_new)
allocate(intvec(mx), dX(mx))
allocate(sigs(mx), xvec(mx), yvec(mx) )
allocate(Ovec(mx), Rvec(mx), Gvec(mx), Bvec(mx))

! initialize arrays
xvec   = 0.0d0
yvec   = 0.0d0
intvec = 0.0d0

xvec = sqrt(gamma_new)*r_new

! differentiate
call vdiff(xvec, dX, cost(1), cost(mx), 0)   ! 0 -> no debugging

! compute integrand
intvec     = 1.0_wp/r2gamma_new - dX**2
intvec(1)  = 10.0_wp
intvec(mx) = 10.0_wp
sigs       = sign(1.0d0, intvec)
intvec     = sqrt(abs(intvec))

if(debug == 1) then
   do i = 1, mx
      write(95, 9995) theta(i), dX(i)
   end do
endif

! integrate
call vintegrate(intvec, yvec, cost(1), cost(mx), 0)
yvec = yvec - yvec((mx-1)/2) ! shift coord origin

! write 3D embedding data

do i = 1, mx 
   ehtest  = 0.35_wp * (1.0_wp - sign(1.0_wp, xvec(i)))
   caustic = 0.47_wp * (1.0_wp - sign(1.0_wp, xvec(i)))

   red   = 0.94_wp - ehtest
   green = 0.70_wp - caustic
   blue  = 0.15_wp + ehtest

   opacity = 1.0_wp
   if( xvec(i) < 0) opacity = 0.0_wp

   Rvec(i) = red
   Gvec(i) = green
   Bvec(i) = blue
   Ovec(i) = opacity
end do

! call axi_outr(mx, index, 'xvec',   xvec)
! call axi_outr(mx, index, 'yvec',   yvec)
! call axi_outr(mx, index, 'intvec',   intvec)
! call axi_outr(mx, index, 'dX',   dX)


!call add2cmesh_2Dto3D(xvec, yvec, Rvec, Gvec, Bvec, Ovec,        &
!        &             mx, nphi, embfile)

call addslice2cmesh(xvec, yvec, Rvec, Gvec, Bvec, Ovec,          &
        &           mx, index, timesteps, 'pairofpants.cmesh')

deallocate(xvec, yvec, sigs, Rvec, Gvec, Bvec, Ovec, intvec, dX)

return
9995 format(2E13.5)
end subroutine embed_constx

subroutine embed_intp(it, index, timesteps, debug)
!===========================================
! purpose: do embedding by interpolation from const x to const theta
use affine,         only: sigma
use gridtranslator, only: uhatofu, omvec, cost
use numservicef90,  only: vdiffx, vintegratex, add2cmesh_2Dto3D, addslice2cmesh
use axisymmetricio, only: axi_outr 

implicit none

integer,       intent(in)                  :: it, index, timesteps, debug

real(kind=wp), allocatable, dimension(:)   :: sigs, xvec, yvec
real(kind=wp), allocatable, dimension(:)   :: Rvec, Gvec, Bvec, Ovec         
real(kind=wp), allocatable, dimension(:)   :: temp, intvec, dX, theta
real(kind=wp)   :: uh, om, phi, dphi
integer         :: i, itime, ierr, embfile=97, iphi, nphi = 8
real            :: red, green, blue, opacity, ehtest, caustic


!! === executable statemets === !!

mx = size(cost)
! allocate arrays
allocate(intvec(mx), dX(mx), theta(mx), temp(mx), xvec(mx), yvec(mx), &
    &    Ovec(mx), Rvec(mx), Gvec(mx), Bvec(mx))

! initialize arrays
yvec   = 0.0d0

theta  = acos(cost) 
xvec   = -omvec*sin(theta)*(uhatofu + 0.5d0*sigma) ! (-) for x > 0 @ EH 

! differentiate
call vdiffx(xvec, dX, theta, 0)   ! 0 -> no debugging

! compute integrand
temp = omvec**2*(uhatofu - 0.5d0*sigma)**2 - dX**2
intvec = sqrt(abs(temp))

do i = 1, mx
   if(debug == 1) then
      write(96, 9995) theta(i), omvec(i)
      write(95, 9995) theta(i), dX(i)
      write(94, 9995) theta(i), -(uhatofu(i) + 0.5_wp*sigma(i))
   endif
end do

! integrate
call vintegratex(intvec, yvec, theta, 0)
yvec = yvec - yvec((mx-1)/2) ! shift coord origin

! write 3D embedding data
do i = 1, mx 
   ehtest  = 0.35_wp * (1.0_wp - sign(1.0_wp, xvec(i)))

   red   = 0.94_wp - ehtest
   green = 0.58_wp - ehtest
   blue  = 0.15_wp + ehtest

   opacity = 1.0_wp
   if (xvec(i) < 0)         opacity = 0.0_wp
!   if (rlambda_new(i) <= 0) green = 0.0

   Rvec(i) = red
   Gvec(i) = green
   Bvec(i) = blue
   Ovec(i) = opacity
end do

call add2cmesh_2Dto3D(xvec, yvec, Rvec, Gvec, Bvec, Ovec,         &
        &             mx, nphi, embfile)

call addslice2cmesh( xvec, yvec, Rvec, Gvec, Bvec, Ovec,          &
        &           mx, index, timesteps, 'pairofpants_XR.cmesh')

call addslice2cmesh(-xvec, yvec, Rvec, Gvec, Bvec, Ovec,          &
        &           mx, index, timesteps, 'pairofpants_XL.cmesh')

deallocate(intvec, dX, theta, temp, xvec, yvec, Ovec, Rvec, Gvec, Bvec)
return
9995 format(2E13.5)
end subroutine embed_intp


end module axihorizon
