module  spheroid_init

contains 
subroutine spheroid_initial (nn, it, time, j, r, rdot, rl, omega, jl)

!-----------------------------------------------------------------------
! we need initial values of: 
!    omega, r_{,\lambda} and j_{,\lambda}
!-----------------------------------------------------------------------

!  use horizon_eth,    only : eth1

   use gridtranslator, only:  spheroid_rj, spheroid_rdot
   use model
   use horizon_eth
   implicit none

   integer,          intent (in) :: nn, it
   double precision, intent (in) :: time
   double precision, dimension (nn,nn,2) :: r, rdot, rl
   double complex,   dimension (nn,nn,2) :: j, jl, omega

   ! initialize r and rdot
   call spheroid_rj  (time, nn, r, j)
   call spheroid_rdot(nn, rdot)

   ! set the rest of the initial data like in the close approximation 
   ! we can set \omega = -1/2 \bar\eth j, assuming linearity
   call eth1(nn, omega, j, 2, -1)
   omega = -0.5 * omega

   !-----------------------------------------------------------------------
   ! initial values of \rho_{,\lambda}
   !-----------------------------------------------------------------------

   call model_rho_lambda (nn, time, rl)

   !-----------------------------------------------------------------------
   ! initial values of j_{,\lambda})
   !-----------------------------------------------------------------------

   call model_j_lambda (nn, time, jl)

return
end subroutine spheroid_initial

subroutine spheroid_initial_axi (nn, time, mass)
!-----------------------------------------------------------------------
! set initial values of: 
!    U,lambda, r_{,\lambda} and gamma_{,\lambda}
!-----------------------------------------------------------------------

! modules used
   use prec,           only: wp
   use axihorizon,     only: allocate_axivars,                              &
        &             mx, r, rdot, Cgamma, r2Cgamma, lngammadot, Cgammadot, &
        &             Ulambda, r2Ulambda, Delta, rlambda, gammalambda, gugg,&
        &             gamma, r2gamma, gammadot
   use gridtranslator, only: spheroid_RRdotCgamma, cos_theta, flipper
   use numservicef90,  only: vdiff
   use axisymmetricio, only: axi_outr

!   use model

   implicit none

! input variables
   real(kind=wp), intent (in) :: time, mass
   integer      , intent (in) :: nn

! output variables: none; indirect outut by resetting module variables

! local varaibles
   real(kind=wp)              :: xlo = -1.0, xhi = 1.0 !  -1 <= x <= 1
   integer                    :: debug = 0
   real(kind=wp), dimension(:), allocatable :: thetaJ, thetaJx

!< === executable statements === >!

!FIXME: this definition coincides with the one in gridtranslator, do only once!

write(*,*)
write(*,*) 'Initializing axisymmetric evolution at time = ', time
write(*,*)

mx = 2*nn + 1
allocate(thetaJ(mx), thetaJx(mx))
call allocate_axivars(mx)

call spheroid_RRdotCgamma(mass, r, rdot, Cgamma, r2Cgamma, lngammadot, &
                        & Cgammadot, gugg)

! set the initial data like in the close approximation 
thetaJ = 0.5_wp*(1.0_wp/Cgamma - Cgamma)*flipper
call vdiff(thetaJ, thetaJx, xlo, xhi, debug)
Ulambda(2:mx-1) = ( thetaJx(2:mx-1) &
    & - 2.0_wp*cos_theta(2:mx-1)*thetaJ(2:mx-1)/(1.0_wp-cos_theta(2:mx-1)**2) )
Ulambda(1)  = 2.0_wp*thetaJx(1)
Ulambda(mx) = 2.0_wp*thetaJx(mx)

r2Ulambda = r**2 * Ulambda

!rlambda = -time/(4.0_wp * mass)
rlambda = -time/(4.0_wp * mass) * ( r/(2.0_wp*mass))
Delta = time * (1.0_wp - ( r/(2.0_wp*mass))) * (1.0_wp + ( r/(2.0_wp*mass)))
gammalambda =  2.0_wp*r*rlambda + time

call axi_outr(mx, 0, 'delta_rlambda_initial', rlambda - gammalambda)
call axi_outr(mx, 0, 'Ulambda_initial', Ulambda)
call axi_outr(mx, 0, 'Delta_initial',   Delta)
call axi_outr(mx, 0, 'thetaJ_initial',  thetaJ)
call axi_outr(mx, 0, 'R_initial',       r)

gammalambda = 0.0_wp

gamma    = (1.0_wp - cos_theta**2) * Cgamma
r2gamma  = (1.0_wp - cos_theta**2) * r2Cgamma
gammadot = (1.0_wp - cos_theta**2) * Cgammadot
 
deallocate(thetaJ, thetaJx)
return
end subroutine spheroid_initial_axi

end module  spheroid_init
