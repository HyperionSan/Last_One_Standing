module horizon

   implicit none

   public horizon_setup, shifter, dot_rho, dot_omega, dot_rho_lambda, &
          dot_j_lambda, htor

   ! global arrays

   double precision, dimension (:,:,:), allocatable, save, public :: &
      rhonew,      rho,    &
      rholnew,     rhol,   &
      rhodotnew,   rhodot

   double precision, dimension (:,:,:), allocatable, save, public :: &
      e2beta, betanew, betarnew, vnew, vrnew, exp_out, beta_red

   double complex, dimension (:,:,:), allocatable, save, public :: &
      jnew,     j,     &
      jlnew,    jl,    &
      omeganew, omega, debug_omega, &
      cBnew,    cBrnew, &
      cKnew,    cKrnew, ethrhol, cnunew, cnurnew

   double complex, dimension (:,:,:), allocatable, save, public :: &
      jrnew, qnew, unew, urnew

   double precision, dimension (:,:,:), allocatable, save, public :: &
      ricci, ricciterm, omegaterm, domegaterm

   double complex, dimension (:,:,:), allocatable, save, public :: &
      shears

   double precision, dimension (:,:,:), allocatable, save, public :: &
      shearqq, shearqp, shearpp, shearplus, twist_theta

   double precision, dimension (:,:,:), allocatable, save, public :: &
      knew, krnew

contains

subroutine horizon_setup (nn)
!-----------------------------------------------------------------------
! allocate all the global horizon arrays.
!-----------------------------------------------------------------------

   implicit none
   integer, intent (in) :: nn

   allocate ( rhonew    (nn,nn,2), rho   (nn,nn,2), &
              rholnew   (nn,nn,2), rhol  (nn,nn,2), &
              rhodotnew (nn,nn,2), rhodot(nn,nn,2), &
              jnew      (nn,nn,2), j     (nn,nn,2), &
              jlnew     (nn,nn,2), jl    (nn,nn,2), &
              omeganew  (nn,nn,2), omega (nn,nn,2), debug_omega(nn,nn,2), &
              e2beta    (nn,nn,2), &
              betanew   (nn,nn,2), &
              betarnew  (nn,nn,2), &
              vnew      (nn,nn,2), &
              jrnew     (nn,nn,2), &
              qnew      (nn,nn,2), &
              unew      (nn,nn,2), &
              urnew     (nn,nn,2), &
              vrnew     (nn,nn,2), &
              exp_out   (nn,nn,2), &
              beta_red  (nn,nn,2), &
              knew      (nn,nn,2), &
              krnew     (nn,nn,2), &
              cBnew     (nn,nn,2), &
              cBrnew    (nn,nn,2), &
              cKnew     (nn,nn,2), &
              cKrnew    (nn,nn,2), &
              cnunew    (nn,nn,2), &
              cnurnew   (nn,nn,2), &
              ethrhol   (nn,nn,2) )

   ! to understand the business of rho,lambda going to zero

   allocate ( ricci (nn,nn,2),      ricciterm (nn,nn,2), &
              omegaterm (nn,nn,2), domegaterm (nn,nn,2) )

   allocate ( shears(nn,nn,2), shearqq(nn,nn,2), shearqp(nn,nn,2), &
              shearpp(nn,nn,2), shearplus(nn,nn,2), twist_theta(nn,nn,2) )

end subroutine horizon_setup

subroutine shifter 
!-----------------------------------------------------------------------
! shift time levels as needed
!-----------------------------------------------------------------------

   implicit none

   rho    = rhonew
   rhol   = rholnew
   rhodot = rhodotnew
   j      = jnew
   jl     = jlnew
   omega  = omeganew

end subroutine shifter

subroutine dot_rho (nn, dt)
!-----------------------------------------------------------------------
! integrate $\ddot \dot \rho = -\rho/4 (\dot j \bar{\dot j} - (\dot k)^2)$
!-----------------------------------------------------------------------

   implicit none

   integer,          intent (in) :: nn
   double precision, intent (in) :: dt

   ! local arrays

   double precision, dimension (nn,nn,2) :: k, kdot, source
   double complex,   dimension (nn,nn,2) :: jdot

   print *, "in dot_rho"

   k = sqrt(1. + j * conjg(j))
   knew = sqrt(1. + jnew * conjg(jnew))
   jdot = (jnew - j) / dt
   kdot = (knew - k) / dt

   source = (jdot * conjg(jdot) - kdot ** 2) / 16.0d0

   rhonew = (rho * (1. - source * dt ** 2) + rhodot * dt ) &
        / (1. + source * dt ** 2)

   ! note: there is a sign error in Phys. Rev. D 64, 024010 (2001)
   ! eqn A4b. the following form is correct.
   rhodotnew = (rhodot * (1. - source * dt ** 2) - 4. * source * rho * dt ) &
        / (1. + source * dt ** 2)

end subroutine dot_rho

subroutine dot_omega (nn, dt)
!-----------------------------------------------------------------------
! integrate $\dot (r^2 \omega) = r^2 \eth (\dot r/r) + \ldots $
!-----------------------------------------------------------------------

   use horizon_eth
   use axisymmetricio, only: axi_outr
   implicit none

   integer,          intent (in) :: nn
   double precision, intent (in) :: dt

   ! local arrays

   double precision, dimension (nn,nn,2) :: kh, kdoth, rhoh, rhodoth

   double complex, dimension (nn,nn,2) :: jh, ckh, crhoh, jdoth, &
      ckdoth, crhodoth, eth_j, ethb_j, eth_k, eth_rho, ethb_rho, eth_jdot, &
      ethb_jdot, eth_kdot, ethb_kdot, eth_rhodot, part0, part1, part2, rhs

   jh = 0.5 * (j + jnew)

   kh = sqrt(1. + jh * conjg(jh))
   rhoh = 0.5 * (rho + rhonew)

   jdoth = (jnew - j) / dt
   kdoth = dble(jdoth * conjg(jh)) / kh
   rhodoth = (rhonew - rho) / dt

   crhoh = rhoh
   ckh = kh
   crhodoth = rhodoth
   ckdoth = kdoth

   call eth1(nn, eth_j     , jh      , 2,  1)
   call eth1(nn, ethb_j    , jh      , 2, -1)
   call eth1(nn, eth_k     , ckh     , 0,  1)
   call eth1(nn, eth_rho   , crhoh   , 0,  1)
   call eth1(nn, ethb_rho  , crhoh   , 0, -1)
   call eth1(nn, eth_jdot  , jdoth   , 2,  1)
   call eth1(nn, ethb_jdot , jdoth   , 2, -1)
   call eth1(nn, eth_kdot  , ckdoth  , 0,  1)
   call eth1(nn, ethb_kdot , ckdoth  , 0, -1)
   call eth1(nn, eth_rhodot, crhodoth, 0,  1)

   part0 = - rhodoth * eth_rho

   part1 = eth_rhodot &
           + (jdoth * conjg(jh) - kdoth * kh) * eth_rho &
           + (jh * kdoth - kh * jdoth) * ethb_rho

   part2 = conjg(jdoth) * eth_j - 4. * kdoth * eth_k &
         + 2. * kdoth * ethb_j + 3. * jdoth * conjg(ethb_j) &
         - 2. * jdoth * conjg(eth_k) + 2. * jh * conjg(eth_kdot) &
         + 2. * conjg(jh) * eth_jdot - 2. * kh * ethb_jdot &
         - 2. * kh * eth_kdot

   rhs = (0.25 * part2 * rhoh + part1) * rhoh + part0
     
   omeganew = (omega * rho ** 2 + rhs * dt) / rhonew ** 2

   debug_omega = -0.5 * ethb_j

!if (modulo(nn,2) == 1) then
! call axi_outr(nn, 0, 'omega_new_abs', abs(omeganew((nn+1)/2,:,1)))
! call axi_outr(nn, 0, 'omega_rhs_abs', abs(rhs((nn+1)/2,:,1)))
! call axi_outr(nn, 0, 'omega_rhs_im', dimag(rhs((nn+1)/2,:,1)))
!else
!   STOP 'we need nn odd for this run!'
!endif

end subroutine dot_omega

subroutine dot_rho_lambda (nn, dt, time, mass)
!-----------------------------------------------------------------------
! integrates \dot r^{2}_{,\lambda} = ... on the horizon
!-----------------------------------------------------------------------

   use horizon_eth
   implicit none

   integer,          intent (in) :: nn
   double precision, intent (in) :: dt, time, mass

   ! local arrays

   double precision, dimension (nn,nn,2) :: rhoh, rhodoth, kh, rhs

   double complex, dimension (nn,nn,2) :: jh, ckh, crhoh, omegah, &
      eth_j, ethb_j, ethb2_j, ethb_k, eth_ethb_k, eth_rho, eth2_rho, &
      ethb_eth_rho, eth_omega, ethb_omega

   ! local variables

   double precision :: th, r_m

   th = time - 0.5 * dt
   r_m = 2. * mass

   rhoh = 0.5 * (rho + rhonew)
   rhodoth = (rhonew - rho) / dt
   crhoh = rhoh
   jh = 0.5 * (j + jnew)
   kh = sqrt(1. + jh * conjg(jh))
   ckh = kh
   omegah = 0.5 * (omega + omeganew)

   call eth1 (nn, eth_j     , jh    , 2,  1)
   call eth1 (nn, ethb_j    , jh    , 2, -1)
   call eth1 (nn, ethb_k    , ckh   , 0, -1)
   call eth1 (nn, eth_omega , omegah, 1,  1)
   call eth1 (nn, ethb_omega, omegah, 1, -1)
   call eth1 (nn, eth_rho   , crhoh , 0,  1)

   call eth2 (nn, ethb2_j     , jh   , 2, -1, -1)
   call eth2 (nn, eth_ethb_k  , ckh  , 0,  1, -1)
   call eth2 (nn, eth2_rho    , crhoh, 0,  1,  1)
   call eth2 (nn, ethb_eth_rho, crhoh, 0, -1,  1)

   ricci = 2. * ckh + ethb2_j - eth_ethb_k &
         + (dble(eth_j) ** 2 + dimag(eth_j) ** 2 &
         - dble(ethb_j) ** 2 - dimag(ethb_j) ** 2) / (4. * kh)

   rhs = dble((ethb_k - conjg(ethb_j)) * (eth_rho / rhoh + omegah) &
   + kh * (  ethb_omega + ethb_eth_rho / rhoh &
              - eth_rho * conjg(eth_rho) / rhoh ** 2 ) &
   + conjg(jh) * (  eth_omega - omegah ** 2 &
                     + eth2_rho / rhoh - (eth_rho / rhoh) ** 2 ) &
   + kh * omegah * conjg(omegah)) - 0.5 * ricci + rhoh ** 2

   rholnew =  (  dt * rhs / r_m ** 2                                       &
               + rhol * (2. - rhodoth * dt + dt * th * rhoh / r_m ** 2)  ) &
             /          (2. + rhodoth * dt - dt * th * rhoh / r_m ** 2)

end subroutine dot_rho_lambda

subroutine dot_rho_lambda_RK2 (nn, dt, time, mass)
!-----------------------------------------------------------------------
! integrates \dot r^{2}_{,\lambda} = ... on the horizon
!  this version uses the standard second order Runge-Kutta rule,
!  see e.g. Numerical Recipes Eq. 16.1.2
!-----------------------------------------------------------------------

   use horizon_eth
   use gridtranslator, only:  spheroid_rdot
   use model, only : model_switch
   implicit none

   integer,          intent (in) :: nn
   double precision, intent (in) :: dt, time, mass

   ! local arrays
   double precision, dimension (nn,nn,2) :: rhoh, rhodoth, kh, rhs,        &
        &                                   rhol_dot_rhs, rhol_mid, k1, k2
   double complex,   dimension (nn,nn,2) :: jh, ckh, crhoh, omegah,        &
        eth_j, ethb_j, ethb2_j, ethb_k, eth_ethb_k, eth_rho, eth2_rho,     & 
        ethb_eth_rho, eth_omega, ethb_omega

   ! local scalar variables
   double precision :: th, r_m
   double precision :: sig     !  sig = sign(cfl) = sign(dt)

   !< executable statements >!

   sig = 1.0d0 ! sign(1.0d0, dt) ! forward or backward in time

   ! compute k1
   th = time - dt     ! time where we already know data
   r_m = 2. * mass

   rhoh = rho                      ! the suffix h was for "halfstep"
   crhoh = rhoh                    ! now let's think of it as "here" maybe?

   rhodoth = rhodot

   jh = j
   kh = sqrt(1. + jh * conjg(jh))
   ckh = kh
   omegah = omega

   call eth1 (nn, eth_j     , jh    , 2,  1)
   call eth1 (nn, ethb_j    , jh    , 2, -1)
   call eth1 (nn, ethb_k    , ckh   , 0, -1)
   call eth1 (nn, eth_omega , omegah, 1,  1)
   call eth1 (nn, ethb_omega, omegah, 1, -1)
   call eth1 (nn, eth_rho   , crhoh , 0,  1)

   call eth2 (nn, ethb2_j     , jh   , 2, -1, -1)
   call eth2 (nn, eth_ethb_k  , ckh  , 0,  1, -1)
   call eth2 (nn, eth2_rho    , crhoh, 0,  1,  1)
   call eth2 (nn, ethb_eth_rho, crhoh, 0, -1,  1)

   ricci = 2. * ckh + ethb2_j - eth_ethb_k &
         + (dble(eth_j) ** 2 + dimag(eth_j) ** 2 &
         - dble(ethb_j) ** 2 - dimag(ethb_j) ** 2) / (4. * kh)

   rhs = dble((ethb_k - conjg(ethb_j)) * (eth_rho / rhoh + omegah) &
   + kh * (  ethb_omega + ethb_eth_rho / rhoh &
              - eth_rho * conjg(eth_rho) / rhoh ** 2 ) &
   - conjg(jh) * (  eth_omega + omegah ** 2 &
                     + eth2_rho / rhoh - (eth_rho / rhoh) ** 2 ) &
   + kh * omegah * conjg(omegah)) - 0.5 * ricci + rhoh ** 2

   rhol_dot_rhs = sig * rhs / (rhoh * 2.*r_m ** 2)               &
                 + rhodoth * (sig * th / r_m ** 2 - rhol/rhoh)
   k1 = dt * rhol_dot_rhs

   ! compute k2
   th = time - 0.5 * dt
   r_m = 2. * mass

   rhoh = 0.5 * (rho + rhonew)
   rhodoth = (rhonew - rho) / dt  ! rhodot at the halfstep
   crhoh = rhoh
   jh = 0.5 * (j + jnew)
   kh = sqrt(1. + jh * conjg(jh))
   ckh = kh
   omegah = 0.5 * (omega + omeganew)

   call eth1 (nn, eth_j     , jh    , 2,  1)
   call eth1 (nn, ethb_j    , jh    , 2, -1)
   call eth1 (nn, ethb_k    , ckh   , 0, -1)
   call eth1 (nn, eth_omega , omegah, 1,  1)
   call eth1 (nn, ethb_omega, omegah, 1, -1)
   call eth1 (nn, eth_rho   , crhoh , 0,  1)

   call eth2 (nn, ethb2_j     , jh   , 2, -1, -1)
   call eth2 (nn, eth_ethb_k  , ckh  , 0,  1, -1)
   call eth2 (nn, eth2_rho    , crhoh, 0,  1,  1)
   call eth2 (nn, ethb_eth_rho, crhoh, 0, -1,  1)

   ricci = 2. * ckh + ethb2_j - eth_ethb_k &
         + (dble(eth_j) ** 2 + dimag(eth_j) ** 2 &
         - dble(ethb_j) ** 2 - dimag(ethb_j) ** 2) / (4. * kh)

   rhs = dble((ethb_k - conjg(ethb_j)) * (eth_rho / rhoh + omegah) &
   + kh * (  ethb_omega + ethb_eth_rho / rhoh &
              - eth_rho * conjg(eth_rho) / rhoh ** 2 ) &
   - conjg(jh) * (  eth_omega + omegah ** 2 &
                     + eth2_rho / rhoh - (eth_rho / rhoh) ** 2 ) &
   + kh * omegah * conjg(omegah)) - 0.5 * ricci + rhoh ** 2

   rhol_mid = rhol + 0.5 * k1 
   rhol_dot_rhs =   sig * rhs / (rhoh * 2.*r_m ** 2)                       &
                 + rhodoth * (sig * th / r_m ** 2 - rhol_mid/rhoh)
   k2 = dt * rhol_dot_rhs

   rholnew =  rhol + k2 ! + O(dt^3)

   ! diagnostics 

   ricciterm = -0.5 * ricci + rhoh ** 2

   omegaterm = dble(-jh * omegah ** 2 + kh * omegah * conjg(omegah))

   domegaterm = dble((ethb_k - conjg(ethb_j)) * omegah &
   + kh * ethb_omega - conjg(jh) * eth_omega )

end subroutine dot_rho_lambda_RK2

subroutine htor (nn, time, mass, debug)
!-----------------------------------------------------------------------
!  convert from horizon to code variables
!-----------------------------------------------------------------------

   use horizon_eth
   implicit none

   integer,          intent (in) :: nn, debug
   double precision, intent (in) :: time, mass

   !-----------------------------------------------------------------------
   ! local variables
   !-----------------------------------------------------------------------

   double precision :: r_m, r_ml, r_mu, r_mlu, min_rl

   !-----------------------------------------------------------------------
   ! local arrays
   !-----------------------------------------------------------------------

   double precision, dimension (nn,nn,2) :: rholu, rlu, rhs, r_l

   double complex, dimension (nn,nn,2) :: ckh, crhoh, crhol,              &
      eth_j, ethb_j, ethb2_j, ethb_k, eth_ethb_k, eth_rho, eth2_rho,      &
      ethb_eth_rho, ethb_rho, eth_rhol, ethb_rhol, eth_omega, ethb_omega, &
      eth_beta, cbetanew, cbetarnew

   double precision :: sig     !  sig = sign(cfl) = sign(dt)

   !< executable statements >!

   sig = 1.0d0 ! sign(1.0d0, dt) ! forward or backward in time - enters eqs!

   if (debug > 0) then
      print *, "time in htor = ", time
   endif
   !-----------------------------------------------------------------------
   ! factor asymptotic dependence of rho, store j, k, and needed eth opers.
   !-----------------------------------------------------------------------

   r_m = 2. * mass
   r_ml = - sig * time / (4. * mass)
   r_mu = 0.
   r_mlu = - sig / (4. * mass)

   knew = sqrt(1. + jnew * conjg(jnew))

   call eth1 (nn, cKnew, dcmplx(knew), 0, 1)

   crhoh = rhonew
   crhol = (sig * rholnew)

   call eth1 (nn, eth_rho ,  crhoh, 0,  1)
   call eth1 (nn, eth_rhol,  crhol, 0,  1)
   
   ethb_rho  = conjg(eth_rho)
   ethb_rhol = conjg(eth_rhol)
 
   ethrhol = eth_rhol
   !-----------------------------------------------------------------------
   ! Compute the outward wavefront expansion on the horizon
   !----------------------------------------------------------------------- 
!XXX
   r_l = r_M * (sig * rholnew) + r_Ml * rhonew    ! R_{,\lambda}
   exp_out = r_l / (r_M * rhonew)

   !-----------------------------------------------------------------------
   ! Compute the betas (full and reduced) on the horizon
   !-----------------------------------------------------------------------

   e2beta   = 1. / r_l                         ! \__  the full beta
   betanew  = 0.5 * sig * log(abs(e2beta))     ! / 

   beta_red = betanew + 0.5*sig*log(abs(r_Ml)) ! this is beta_R = reduced beta 

   cbetanew = betanew 
   call eth1 (nn, cBnew, cbetanew, 0, 1)

   min_rl = min( minval(r_l(:,:,1)), minval(r_l(:,:,2)) )
   if (min_rl <= 0.0) then
      write(10,*) 'BH ?! min(R,lambda) =', min_rl, '@ time', time
!      call flush(10)
   endif
!XXX
   !-----------------------------------------------------------------------
   ! compute j_{,r} and k_{,r} on the horizon
   !-----------------------------------------------------------------------

   jrnew = e2beta * jlnew

   krnew = (jrnew * conjg(jnew) + jnew * conjg(jrnew)) / (2. * knew)

   call eth1 (nn, cKrnew, dcmplx(krnew), 0, 1)

!! my 'nu' stuff

   call eth1(nn, cnunew, jnew, 2, -1)
   call eth1(nn, cnurnew, jrnew, 2, -1)

!! end my 'nu' stuff 
   !-----------------------------------------------------------------------
   ! compute beta_{,r} on the horizon
   !-----------------------------------------------------------------------

   betarnew = (r_m * rhonew) / 8.0d0 * (jrnew * conjg(jrnew) - krnew ** 2)

   cbetarnew = betarnew
   call eth1 (nn, cBrnew, cbetarnew, 0, 1)

   !-----------------------------------------------------------------------
   ! compute u on the horizon
   !-----------------------------------------------------------------------

   unew = -e2beta * (- jnew * ethb_rho + knew * eth_rho) / (r_m * rhonew ** 2)

   !-----------------------------------------------------------------------
   ! compute q on the horizon
   !-----------------------------------------------------------------------

   eth_beta = - 0.5 * e2beta * (r_m * eth_rhol + r_ml * eth_rho)
   qnew = - 2. * omeganew &
          + 2. * eth_beta &
          + 2. * eth_rho / rhonew &
          + r_m * (   jrnew * (knew * ethb_rho - conjg(jnew) * eth_rho) &
                    + krnew * (knew * eth_rho - jnew * ethb_rho) )

   !-----------------------------------------------------------------------
   ! compute u_{,r} on the horizon
   !-----------------------------------------------------------------------

   urnew = e2beta * (knew * qnew - jnew * conjg(qnew)) &
         / (r_m ** 2 * rhonew ** 2)

   !-----------------------------------------------------------------------
   ! compute v on the horizon
   !-----------------------------------------------------------------------

   vnew = - 2. * r_m ** 2 * rhonew * rhodotnew &
        + e2beta * r_m / (2. * rhonew) &
        * ( - conjg(jnew) * eth_rho ** 2 - jnew * ethb_rho ** 2 &
            + 2. * knew * eth_rho * ethb_rho )

   !-----------------------------------------------------------------------
   ! compute \rho_{,\lambda u} on the horizon
   !-----------------------------------------------------------------------

   ckh = knew

   call eth1 (nn, eth_j     , jnew    , 2,  1)
   call eth1 (nn, ethb_j    , jnew    , 2, -1)
   call eth1 (nn, ethb_k    , ckh     , 0, -1)
   call eth1 (nn, eth_omega , omeganew, 1,  1)
   call eth1 (nn, ethb_omega, omeganew, 1, -1)
   call eth1 (nn, eth_rho   , crhoh   , 0,  1)

   call eth2 (nn, ethb2_j     , jnew , 2, -1, -1)
   call eth2 (nn, eth_ethb_k  , ckh  , 0,  1, -1)
   call eth2 (nn, eth2_rho    , crhoh, 0,  1,  1)
   call eth2 (nn, ethb_eth_rho, crhoh, 0, -1,  1)

   ricci = 2. * knew + ethb2_j - eth_ethb_k &
         + (dble(eth_j) ** 2 + dimag(eth_j) ** 2 &
            - dble(ethb_j) ** 2 - dimag(ethb_j) ** 2) / (4. * knew)

   rhs = dble((ethb_k - conjg(ethb_j)) * (eth_rho / rhonew + omeganew) &
   + knew * (  ethb_omega + ethb_eth_rho / rhonew &
              - eth_rho * conjg(eth_rho) / rhonew ** 2 ) &
   + conjg(jnew) * (  eth_omega - omeganew ** 2 &
                     + eth2_rho / rhonew - (eth_rho / rhonew) ** 2 ) &
   + knew * omeganew * conjg(omeganew)) - 0.5 * ricci + rhonew ** 2

   rholu = sig * rhs / (8. * mass ** 2) &
         - (rhodotnew - sig * time * rhonew / (4. * mass ** 2)) * rholnew

   rlu = r_mlu * rhonew + r_ml * rhodotnew + r_mu * (sig*rholnew)+r_m * rholu
   rlu = sig * rlu
    ! another overall sign ???
   !-----------------------------------------------------------------------
   ! compute v,r on the horizon
   !-----------------------------------------------------------------------

   vrnew = vnew * (2. * betarnew + 1. / (r_m * rhonew)) &
         + 4. * r_m ** 2 * rhonew * rhodotnew * betarnew &
         - 2. * r_m * rhonew * rlu * e2beta &
         + 2. * e2beta / rhonew * ( - conjg(jnew) * omeganew * eth_rho &
                                    - jnew * conjg(omeganew) * ethb_rho &
                                    + knew * (conjg(omeganew) * eth_rho &
                                            + omeganew * ethb_rho) ) &
         + e2beta ** 2 * r_m / rhonew * ( - conjg(jnew) * eth_rhol * eth_rho &
                                          - jnew * ethb_rhol * ethb_rho &
                                          + knew * ethb_rhol * eth_rho &
                                          + knew * eth_rhol * ethb_rho ) &
         + e2beta * ( conjg(jnew) * eth_rho ** 2 + jnew * ethb_rho ** 2 &
                      - 2. * knew * eth_rho * ethb_rho ) &
         + e2beta * r_m / (2. * rhonew) * ( - conjg(jrnew) * eth_rho ** 2 &
                                            - jrnew * ethb_rho ** 2 &
                                            + 2. * krnew * eth_rho * ethb_rho)

   call shearntwist (nn, mass)


  e2beta = e2beta * sig

end subroutine htor

subroutine schwarzschild (nn, time, mass, debug)
!-----------------------------------------------------------------------
!  supply Bondi variables when the horizon is Schwarzschild
!-----------------------------------------------------------------------

   use horizon_eth
   use particle, only : U_2, bh_mass
   use null_grid, only: z, pp
   implicit none

   integer,          intent (in) :: nn, debug
   double precision, intent (in) :: time, mass

   !-----------------------------------------------------------------------
   ! local arrays
   !-----------------------------------------------------------------------

   double complex, dimension (nn,nn,2) :: cbetanew, cbetarnew
   !< executable statements >!

   if (debug > 0) then
      print *, "time in schwarzschild = ", time
   endif
   !-----------------------------------------------------------------------
   ! factor asymptotic dependence of rho, store j, k, and needed eth opers.
   !-----------------------------------------------------------------------

   knew = sqrt(1. + jnew * conjg(jnew))

   call eth1 (nn, cKnew, dcmplx(knew), 0, 1)

   !-----------------------------------------------------------------------
   ! Compute the betas (full and reduced) on the horizon
   !-----------------------------------------------------------------------

   e2beta   = 1.0d0
   betanew  = 0.0d0

! beta_red = betanew + 0.5*sig*log(abs(r_Ml)) ! this is beta_R = reduced beta

   cbetanew = betanew
   call eth1 (nn, cBnew, cbetanew, 0, 1)

   !-----------------------------------------------------------------------
   ! compute j_{,r} and k_{,r} on the horizon
   !-----------------------------------------------------------------------

   jrnew = (0.0d0, 0.0d0)

   krnew = 0.0d0

   call eth1 (nn, cKrnew, dcmplx(krnew), 0, 1)


   !-----------------------------------------------------------------------
   ! compute beta_{,r} on the horizon
   !-----------------------------------------------------------------------

   betarnew = 0.0d0

   cbetarnew = betarnew
   call eth1 (nn, cBrnew, cbetarnew, 0, 1)

   !-----------------------------------------------------------------------
   ! compute u on the horizon
   !-----------------------------------------------------------------------

!   unew = (0.0d0, 0.0d0)
   unew(:,:,1) = U_2*(1.+z**2)/pp
   unew(:,:,2) = - unew(:,:,1)

   !-----------------------------------------------------------------------
   ! compute q on the horizon
   !-----------------------------------------------------------------------

   qnew = (0.0d0, 0.0d0)

   !-----------------------------------------------------------------------
   ! compute u_{,r} on the horizon
   !-----------------------------------------------------------------------

   urnew = (0.0d0, 0.0d0)

   !-----------------------------------------------------------------------
   ! compute v on the horizon
   !-----------------------------------------------------------------------

! Here, V is the Bondi variable, with value r in the Minkowski case and
! r-2M in the Schwarzschild case. N.B. - the boundary is hard-coded at r=2.

   vnew = 2.0d0*(1.0d0-bh_mass)

   !-----------------------------------------------------------------------
   ! compute v,r on the horizon
   !-----------------------------------------------------------------------

   vrnew = 1.0d0

   call shearntwist (nn, mass)

end subroutine schwarzschild

subroutine dot_j_lambda (nn, dt, time, mass)
!-----------------------------------------------------------------------
! integrates \dot j_{,\lambda} = ... on the horizon
!-----------------------------------------------------------------------

   use horizon_eth
   implicit none

   integer,          intent (in) :: nn
   double precision, intent (in) :: dt, time, mass

   ! local arrays and variables

   double precision, dimension (nn,nn,2) :: rhoh, rholh, dotrhoh, kh, klh, kdoth
   double complex,   dimension (nn,nn,2) :: jh, jlh, jdoth, ckh, crhoh, &
                                            eth_rhoh, ethb_rhoh, eth_jh, &
                                            ethb_jh, eth_kh, ethb_kh, &
                                            eth_omegah, ethb_omegah, &
                                            omegah, delta1, delta2

   double precision :: th, r_m ! , r_ml  FIXME: r_ml not needed here?
   double precision :: sig     !  sig = sign(cfl) = sign(dt)

   !< executable statements >!

   sig = 1.0d0 ! sign(1.0d0, dt) ! forward or backward in time

   ! rk-2, first step:

   th = time - dt
   r_m = 2. * mass
   ! r_ml = - th / (4. * mass)
   rhoh = rho
   rholh = rhol
   dotrhoh = rhodot
   jh = j
   jlh = jl
   jdoth = (jnew - j) / dt
   kh = sqrt(1. + jh * conjg(jh))
   ckh = kh
   crhoh = rhoh
   klh = dble(jlh * conjg(jh)) / kh
   kdoth = dble(jdoth * conjg(jh)) / kh
   omegah = omega

   call eth1 (nn, eth_rhoh,    crhoh , 0,  1)
   call eth1 (nn, ethb_rhoh,   crhoh , 0, -1)
   call eth1 (nn, eth_jh,      jh    , 2,  1)
   call eth1 (nn, ethb_jh,     jh    , 2, -1)
   call eth1 (nn, eth_kh,      ckh   , 0,  1)
   call eth1 (nn, ethb_kh,     ckh   , 0, -1)
   call eth1 (nn, eth_omegah,  omegah, 1,  1)
   call eth1 (nn, ethb_omegah, omegah, 1, -1)

   delta1 = - dotrhoh / rhoh * jlh &
       - (rholh / rhoh - sig * th / (2. * r_m ** 2)) * jdoth &
       + 0.5 * jh * (conjg(jdoth) * jlh + jdoth * conjg(jlh) & 
                         - 2. * kdoth * klh) &
       + sig * ((1. + kh ** 2) &
       * (eth_omegah + omegah ** 2 - 2. * omegah * eth_rhoh / rhoh) &
       - omegah * (jh * ethb_kh - kh * ethb_jh) &
       - conjg(omegah) * (kh * eth_jh - jh * eth_kh) &
       + jh * (jh * conjg(eth_omegah) &
                  - kh * (conjg(ethb_omegah) + ethb_omegah) &
                  + jh * conjg(omegah) ** 2 &
                  - 2. * kh * omegah * conjg(omegah) &
                  + 2. * (omegah * ethb_rhoh + conjg(omegah) * eth_rhoh) &
                        / rhoh &
                  - 2. * jh * conjg(omegah) * ethb_rhoh / rhoh)) &
        / (2. * (r_m * rhoh) ** 2)

   jlnew = jl + 0.5 * dt * delta1

   ! rk-2, second step:

   th = time - 0.5 * dt
   r_m = 2. * mass
   ! r_ml = - th / (4. * mass)
   rhoh = 0.5 * (rhonew + rho)
   rholh = 0.5 * (rholnew + rhol)
   dotrhoh = 0.5 * (rhodotnew + rhodot)
   jh = 0.5 * (jnew + j)
   jlh = 0.5 * (jlnew + jl)
   jdoth = (jnew - j) / dt
   kh = sqrt(1. + jh * conjg(jh))
   ckh = kh
   crhoh = rhoh
   klh = dble(jlh * conjg(jh)) / kh
   kdoth = dble(jdoth * conjg(jh)) / kh
   omegah = 0.5 * (omeganew + omega)

   call eth1 (nn, eth_rhoh,    crhoh , 0,  1)
   call eth1 (nn, ethb_rhoh,   crhoh , 0, -1)
   call eth1 (nn, eth_jh,      jh    , 2,  1)
   call eth1 (nn, ethb_jh,     jh    , 2, -1)
   call eth1 (nn, eth_kh,      ckh   , 0,  1)
   call eth1 (nn, ethb_kh,     ckh   , 0, -1)
   call eth1 (nn, eth_omegah,  omegah, 1,  1)
   call eth1 (nn, ethb_omegah, omegah, 1, -1)

   delta2 = - dotrhoh / rhoh * jlh &
       - (rholh / rhoh - sig * th / (2. * r_m ** 2)) * jdoth &
       + 0.5 * jh * (conjg(jdoth) * jlh + jdoth * conjg(jlh) &
                         - 2. * kdoth * klh) &
       + sig * ((1. + kh ** 2) &
       * (eth_omegah + omegah ** 2 - 2. * omegah * eth_rhoh / rhoh) &
       - omegah * (jh * ethb_kh - kh * ethb_jh) &
       - conjg(omegah) * (kh * eth_jh - jh * eth_kh) &
       + jh * (jh * conjg(eth_omegah) &
                  - kh * (conjg(ethb_omegah) + ethb_omegah) &
                  + jh * conjg(omegah) ** 2 &
                  - 2. * kh * omegah * conjg(omegah) &
                  + 2. * (omegah * ethb_rhoh + conjg(omegah) * eth_rhoh) &
                        / rhoh &
                  - 2. * jh * conjg(omegah) * ethb_rhoh / rhoh)) &
        / (2. * (r_m * rhoh) ** 2)

   jlnew = jl + dt * delta2

end subroutine dot_j_lambda

subroutine shearntwist (nn, mass)

   use null_grid, only : qs, ps
   implicit none

   integer,          intent (in) :: nn
   double precision, intent (in) :: mass

   double precision, dimension (nn,nn,2) :: bigp, k, kl
   double complex,   dimension (nn,nn,2) :: mq, mp
   double precision :: r_m

   bigp(:,:,1) = 1. + qs * qs + ps * ps
   bigp(:,:,2) = 1. + qs * qs + ps * ps

   r_m = 2. * mass

   k  = sqrt(1. + j  * conjg(j))
   kl = dble(     jl * conjg(j)) / k

   ! spin-weight 2 field

   shears = 0.5 * (r_m * rho) ** 2 * ( (1. + k) * jl &
          - 2. * j * kl + j ** 2 * conjg(jl) / (1. + k) )

   Print *, "In shear : maxval(abs(shears))", maxval(abs(shears))

   ! the dyad we choose

   mq = bigp / (2. * sqrt(2. * (1. + k))) * (1. + k - j) 
   mp = (0.,1.) * bigp / (2. * sqrt(2. * (1. + k))) * (1. + k + j) 

   ! spin-weight 0 tensor

   shearqq = 2. * dble(conjg(shears) * mq * mq)
   shearqp = 2. * dble(conjg(shears) * mq * mp)
   shearpp = 2. * dble(conjg(shears) * mp * mp)

   ! the *short* way about it (valid only at q=0 or p=0)

   shearplus = 0.5 * (k * dble(jl) - dble(j) * kl) &
             / ((k - dble(j)) * (k + dble(j)))

   ! theta component of the twist

   twist_theta = dble(omega) * bigp**2 / (4. * r_m * rho* (k+ dble(j)))

!write (*,*) 'In shearntwist:'
!write (*,*) 'max of bigp ', maxval(bigp)
!write (*,*) 'max of dble(omega) ', maxval(dble(omega))
!write (*,*) 'max of dble(j) ', maxval(dble(j))
!write (*,*) 'max of k ', maxval(k)
!write (*,*) 'max of rho ', maxval(rho)

end subroutine shearntwist

end module horizon
