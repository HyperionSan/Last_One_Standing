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
      omeganew, omega

   double complex, dimension (:,:,:), allocatable, save, public :: &
      jrnew, qnew, unew, urnew

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
              omeganew  (nn,nn,2), omega (nn,nn,2), &
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
              beta_red  (nn,nn,2)   )

end subroutine horizon_setup

subroutine shifter 
!-----------------------------------------------------------------------
! shift down time levels as needed
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
! integrate $\dot \rho = -\rho/4 (\dot j \bar{\dot j} - (\dot k)^2)$
!-----------------------------------------------------------------------

   implicit none

   integer,          intent (in) :: nn
   double precision, intent (in) :: dt

   ! local arrays

   double precision, dimension (nn,nn,2) :: k, knew, kdot, source
   double complex,   dimension (nn,nn,2) :: jdot

   k = sqrt(1. + j * conjg(j))
   knew = sqrt(1. + jnew * conjg(jnew))
   jdot = (jnew - j) / dt
   kdot = (knew - k) / dt

   source = (jdot * conjg(jdot) - kdot ** 2) / 16.0d0

   rhonew = (rho * (1. - source * dt ** 2) + rhodot * dt ) &
        / (1. + source * dt ** 2)

   rhodotnew = (rhodot * (1. - source * dt ** 2) + 4. * source * rho * dt ) &
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

!if (modulo(nn,2) == 1) then
! call axi_outr(nn, 0, 'omega_new_abs', abs(omeganew((nn+1)/2,:,1)))
! call axi_outr(nn, 0, 'omega_rhs_abs', abs(rhs((nn+1)/2,:,1)))
! call axi_outr(nn, 0, 'omega_rhs_im', aimag(rhs((nn+1)/2,:,1)))
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

   double precision, dimension (nn,nn,2) :: rhoh, rhodoth, kh, ricci, rhs

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
   implicit none

   integer,          intent (in) :: nn
   double precision, intent (in) :: dt, time, mass

   ! local arrays

   double precision, dimension (nn,nn,2) :: rhoh, rhodoth, kh, ricci, rhs, &
      rdot, rhol_dot_rhs ,rhol_mid, k1, k2
   double complex,   dimension (nn,nn,2) :: jh, ckh, crhoh, omegah,        &
      eth_j, ethb_j, ethb2_j, ethb_k, eth_ethb_k, eth_rho, eth2_rho,       & 
      ethb_eth_rho, eth_omega, ethb_omega

   ! local scalar variables

   double precision :: th, r_m

   !< executable statements >!

   ! compute k1
   th = time - dt
   r_m = 2. * mass

   call spheroid_rdot(nn, rdot)

   rhoh = rho                      ! the suffix h was for "halfstep"
   rhodoth = rdot                  ! before, i kept it for minimal code
   crhoh = rhoh                    ! change
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
   + conjg(jh) * (  eth_omega - omegah ** 2 &
                     + eth2_rho / rhoh - (eth_rho / rhoh) ** 2 ) &
   + kh * omegah * conjg(omegah)) - 0.5 * ricci + rhoh ** 2

   rhol_dot_rhs =   rhs / (rhoh * 2.*r_m ** 2)                           &
                 + rhodoth * (th / r_m ** 2 - rhol/rhoh)
   k1 = dt * rhol_dot_rhs

   ! compute k2
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

   rhol_mid = rhol + 0.5 * k1 
   rhol_dot_rhs =   rhs / (rhoh * 2.*r_m ** 2)                           &
                 + rhodoth * (th / r_m ** 2 - rhol_mid/rhoh)
   k2 = dt * rhol_dot_rhs

   rholnew =  rhol + k2 ! + O(dt^3)


end subroutine dot_rho_lambda_RK2


subroutine htor (nn, time, dt, mass)
!-----------------------------------------------------------------------
!  convert from horizon to code variables
!-----------------------------------------------------------------------

   use horizon_eth
   implicit none

   integer,          intent (in) :: nn
   double precision, intent (in) :: time, dt, mass

   !-----------------------------------------------------------------------
   ! local variables
   !-----------------------------------------------------------------------

   double precision :: r_m, r_ml, r_mu, r_mlu, min_rl

   !-----------------------------------------------------------------------
   ! local arrays
   !-----------------------------------------------------------------------

   double precision, dimension (nn,nn,2) :: rholu, knew, krnew, rlu,      &
                                          & ricci, rhs, r_l

   double complex, dimension (nn,nn,2) :: ckh, crhoh, crhol,              &
      eth_j, ethb_j, ethb2_j, ethb_k, eth_ethb_k, eth_rho, eth2_rho,      &
      ethb_eth_rho, ethb_rho, eth_rhol, ethb_rhol, eth_omega, ethb_omega, &
      eth_beta

   !-----------------------------------------------------------------------
   ! factor asymptotic dependence of rho, store j, k, and needed eth opers.
   !-----------------------------------------------------------------------

   r_m = 2. * mass
   r_ml = - time / (4. * mass)
   r_mu = 0.
   r_mlu = - 1. / (4. * mass)

   knew = sqrt(1. + jnew * conjg(jnew))

   crhoh = rhonew
   crhol = rholnew

   call eth1 (nn, eth_rho ,  crhoh, 0,  1)
   call eth1 (nn, eth_rhol,  crhol, 0,  1)
   
   ethb_rho  = conjg(eth_rho)
   ethb_rhol = conjg(eth_rhol)

   !-----------------------------------------------------------------------
   ! Compute the betas (full and reduced) on the horizon
   !-----------------------------------------------------------------------

   r_l = r_M * rholnew + r_Ml * rhonew    ! R_{,\lambda}

   e2beta   = 1. / r_l                    ! \__  the full beta
   betanew  = 0.5 * log(e2beta)           ! / 

   beta_red = betanew + 0.5*log(r_Ml)     ! this is beta_R = reduced beta 

   min_rl = min( minval(r_l(:,:,1)), minval(r_l(:,:,2)) )
   if (min_rl <= 0.0) then
      write(10,*) 'BH ?! min(R,lambda) =', min_rl, '@ time', time
      call flush(10)
   endif

   !-----------------------------------------------------------------------
   ! Compute the outward wavefront expansion on the horizon
   !----------------------------------------------------------------------- 

   exp_out = r_l/(r_M * rhonew)

   !-----------------------------------------------------------------------
   ! compute j_{,r} and k_{,r} on the horizon
   !-----------------------------------------------------------------------

   jrnew = e2beta * jlnew
   krnew = (jrnew * conjg(jnew) + jnew * conjg(jrnew)) / (2. * knew)

   !-----------------------------------------------------------------------
   ! compute beta_{,r} on the horizon
   !-----------------------------------------------------------------------

   betarnew = (r_m * rhonew) / 8.0d0 * (jrnew * conjg(jrnew) - krnew ** 2)

   !-----------------------------------------------------------------------
   ! compute u on the horizon
   !-----------------------------------------------------------------------

   unew = e2beta * (- jnew * ethb_rho + knew * eth_rho) / (r_m * rhonew ** 2)

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

   rholu = rhs / (8. * mass ** 2) &
         - (rhodotnew - time * rhonew / (4. * mass ** 2)) * rholnew

   rlu = r_mlu * rhonew + r_ml * rhodotnew + r_mu * rholnew + r_m * rholu

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

end subroutine htor


subroutine htor_dummy (nn, time, dt, mass)
!-----------------------------------------------------------------------
!  fill the evolution code variables with dummy values
!-----------------------------------------------------------------------

   implicit none

   integer,          intent (in) :: nn
   double precision, intent (in) :: time, dt, mass

   !-----------------------------------------------------------------------
   ! local variables
   !-----------------------------------------------------------------------

   double precision :: r_m, r_ml, r_mu, r_mlu, min_rl

   !-----------------------------------------------------------------------
   ! local arrays
   !-----------------------------------------------------------------------

   double precision, dimension (nn,nn,2) :: rholu, knew, krnew, rlu, r_l

   !-----------------------------------------------------------------------
   ! factor asymptotic dependence of rho, store j, k, and needed eth opers.
   !-----------------------------------------------------------------------

   r_m = 2. * mass
   r_ml = - time / (4. * mass)
   r_mu = 0.
   r_mlu = - 1. / (4. * mass)

   knew = 1.
  
   !-----------------------------------------------------------------------
   ! Compute the betas (full and reduced) on the horizon
   !-----------------------------------------------------------------------

   rhonew  = 1.0
   rholnew = 0.0
   r_l = 1.0   ! R_{,\lambda}

   e2beta   = 1. / r_l                    ! \__  the full beta
   betanew  = 0.5 * log(e2beta)           ! / 

   beta_red = betanew + 0.5*log(abs(r_ml))     ! this is beta_R = reduced beta 

   !-----------------------------------------------------------------------
   ! Compute the outward wavefront expansion on the horizon
   !----------------------------------------------------------------------- 

   exp_out = r_l/r_M

   !-----------------------------------------------------------------------
   ! compute beta_{,r} on the horizon
   !-----------------------------------------------------------------------

   betarnew = r_m / 8.0d0 

   !-----------------------------------------------------------------------
   ! compute u on the horizon
   !-----------------------------------------------------------------------

   unew = 0.0

   !-----------------------------------------------------------------------
   ! compute q on the horizon
   !-----------------------------------------------------------------------

   qnew = 0.0

   !-----------------------------------------------------------------------
   ! compute u_{,r} on the horizon
   !-----------------------------------------------------------------------

   urnew = 0.0

   !-----------------------------------------------------------------------
   ! compute v on the horizon
   !-----------------------------------------------------------------------

   vnew = 0.0

   !-----------------------------------------------------------------------
   ! compute \rho_{,\lambda u} on the horizon
   !-----------------------------------------------------------------------

 
   rlu = 0.0

   !-----------------------------------------------------------------------
   ! compute v,r on the horizon
   !-----------------------------------------------------------------------

   vrnew = 0.0

end subroutine htor_dummy


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

   double precision :: th, r_m, r_ml


   ! rk-2, first step:

   th = time - dt
   r_m = 2. * mass
   r_ml = - th / (4. * mass)
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
       - (rholh / rhoh - th / (2. * r_m ** 2)) * jdoth &
       + 0.5 * jh * (conjg(jdoth) * jlh + jdoth * conjg(jlh) &
                         - 2. * kdoth * klh) / rhoh ** 2 &
       + ((1. + kh ** 2) &
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
   r_ml = - th / (4. * mass)
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
       - (rholh / rhoh - th / (2. * r_m ** 2)) * jdoth &
       + 0.5 * jh * (conjg(jdoth) * jlh + jdoth * conjg(jlh) &
                         - 2. * kdoth * klh) / rhoh ** 2 &
       + ((1. + kh ** 2) &
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

end module horizon
