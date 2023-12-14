module news2

   use ascii_io
   use null_params, only : iot2, it_start
   implicit none

   double precision, dimension (:,:,:), allocatable, save, public :: uBondi
   double complex,   dimension (:,:,:), allocatable, save, public :: NewsB
   double precision, dimension (:,:,:), allocatable, save, public :: dMdOmega

   double complex,   dimension (:,:,:), allocatable, save, private :: &
      zeta, zetabar, Jo, Jn, Jo_l, Jn_l, cBo, cBn, Uo, Un, comegao, comegan, &
      News, zBondio, zBondin, zBondih, Uyo, Uyn
 
   double precision, dimension (:,:,:), allocatable, save, private :: &
      P, betao, betan, omegao, omegan, redshiftB, &
      uBondio, uBondin, deltao, deltan

   integer, dimension (:,:,:), allocatable, save, private :: &
       nsmask
   integer, private :: ret
   double precision, dimension(:,:,:), allocatable :: circle
   double precision :: time_of_news
   logical, save, private :: initial = .true.
contains

subroutine news_allocate

   use null_grid, only : nn, dd, z, zb, pp
   use null_params, only: time
   implicit none

   integer ip

   allocate ( zeta(nn,nn,2), zetabar(nn,nn,2), &
      Jo(nn,nn,2), Jn(nn,nn,2), Jo_l(nn,nn,2), Jn_l(nn,nn,2), &
      cBo(nn,nn,2), cBn(nn,nn,2), &
      Uo(nn,nn,2), Un(nn,nn,2), &
      comegao(nn,nn,2), comegan(nn,nn,2), &
      deltao(nn,nn,2), deltan(nn,nn,2), &
      News(nn,nn,2), &
      P(nn,nn,2), &
      betao(nn,nn,2), betan(nn,nn,2), &
      omegao(nn,nn,2), omegan(nn,nn,2), &
      redshiftB(nn,nn,2), &
      uBondi(nn,nn,2) , uBondio(nn,nn,2), uBondin(nn,nn,2), &
      zBondio(nn,nn,2), zBondin(nn,nn,2), zBondih(nn,nn,2), &
      Uyo(nn,nn,2), Uyn(nn,nn,2), &
      NewsB(nn,nn,2), &
      dMdOmega(nn,nn,2), nsmask(nn,nn,2) )

   do ip = 1, 2
      zeta(:,:,ip) = z
      zetabar(:,:,ip) = zb
   end do
   P = 1.0d0 + zeta * zetabar

   omegao = (1.0d0, 0.0d0)
   comegao = (0.0d0, 0.0d0)
   deltao = 0.0d0      ! deltao and mask will be changes 
   uBondio = time      ! prior to the first run of newszbondi
   nsmask = 0          ! see newszbondi

   do ip = 1, 2
      zBondio(:,:,ip) = z
   end do


   allocate (circle(nn,nn,2))
   do ip = 1, 2
!     where (pp - 1.0d0 < (1.0d0 - 3.0d0*dd)**2)
      where (pp < 2.0d0)
         circle(:,:,ip) = 1.0d0
      elsewhere
         circle(:,:,ip) = 0.0d0
      end where
   end do

! ...later...
!
!  open (unit = 12, file = 'news.in', status = 'unknown')
!  write (unit = 12, nml = news_input)
!  close (unit = 12)

end subroutine news_allocate

subroutine news_omega (Uo, Un, omegao, omegan)

   use null_grid, only : nn, it, dt
   use null_params, only: time

   implicit none

   double complex,   dimension (nn,nn,2), intent (in)  :: Uo, Un
   double precision, dimension (nn,nn,2), intent (in)  :: omegao
   double precision, dimension (nn,nn,2), intent (out) :: omegan

   double complex,   dimension (nn,nn,2) :: U 
   double precision, dimension (nn,nn,2) :: omega, omega_u

   ! Runge-Kutta for \omega_{,u}

   call news_dot_omega (Uo, omegao, omega_u)
   omegan = omegao + dt * omega_u

   U = 0.5d0 * (Uo + Un)
   omega = 0.5d0 * (omegao + omegan)
   call news_dot_omega (U, omega, omega_u)
   omegan = omegao + dt * omega_u

   !if (mod(it,iot2) .eq. 0 .and. it .gt. 1) then
   if (mod(it,iot2) .eq. 0 .and. it .gt. it_start) then
     ret = gft_write ('omega', time, omegan(:,:,1))
   end if

end subroutine news_omega 

subroutine news_dot_omega (U, omega, omega_u)

   use null_grid, only : nn
   use horizon_eth
   implicit none

   double complex,   dimension (nn,nn,2), intent (in)  :: U
   double precision, dimension (nn,nn,2), intent (in)  :: omega
   double precision, dimension (nn,nn,2), intent (out) :: omega_u

   double complex,   dimension (nn,nn,2) :: eth_omega, eth_Ub

   call eth1 (nn, eth_omega, dcmplx(omega), 0, 1)
   call eth1 (nn, eth_Ub, conjg(U), -1, 1)

   omega_u = - dble(eth_omega * conjg(U) + 0.5d0 * omega * eth_Ub)

end subroutine news_dot_omega

subroutine news_comega (omegao, omegan, comegao, comegan)

   use null_grid, only : nn, it, dt
   use null_params, only: time
   use horizon_eth
   implicit none

   double precision, dimension (nn,nn,2), intent (in)  :: omegao, omegan
   double complex,   dimension (nn,nn,2), intent (in)  :: comegao
   double complex,   dimension (nn,nn,2), intent (out) :: comegan

   double complex, dimension (nn,nn,2) :: eth_omega_u

   ! Mid-point rule for \eth\omega_{,u}

   call eth1 (nn, eth_omega_u, dcmplx((omegan - omegao) / dt), 0, 1)
   comegan = comegao + dt * eth_omega_u

   if (mod(it - (it_start+1), iot2) == 0 .and. it > it_start) then
!     ret = gft_write ('rcomega', time,  dble(comegan(:,:,1)))
!     ret = gft_write ('icomega', time, dimag(comegan(:,:,1)))
   end if

end subroutine news_comega


subroutine news_uframe (J, J_u, J_l, J_l_u, beta, cB, U, omega, &
                        comega, News)

   use null_grid, only : nn, it, dt
   use null_params, only: time
   use horizon_eth
   use null_interp
   implicit none

   double precision, dimension (nn,nn,2) :: beta, omega
   double complex,   dimension (nn,nn,2) :: J, J_u, J_l, J_l_u, &
                                            cB, U, comega, News

   ! temporary arrays

   double precision, dimension (nn,nn,2) :: a, K, K_l, K_u, K_l_u
   double complex,   dimension (nn,nn,2) :: Jb, Ub, s1, s2, s3, s4, s5

   double complex,   dimension (nn,nn,2) :: eth_J, ethb_J, eth_J_l, ethb_J_l, &
      eth_K, eth_K_l, eth_beta, eth2_beta, eth_ethb_beta, eth_U, ethb_U, &
      eth_omega, eth2_omega, eth_ethb_omega, &
      eth_a, eth2_a, eth_ethb_a

   Jb = conjg(J)
   Ub = conjg(U)
   K = sqrt(1.0d0 + J * Jb)

   K_u = dble( J_u * Jb ) / K
   K_l = dble( J_l * Jb ) / K
   K_l_u = dble( J_u * conjg(J_l) + J_l_u * Jb )/ K - K_l * K_u / K 

   call eth1 (nn, eth_K,   dcmplx(K),   0, 1)
   call eth1 (nn, eth_U,  U, 1,  1)
   call eth1 (nn, ethb_U, U, 1, -1)
   call eth1 (nn, eth_J,  J, 2,  1)
   call eth1 (nn, ethb_J, J, 2, -1)
   call eth1 (nn, eth_J_l, J_l, 2,  1)
   call eth1 (nn, ethb_J_l, J_l,  2, -1)
   call eth1 (nn, eth_K_l, dcmplx(K_l), 0, 1)
   eth_omega = comega
   call eth1 (nn, eth2_omega,     comega, 1,  1)
   call eth1 (nn, eth_ethb_omega, comega, 1, -1)
   eth_beta = cB
   call eth1 (nn, eth2_beta,     cB, 1,  1)
   call eth1 (nn, eth_ethb_beta, cB, 1, -1)

   a = omega * exp(2.0d0 * beta)

! debug
!  eth2_omega = (0.0d0, 0.0d0)
!  eth_ethb_omega = (0.0d0, 0.0d0)
! debug

! debug
!  eth2_beta = (0.0d0, 0.0d0)
!  eth_ethb_beta = (0.0d0, 0.0d0)
! debug

   eth_a = exp(2.0d0 * beta) * ( eth_omega + 2.0d0 * omega * eth_beta )
	
   eth2_a = exp(2.0d0 * beta) * ( 4.0d0 * eth_beta * eth_omega &
                             + 4.0d0 * omega * eth_beta**2 &
                             + eth2_omega + 2.0d0 * omega * eth2_beta )

   eth_ethb_a = exp(2.0d0 * beta) * ( 4.0d0 * dble(eth_beta * conjg(eth_omega)) &
                                 + 4.0d0 * omega * eth_beta * conjg(eth_beta) &
                                 + eth_ethb_omega + 2.0d0 * omega * eth_ethb_beta )

! debug
!  eth2_a = (0.0d0, 0.0d0)
!  eth_ethb_a = (0.0d0, 0.0d0)
! debug

   s1 = ( -2.0d0 * K_l_u * J * (K + 1.0d0) + J_l_u * (K + 1.0d0)**2 &
          + conjg(J_l_u) * J**2 ) / (K + 1.0d0)

   s2 = 0.5d0 / ( K + 1.0d0) * ( &
      	     (K + 1.0d0)* (eth_J_l *Ub * (K+1.0d0) - 2.0d0* eth_K_l * J *Ub ) &
             + eth_U * (K+1.0d0)* ( -2.0d0 * J * conjg(J_l) + K_l * 2.0d0 * (K+1.0d0) ) &
             + conjg(ethb_U) * (K+1.0d0) * ( -2.0d0* J * K_l + J_l * 2.0d0 * (K+1.0d0) ) &
             + ethb_J_l * U * (K+1.0d0)**2 - conjg(eth_K_l) * 2.0d0 * U * J * (K+1.0d0) &
             + ethb_U * 2.0d0 * J * ( J * conjg(J_l) - (K+1.0d0) * K_l) &
             + J**2 * ( U * conjg(eth_J_l) + conjg(ethb_J_l * U) ) &
             + J * 2.0d0 * conjg(eth_U) * ( J * K_l - J_l * (K+1.0d0) ) )

   s3 = ( J_l * (K + 1.0d0)**2 -2.0d0 * K_l * J * (K + 1.0d0) &
          + conjg(J_l) * J**2) / (K + 1.0d0)

   s4 = 0.5d0 / ( K + 1.0d0) * ( eth_a * eth_omega * (K + 1.0d0)**2 &
      		- (K+1.0d0) * J * 2.0d0* dble( eth_a * conjg(eth_omega) ) &
                + J**2 * conjg(eth_a * eth_omega) )

   s5 = 0.25d0 / ( K + 1.0d0) * ( 2.0d0 * eth2_a * (K+1.0d0)**2 &
             + 2.0d0 * J**2 * conjg(eth2_a) &
             - 4.0d0 * eth_ethb_a * J * (K+1.0d0) &
             + Jb * eth_a * eth_J* (K+1.0d0)**2 &
             + J * eth_a * conjg(ethb_J) * (K+1.0d0)**2 &
             - eth_a * eth_K * 2.0d0 * (K+1.0d0) * ( J*Jb + (K+1.0d0) ) & 
             + eth_a * ethb_J * (K+1.0d0) * ( -J*Jb + (K+1.0d0) ) &
             - J**2 * eth_a * conjg(eth_J) * K &
             +  J**2 * Jb * 2.0d0* eth_a * conjg(eth_K) &
             - conjg(eth_a) * eth_J * (K+1.0d0) * ( J*Jb + K+1.0d0 ) &
             - conjg(ethb_J) * conjg(eth_a) * J**2 * ( K + 2.0d0) &
             + J * 2.0d0 * (K+1.0d0)**2 * eth_K * conjg(eth_a)  &
             + J**2 * Jb * ethb_J * conjg(eth_a) &
             + J**3 * conjg(eth_a * eth_J) &
             - 2.0d0* J**2 *K*conjg(eth_K * eth_a) )

!   News = 0.25d0 * ( s1 + s2 + 0.5d0 * dble(ethb_U) * s3 &
!   - 4.0d0 * s4 / omega**2 + 2.0d0 * s5 / omega ) / ( omega**2 * exp(2.0d0 * beta) )

! change sign of s3 to compensate for a bug in Eqs. 30, 37, and 38 of
! HPN
   News = 0.25d0 * ( s1 + s2 - 0.5d0 * dble(ethb_U) * s3 &
   - 4.0d0 * s4 / omega**2 + 2.0d0 * s5 / omega ) / ( omega**2 * exp(2.0d0 * beta) )

!debug:
!  s2 = (0.0d0, 0.0d0)
!  s4 = (0.0d0, 0.0d0)
!  s5 = (0.0d0, 0.0d0)
!  News = 0.25d0 * Delta * ( s1 + s2 + 0.5d0 * dble(ethb_U) * s3 &
!  - 4.0d0 * s4 / omega**2 + 2.0d0 * s5 / omega ) / ( omega**2 * exp(2.0d0 * beta) )
!debug:

!   call null_cnsint2 (news(:,:,1), news(:,:,2), 2)
   ! your mileage may vary...
   call null_cnsint (news(:,:,1), news(:,:,2), 2)
   call null_cnsint (news(:,:,2), news(:,:,1), 2)
!  call null_ccircint (news(:,:,1), news(:,:,2), 2)

   !if (mod(it,iot2) .eq. 0 .and. it .gt. 1) then
   if (mod(it,iot2) .eq. 0 .and. it .gt. it_start) then
!     ret = gft_write ('rN', time,  dble(News(:,:,1)))
!     ret = gft_write ('iN', time, dimag(News(:,:,1)))
   end if

end subroutine news_uframe


subroutine cscrival (jnn, jns, J)

   use null_grid, only : nn, nx
   implicit none

   double complex, dimension (nn,nn,nx) :: jnn, jns
   double complex, dimension (nn,nn,2) :: J

   J(:,:,1) = jnn(:,:,nx)
   J(:,:,2) = jns(:,:,nx)
   
end subroutine cscrival

subroutine rscrival (bnn, bns, beta)

   use null_grid, only : nn, nx
   implicit none

   double precision, dimension (nn,nn,nx) :: bnn, bns
   double precision, dimension (nn,nn,2) :: beta

   beta(:,:,1) = bnn(:,:,nx)
   beta(:,:,2) = bns(:,:,nx)
   
end subroutine rscrival

subroutine cscrivalh (unn, uns, U)

   use null_grid, only : nn, nx
   implicit none

   double complex, dimension (nn,nn,nx) :: unn, uns
   double complex, dimension (nn,nn,2) :: U

   U(:,:,1) = 0.5d0 * (unn(:,:,nx) + unn(:,:,nx-1))
   U(:,:,2) = 0.5d0 * (uns(:,:,nx) + uns(:,:,nx-1))
   
end subroutine cscrivalh

subroutine cscridbydl (jnn, jns, J_l)

   use null_grid, only : nn, nx, rwt, dx
   implicit none

   double complex, dimension (nn,nn,nx) :: jnn, jns
   double complex, dimension (nn,nn,2) :: J_l

   J_l(:,:,1) = - 0.5d0 * ( 3.0d0 * jnn(:,:,nx) - 4.0d0 * jnn(:,:,nx-1) &
                             + jnn(:,:,nx-2) ) / (dx ) * rwt
   
   J_l(:,:,2) = - 0.5d0 * ( 3.0d0 * jns(:,:,nx) - 4.0d0 * jns(:,:,nx-1) &
                             + jns(:,:,nx-2) ) / (dx ) * rwt
   
end subroutine cscridbydl

subroutine scrivals (Jo, Jn, Jo_l, Jn_l, betao, betan, cBo, cBn, Uo, Un)

   use null_grid, only : nn, it, dt
   use null_params, only: time
   use null_vars
   implicit none
   
   double complex,   dimension (nn,nn,2), intent (out) :: Jo, Jn, Jo_l, Jn_l, &
                                                          cBo, cBn, Uo, Un
   double precision, dimension (nn,nn,2), intent (out) :: betao, betan

   ! get scri quantities from null code variables
   ! AXISYMETRY REQUIRED 

   call cscrival (jon, jon, Jo)
   call cscrival (jnn, jnn, Jn)

   call cscridbydl (jon, jon, Jo_l)
   call cscridbydl (jnn, jnn, Jn_l)

   call rscrival (bon, bon, betao)
   call rscrival (bnn, bnn, betan)

   call cscrival (cbon, cbon, cBo)
   call cscrival (cbnn, cbnn, cBn)

   call cscrivalh (uon, uon, Uo)
   call cscrivalh (unn, unn, Un)


   !if (mod(it,iot2) .eq. 0 .and. it .gt. 1) then
   if (mod(it,iot2) .eq. 0 .and. it .gt. it_start) then
!     ret = gft_write ('rJ',   time,    dble(Jn(:,:,1)))
!     ret = gft_write ('iJ',   time,   dimag(Jn(:,:,1)))
!     ret = gft_write ('rJl',  time,  dble(Jn_l(:,:,1)))
!     ret = gft_write ('iJl',  time, dimag(Jn_l(:,:,1)))
!     ret = gft_write ('beta', time,      betan(:,:,1))
!     ret = gft_write ('rU',   time,    dble(Un(:,:,1)))
!     ret = gft_write ('iU',   time,   dimag(Un(:,:,1)))
   end if

end subroutine scrivals

subroutine news_zBondi (Uo, Un, Uyo, Uyn, zBondio, zBondin, zBondih, &
                     News, NewsB, uBondio,uBondi, uBondin, redshiftB, &
                    dMdOmega, omega, Jh, eth_J, beth_J, eth_U, beth_U, du_J, beta, nsmask)
   use null_grid, only : nn, it, dt, ii, qs, ps, qsize
   use null_params, only : time, ZERO, EDGE
   use null_interp2
   use null_interp
   implicit none

   double complex,   dimension (nn,nn,2), intent (in)    :: Uo, Un, Jh, eth_J, beth_J, eth_U, beth_u, du_J
   double complex,   dimension (nn,nn,2), intent (inout) :: Uyo, Uyn
   double complex,   dimension (nn,nn,2), intent (inout) :: zBondio, zBondin, &
                                                            zBondih
   double complex,   dimension (nn,nn,2), intent (in)    :: News
   double complex,   dimension (nn,nn,2), intent (out)   :: NewsB
   double precision, dimension (nn,nn,2), intent (in)    :: omega, beta
   double precision, dimension (nn,nn,2), intent (inout) :: uBondi, uBondin, uBondio
   double precision, dimension (nn,nn,2), intent (inout) :: redshiftB
   double precision, dimension (nn,nn,2), intent (out)   :: dMdOmega
   integer,          dimension (nn,nn,2), intent (inout) :: nsmask

   double complex,   dimension (nn,nn,2) :: Uh, Uyh, Po, Pm
   double complex,   dimension (nn,nn,2) :: delta_RHS, delta_RHS_B
   double precision,   dimension (nn,nn,2) :: delta, K, betaB, omegaB
   double precision, save :: pii = 3.1415926535897932385d0
   integer ::  i, j, ip, ipo, ipm
   double precision qh, ph, qho, pho
   double complex zho, zbho
 
   logical, save :: ThisIsTheFirstTime = .true.
  
   if (initial) then
       Uyo = Uo
      initial = .false.
      where (dble(zbondio*conjg(zbondio))>1.0d0)  ! correctly mask all point initialy outside equator
        nsmask = 1
        deltao = pii -2.0d0 *atan2(dimag(zbondio), dble(zbondio))
        zbondio = 1.0d0 / zbondio
        Uyo = Uyo * (-zbondio /conjg(zbondio))
      end where
   end if

   ! "old" level of U(y^A) is known, step forward half a dt


   Po = 1.0d0+ zbondio*conjg(zbondio)
   zBondin = zBondio + 0.25d0 * Po * Uyo * dt
   Uh = 0.5d0 * (Un + Uo)
   Pm = 1.0d0+ zbondin*conjg(zbondin)
   ! interpolation of Uh (at mid level) to get Uyh. 

   do ip = 1, 2
      do j = 1, nn
         do i = 1, nn
            ipm = mod( ip + nsmask(i,j,ip) - 1 ,2 ) +1
            qh =  dble(zBondin(i,j,ip))
            ph = dimag(zBondin(i,j,ip))
            if ((abs(qh) .le. qsize) .and. (abs(ph) .le. qsize)) then 
               call cinterp(Uh(:,:,ipm), Uyh(i,j,ip), qh, ph)
            else
               if (ipm .eq. 1) then
                  ipo = 2
               else
                  ipo = 1
               end if
               qho =  qh / (qh * qh + ph * ph)
               pho = -ph / (qh * qh + ph * ph)
               zho  = qho + ii * pho
               zbho = qho - ii * pho
               call cinterp(Uh(:,:,ipo), Uyh(i,j,ip), qho, pho)
               Uyh(i,j,ip) = Uyh(i,j,ip) * (-zbho / zho)
            end if
         end do
      end do
   end do

   ! now step forward a full dt, with the values of Uyh (at mid level)
   ! to get a second-order accurate (RK2) value for z(Bondi) at the new level

   zBondin = zBondio + .5*Pm * Uyh * dt
   zBondih = 0.5d0 * (zBondin + zBondio)

   ! interpolation of Un (at the new level) to get Uyn. these values 
   ! will be needed the next time around.
   ! shamelessly reuse the scalar variables [qh, ph, qho, pho, ipo]
   ! to mean the coordinates and index at the new level in the loops below.

   ! the uBn, uBo live on the computational grid at scri, and *not* in the
   ! inertial observers grid, so we interpolate to get the change in Bondi
   ! time for those observers.

   do ip = 1, 2
      do j = 1, nn
         do i = 1, nn
            ipm = mod( ip + nsmask(i,j,ip) - 1 ,2 ) +1
            qh =  dble(zBondin(i,j,ip))
            ph = dimag(zBondin(i,j,ip))
            if ((abs(qh) .le. qsize) .and. (abs(ph) .le. qsize)) then 
               call cinterp(Un(:,:,ipm), Uyn(i,j,ip), qh, ph)
            else
               if (ipm .eq. 1) then
                  ipo = 2
               else
                  ipo = 1
               end if
               qho =  qh / (qh * qh + ph * ph)
               pho = -ph / (qh * qh + ph * ph)
               zho  = qho + ii * pho
               zbho = qho - ii * pho
               call cinterp(Un(:,:,ipo), Uyn(i,j,ip), qho, pho)
               Uyn(i,j,ip) = Uyn(i,j,ip) * (-zbho / zho)
            end if
         end do
      end do
   end do

   K = sqrt(1.0d0 + dble(Jh*conjg(Jh)))
   delta_RHS =  conjg(du_J)*Jh / (K + 1.0d0) + &
                 .5d0 /(K + 1.0d0) * (Jh*(Uh*conjg(eth_J)+conjg(Uh)*conjg(beth_J))) +&
                  Jh*conjg(eth_U) + K*beth_U


   ! interpolation of the News, uBondi and redshift (at mid level) 
   ! to the points of the grid of inertial observers at scri.

   do ip = 1, 2
      do j = 1, nn
         do i = 1, nn
            ipm = mod( ip + nsmask(i,j,ip) - 1 ,2 ) +1
            qh =  dble(zBondih(i,j,ip))
            ph = dimag(zBondih(i,j,ip))
            if ((abs(qh) .le. qsize) .and. (abs(ph) .le. qsize)) then 
               call cinterp(Uh(:,:,ipm), Uyh(i,j,ip), qh, ph)
               call cinterp(delta_RHS(:,:,ipm), delta_RHS_B(i,j,ip), qh, ph)
               call rinterp(beta(:,:,ipm), betaB(i,j,ip), qh, ph)
               call cinterp(News(:,:,ipm), NewsB(i,j,ip), qh, ph)
               call rinterp(omega(:,:,ipm), omegaB(i,j,ip), qh, ph)
            else
               if (ipm .eq. 1) then
                  ipo = 2
               else
                  ipo = 1
               end if
               qho =  qh / (qh * qh + ph * ph)
               pho = -ph / (qh * qh + ph * ph)
               zho  = qho + ii * pho
               zbho = qho - ii * pho
               call cinterp(Uh(:,:,ipo), Uyh(i,j,ip), qho, pho)
               Uyh(i,j,ip) = Uyh(i,j,ip) * (-zbho / zho)
               call cinterp(News(:,:,ipo), NewsB(i,j,ip), qho, pho)
               NewsB(i,j,ip) = NewsB(i,j,ip) * (-zbho / zho)**2
               call cinterp(delta_RHS(:,:,ipo), delta_RHS_B(i,j,ip), qho, pho)
               call rinterp(beta(:,:,ipo), betaB(i,j,ip), qho, pho)
               call rinterp(omega(:,:,ipo), omegaB(i,j,ip), qho, pho)
            end if
         end do
      end do
   end do
   delta_RHS_B = .5d0 * delta_RHS_B  +  uyh*conjg(zbondih) 
   deltan = deltao + dt * dimag(delta_RHS_B)

   redshiftB = omegaB * exp(2*betaB) 
   uBondin = uBondio + dt * omegaB * exp(2*betaB)
   uBondi = .5d0 * (uBondin + uBondio)

   delta = .5d0 * (deltan + deltao)
   NewsB = NewsB * ( cos(2.0d0*delta) - ii*sin(2.0d0*delta) )
   call null_cnsint (NewsB(:,:,2), NewsB(:,:,1), 2)
   call null_cnsint (NewsB(:,:,1), NewsB(:,:,2), 2)

   ! La definicion es z=\tan(\theta/2) e^{i\phi}, asi que aca le 
   ! sacamos el z/zb que le pusimos en Spheroid_RJ.

   ! Translation: the definition is z=\tan(\theta/2) e^{i\phi}, thus here
   ! we take out the factor of zb/z that we stuck J with in Spheroid_RJ

   where (zeta .ne. (0.0, 0.0))
      NewsB = NewsB * zetabar / zeta
   end where

   dMdOmega = NewsB * conjg(NewsB) * redshiftB * circle

   ! update the time levels of Uy (the "shifted" array of U)

   zBondio = zBondin
   Uyo = Uyn
   deltao = deltan
   uBondio = uBondin

   where ( dble(zbondio * conjg(zbondio)) > 1.0d0 )
     nsmask = mod(nsmask +1, 2)  ! flip patch (or back, as the case may be)
     deltao = deltao + pii -2.0d0 *atan2(dimag(zbondio), dble(zbondio))  ! e^{i delta} does
                                       !have spinweight with respect to this n/s flipping
     zbondio = 1.0d0 / zbondio
     Uyo = Uyo * (-zbondio /conjg(zbondio))
   end where

   

   !if (mod(it,iot2) .eq. 0 .and. it .gt. 1) then
   if (mod(it,iot2) .eq. 0 .and. it .gt. it_start) then
!     ret = gft_write ('q', time,  dble(zBondin(:,:,1)) - qs)
!     ret = gft_write ('p', time, dimag(zBondin(:,:,1)) - ps)

!     ret = gft_write ('rUh', time,  dble(Uh(:,:,1)))
!     ret = gft_write ('iUh', time, dimag(Uh(:,:,1)))

!     ret = gft_write ('rUyh', time,  dble(Uyh(:,:,1)))
!     ret = gft_write ('iUyh', time, dimag(Uyh(:,:,1)))

!     ret = gft_write ('rNB', time,  dble(NewsB(:,:,1)))
!     ret = gft_write ('iNB', time, dimag(NewsB(:,:,1)))
!     ret = gft_write ('NB', time, dble(NewsB(:,:,1)*conjg(NewsB(:,:,1))))

      ret = gft_write ('uB',  time,      uBondi(:,:,1))
      ret = gft_write ('red', time,   redshiftB(:,:,1))
  !    ret = gft_write_mask ('rNc', time,  dble(NewsB(:,:,1)), circle(:,:,1))
  !    ret = gft_write_mask ('iNc', time, dimag(NewsB(:,:,1)), circle(:,:,1))
      ret = gft_write ('rNc', time,  dble(NewsB(:,:,1)))
      ret = gft_write ('NNB', time,  dble(NewsB(:,:,1)* conjg(NewsB(:,:,1))))
      ret = gft_write ('rNcs', time,  dble(NewsB(:,:,2)))
      ret = gft_write ('iNc', time, dimag(NewsB(:,:,1)))
      ret = gft_write ('rNh', time,  dble(News(:,:,1)))
      ret = gft_write ('NNBh', time,  dble(News(:,:,1)*conjg(News(:,:,1))))
      ret = gft_write ('delta', time,  delta(:,:,1))

   end if
   write(404,*) time_of_news, maxval(abs(NewsB(:,:,1)- NewsB(:,:,2)))
   call flush(404)
   call axichecker(time_of_news, uBondi(:,:,1), 'uB')
   call axichecker(time_of_news, redshiftB(:,:,1), 'red')
   call axichecker(time_of_news, dble(NewsB(:,:,1)), 'Nr')
   call axichecker(time_of_news, dimag(NewsB(:,:,1)), 'Nc')

   if (ThisIsTheFirstTime) then
     ThisIsTheFirstTime = .false.
     open(unit=321, file="NewsLine", status='replace', FORM='unformatted')
     open(unit=322, file="TimeLine", status='replace', FORM='unformatted')
     open(unit=323, file="QLine", status='replace', FORM='unformatted')
     open(unit=324, file="DataLine", status='replace', FORM='unformatted')
     open(unit=325, file="ADataLine", status='replace')

     write(321) NewsB(ZERO:EDGE,ZERO,1)
     write(322) UBondi(ZERO:EDGE,ZERO,1)

     write(323) qs(ZERO:EDGE,ZERO)
     write(324) int(EDGE - ZERO)
     write(325,*) int(EDGE - ZERO)
     
     close(321)
     close(322)
     close(323)
     close(324)
     close(325)
  else
     open(unit=321, file="NewsLine", status='old', FORM='unformatted', position='append')
     open(unit=322, file="TimeLine", status='old', FORM='unformatted', position='append')
     write(321) NewsB(ZERO:EDGE,ZERO,1)
     write(322) UBondi(ZERO:EDGE,ZERO,1)
     close(321)
     close(322)
  endif

end subroutine news_zBondi

subroutine getnews !(zetabar, &
                   ! Jo, Jn, Jo_l, Jn_l, betao, betan, cBo, cBn, Uo, Un, &
                   ! omegao, omegan, comegao, comegan, Deltao, Deltan, &
                   ! News, redshift, uBo, uBn, uBh, zBondio, zBondin)

   use null_grid, only : nn, it, dt
   use null_params, only: time
   use horizon_eth
   implicit none

!  double complex,   dimension (nn,nn,2), intent (in) :: zetabar
!  double complex,   dimension (nn,nn,2), intent (out) :: Jo, Jn, Jo_l, Jn_l, &
!                                                         cBo, cBn, Uo, Un
!  double precision, dimension (nn,nn,2), intent (out) :: betao, betan
!  double precision, dimension (nn,nn,2), intent (inout) :: omegao, omegan
!  double complex,   dimension (nn,nn,2), intent (inout) :: comegao, comegan
!  double complex,   dimension (nn,nn,2), intent (inout) :: Deltao, Deltan
!  double complex,   dimension (nn,nn,2), intent (inout) :: News
!  double precision, dimension (nn,nn,2), intent (inout) :: redshift
!  double precision, dimension (nn,nn,2), intent (inout) :: uBo, uBn, uBh
!  double complex,   dimension (nn,nn,2), intent (inout) :: zBondio, zBondin

   ! local arrays

   double complex,   dimension (nn,nn,2) :: J, J_u, J_l, J_l_u, &
                                            cB, U, comega, eth_J, beth_J, eth_U, beth_U
   double precision, dimension (nn,nn,2) :: beta, omega

   call scrivals (Jo, Jn, Jo_l, Jn_l, betao, betan, cBo, cBn, Uo, Un)
   call news_omega (Uo, Un, omegao, omegan)
   call news_comega (omegao, omegan, comegao, comegan)
   ! mid level values, from the characteristic evolution

   J = 0.5d0 * (Jn + Jo)
   J_u = (Jn - Jo) / dt
   J_l = 0.5d0 * (Jn_l + Jo_l)
   J_l_u = (Jn_l - Jo_l) / dt
   beta = 0.5d0 * (betan + betao)
   cB = 0.5d0 * (cBn + cBo)
   U = 0.5d0 * (Un + Uo)

   call eth1 (nn, eth_U, U, 1, 1)
   call eth1 (nn, beth_U, U, 1, -1)
   call eth1 (nn, eth_J, J, 2, 1)
   call eth1 (nn, beth_J, J, 2, -1)

   time_of_news = time - .5d0 * dt    ! assumes 'time' is time on the new level

   ! mid level values, from quantities computed above

   omega = 0.5d0 * (omegan + omegao)
   comega = 0.5d0 * (comegan + comegao)


   call news_uframe (J, J_u, J_l, J_l_u, beta, cB, U, omega, &
                     comega, News)

   call news_zBondi (Uo, Un, Uyo, Uyn, zBondio, zBondin, zBondih, &
                     News, NewsB, uBondio, uBondi, uBondin, redshiftB, dMdOmega, &
                     omega, J, eth_J, beth_J, eth_U, beth_U, J_u, beta, nsmask)

   omegao = omegan
   comegao = comegan

   !if (mod(it,iot2) .eq. 0 .and. it .gt. 1) then
   if (mod(it,iot2) .eq. 0 .and. it .gt. it_start) then
!     ret = gft_write ('redshift', time, redshift(:,:,1))
!     ret = gft_write ('uB', time, uBh(:,:,1))
   end if

end subroutine getnews
subroutine news_check
  use checkpoint_defs
  implicit none
! news
  
  open(unit=NdeltaFile, file=folderu // 'Newsdelta', &
                       status='replace', FORM='unformatted')
  open(unit=NmaskFile, file=folderu // 'Newsmask', &
                       status='replace', FORM='unformatted')
  open(unit=NcomegaFile, file=folderu // 'Newscomega', &
                       status='replace', FORM='unformatted')
  open(unit=NomegaFile, file=folderu // 'Newsomega', &
                       status='replace', FORM='unformatted')
  open(unit=NuBFile, file=folderu // 'NewsuB', &
                       status='replace', FORM='unformatted')
  open(unit=NzBondiFile, file=folderu // 'NewszBondi', &
                       status='replace', FORM='unformatted')
  open(unit=NUyFile, file=folderu // 'NewsUy', &
                       status='replace', FORM='unformatted')

   write( NdeltaFile     )  deltao
   write( NmaskFile      )  nsmask
   write( NcomegaFile    )  comegao
   write( NomegaFile     )  omegao
   write( NuBFile        )  uBondio
   write( NzBondiFile    )  zBondio
   write( NUyFile        )  Uyo
 
  
  close( NdeltaFile     )
  close( NmaskFile      )
  close( NcomegaFile    )
  close( NomegaFile     )
  close( NuBFile        )
  close( NzBondiFile    )
  close( NUyFile        )

 end subroutine news_check

  subroutine recover_news
    use checkpoint_defs
    implicit none 

  open(unit=NdeltaFile, file=folderr // 'Newsdelta', &
                       status='old', FORM='unformatted')
  open(unit=NmaskFile, file=folderr // 'Newsmask', &
                       status='old', FORM='unformatted')
  open(unit=NcomegaFile, file=folderr // 'Newscomega', &
                       status='old', FORM='unformatted')
  open(unit=NomegaFile, file=folderr // 'Newsomega', &
                       status='old', FORM='unformatted')
  open(unit=NuBFile, file=folderr // 'NewsuB', &
                       status='old', FORM='unformatted')
  open(unit=NzBondiFile, file=folderr // 'NewszBondi', &
                       status='old', FORM='unformatted')
  open(unit=NUyFile, file=folderr // 'NewsUy', &
                       status='old', FORM='unformatted')

   read( NdeltaFile     )  deltao
   read( NmaskFile      )  nsmask
   read( NcomegaFile    )  comegao
   read( NomegaFile     )  omegao
   read( NuBFile        )  uBondio
   read( NzBondiFile    )  zBondio
   read( NUyFile        )  Uyo
 
  
  close( NdeltaFile     )
  close( NmaskFile      )
  close( NcomegaFile    )
  close( NomegaFile     )
  close( NuBFile        )
  close( NzBondiFile    )
  close( NUyFile        )

  initial = .false.   ! we will assume that no one will attemp a recovery before
                      ! this is set to false

  end subroutine recover_news

  subroutine axichecker(t, bnn, name)
    use null_grid
    use null_params
    implicit none
    double precision, intent(in) :: t
    double precision, dimension(nn,nn), intent(in) ::  bnn
    character*(*) :: name

    integer :: iunit1 = 991
    integer :: iunit2 = 992
    integer :: iunit3 = 993
    integer :: iunit4 = 994
    integer :: iunit5 = 995
    integer :: iunit6 = 996
    logical, save :: neverbeencalled = .true.

    integer, save :: i1, i2, i3, i4, i5, i6, id1, id2, id3, id4, id5, id6    
    double precision, save :: q1, q2, q3, q4, q5, q6
    double precision, save :: c11, c12, c13, c14, c21, c22, c23, c24, c31, c32, c33, c34    
    double precision, save :: c41, c42, c43, c44, c51, c52, c53, c54, c61, c62, c63, c64

    double precision :: int1, int2, int3, int4, int5, int6
    if (neverbeencalled) then
      neverbeencalled = .false.
      i1 = nint((.2d0 + qsize) / dd) + 3
      i2 = nint((.4d0 + qsize) / dd) + 3
      i3 = nint((.6d0 + qsize) / dd) + 3
      i4 = nint((.8d0 + qsize) / dd) + 3
      i5 = nint((1.0d0 + qsize) / dd) + 3
      i6 = nint((1.2d0 + qsize) / dd) + 3
    
      q1 = qs(i1,i1)    ! should be .2, .4, .6, .8, 1, 1.2
      q2 = qs(i2,i2)
      q3 = qs(i3,i3)
      q4 = qs(i4,i4)
      q5 = qs(i5,i5)
      q6 = qs(i6,i6)

      call get_coefs(q1, id1, c11, c12, c13, c14)
      call get_coefs(q2, id2, c21, c22, c23, c24)
      call get_coefs(q3, id3, c31, c32, c33, c34)
      call get_coefs(q4, id4, c41, c42, c43, c44)
      call get_coefs(q5, id5, c51, c52, c53, c54)
      call get_coefs(q6, id6, c61, c62, c63, c64)

    endif 

    open(unit = iunit1,  file = name // '..2.axi', position = 'append')
    open(unit = iunit2,  file = name // '..4.axi', position = 'append')
    open(unit = iunit3,  file = name // '..6.axi', position = 'append')
    open(unit = iunit4,  file = name // '..8.axi', position = 'append')
    open(unit = iunit5,  file = name // '.1.axi', position = 'append')
    open(unit = iunit6,  file = name // '.1.2.axi', position = 'append')

    int1 = bnn(id1 -1, id1 -1) * c11 + bnn (id1, id1) * c12 + &
           bnn(id1+1, id1+1) * c13 + bnn(id1+2, id1+2) * c14 
 
    int2 = bnn(id2 -1, id2 -1) * c21 + bnn (id2, id2) * c22 + &
           bnn(id2+1, id2+1) * c23 + bnn(id2+2, id2+2) * c24 

    int3 = bnn(id3 -1, id3 -1) * c31 + bnn (id3, id3) * c32 + &
           bnn(id3+1, id3+1) * c33 + bnn(id3+2, id3+2) * c34 

    int4 = bnn(id4 -1, id4 -1) * c41 + bnn (id4, id4) * c42 + &
           bnn(id4+1, id4+1) * c43 + bnn(id4+2, id4+2) * c44 

    int5 = bnn(id5 -1, id5 -1) * c51 + bnn (id5, id5) * c52 + &
           bnn(id5+1, id5+1) * c53 + bnn(id5+2, id5+2) * c54 

    int6 = bnn(id6 -1, id6 -1) * c61 + bnn (id6, id6) * c62 + &
           bnn(id6+1, id6+1) * c63 + bnn(id6+2, id6+2) * c64 

    if (abs(int1) > 0.0) then
      write( unit = iunit1, fmt = * )  t, (bnn(i1, ZERO) - int1) / int1
    else
      write( unit = iunit1, fmt = * )  t, (bnn(i1, ZERO) - int1)
    endif
    if (abs(int2) > 0.0) then
      write( unit = iunit2, fmt = * )  t, (bnn(i2, ZERO) - int2) / int2 
    else
      write( unit = iunit2, fmt = * )  t, (bnn(i2, ZERO) - int2)
    endif
    if (abs(int3) > 0.0) then
      write( unit = iunit3, fmt = * )  t, (bnn(i3, ZERO) - int3) / int3
    else
      write( unit = iunit3, fmt = * )  t, (bnn(i3, ZERO) - int3)
    endif
    if (abs(int4) > 0.0) then
      write( unit = iunit4, fmt = * )  t, (bnn(i4, ZERO) - int4) / int4
    else
      write( unit = iunit4, fmt = * )  t, (bnn(i4, ZERO) - int4)
    endif
    if (abs(int5) > 0.0) then
      write( unit = iunit5, fmt = * )  t, (bnn(i5, ZERO) - int5) / int5
    else
      write( unit = iunit5, fmt = * )  t, (bnn(i5, ZERO) - int5)
    endif
    if (abs(int6) > 0.0) then
      write( unit = iunit6, fmt = * )  t, (bnn(i6, ZERO) - int6) / int6
    else
      write( unit = iunit6, fmt = * )  t, (bnn(i6, ZERO) - int6)
    endif

    close(unit = iunit1)
    close(unit = iunit2)
    close(unit = iunit3)
    close(unit = iunit4)
    close(unit = iunit5)
    close(unit = iunit6)

  end subroutine
  subroutine get_coefs(qwant, igot, c1, c2, c3, c4)
    use null_grid
    use null_params
    implicit none
    double precision, intent(in) :: qwant
    integer, intent(out) :: igot
    double precision, intent(out) :: c1, c2, c3, c4

    double precision :: xi, xim1, xip1, xip2
    igot = (qwant/sqrt(2.0d0) + qsize)/ dd  + 3
    xi   = sqrt(qs(igot  ,igot  )**2 + ps(igot  ,igot  )**2)
    xim1 = sqrt(qs(igot-1,igot-1)**2 + ps(igot-1,igot-1)**2)
    xip1 = sqrt(qs(igot+1,igot+1)**2 + ps(igot+1,igot+1)**2)
    xip2 = sqrt(qs(igot+2,igot+2)**2 + ps(igot+2,igot+2)**2)

    c1 = (qwant - xi)*(qwant - xip1)*(qwant - xip2) / ( &
           (xim1    - xi)*(xim1    - xip1)*(xim1    - xip2) )

    c2 = (qwant - xim1)*(qwant - xip1)*(qwant - xip2) / ( &
           (xi    - xim1)*(xi    - xip1)*(xi    - xip2) )

    c3 = (qwant - xim1)*(qwant - xi)*(qwant - xip2) / ( &
           (xip1    - xim1)*(xip1    - xi)*(xip1    - xip2) )

    c4 = (qwant - xim1)*(qwant - xi)*(qwant - xip1) / ( &
           (xip2    - xim1)*(xip2    - xi)*(xip2    - xip1) )
  end subroutine get_coefs
 
end module news2
