      MODULE CALC_NEWSNORTH
	implicit none
	contains

	subroutine calcnewsnorth(news, un, uo, betan, betao)

      use null_grid
      use null_eth
      use null_newsvarnorth
      use null_vars

! news <-- ::wn, wo , wom1, alphan, alphao, alphaom1
! jscri <-- :: jn, jo, jom1, j_ln, j_lo, j_lom1

! actually alpha stands for delta= alpha - gamma

      implicit none

	double complex, dimension(nn,nn) :: news
	double complex, dimension(nn,nn) :: j, un, uo, u
	double precision, dimension(nn,nn) :: k,k_u, k_l, k_l_u,alpha, betan, betao, beta

!temporary arrays

	double precision, dimension(nn,nn) :: a, w

	double complex, dimension(nn,nn) :: eth_j, ethb_j, eth_k, eth_u, eth_ub,&
     			 eth_k_l, eth_ethb_w,eth_beta,&
     			 eth_a, eth2_a, eth_ethb_a, eth_w, ethb_w, ethb_u,&
                         eth_j_l, ethb_j_l, eth2_w,eth2_beta,eth_ethb_beta,&
                         eth_alpha

	double complex, dimension(nn,nn) :: t1, t2
	double complex, dimension(nn,nn) :: j_u, zeta

	double complex, dimension(nn,nn) :: s1, s2, s3, s4, s5
	double complex, dimension(nn,nn):: j_l_u,  j_l

	double complex, dimension(nn,nn,nx):: temp

!************** eval de derivatives of j!

      w = .5*(nwo+nwn)

      u = .5*(uo+un)
      beta = .5*(betao+betan)
      j = .5*(jon(:,:,nx)+jnn(:,:,nx))
      k = sqrt(1.+j*conjg(j) )
 
      alpha = .5*(alphaon+alphann)

!    call wt_rsmooth(nn, w)
!    call wt_rsmooth(nn,alpha)   

      j_u = (jnn(:,:,nx) - jon(:,:,nx)) / dt
      j_l = .5*(j_lon+j_lnn)

!      j_l_u = (j_lnn - j_lom1n)/ (2.* dt)

	temp =  (jnn - jon)/ dt

	j_l_u =  - .5 * (3. * temp(:,:,nx) - 4. * temp(:,:,nx-1) &
     &                           + temp(:,:,nx-2) ) / dx

!        call wt_csmooth(nn,j_l_u)

      k_u = real( j_u * conjg(j) ) / k
      k_l = real( j_l * conjg(j) ) / k
      k_l_u = real( j_u * conjg(j_l) + j_l_u * conjg(j) )/ k - &
     &        k_l * k_u / k 


      call null_d1 (eth_k, dcmplx(k), 0, 1)
      call null_d1 (eth_u, u, 1, 1)
      call null_d1 (eth_ub, conjg(u), -1, 1)
      call null_d1 (eth_j, j, 2, 1)
      call null_d1 (ethb_j, j, 2, -1)
      call null_d1 (eth_j_l, j_l, 2, 1)

      call null_d1 (ethb_j_l, j_l, 2, -1)
      call null_d1 (eth_k_l, dcmplx(k_l), 0, 1)

      call null_d1 (eth_w, dcmplx(w), 0, 1)
      call null_d2 (eth2_w, dcmplx(w), 0, 1, 1)
      call null_d2 (eth_ethb_w, dcmplx(w), 0, 1, -1)
      call null_d1 (eth_beta, dcmplx(beta), 0, 1)
      call null_d2 (eth2_beta, dcmplx(beta), 0, 1, 1)
      call null_d2 (eth_ethb_beta, dcmplx(beta), 0, 1, -1)
      call null_d1 (eth_alpha, dcmplx(alpha), 0, 1)

	 ethb_u = conjg(eth_ub)


! news *****************************************
!
	
	a = w * exp(2. * beta)

	eth_a = exp(2. * beta) * ( eth_w + 2. * w * eth_beta )
	
	eth2_a = exp(2. * beta) * ( 4. * eth_beta * eth_w +  &
     &                              4. * w * eth_beta**2 + eth2_w + &
     &				    2. * w * eth2_beta )


	eth_ethb_a = exp(2. * beta) * ( 4. * real(eth_beta * conjg(eth_w)) +  &
     &                              4. * w * eth_beta * conjg(eth_beta) + &
     &				    eth_ethb_w + &
     &				    2. * w * eth_ethb_beta )





	s1 =   ( -2. * k_l_u * j * (k+1.) +  &
     &	          j_l_u * (k+1.)**2 +  &
     &            conjg( j_l_u) * j**2 ) /  ( k + 1.)



	s3 = ( j_l * (k+1.)**2    &
     &	       -2. * k_l * j * (k+1.) +  &
     &         conjg( j_l ) * j**2 ) / (k + 1.)


	s4 = 0.5 / ( k + 1.) * ( &
     &          eth_a * eth_w * (k+1.)**2  &
     &		- (k+1.) * j * 2.* real( eth_a * conjg(eth_w) )  &
     &          + j**2 * conjg(eth_a * eth_w) )


	s2 = 0.5 / ( k + 1.) * ( &
     &	     (k+1.)* (eth_j_l *conjg(u) * (k+1.) - 2.* eth_k_l * j *conjg(u) ) &
     &       + eth_u * (k+1.)* ( -2. * j * conjg(j_l) + k_l * 2. * (k+1.) ) &
     &       + conjg(ethb_u) * (k+1.) * ( -2.* j * k_l + j_l * 2. * (k+1.) ) &
     &       + ethb_j_l * u * (k+1.)**2 - conjg(eth_k_l) * 2. * u * j * (k+1.) &
     &       + ethb_u * 2. * j * ( j * conjg(j_l) - (k+1.) * k_l)  &
     &       + j**2 * ( u * conjg(eth_j_l) + conjg(ethb_j_l * u) ) &
     &       + j * 2. * conjg(eth_u) * ( j * k_l - j_l * (k+1.) )  )


	s5 = 0.25 / ( k + 1.) * ( &
     &          2. * eth2_a * (k+1.)**2 +  2. * j**2 * conjg(eth2_a) &
     &        - 4. * eth_ethb_a * j * (k+1.) + &
     &                            conjg(j) * eth_a * eth_j* (k+1.)**2 &
     &       + j * eth_a * conjg(ethb_j) * (k+1.)**2 - &
     &             eth_a * eth_k * 2. * (k+1.) * ( j*conjg(j) + (k+1.) ) & 
     &            + eth_a * ethb_j * (k+1.) * ( -j*conjg(j) + (k+1.) ) &
     &       - j**2 * eth_a * conjg(eth_j) * k +  &
     &                            j**2 * conjg(j) * 2.* eth_a * conjg(eth_k) &
     &       - conjg(eth_a) * eth_j * (k+1.) * ( j*conjg(j) + k+1. ) &
     &       - conjg(ethb_j) * conjg(eth_a) * j**2 * ( k + 2.) &
     &       + j * 2. * (k+1.)**2 * eth_k * conjg(eth_a)  &
     &       + j**2 * conjg(j) * ethb_j * conjg(eth_a) &
     &       + j**3 * conjg(eth_a * eth_j) - 2.* j**2 *k*conjg(eth_k * eth_a) )

	
	news = 0.25 * (cos(-2.*alpha)+ii*sin(-2.*alpha)) * ( s1 + s2 + &
     &                              0.5 * real(ethb_u) * s3  &
     &                             - 4. * s4 /w**2 + 2. * s5/w ) / &
     &                  ( w**2 * exp(2. * beta) )


!        call wt_csmooth(nn,news)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!_

	return
	end subroutine calcnewsnorth

      end module calc_newsnorth
