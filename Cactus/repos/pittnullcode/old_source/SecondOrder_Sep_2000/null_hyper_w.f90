module null_hyper_w

   implicit none

   double complex, dimension (:,:), allocatable, save, private :: &
       & j, eth_j, ethb_j, ethb2_j, ethb_k, eth_ethb_k, eth_beta, &
       & eth2_beta, eth_ethb_beta, u, dx_u, ethb_u, ethb_dx_u, &
       & van

   double precision, dimension (:,:), allocatable, save, private :: &
       & k, beta, ricci, rhs, n_w 

   integer, dimension (:,:), allocatable, save, private :: mask

contains

subroutine null_hyper_w_allocate

   use null_grid, only : nn

   allocate(j(nn,nn), eth_j(nn,nn), ethb_j(nn,nn), ethb2_j(nn,nn),&
         & ethb_k(nn,nn), eth_ethb_k(nn,nn), eth_beta(nn,nn),&
         & eth2_beta(nn,nn), eth_ethb_beta(nn,nn), u(nn,nn), dx_u(nn&
         & ,nn), ethb_u(nn,nn), ethb_dx_u(nn,nn), k(nn,nn), beta(nn&
         & ,nn), ricci(nn,nn), rhs(nn,nn), n_w(nn,nn), mask(nn,nn) )

end subroutine null_hyper_w_allocate

subroutine null_w (i, jns, bns, uns, wns, masks)

   use null_grid
   use null_mask
   use null_eth
 
   implicit none

   integer,                             intent (in)    :: i
   double complex,   dimension (:,:,:), intent (in)    :: jns
   double precision, dimension (:,:,:), intent (in)    :: bns 
   double complex,   dimension (:,:,:), intent (in)    :: uns 
   double precision, dimension (:,:,:), intent (inout) :: wns
   integer,          dimension (:,:),   intent (in)    :: masks

   double precision xhere

   call null_remask (i, masks+1, mask)

   xhere = xh(i-1)

   beta = 0.5 * (bns(:,:,i) + bns(:,:,i-1))

   call null_d1 (eth_beta, dcmplx(beta), 0, 1)
   call null_d2 (eth_ethb_beta, dcmplx(beta), 0, 1, -1)
   call null_d2 (eth2_beta, dcmplx(beta), 0, 1, 1)

   u = uns(:,:,i-1)
   dx_u = (uns(:,:,i) - uns(:,:,i-2)) / (2. * dx)

   call null_d1 (ethb_u, u, 1, -1)
   call null_d1 (ethb_dx_u, dx_u, 1, -1)

   j = 0.5 * (jns(:,:,i) + jns(:,:,i-1))
   k = sqrt (1. + dble(j) ** 2 + dimag(j) ** 2)

   call null_d1 (eth_j, j, 2, 1)
   call null_d1 (ethb_j, j, 2, -1)

   call null_d2 (ethb2_j, j, 2, -1, -1)
   call null_d1 (ethb_k, dcmplx(k), 0, -1)
   call null_d2 (eth_ethb_k, dcmplx(k), 0, 1, -1)

   ricci = 2. * k + dble(ethb2_j) - dble(eth_ethb_k) + (&
            & dble(eth_j)  ** 2 + dimag(eth_j)  ** 2 - dble(ethb_j)&
            & ** 2 - dimag(ethb_j) ** 2 ) / (4. * k)

   n_w = exp(2. * beta) * ( (1. - k) * ( dble(eth_ethb_beta) +&
            & dble(eth_beta) ** 2 + dimag(eth_beta) ** 2 ) + dble( j&
            & * ( conjg(eth_beta) **2 + conjg(eth2_beta) ) ) + dble(&
            & eth_beta * ( conjg(ethb_j) - ethb_k ) )  )

   n_w = - xhere ** 4 * rwt ** 2 * exp(-2. * beta) * ( dble(&
            & conjg(j) * dx_u ** 2 ) + k * ( dble(dx_u) ** 2 +&
            & dimag(dx_u) ** 2 ) ) / 4. + n_w

   rhs = (1. - xhere) / (rwt * xhere) * (0.5 * exp(2. * beta) *&
            & ricci - 1. - exp(2.*beta) * (dble(eth_beta) ** 2 +&
            & dimag(eth_beta) ** 2 + eth_ethb_beta ) ) +  0.25 *&
            & xhere * (1. - xhere) * 2. * dble(ethb_dx_u) +  2. *&
            & dble(ethb_u) + (1. - xhere) / (rwt * xhere) *  n_w

   wns(:,:,i) = wns(:,:,i) * (1 - mask) + mask * ((xhere * (1. -&
            & xhere) - dx) * wns(:,:,i-1) + rhs * dx) / (xhere * (1.&
            & - xhere) + dx)

end subroutine null_w

end module null_hyper_w
