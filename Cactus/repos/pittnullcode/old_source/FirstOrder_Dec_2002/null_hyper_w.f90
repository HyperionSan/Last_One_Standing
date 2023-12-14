module null_hyper_w

   double complex, dimension (:,:), allocatable, save, private :: &
      U, dx_U, ck, cB, nu, J, ethb_U, ethb_dx_U, eth_J, ethb_nu, ethb_ck, &
      eth_cb, ethb_cb

   double precision, dimension (:,:), allocatable, save, private :: &
      beta, K, Ricci, rhs

   integer, dimension (:,:), allocatable, save, private :: mask

contains

subroutine null_hyper_w_allocate

   use null_grid

   allocate (U(nn,nn), dx_U(nn,nn), ck(nn,nn), cB(nn,nn), nu(nn,nn), &
             J(nn,nn), ethb_U(nn,nn), ethb_dx_U(nn,nn), eth_J(nn,nn), &
             ethb_nu(nn,nn), ethb_ck(nn,nn), eth_cb(nn,nn), ethb_cb(nn,nn), &
             beta(nn,nn), K(nn,nn), Ricci(nn,nn), rhs(nn,nn), mask(nn,nn))

end subroutine null_hyper_w_allocate

subroutine null_w (i, jns, nuns, ckns, bns, cbns, uns, wns, masks,switch)

   use null_grid
   use null_mask
   use null_eth
   use particle
 
   implicit none

   integer,                             intent (in)    :: i
   double complex,   dimension (:,:,:), intent (in)    :: jns, nuns, ckns
   double precision, dimension (:,:,:), intent (in)    :: bns 
   double complex,   dimension (:,:,:), intent (in)    :: cbns, uns 
   double precision, dimension (:,:,:), intent (inout) :: wns
   integer,          dimension (:,:),   intent (in)    :: masks
   logical,                             intent (in)    :: switch

   double precision xhere,pr

   call null_remask (i, masks+1, mask)

   xhere = xh(i-1)

   U = uns(:,:,i-1)
   dx_U = (uns(:,:,i) - uns(:,:,i-2)) / (2. * dx)

   ck   = 0.5 * (ckns(:,:,i) + ckns(:,:,i-1))
   cB   = 0.5 * (cbns(:,:,i) + cbns(:,:,i-1))
   nu   = 0.5 * (nuns(:,:,i) + nuns(:,:,i-1))
   beta = 0.5 * ( bns(:,:,i) +  bns(:,:,i-1))
   J    = 0.5 * ( jns(:,:,i) +  jns(:,:,i-1))
   K    = sqrt(1. + J * conjg(J))

   call null_d1 (ethb_U,    U ,   1, -1)
   call null_d1 (ethb_dx_U, dx_U, 1, -1)

   call null_d1 (eth_J,   J,  2,  1)
   call null_d1 (ethb_nu, nu, 1, -1)
   call null_d1 (ethb_ck, ck, 1, -1)
   call null_d1 (ethb_cb, cb, 1, -1)
   call null_d1 (eth_cb,  cb, 1,  1)

   Ricci = dble( 2. * k + ethb_nu - ethb_ck &
         + ( eth_J * conjg(eth_J) - nu * conjg(nu) ) / (4. * k) )

   rhs = (1. - xhere) / (rwt * xhere) * ( exp(2. * beta) * (0.5 * Ricci &
      - K * (ethb_cb + cb * conjg(cb)) + conjg(J) * (eth_cb + cb**2) &
      + conjg(cb) * (nu - ck) ) - 1.) &
      + 2. * ethb_U + 0.5 * xhere * (1. - xhere) * ethb_dx_U &
      - 0.25 * exp(-2. * beta) * rwt * xhere**3 * (1. - xhere) &
      * conjg(dx_U) * (K * dx_U + J * conjg(dx_U))

   if(switch) then
   pr=0.0d0 !The pressure - adjust if necessary
   rhs = rhs - 4. * PI * ( (1.-xhere)/xhere * &
        (density(:,:,i)+pr + density(:,:,i-1)+pr)/2 * ( K *  &
        (Real(p_van)**2 + Aimag(p_van)**2) - real( p_van**2 * conjg(J) ) ) &
         + xhere/(1.-xhere) * &
         (density(:,:,i)-pr + density(:,:,i-1)-pr)/2  ) * exp(2.*beta)
   endif


   wns(:,:,i) = wns(:,:,i) * (1 - mask) &
      + mask * ((xhere * (1. - xhere) - dx) * wns(:,:,i-1) + rhs * dx) &
              / (xhere * (1. - xhere) + dx)

end subroutine null_w

end module null_hyper_w
