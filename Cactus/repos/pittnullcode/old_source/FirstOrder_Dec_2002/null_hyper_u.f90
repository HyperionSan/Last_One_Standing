module null_hyper_u

   implicit none

   double complex, dimension (:,:), allocatable, save, private :: &
      J, nu, ck, cb, dx_J, dx_nu, dx_ck, dx_cb, dx_U, eth_J, eth_dx_J, &
      n_Q, rhs, Q

   double precision, dimension (:,:), allocatable, save, private :: &
      K, dx_K, beta

   integer, dimension (:,:), allocatable, save, private :: mask

contains

subroutine null_hyper_u_allocate

   use null_grid

   allocate (J(nn,nn), nu(nn,nn), ck(nn,nn), cb(nn,nn), dx_J(nn,nn), &
             dx_nu(nn,nn), dx_ck(nn,nn), dx_cb(nn,nn), dx_U(nn,nn), &
             eth_J(nn,nn), eth_dx_J(nn,nn), n_Q(nn,nn), rhs(nn,nn), &
             Q(nn,nn), K(nn,nn), dx_K(nn,nn), beta(nn,nn), mask(nn,nn))

end subroutine null_hyper_u_allocate

subroutine null_u (i, jns, nuns, ckns, bns, cbns, uns, masks,switch)

   use null_grid
   use null_mask
   use null_eth
   use particle

   implicit none

   integer,                             intent (in)    :: i
   double complex,   dimension (:,:,:), intent (in)    :: jns, nuns, ckns, cbns
   double precision, dimension (:,:,:), intent (in)    :: bns 
   double complex,   dimension (:,:,:), intent (inout) :: uns
   integer,          dimension (:,:),   intent (in)    :: masks
   logical,                             intent (in)    :: switch

   double precision :: xhere,pr

   call null_remask (i, masks+1, mask)

   if (i .lt. nx) then

      xhere = xh(i-1)

      J = 0.5 * ( jns(:,:,i-1) + jns(:,:,i) )
      dx_J = (jns(:,:,i) - jns(:,:,i-1)) / dx

      K = sqrt(1. + J * conjg(J))
      dx_K = dble(dx_J * conjg(J)) / K

      nu = 0.5 * ( nuns(:,:,i-1) + nuns(:,:,i) )
      dx_nu = (nuns(:,:,i) - nuns(:,:,i-1)) / dx

      ck = 0.5 * ( ckns(:,:,i-1) + ckns(:,:,i) )
      dx_ck = (ckns(:,:,i) - ckns(:,:,i-1)) / dx

      cb = 0.5 * ( cbns(:,:,i-1) + cbns(:,:,i) )
      dx_cb = (cbns(:,:,i) - cbns(:,:,i-1)) / dx

      call null_d1 (eth_J, J, 2, 1)
      call null_d1 (eth_dx_J, dx_J, 2, 1)

      n_Q = 2. * conjg(J) * conjg(ck) * dx_J * J &
          + 3. * K * conjg(nu) * (dx_J * K - dx_K * J) &
          + 4. * J * ck * dx_K * conjg(J) &
          - conjg(J) * conjg(nu) * dx_J * J &
          - 2. * K * conjg(ck) * dx_K * J &
          - 4. * K * ck * dx_J * conjg(J) &
          - 2. * K * ck * conjg(dx_J) * J &
          + J**2 * conjg(eth_J) * dx_K &
          + J**2 * conjg(nu) * conjg(dx_J) &
          + 2. * conjg(J)**2 * eth_J * dx_J &
          + 2. * K**2 * dx_K * (ck + nu) &
          + K**2 * eth_J * conjg(dx_J) &
          - K * conjg(eth_J) * dx_J * J &
          - conjg(J) * nu * (dx_K * J + dx_J * K) &
          - 3. * conjg(J) * eth_J * dx_K * K &
          + 2. * dx_ck + 2. * dx_nu &
          - 2. * K * (dx_nu + dx_ck) &
          + 2. * J * conjg(dx_ck) &
          + 2. * conjg(J) * eth_dx_J

      rhs = - xhere * (1. - xhere) * (dx_nu + dx_ck - 2. * dx_cb) &
          - 4. * cb + 0.5 * (1. - xhere) * (xhere) * n_Q  

   if(switch) then
      pr=0.0d0 !The pressure - adjust if necessary
      rhs = rhs+ 16. * PI * xhere/(1.-xhere) * (density(:,:,i)+pr &
          +density(:,:,i-1)+pr)/2 * p_vdo(1) * p_van
   endif

      xhere = x(i-1)

      J = jns(:,:,i-1)
      K = sqrt(1. + J * conjg(J))
      beta = bns(:,:,i-1)
      dx_U = (uns(:,:,i-1) - uns(:,:,i-2)) / dx
      Q = rwt * xhere ** 2 * exp(-2. * beta) * (J * conjg(dx_U) + K * dx_U)

      xhere = xh(i-1)
      Q = ((xhere * (1. - xhere) - dx) * Q +  dx * rhs) &
         / (xhere * (1. - xhere) + dx)

   else if (i .eq. nx) then

      cb = cbns(:,:,nx)
      Q = -2. * cb

   end if

   xhere = x(i)

   J = jns(:,:,i)
   K = sqrt(1. + J * conjg(J))
   beta = bns(:,:,i)

   uns(:,:,i) = uns(:,:,i) * (1 - mask) + mask * (uns(:,:,i-1) &
   + dx * exp(2. * beta) / (rwt * xhere ** 2) * (- J * conjg(Q) + K * Q) )

end subroutine null_u

end module null_hyper_u
