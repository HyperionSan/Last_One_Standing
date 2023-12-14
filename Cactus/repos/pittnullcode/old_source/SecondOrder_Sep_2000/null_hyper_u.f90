module null_hyper_u

  implicit none

  integer, dimension (:,:), allocatable, save, private :: mask

  double complex,   dimension (:,:), allocatable, save, private :: &
       & j, jb, dx_j, dx_jb, eth_j, ethb_j, ethb_jb, eth_dx_j, &
       & ethb_dx_j, eth_k, ethb_k, eth_dx_k, ethb_dx_k, eth_beta, &
       & eth_dx_beta, rhs, n_q, q, dx_u

  double precision, dimension (:,:), allocatable, save, private :: &
       & k, dx_k, beta, dx_beta

contains

subroutine null_hyper_u_allocate

   use null_grid, only : nn
    
   allocate (j(nn,nn), jb(nn,nn), dx_j(nn,nn), dx_jb(nn,nn),&
         & eth_j(nn,nn), ethb_j(nn,nn), ethb_jb(nn,nn), eth_dx_j(nn&
         & ,nn), ethb_dx_j(nn,nn), eth_k(nn,nn), ethb_k(nn,nn),&
         & eth_dx_k(nn,nn), ethb_dx_k(nn,nn), eth_beta(nn,nn),&
         & eth_dx_beta(nn,nn), rhs(nn,nn), n_q(nn,nn), q(nn,nn),&
         & dx_u(nn,nn), k(nn,nn), dx_k(nn,nn), beta(nn,nn),&
         & dx_beta(nn,nn), mask(nn,nn) )

end subroutine null_hyper_u_allocate

subroutine null_u (i, jns, bns, uns, masks)

   use null_grid
   use null_mask
   use null_eth

   implicit none

   integer,                             intent (in)    :: i
   double complex,   dimension (:,:,:), intent (in)    :: jns 
   double precision, dimension (:,:,:), intent (in)    :: bns 
   double complex,   dimension (:,:,:), intent (inout) :: uns
   integer,          dimension (:,:),   intent (in)    :: masks

   double precision :: xhere

   call null_remask (i, masks+1, mask)

   if (i .eq. nx) then

      beta = bns(:,:,nx)
      call null_d1 (eth_beta, dcmplx(beta), 0, 1)

      q = -2. * eth_beta

   else if (i .lt. nx) then

      xhere = xh(i-1)

      j = 0.5 * ( jns(:,:,i-1) + jns(:,:,i) )
      jb = conjg(j)
      dx_j = (jns(:,:,i) - jns(:,:,i-1)) / dx
      dx_jb = conjg(dx_j)
      k = sqrt(1. + dble(j) ** 2 + dimag(j) ** 2)
      dx_k = ( sqrt(1. + dble(jns(:,:,i))   ** 2 + dimag(jns(:,:,i))   ** 2) &
           - sqrt(1. + dble(jns(:,:,i-1)) ** 2 + dimag(jns(:,:,i-1)) ** 2) &
           ) / dx
      beta = 0.5 * ( bns(:,:,i-1) + bns(:,:,i) )
      dx_beta = (bns(:,:,i) - bns(:,:,i-1)) / dx

      call null_d1 (eth_j, j, 2, 1)
      call null_d1 (ethb_j, j, 2, -1)
      call null_d1 (eth_dx_j, dx_j, 2, 1)
      call null_d1 (ethb_dx_j, dx_j, 2, -1)
      call null_d1 (eth_k, dcmplx(k), 0, 1)
      call null_d1 (eth_dx_k, dcmplx(dx_k), 0, 1)
      call null_d1 (eth_beta, dcmplx(beta), 0, 1)
      call null_d1 (eth_dx_beta, dcmplx(dx_beta), 0, 1)

      ethb_jb = conjg(eth_j)
      ethb_k = conjg(eth_k)
      ethb_dx_k = conjg(eth_dx_k)

      n_q =(2.*conjg(j)*conjg(eth_k)*dx_j*j-3.*k*conjg(ethb_j)*dx_k&
           & *j+4.*j* eth_k*dx_k*conjg(j)-conjg(j)*conjg(ethb_j)&
           & *dx_j*j-2.*k*conjg(eth_k) *dx_k*j-4.*k*eth_k*dx_j&
           & *conjg(j)-2.*k*eth_k*conjg(dx_j)*j+3.*k**2*&
           & conjg(ethb_j)*dx_j+j**2*conjg(eth_j)*dx_k+j**2&
           & *conjg(ethb_j)*conjg( dx_j)+2.*conjg(j)**2*eth_j*dx_j&
           & +2.*k**2*eth_k*dx_k+2.*k**2*ethb_j*dx_k+ k**2*eth_j&
           & *conjg(dx_j)-k*conjg(eth_j)*dx_j*j-conjg(j)*ethb_j*dx_k&
           & * j-conjg(j)*ethb_j*dx_j*k-3.*conjg(j)*eth_j*dx_k*k+2.&
           & *eth_dx_k+2.* ethb_dx_j-2.*k*ethb_dx_j-2.*k*eth_dx_k+2.&
           & *j*conjg(eth_dx_k)+2.* conjg(j)*eth_dx_j)*rwt*xhere**2&
           & * 0.5


      rhs = - xhere * (1. - xhere) * (ethb_dx_j + eth_dx_k - 2. *&
           & eth_dx_beta) - 4. * eth_beta + (1. - xhere) / (rwt *&
           & xhere) * n_q  

      xhere = x(i-1)

      j = jns(:,:,i-1)
      k = sqrt(1. + dble(j) ** 2 + dimag(j) ** 2)
      beta = bns(:,:,i-1)
      dx_u = (uns(:,:,i-1) - uns(:,:,i-2)) / dx
      q = rwt * xhere ** 2 * exp(-2. * beta) * ( j * conjg(dx_u) + k * dx_u)

      xhere = xh(i-1)
      q = ((xhere * (1. - xhere) - dx) * q +  dx * rhs) / (xhere *&
           & (1. - xhere) + dx)

   end if

   xhere = x(i)

   j = jns(:,:,i)
   k = sqrt(1. + dble(j) ** 2 + dimag(j) ** 2)
   beta = bns(:,:,i)

   uns(:,:,i) = uns(:,:,i) * (1 - mask) + mask * (uns(:,:,i-1) + dx&
        & * exp(2. * beta) / (rwt * xhere ** 2) * (- j * conjg(q) + k * q) )

end subroutine null_u

end module null_hyper_u
