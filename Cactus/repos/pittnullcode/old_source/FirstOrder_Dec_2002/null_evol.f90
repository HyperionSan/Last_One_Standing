module null_evol 

   implicit none

   double precision, dimension (:,:), allocatable, save, private :: &
      xc, beta, dx_beta, W, K, dx_K, dx_W, f1, f2, P_u

   double complex, dimension (:,:), allocatable, save, private :: &
      J, U, cB, cK, d2x_jave, dx_U, dx_j, dx_jave, dx_mu, dx_nu, eth_J, eth_U, &
      ethb_U, eth_cB, eth_dx_J, eth_dx_U, ethb_dx_U, ethb_cB,fac, h_c, jave, &
      mu, nu, J_H

   integer, dimension (:,:), allocatable, save, private :: mask

contains

subroutine null_evol_allocate

   use null_grid, only : nn

   allocate ( J(nn,nn), J_H(nn,nn), Jave(nn,nn), K(nn,nn), P_u(nn,nn), &
              U(nn,nn), W(nn,nn), beta(nn,nn), cB(nn,nn), cK(nn,nn), &
              d2x_Jave(nn,nn), dx_J(nn,nn), dx_Jave(nn,nn), dx_K(nn,nn), &
              dx_U(nn,nn), dx_W(nn,nn), dx_beta(nn,nn), dx_mu(nn,nn), &
              dx_nu(nn,nn), eth_J(nn,nn), eth_U(nn,nn), ethb_U(nn,nn), &
              eth_cB(nn,nn), eth_dx_J(nn,nn), &
              eth_dx_U(nn,nn), ethb_dx_U(nn,nn), ethb_cB(nn,nn), f1(nn,nn), &
              f2(nn,nn), fac(nn,nn), h_c(nn,nn), mask(nn,nn), mu(nn,nn), &
              nu(nn,nn), xc(nn,nn) )

end subroutine null_evol_allocate

subroutine null_j (i, jns, nuns, ckns, bns, cbns, uns, wns, &
                      jos, nuos, ckos, bos, cbos, uos, wos, &
                   masks,switch)

   use null_grid, only : nn, nx, dd, dx, dt, rwt, x, xh, pi, ii
   use null_mask
   use null_eth
   use null_params, only : disip
   use particle

   implicit none

   integer,                             intent (in)    :: i
   double complex,   dimension (:,:,:), intent (inout) :: jns
   double complex,   dimension (:,:,:), intent (in)    :: jos, &
      nuns, nuos, ckns, ckos, cbns, cbos, uns, uos
   double precision, dimension (:,:,:), intent (in)    :: bns, bos, wns, wos
   integer,          dimension (:,:),   intent (in)    :: masks
   double precision :: rh,pr
   logical :: switch
   double precision :: disip_rescaled

   integer cran
   integer, parameter :: cran_max = 2

   call null_remask (i, masks+1, mask)

   disip_rescaled = disip * dx
   xc = xh(i-1)

   U    = 0.25 * (uns(:,:,i-1) + uos(:,:,i-1) + uns(:,:,i-2) + uos(:,:,i))
   dx_U = 0.5 * (uns(:,:,i-1) - uns(:,:,i-2) + uos(:,:,i) - uos(:,:,i-1)) / dx

   J    = 0.5 * ( jns(:,:,i-1) +  jos(:,:,i))
   beta = 0.5 * ( bns(:,:,i-1) +  bos(:,:,i))
   cB   = 0.5 * (cbns(:,:,i-1) + cbos(:,:,i))
   cK   = 0.5 * (ckns(:,:,i-1) + ckos(:,:,i))
   nu   = 0.5 * (nuns(:,:,i-1) + nuos(:,:,i))
   W    = 0.5 * ( wns(:,:,i-1) +  wos(:,:,i))

   K    = sqrt (1. + J * conjg(J))

   if (i .ne. nx) then

      dx_J = 0.5 * (jns(:,:,i-1) - jns(:,:,i-2) &
                  + jos(:,:,i+1) - jos(:,:,i)) / dx

      dx_beta = 0.5 * (bns(:,:,i-1) - bns(:,:,i-2) &
                     + bos(:,:,i+1) - bos(:,:,i)) / dx

      dx_W = 0.5 * (wns(:,:,i-1) - wns(:,:,i-2) &
                  + wos(:,:,i+1) - wos(:,:,i)) / dx

      dx_nu = 0.5 * (nuns(:,:,i-1) - nuns(:,:,i-2) &
                   + nuos(:,:,i+1) - nuos(:,:,i)) / dx

   else

      dx_J = 0.5 * (jns(:,:,i-1) - jns(:,:,i-2) &
      + 0.5 * (3. * jos(:,:,i) - 4. * jos(:,:,i-1) + jos(:,:,i-2))) / dx

      dx_beta = 0.5 * (bns(:,:,i-1) - bns(:,:,i-2) &
         + 0.5 * (3. * bos(:,:,i) - 4. * bos(:,:,i-1) + bos(:,:,i-2))) / dx

      dx_W = 0.5 * (wns(:,:,i-1) - wns(:,:,i-2) &
      + 0.5 * (3. * wos(:,:,i) - 4. * wos(:,:,i-1) + wos(:,:,i-2))) / dx

      dx_nu = 0.5 * (nuns(:,:,i-1) - nuns(:,:,i-2) &
       + 0.5 * (3. * nuos(:,:,i) - 4. * nuos(:,:,i-1) + nuos(:,:,i-2))) / dx

   end if
    
   call null_d1 (eth_J,  J, 2,  1)
   call null_d1 (eth_U,  U, 1,  1)
   call null_d1 (ethb_U, U, 1, -1)
   call null_d1 (eth_dx_J,  dx_J, 2,  1)
   call null_d1 (eth_dx_U,  dx_U, 1,  1)
   call null_d1 (ethb_dx_U, dx_U, 1, -1)
   call null_d1 (eth_cB,  cB, 1,  1)
   call null_d1 (ethb_cB, cB, 1, -1)
   call null_d1 (mu, J, 2, 1)
   call null_d1 (dx_mu, dx_J, 2, 1)

   dx_K = dble( dx_J * conjg(J) ) / K 
   fac = ( - dx_K * J + K * dx_J )  / K

   ! non-linear terms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   J_H = (1. - xc) / (rwt * xc) * exp(2. * beta) * ( -K * eth_J * conjg(cB) &
       + (K * nu + (K**2 - 1.) * eth_J - 2. * K * ck) * cB &
       + J * ((2. * ck - nu) * conjg(cB) - 2. * K * (ethb_cB + cB * conjg(cB)) &
              + 2. * dble((nu - ck) * conjg(cB) + conjg(J) * (eth_cB + cB * cB)))) &
       + 0.5 * rwt * (1. - xc) * xc**3 * exp(-2. * beta) &
       * ( (K * dx_U + J * conjg(dx_U))**2 &
          - J * dble(conjg(dx_U) * (K * dx_U + J * conjg(dx_U)))) &
       - 0.5 * (  nu * (xc * (1. - xc) * dx_U + 2. * U) &
                + eth_J * conjg(xc * (1. - xc) * dx_U + 2. * U) ) &
       + J * ii * dimag(xc * (1. - xc) * ethb_dx_U + 2. * ethb_U) &
       - xc * (1. - xc) * dx_J * dble(ethb_U) &
       + xc * (1. - xc) * (conjg(U) * eth_J + U * nu) * ii * dimag(J * conjg(dx_J)) &
       - xc * (1. - xc) * (conjg(U) * eth_dx_J + U * dx_nu) &
       - 2. * xc * (1. - xc) * (J * dx_K - K * dx_J) &
            * ( dble(conjg(U) * ck) + ii * dimag(K * ethb_U - conjg(J) * eth_U) )
       
   ! the rhs ...

   h_c = - K * (xc * (1. - xc) * eth_dx_U + 2. * eth_U ) &
         + 2. * (1. - xc) / (rwt * xc) * exp (2. * beta) * (eth_cB + cB * cB) &
         - (xc * (1. - xc) * dx_W + W) * J + J_H



   if(switch) then
   rh = rwt*xh(i-1)*(1. - xh(i-1))
   pr=0.0d0 !The pressure - adjust if necessary
   h_c = h_c + 4*PI*exp(2*beta)*(density(:,:,i-1)+pr)/2* &
               ((J*conjg(p_van)-K*p_van)**2 +p_van**2)/rh
   h_c = h_c + 4*PI*exp(2*beta)*(density_o(:,:,i)+pr)/2* &
               ((J*conjg(p_vano)-K*p_vano)**2 +p_vano**2)/rh
   endif


   ! interpolate the $f = x j$ values
   !crank initial

   f1 = dt * (xc * (1. - xc) * dx_W + W)
   f2 = dt * (1. - xc)**2 * ((1. - xc)/rwt + xc * W)

   jns(:,:,i) = (1-mask)* jns(:,:,i) + jos(:,:,i) * mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (i .ne. nx) then

      jave = 0.5 * (x(i) * jos(:,:,i) + x(i-1) * jos(:,:,i-1)) &
           + (1. - xc) * ((1. - 0.75 * disip_rescaled) &
                 * (x(i) * jos(:,:,i) - x(i-1) * jos(:,:,i-1)) &
                          + 0.25 * disip_rescaled &
                 * (x(i+1) * jos(:,:,i+1) - x(i-2) * jos(:,:,i-2))) / dx 

      d2x_jave = x(i+1) * jos(:,:,i+1) - 2. * x(i) * jos(:,:,i) &
               + x(i-1) * jos(:,:,i-1)

   else

      jave  = 0.5 * (x(i) * jos(:,:,i) + x(i-1) * jos(:,:,i-1) ) &
      + (1. - xc) * (x(i) * jos(:,:,i) - x(i-1) * jos(:,:,i-1) ) / dx

      d2x_jave = ( 2. * x(i)   * jos(:,:,i)   - 5. * x(i-1) * jos(:,:,i-1) &
                 + 4. * x(i-2) * jos(:,:,i-2) -      x(i-3) * jos(:,:,i-3) )

   end if

   do cran = 1, cran_max

    if(cran.eq.1) then
      P_u = J * (1. - xc) * (xc * 2. * dble(( jns(:,:,i-1) &
                                             - jos(:,:,i-1)) &
                                            / (dt) * conjg(fac)) &
                             - 8. * dx_beta * ((1. - xc)/rwt + xc * W))
    else
      P_u = J * (1. - xc) * (xc * 2. * dble((jns(:,:,i) + jns(:,:,i-1) &
                                           - jos(:,:,i) - jos(:,:,i-1)) &
                                            / (2. * dt) * conjg(fac)) &
                             - 8. * dx_beta * ((1. - xc)/rwt + xc * W))
    end if

      jns(:,:,i) = jns(:,:,i) * (1 - mask) &
      + mask * (- x(i-1) * jns(:,:,i-1) * ( (1. - f1 * 0.25) &
                * (0.5 - (1. - xc) / dx) + 0.5 * f2 / dx**2) &
                + x(i-2) * jns(:,:,i-2) * f2 * 0.25 / dx**2 &
                + jave * (1. + f1 * 0.25) + f2/dx**2 * d2x_jave*0.25 &
                + 0.5 * (h_c + P_u) * dt) &
                / (x(i) * ((1. - f1 * 0.25) &
                           * (0.5 + (1. - xc) / dx) - 0.25 * f2 / dx**2))

   end do

end subroutine null_j

end module null_evol
