module null_hyper_beta

   integer, dimension (:,:), allocatable, save, private :: mask
   double precision, dimension (:,:), allocatable, save, private :: &
      K, dx_beta, dx_K
   double complex, dimension (:,:), allocatable, save, private :: &
      J, dx_J, ethb_dx_J, eth_dx_beta, eth_dx_K

contains

subroutine null_hyper_beta_allocate

   use null_grid

   allocate (mask(nn,nn), K(nn,nn), dx_beta(nn,nn), dx_K(nn,nn), J(nn,nn), &
             dx_J(nn,nn), ethb_dx_J(nn,nn), eth_dx_beta(nn,nn), eth_dx_K(nn,nn))

end subroutine null_hyper_beta_allocate

subroutine null_beta (i, jns, bns, masks,switch)

   use null_grid
   use null_mask
   use particle

   implicit none

   integer,                             intent (in)    :: i
   double complex,   dimension (:,:,:), intent (in)    :: jns
   double precision, dimension (:,:,:), intent (inout) :: bns 
   integer,          dimension (:,:),   intent (in)    :: masks
   logical,                             intent (in)    :: switch

   double precision xhere, pr

   call null_remask (i, masks+1, mask)

   xhere = xh(i-1)

   J = 0.5 * (jns(:,:,i) + jns(:,:,i-1))
   K = sqrt(1. + J * conjg(J))

   dx_J = (jns(:,:,i) - jns(:,:,i-1)) / dx
   dx_K = dble( dx_J * conjg(J) ) / K

   bns(:,:,i) = bns(:,:,i) * (1 - mask) + mask * ( bns(:,:,i-1) &
   + dx * xhere * (1. - xhere) / 8. * (dx_J * conjg(dx_J) - dx_K * dx_K) )
   

   if(switch) then
   pr=0.0d0 !The pressure - adjust if necessary
   bns(:,:,i)=bns(:,:,i) &
       +2.*dx*Pi*xhere/(1.-xhere)**3*(density(:,:,i)+pr &
                                     +density(:,:,i-1)+pr)/2 *p_vdo(1)**2 
   endif
   
end subroutine null_beta

subroutine null_cb (i, bns, cbns, masks)

   use null_grid
   use null_mask
   use null_eth

   implicit none

   integer,                             intent (in)    :: i
   double precision, dimension (:,:,:), intent (in)    :: bns
   double complex,   dimension (:,:,:), intent (inout) :: cbns 
   integer,          dimension (:,:),   intent (in)    :: masks

   call null_remask (i, masks+1, mask)

   dx_beta = (bns(:,:,i) - bns(:,:,i-1)) / dx

   call null_d1 (eth_dx_beta, dcmplx(dx_beta), 0, 1)

   cbns(:,:,i) = cbns(:,:,i) * (1 - mask) &
   + mask * (cbns(:,:,i-1) + dx * eth_dx_beta)

end subroutine null_cb

subroutine null_ck (i, jns, ckns, masks)

   use null_grid
   use null_mask
   use null_eth

   implicit none

   integer,                           intent (in)    :: i
   double complex, dimension (:,:,:), intent (in)    :: jns
   double complex, dimension (:,:,:), intent (inout) :: ckns 
   integer,        dimension (:,:),   intent (in)    :: masks

   call null_remask (i, masks+1, mask)

   J = 0.5 * (jns(:,:,i) + jns(:,:,i-1))
   K = sqrt(1. + J * conjg(J))

   dx_J = (jns(:,:,i) - jns(:,:,i-1)) / dx
   dx_K = dble( dx_J * conjg(J) ) / K

   call null_d1 (eth_dx_K, dcmplx(dx_K), 0, 1)

   ckns(:,:,i) = ckns(:,:,i) * (1 - mask) &
   + mask * (ckns(:,:,i-1) + dx * eth_dx_K)

end subroutine null_ck

subroutine null_nu (i, jns, nuns, masks)

   use null_grid
   use null_mask
   use null_eth

   implicit none

   integer,                           intent (in)    :: i
   double complex, dimension (:,:,:), intent (in)    :: jns
   double complex, dimension (:,:,:), intent (inout) :: nuns 
   integer,        dimension (:,:),   intent (in)    :: masks

   call null_remask (i, masks+1, mask)

   dx_J = (jns(:,:,i) - jns(:,:,i-1)) / dx

   call null_d1 (ethb_dx_J, dx_J, 2, -1)

   nuns(:,:,i) = nuns(:,:,i) * (1 - mask) &
   + mask * (nuns(:,:,i-1) + dx * ethb_dx_J)

end subroutine null_nu

end module null_hyper_beta
