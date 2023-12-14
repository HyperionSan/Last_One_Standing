module null_code

contains

subroutine null_read_params

   use null_grid, only : nn, nx, nt, rwt, time, cfl
   use null_params, only : disip, mass, it_skip, null_output

   implicit none

   namelist /null_input/ nn, nx, nt, rwt, time, cfl, &
         mass, it_skip, null_output, disip

   open (unit = 10, file = "null.in", status = "old" )
   read (unit = 10, nml = null_input)
   close (unit = 10)

end subroutine null_read_params

subroutine null_allocate	

   use null_grid
   use null_vars
   use null_eth
   use null_hyper_beta
   use null_hyper_u
   use null_hyper_w
   use null_evol
!  use calc_newsnorth
!  use calc_newssouth 
!  use null_coortranvars
!  use null_newsvarnorth
!  use null_newsvarsouth

   implicit none

   call null_grid_allocate
   call null_vars_allocate
   call null_eth_allocate
   call null_hyper_beta_allocate
   call null_hyper_u_allocate
   call null_hyper_w_allocate
   call null_evol_allocate
!  call null_coordtran_allocate
!  call null_newsnorth_allocate
!  call null_newssouth_allocate

end subroutine null_allocate

subroutine null_evolve

   use null_grid
   use null_vars
   use null_hyper_beta
   use null_hyper_u
   use null_hyper_w
   use null_evol
   use null_interp

   implicit none

   integer i

   do i = min(minval(masks),minval(maskn))+1, nx	 

      call null_j (i, jnn, bnn, unn, wnn, jon, bon, uon, won, maskn, e2bmn)
      call null_j (i, jns, bns, uns, wns, jos, bos, uos, wos, masks, e2bms)

      call null_cnsint (jns(:,:,i), jnn(:,:,i), 2)
      call null_cnsint (jnn(:,:,i), jns(:,:,i), 2)

      call null_beta (i, jns, bns, masks)
      call null_beta (i, jnn, bnn, maskn)

      call null_u (i, jnn, bnn, unn, maskn)
      call null_u (i, jns, bns, uns, masks)

      call null_cnsint (uns(:,:,i), unn(:,:,i), 1)
      call null_cnsint (unn(:,:,i), uns(:,:,i), 1)

      call null_w (i, jns, bns, uns, wns, masks)
      call null_w (i, jnn, bnn, unn, wnn, maskn)

      call null_rnsint (wns(:,:,i), wnn(:,:,i))
      call null_rnsint (wnn(:,:,i), wns(:,:,i))

   end do

end subroutine null_evolve
  
subroutine null_copy

   use null_vars

   implicit none

   jon = jnn
   jos = jns
   bon = bnn
   bos = bns
   uon = unn
   uos = uns
   won = wnn
   wos = wns

end subroutine null_copy

end module null_code
