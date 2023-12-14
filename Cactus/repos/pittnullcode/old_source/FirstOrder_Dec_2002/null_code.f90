module null_code

contains

subroutine null_read_params

   use null_params
   implicit none

   namelist /null_input/ nn, nx, nt, iot0, iot1, iot2, iot3, rwt, cfl, disip, &
                         time, ID_CASE, R_left, R_right, ID_AMP,      &
                         horizon_output, scri_output, boundary_output, &
                         it_cfl, dt_fix, cfl_compute, circ_rad, NIntPts, &
                         qsize, pure_schwarzs, do_checkpoint, &
                         it_checkpoint, chkrcvr, calc_news, horizon_on_grid, &
                         cdumpj, cdumpit, docdump, stdang, stdrad

   open (unit = 10, file = "null.in", status = "old" )
   read (unit = 10, nml = null_input)
   close(unit = 10)

  write(*,*) do_checkpoint, it_checkpoint, chkrcvr

   if (pure_schwarzs .AND. ( NIntPts .ne. 0)) then 
     write(*,*) "!!!!!!!!! WARNING !!!!!!!!!!!!"
     write(*,*) "pure_schwarzs is set, but NIntPts .ne. 0"
     write(*,*) "NIntPts will be reset to 0"
     NIntPts = 0
  end if
end subroutine null_read_params

subroutine null_allocate 

   use null_grid
   use null_vars
   use null_hyper_beta
   use null_hyper_u
   use null_hyper_w
   use null_evol

   implicit none

   call null_grid_allocate
   call null_vars_allocate
   call null_hyper_beta_allocate
   call null_hyper_u_allocate
   call null_hyper_w_allocate
   call null_evol_allocate

end subroutine null_allocate

subroutine null_evolve

   use null_params, only: it_start, it_cfl
   use null_grid
   use null_vars
   use null_hyper_beta
   use null_hyper_u
   use null_hyper_w
   use null_evol
   !use null_interp_bela
   use null_interp
   use null_cfl_test
   use particle
   use boundary, only : xbdry

   implicit none

   integer i

   do i = minval(maskn)+1, nx  
      IF(it > it_start) THEN

        call null_j (i, jnn, nunn, cknn, bnn, cbnn, unn, wnn, &
                      jon, nuon, ckon, bon, cbon, uon, won, &
                   maskn,p_patchN)

        call null_cnsint (jnn(:,:,i), jnn(:,:,i), 2)
      END IF

      call null_nu (i, jnn, nunn, maskn)

      call null_cnsint (nunn(:,:,i), nunn(:,:,i), 1)

      call null_ck (i, jnn, cknn, maskn)

      call null_cnsint (cknn(:,:,i), cknn(:,:,i), 1)

      call null_beta (i, jnn, bnn, maskn,p_patchN)

      call null_rnsint (bnn(:,:,i), bnn(:,:,i))

      call null_cb (i, bnn, cbnn, maskn)

      call null_cnsint (cbnn(:,:,i), cbnn(:,:,i), 1)

      call null_u (i, jnn, nunn, cknn, bnn, cbnn, unn, maskn,p_patchN)

      call null_cnsint (unn(:,:,i), unn(:,:,i), 1)

      call null_w (i, jnn, nunn, cknn, bnn, cbnn, unn, wnn, maskn,p_patchN)

      call null_rnsint (wnn(:,:,i), wnn(:,:,i))
   end do

end subroutine null_evolve
  
subroutine null_copy

   use null_vars

   implicit none

   jon = jnn

   nuon = nunn

   ckon = cknn

   bon = bnn

   cbon = cbnn

   uon = unn

   won = wnn

end subroutine null_copy

end module null_code
