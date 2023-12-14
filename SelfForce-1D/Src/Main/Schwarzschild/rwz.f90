program rwz

  use parameters
  use time_info
  use DG_structures
  use grid
  use rwz_pert_schw
  use all_integrators
  use observers
  use geodesic_schwarzschild
  use osculating_schwarzschild
  use world_tube
  use rwz_pert_schw_eff
!  use self_force_observer
  use self_force_base
  use numerics, only : time_window
  use analytic_circular_orbit
!  use singular_observer
  use flux_observer
  use metric_observer
  use strain_observer
  use acceleration_history

  implicit none

  class(equation), pointer  :: eqp
  type(rwz_schw), allocatable, target :: my_eq
  integer(ip) :: i, j, k
  type(ref_element) :: relem
  class(equation_pointer), dimension(:), allocatable :: eqs
  type(cobserver), dimension(:,:), allocatable :: cobs
!  type(sing_observer) :: sobs
  type(flx_observer) :: flxobs
  type(met_observer) :: metobs
  type(str_observer) :: strobs
  class(equation), pointer :: orbit
  complex(wp), dimension(:,:), allocatable :: psi, dpsidt, dpsidr
  real(wp) :: r, phi, ur, En, Lz, tfac, dtfac_dt, dt2fac_dt2

  call read_parameters ()

  print*,'p_orb = ', p_orb
  print*,'ecc = ', ecc
  relem = ref_element( order )

  call init_grid_coordinates ( relem )

  call init_time ( t_initial )

! Construct the acceleration history type if necessary.
  if ( a_timelevels > 0 ) then
    if ( a_timelevels >= a_smooth_derivative_order+1 ) then
      ah = accel_history ( a_timelevels, a_smooth_derivative_order )
    else
      print*,'a_timelevels has to be larger than or equal to a_smooth_derivative_order+1'
      stop
    end if
  end if

  if (use_world_tube) then
    wt = wtube()
  end if

  if ( use_particle ) then
    if ( use_generic_orbit ) then
      if ( use_osculating_orbit ) then
        allocate ( osc_schw :: orbit )
      else
        allocate ( geod_schw :: orbit )
      end if
    else
      allocate ( circular_orbit :: orbit )
    end if 
    call orbit%init ( )
  end if

  print*,'rho%vname = ', rho%vname
  print*,'equation_name = ', equation_name
  print*,'n_elems = ', n_elems
  print*,'order = ', order
  allocate(my_eq)
  call my_eq%init ( )

  if ( output_coords_for_exact ) then
    call my_eq%output_coords ( )
  end if

!  if ( use_particle ) then
!    call sobs%init ( (/ p_orb /), my_eq%r_star, my_eq%effs )
!  end if

  call flxobs%init ( (/ rho%elems(1)%var(1), rho%elems(n_elems)%var(order+1) /), rho, my_eq )
  call metobs%init ( (/ rho_particle /), rho, my_eq )
  call strobs%init ( (/ rho%elems(n_elems)%var(order+1) /), rho, my_eq )

  if ( use_particle ) then
    allocate(eqs(2))
    eqs(2)%p => orbit
  else
    allocate(eqs(1))
  end if

  eqs(1)%p => my_eq
  print*,'my_eq%eq_data(1,1)%vname = ', my_eq%eq_data(1,1)%vname
  print*,'my_eq%eq_data(2,1)%vname = ', my_eq%eq_data(2,1)%vname
  print*,'my_eq%eq_data(3,1)%vname = ', my_eq%eq_data(3,1)%vname

  if (use_particle) then
    call time_window ( get_current_time (), tsigma, torder, tfac, &
                                     dtfac_dt, dt2fac_dt2 )
! Tell the effective source about the time 
    if ( turn_on_source_smoothly ) then
      call my_eq%effs%set_time_window ( tfac, dtfac_dt, dt2fac_dt2, -1 )
    else
      call my_eq%effs%set_time_window ( tfac, dtfac_dt, dt2fac_dt2, &
                                      exact_initial_data_lmax )
    end if
    print*,'Setting particle position'
    call orbit_info%get_orbit ( r, phi, ur, En, Lz )
    call my_eq%effs%set_particle_pos ( r, phi, ur, En, Lz,  &
                               0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp )
    print*,'Particle position set'
    if ( wt%wsize > 0 ) then
      print*,'Evaluating effective source'
      call my_eq%effs%evaluate_source ( my_eq%r_schw, wt )
      print*,'Effective source evaluated'
    end if
  end if

  if (use_exact_initial_data) then
    call my_eq%read_all_modes ( )
  end if

!  if (out1d_every>0) then
!    if ( use_particle ) then
!      do i = 1, my_eq%effs%nmodes
!        call my_eq%effs%source(1,i)%output(rho)
!      end do
!    end if
!  end if
!  call my_eq%r_schw%output(rho)
!  do i = 1, 4
!    call my_eq%eq_coeffs(i)%output(rho)
!  end do
!  do i = 1, my_eq%nmodes
!    call my_eq%eq_lcoeffs(i)%output(rho)
!  end do
!  do i = 1, 2
!    call my_eq%eq_lambda(i)%output(rho)
!    do j = 1, 2
!      call my_eq%eq_s(i,j)%output(rho)
!      call my_eq%eq_sinv(i,j)%output(rho)
!    end do
!  end do 
!  call my_eq%rhs ( )

  if (out1d_every>0) then
    do i = 1, my_eq%nmodes
      do j = 1, 3
        call my_eq%eq_data(j,i)%output(rho)
!        call my_eq%eq_rhs_data(j,i)%output(rho)
      end do
!      do j = 1, 2
!        call my_eq%eq_flux_data(j,i)%output(rho)
!      end do
    end do
  end if
!  if ( use_particle ) then
!    print*,'Calling orbit output'
!    call orbit%output()
!  end if

!  if ( use_field_observer ) then
!    allocate ( cobs(3,my_eq%nmodes) )
!    do i = 1, my_eq%nmodes
!      do j = 1, 3
!        call cobs(j,i)%init( (/ rho_min, rho_max /), rho, my_eq%eq_data(j,i) )
!        call cobs(j,i)%extract ( )
!        call cobs(j,i)%output ( )
!      end do
!    end do
!  end if
  call my_eq%init_coords_output ( )

  if (use_particle) then
!    call sfobs%init( (/ rho_particle /), rho, my_eq )
! During this call the self_force is extracted.
    call my_eq%save_globals_2 ( )
!    call sfobs%output ( )
!    call sf%output ( )
  end if

  call flxobs%extract ( )
!  call metobs%extract ( )
!  call strobs%extract ( )

  call flxobs%output ( )
!  call metobs%output ( )
!  call strobs%output ( )

  call choose_integrator ( )

  if (use_particle .and. use_generic_orbit) then
    call set_dtime ( coufac*delta_rho_min/my_eq%time_dep_coord%maxspeed )
  else
    call set_dtime ( coufac*delta_rho_min )
  end if
  print*,'dtime = ', get_current_dtime ( )

  call mol_int%init(eqs)
  print*,'Integrator is ', mol_int%iname
  print*,'ntemp for integrator ', mol_integrator, ' = ', mol_int%ntemp()

  if ( use_particle ) then
    print*,'turn_on_source_smoothly = ', turn_on_source_smoothly
    print*,'exact_initial_data_lmax = ', exact_initial_data_lmax, lmax
    if ( turn_on_source_smoothly .or. &
         (use_exact_initial_data .and. exact_initial_data_lmax < lmax ) ) then
      if (get_current_dtime()*100.0_wp>4.0_wp*tsigma) then
        print*,'In order to turn on the source smoothly, using dtime = ', 4.0_wp*tsigma/100
        call set_dtime ( 4.0_wp*tsigma/100 )
        short_timesteps_active = .true.
      end if
    end if
  end if

  print*,'Starting evolution at time = ', get_current_time ( )
  k = 0
  do 
    if ( get_current_time () > t_final ) exit
    k = k+1
    call mol_int%step( )

    if (use_filter) then
      call my_eq%apply_filter ( )
    end if

    if (out1d_every>0 .and. mod(k,out1d_every)==0) then
      do i = 1, my_eq%nmodes
        do j = 1, 3
          call my_eq%eq_data(j,i)%output(rho)
!          call my_eq%eq_rhs_data(j,i)%output(rho)
        end do
!        do j = 1, 2
!          call my_eq%eq_flux_data(j,i)%output(rho)
!        end do
!        if (use_particle) then
!          call my_eq%effs%source(1,i)%output(rho)
!        end if
!      end do
!      do j = 1, 4
!        call my_eq%eq_coeffs(j)%output(rho)
      end do 
!      call my_eq%r_star%output(rho)
!      call my_eq%r_schw%output(rho)
    end if
    if (out0d_every>0 .and. mod(k,out0d_every)==0) then
!      if (use_field_observer) then
!        do i = 1, my_eq%nmodes
!          do j = 1, 3
!            call cobs(j,i)%extract ( )
!            call cobs(j,i)%output ( )
!          end do
!        end do
!      end if
      if (use_particle) then
! During this call the self_force is extracted.
        call my_eq%save_globals_2 ( )
!        call sfobs%output ( )
!        call orbit%output ( )
!        call sf%output ( )
      end if
      call flxobs%extract ( )
!      call metobs%extract ( )
!      call strobs%extract ( )
      call flxobs%output ( )
!      call metobs%output ( )
!      call strobs%output ( )
    end if
  end do

  print*,'Program ended after ', k, ' iterations.'
end program rwz
