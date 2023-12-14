module hdrive

contains

  subroutine hdriver(timefrom3dcode, dt, nn, it, nt)

    use horizon
    use model
    use hio
    use gridtranslator
    use spheroid_init
    use null_params, only : time_real_start, real_start_dt, it_start
    use checkpoint,  only: recover_horizon, checkpointrecover
    use checkpoint_defs,  only: checkTimeLvls
    use axihorizon, only: step_ext_curv_axi, interpolate_axi2stereo, &
         axishifter, axishifter_up, axio
    implicit none
    double precision, intent (in) :: timefrom3dcode, dt
    integer,          intent (in) :: nn, nt, it

    double precision              :: time

    integer                       :: debug_lvl = 0
    integer                       :: jt ! the initial data time level
    integer, save                 :: times_called = 0
    integer, save                 :: it_initial = - 99999   
    integer                       :: recover_nt, recover_it
    double precision              :: recover_time
    time = timefrom3dcode   ! this is the time of the new level,
    ! which we have to fill with data
    ! at first call, this is the initial time


    ! the logic of the time loop is as follows:
    ! our intial slice is @ time level it = 1, which is the first time
    ! this subroutine is called. 
    ! then we first need to do some "extra" initialization calls with it set to 0
    !
    ! then we continue to calculate the new fields at time level it + 1


!!!! WARNING WARNING WARNING  !!!!
    ! I (yosef) changed the routine so that the first 'new' level is filled when it = it_start
    ! this may have broken (has broken ?) the axisymetric stuff 
    ! dont use it with this code -- also the checkpoint/recovery system does not work with 
    ! the axi stuff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (debug_lvl > 0) then
       print *, 'in hdrive(hdriver): it = ', it, 'mass = ', mass
    endif

    times_called = times_called + 1
    if ( times_called == 1) then
       write(*,*) "hdriver called initially with it = ", it
       it_initial = it
    endif
    if (checkpointrecover) then
       call model_setup (nn, time_real_start)   ! model_setup sets t_ = initial_time so we need to give it the real initial time
       call horizon_setup (nn)
       jt = 0

       it_initial = it_initial - 1  ! don't 're-initialize', and use std routine to fill in new vals when we are done
       checkpointrecover=.false.
       if ( model_switch == 1 ) then                  ! For conformal model we need to initialize the integrator 
          call solve_focuseq(time_real_start, jt, nt, real_start_dt, debug_lvl) 

          ! the point of this little excersize is to
          ! get the adaptive rungge kutta integrator
          ! to give the same values it did on
          ! the original run. To this end we recall it
          ! with the same time levels as the original
          do recover_it = it_start+1, it -1
             recover_time = checkTimeLvls(recover_it)
             call solve_focuseq(recover_time, recover_it, nt, dt, debug_lvl)
          end do
       end if
       call recover_horizon   ! this will read in all the 'old' values
    endif


    ! the initialization calls
    if (it == it_initial) then    ! wont be called if checkpointrecover was set
       print*
       print*, '<<< initializing hdriver >>>'
       print*

       jt = 0

       call model_setup (nn, time)
       call horizon_setup (nn)

       select case (model_switch)

       case (-1)                 !Schwarzschild case
          print*, 'Horizon will be Schwarzschild (fixed)'
          rho = 1.0;  rhodot = 0.0; j = 0.0
          rhonew = 1.0;  rhodotnew = 0.0; jnew = 0.0

       case (0)                  ! use perturbation equations
          !   call model_initial (nn, time, dt, j, rho, rhodot, rhol, omega, jl)
          call model_initial (nn, time, dt, jnew, rhonew, rhodotnew, rholnew, omeganew, jlnew)
          j      = jnew
          rho    = rhonew                 ! shearntwist needs these to be filled .. but they are never used
          rhodot = rhodotnew              ! in any evolution
          rhol   = rholnew
          omega  = omeganew
          jl     = jlnew

       case (1)                  ! use spheroidal model
          call solve_focuseq(time_real_start, jt, nt, dt, debug_lvl)
          !   call spheroid_initial (nn, jt, time, j, rho, rhodot, rhol, omega, jl)
          call spheroid_initial (nn, jt, time, jnew, rhonew, rhodotnew, rholnew, omeganew, jlnew)
          j      = jnew
          rho    = rhonew                 ! shearntwist needs these to be filled .. but they are never used
          rhodot = rhodotnew              ! in any evolution
          rhol   = rholnew
          omega  = omeganew
          jl     = jlnew

       case(2)                   ! use spheroidal model & manifest symmetry 
          stop 'this is not known to work with the current code, see warning'
          call solve_focuseq(time, jt, nt, dt, debug_lvl)
          call spheroid_initial_axi(nn, time, mass)
          call axishifter_up ! to have next routine work correctly
          if (null_evolution_switch == 1) then
             call interpolate_axi2stereo(jt, time, mass) ! the *new vals
          endif

       case DEFAULT
          write(*,*) 'this case of model_switch is not supported'
          STOP    'STOPPING EXECUTION'

       end select

       if (model_switch > 0) then
          call tcross(10000) ! 'estimate crossing time @ equator
       end if
       print*, '<<< DONE initializing hdriver >>>'
       print*

    else


       select case (model_switch)
       case(-1)                  !Schwarzschild - do nothing
          rhonew = 1.0;  rhodotnew = 0.0; jnew = 0.0
       case (0)                  ! use perturbation equations
          call model_sphmetric (nn, time, jnew)
          print *, "calling dot_rho"
          call dot_rho (nn, dt)
       case (1)  ! obtain r and j from spheroidal model
          call solve_focuseq(time, it, nt, dt, debug_lvl) 
          call spheroid_rj(time, nn, rhonew, jnew)
          call spheroid_rdot(nn, rhodotnew)
       case (2)
          call solve_focuseq(time, it, nt, dt, debug_lvl)
          call step_ext_curv_axi (it, dt, time)
          call axio (it, time, nt, mass) 
          if (null_evolution_switch == 1) then
             call interpolate_axi2stereo (it, time, mass) 
          endif
          call axishifter
       end select

       if (model_switch /= 2.and.model_switch /=-1) then
          call dot_omega (nn, dt)
          call dot_rho_lambda_RK2 (nn, dt, time, mass)
          call dot_j_lambda (nn, dt, time, mass)
       endif
    end if

    if (debug_lvl > 0) then
        write(*,*) 'min & max rhonew in Hdriver =', minval(rhonew),maxval(rhonew)
    endif

    if (null_evolution_switch == 1) then
       if (model_switch == -1) then
          call schwarzschild  (nn, time, mass, debug_lvl)
       else
          call htor (nn, time, mass, debug_lvl)
       endif
    endif

    !  call io (nn, it, time)				! Oh, be quiet !!!
    if (model_switch == 1 .AND. geomview_output) then
       call gio (nn, it, time)
    end if
    call shifter

  end subroutine hdriver

end module hdrive
