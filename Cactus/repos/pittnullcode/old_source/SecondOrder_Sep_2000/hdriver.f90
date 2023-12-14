module hdrive

contains

subroutine hdriver(timefrom3dcode, dt, dd, nn, it, nt)

   use horizon
   use model
   use hio
   use gridtranslator
   use spheroid_init
   use axihorizon, only: step_ext_curv_axi_beta, interpolate_axi2stereo, &
                         axishifter, axishifter_up, axio
   implicit none
   double precision, intent (in) :: timefrom3dcode, dt, dd
   integer,          intent (in) :: nn, nt, it

   double precision              :: time

   integer                       :: debug_lvl = 0
   integer                       :: jt ! the initial data time level
   time = timefrom3dcode
 


! the logic of the time loop is as follows:
! our intial slice is @ time level it = 1, which is the first time
! this subroutine is called. 
! then we first need to do some "extra" initialization calls with it set to 0
!
! then we continue to calculate the new fields at time level it + 1

! the initialization calls
   if (it == 1) then
      jt = 0
      call model_setup (nn)

      call horizon_setup (nn)

      select case (model_switch)
      case (0)                  ! use perturbation equations
         call model_initial (nn, time, dt, j, rho, rhodot, rhol, omega, jl)
      case (1)                  ! use spheroidal model
         call solve_focuseq(time, jt, nt, dt, debug_lvl)
         call spheroid_initial(nn, jt, time, j, rho, rhodot, rhol, omega, jl)
      case(2)                   ! use spheroidal model & manifest symmetry 
         call solve_focuseq(time, jt, nt, dt, debug_lvl)
         call spheroid_initial_axi(nn, time, mass)
         call axishifter_up ! to have next routine work correctly
         if (ext_curv_switch == 1) then
            call interpolate_axi2stereo(jt, dt, time, mass) ! the *new vals
         endif
      case DEFAULT
         write(*,*) 'this case of model_switch is not supported'
         STOP    'STOPPING EXECUTION OF FULL CODE'
      end select

      call tcross(10000) ! 'estimate crossing time @ equator
   endif

! now calculate the new field values at time level (it + 1)
   select case (model_switch)
   case (0)                  ! use perturbation equations
      call model_sphmetric (nn, time, jnew)
      call dot_rho (nn, dt)
   case (1)  ! obtain r and j from spheroidal model
      call solve_focuseq(time, it, nt, dt, debug_lvl) 
      call spheroid_rj(time, nn, rhonew, jnew)
   case (2)
      call solve_focuseq(time, it, nt, dt, debug_lvl)
      call step_ext_curv_axi_beta (it, dt, time, mass)
      call axio (it, time, nt, mass) 
      if (ext_curv_switch == 1) then
         call interpolate_axi2stereo (it, dt, time, mass) 
      endif
      call axishifter
   end select

   if ((model_switch < 2).AND.(ext_curv_switch == 1)) then     
      call dot_omega (nn, dt)

      if (rhol_integrator_switch == 0) then
         call dot_rho_lambda     (nn, dt, time, mass)
      else
         call dot_rho_lambda_RK2 (nn, dt, time, mass)
      endif

      call dot_j_lambda (nn, dt, time, mass)
   endif

   if (ext_curv_switch == 1) then 
      call htor (nn, time, dt, mass)
      
      call io (nn, it, time)
      !   call gio (nn, it, time)
   else
      call htor_dummy(nn, time, dt, mass)
   endif
   
   call shifter

end subroutine hdriver

end module
