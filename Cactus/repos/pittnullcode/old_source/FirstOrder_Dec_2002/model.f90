Module Model

   use null_grid, only : dd, q=>qs, p=>ps, pp
   Implicit None

   ! Declare public interfaces. All others should be private

   Public Model_setup, Model_Sphmetric, Model_Twist, &
          Model_J_r, Model_Initial

   ! Any and all parameters for the axysimmetric model.
   ! Which ones need to be public? I think only the mass (for now...)

   Double Precision, Public :: Mass    =    1.0,     &
                               R_inf           ,     &
                               R_a     =    1.0,     &
                               eps     =    1.0D-04, &
                               t_0     =  -10.0,     &
                               t_      = -100.0,     &
                               p_0     =    1.0,     &
                               a_inf   =    1.0,     &
                               p_model,              &
                               red_t_min = -.99,     &
                               red_t_max = -.90        

   Integer, Public ::          io_conv_skip          =  1,  &
                               model_switch          =  1,  &
                               L_inf_switch          =  1,  &
                               null_evolution_switch =  1

   Logical, Public ::          p_switch        = .false.,   &
                               geomview_output = .false.,   &
                               pure_schwarzs   = .false.

   ! Local variables, whatever you need to compute J

   Double Precision, Dimension (:,:,:), Allocatable, Save, Public :: t_hat

   Private

Contains

Subroutine Model_Setup (nn, time)
!-----------------------------------------------------------------------
! Allocationf of variables, grid, etc.
!-----------------------------------------------------------------------

   Implicit None

   Integer,          Intent (in) :: nn
   Double Precision, Intent (in) :: time

   ! Read your parameters. I think we want to keep them separate
   ! from the grid, etc.

   Namelist /Model_Input/ Mass, R_a, t_0, eps, p_0, &
                          a_inf, model_switch, L_inf_switch, &
                          null_evolution_switch, io_conv_skip, &
                          red_t_min, red_t_max, p_switch, geomview_output

   Open (unit = 10, file = "model.in", status = "old" )
   Read (unit = 10, nml = Model_Input)
   Close(unit = 10)

   ! Set the initial time to that read from the "null.in" file

   t_ = time

   allocate (t_hat(nn,nn,2))

   ! Whatever else you may need
   R_inf   = 2.0 * Mass
  
   p_model = p_0 * R_a / sqrt(13.0d0)

End Subroutine Model_Setup

Subroutine Model_Sphmetric (nn, time, J)
!-----------------------------------------------------------------------
! J(u) on the horizon
!-----------------------------------------------------------------------

   use null_grid, only : z, zb
   Implicit None

   Integer, Intent (in) :: nn
   Double Precision, Intent (in) :: time
   Double Complex, Dimension (nn,nn,2), Intent (out) :: J

   Double Precision :: PU, red_time

! Take J ~ sin^2(theta), taking care of transforming the function properly.

   red_time = time / mass
   if ( red_time > red_t_min .and. red_time < red_t_max ) then
       PU = 4096.0d0 * ((red_time - red_t_min) * (red_t_max - red_time))**6 &
           / ( red_t_max - red_t_min )**12
   else
       PU = 0.0d0
   end if 
   J(:,:,1) =  eps * PU * 4.0*(z*z)/(1+z*zb)**2 

   J(:,:,2) =  eps * PU * 4.0*(z*z)/(1+z*zb)**2   

!   t_hat = time
!
!   J(:,:,1) = - eps * R_a / (t_hat(:,:,1) - t_0 - R_a) &
!   *4.0*(zb*zb)/(1+z*zb)**2 
!
!   J(:,:,2) = - eps * R_a / (t_hat(:,:,2) - t_0 - R_a) &
!   *4.0*(zb*zb)/(1+z*zb)**2   

End Subroutine Model_Sphmetric

Subroutine Model_Twist (nn, time, J, omega)
!-----------------------------------------------------------------------
! Set \omega = -\bar\eth J/2
!-----------------------------------------------------------------------

   use null_grid, only : z, zb
   Implicit None

   Integer,                             Intent (in) :: nn
   Double Precision,                    Intent (in) :: time
   Double Complex, Dimension (nn,nn,2), Intent (in) :: J
   Double Complex, Dimension (nn,nn,2), Intent (out) :: omega


   omega(:,:,1) = (0.0d0, 0.0d0)
   omega(:,:,2) = (0.0d0, 0.0d0)

!   omega(:,:,1) = - eps * R_a / (t_hat(:,:,1) - t_0 - R_a) &
!                   * 4.0*z*(z*zb - 1.0)/(z*zb + 1.0)**2
!
!   omega(:,:,2) = - eps * R_a / (t_hat(:,:,2) - t_0 - R_a) &
!                   * 4.0*z*(z*zb - 1.0)/(z*zb + 1.0)**2

End Subroutine Model_Twist

Subroutine Model_J_r (nn, time, Jr)
!-----------------------------------------------------------------------
! \J_{,r} on the horizon
!-----------------------------------------------------------------------

   use null_grid, only : z, zb
   Implicit None

   Integer, Intent (in) :: nn
   Double Precision, Intent (in) :: time
   Double Complex, Dimension (nn,nn,2), Intent (out) :: Jr

   Double Precision :: u, u_, lograt, r_M

   r_M = 2. * Mass
   u = time - t_0 - R_a
   u_ = t_ - t_0 - R_a
   lograt = log(u / u_)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   Jr(:,:,2) = 3 / r_M * (8 * z ** 3 * zb ** 3 + 12 * z ** 2 * zb ** 2 &
   - 8 * z * zb - 11) * (1 + z * zb) * z ** 2 * eps * R_a * lograt &
   / (( - 90 * z ** 2 * zb ** 2 - 30 * z * zb + 48 * z ** 5 * zb ** 5 &
   + 8 * z ** 6 * zb ** 6 + 4 + 54 * z ** 4 * zb ** 4 &
   - 42 * z ** 3 * zb ** 3) * lograt * R_a * eps - u &
   - 240 * u * z ** 4 * zb ** 4 - 160 * u * z ** 3 * zb ** 3 &
   - 64 * u * z ** 6 * zb ** 6 - 192 * u * z ** 5 * zb ** 5 &
   - 60 * u * z ** 2 * zb ** 2 - 12 * u * z * zb)

   Jr(:,:,1) = - 3 / r_M * (11 * z ** 3 * zb ** 3 + 8 * z ** 2 * zb ** 2 &
   - 12 * z * zb - 8) * (1 + z * zb) * z ** 2 * eps * R_a * lograt &
   / ((54 * z ** 2 * zb ** 2 + 48 * z * zb - 30 * z ** 5 * zb ** 5 &
   + 4 * z ** 6 * zb ** 6 + 8 - 90 * z ** 4 * zb ** 4 &
   - 42 * z ** 3 * zb ** 3) * lograt * R_a * eps - u * z ** 6 * zb ** 6 &
   - 240 * u * z ** 2 * zb ** 2 - 160 * u * z ** 3 * zb ** 3 - 64 * u &
   - 192 * u * z * zb - 60 * u * z ** 4 * zb ** 4 - 12 * u * z ** 5 * zb ** 5)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Subroutine Model_J_r

Subroutine Model_Initial (nn, time, dt, J, R, Rdot, Rl, omega, Jl)
!-----------------------------------------------------------------------
! This is an initial value problem, after all: we need values of
! J, R, dot(R), R_{,\lambda}, \omega and J_{,\lambda}
!-----------------------------------------------------------------------

   Implicit None

   Integer,          Intent (in) :: nn
   Double Precision, Intent (in) :: time, dt
   Double Precision, Dimension (nn,nn,2) :: R, Rdot, Rl
   Double Complex,   Dimension (nn,nn,2) :: J, Jl, omega

   ! Initial values of J

   Call Model_Sphmetric (nn, time, J)

   !-----------------------------------------------------------------------
   ! Initial values of \rho and \dot \rho
   !-----------------------------------------------------------------------

   R = 1.0
   Rdot = 0.0

   !-----------------------------------------------------------------------
   ! Initial values of \rho_{,\lambda} and J_{,\lambda}
   !-----------------------------------------------------------------------

   rl(:,:,:) = 0.0
   Jl(:,:,:) = cmplx(0.0, 0.0)

   !-----------------------------------------------------------------------
   ! Set \omega = -\bar\eth J/2
   !-----------------------------------------------------------------------

   Call Model_Twist (nn, time, J, omega)

End Subroutine Model_Initial

End Module Model
