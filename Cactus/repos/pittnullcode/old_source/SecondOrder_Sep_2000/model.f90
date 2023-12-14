Module Model

   Implicit None

   ! Declare public interfaces. All others should be private

   Public Model_setup, Model_Sphmetric, Model_Twist, Model_rho_lambda, &
          Model_J_lambda, Model_J_r, Model_Initial

   ! Any and all parameters for the axysimmetric model.
   ! Which ones need to be public? I think only the mass (for now...)

   Double Precision, Public :: Mass    =    1.0,     &
                               R_inf   =    2.0,     &
                               R_a     =    1.0,     &
                               eps     =    1.0D-04, &
                               t_0     =  -10.0,     &
                               t_      = -100.0,     &
                               p_0     =    1.0,     &
                               a_inf   =    1.0,     &
                               p_model

   Integer, Public :: iot0 = 1, iot2 = 1, model_switch, L_inf_switch,  &
                    & rhol_integrator_switch, ext_curv_switch, subsample = 1

   ! Local variables, whatever you need to compute J

   Double Precision, Dimension (:,:),   Allocatable, Save, Public  :: q, p, pp
   Double Complex,   Dimension (:,:),   Allocatable, Save, Private :: z, zb
   Double Precision, Dimension (:,:,:), Allocatable, Save, Public  :: u_0, t_hat

   Private

Contains

Subroutine Model_Setup (nn)
!-----------------------------------------------------------------------
! Allocationf of variables, grid, etc.
!-----------------------------------------------------------------------

   Implicit None

   Integer, Intent (in) :: nn

   Double Precision :: dd

   Double Precision, Parameter :: pi = 3.1415926535897932385
   Double Complex,   Parameter :: ii = (0.0, 1.0)
   Integer :: i, j

   ! Read your parameters. I think we want to keep them separate
   ! from the grid, etc.

   Namelist /Model_Input/ Mass, R_inf, R_a, eps, t_, t_0, p_0, &
                          iot0, iot2, model_switch, L_inf_switch, a_inf, &
                          rhol_integrator_switch, ext_curv_switch, subsample

   Open (unit = 10, file = "model.in", status = "old" )
   Read (unit = 10, nml = Model_Input)
   Close(unit = 10)

   ! Array, arrays... the Null Grid module deals with *one* patch.
   ! Life is simpler if you allocate and deal with the two patches
   ! at once, so declare your own variables. 
   ! Pass the size of the array to the allocation routine, and 
   ! the time and the output array to the routine that computes "J".

   Allocate (q(nn,nn), p(nn,nn), pp(nn,nn))
   Allocate (z(nn,nn), zb(nn,nn))
   Allocate (u_0(nn,nn,2), t_hat(nn,nn,2))

   ! Grid spacing - don't change this

   dd = 2.0 / Dble(nn - 5)

   ! real and imaginary part of the sterographic coordinates

   Do j = 1, nn
      Do i = 1, nn
         q(i,j) = -1.0 + (i-3) * dd
         p(i,j) = -1.0 + (j-3) * dd
      End Do
   End Do

   ! This pops up all the time

   pp = 1.0 + q ** 2 + p ** 2

   ! The stereographic coordinate and its complex conjugate...

   z  = q + ii * p
   zb = q - ii * p

   ! Whatever else you may need


   ! ... !

   !u_0 on the initial slice of the horizon
   ! check these expressions!!
   u_0(:,:,1) = - R_a * (1. + eps * (0.5 &
      - (z * zb / (2. + z * zb)) ** 2))

   u_0(:,:,2) = - R_a * (1. + eps * (0.5 &
      - 1. / (2. * z * zb + 1.) ** 2))

   p_model = p_0 * R_a / sqrt(13.0d0)

End Subroutine Model_Setup


Subroutine Model_Sphmetric (nn, time, J)
!-----------------------------------------------------------------------
! J(u) on the horizon
!-----------------------------------------------------------------------

   Implicit None

   Integer, Intent (in) :: nn
   Double Precision, Intent (in) :: time
   Double Complex, Dimension (nn,nn,2), Intent (out) :: J

! Take J ~ sin^2(theta), taking care of transforming the function properly.

   t_hat = time ! + 5. * p_model &
         ! * log((u_0 + time - t_0) / (u_0 + t_ - R_a))

   J(:,:,1) = - eps * R_a / (t_hat(:,:,1) - t_0 - R_a) &
   *4.0*(z*z)/(1+z*zb)**2 

   J(:,:,2) = - eps * R_a / (t_hat(:,:,2) - t_0 - R_a) &
   *4.0*(z*z)/(1+z*zb)**2   

!   J(:,:,1) = - eps * R_a / (t_hat(:,:,1) - t_0 - R_a) &
!   *(1+z*zb)*z**2*(2*z*zb+1)/(2+z*zb)**4 

!   J(:,:,2) = - eps * R_a / (t_hat(:,:,2) - t_0 - R_a) &
!   *(z*zb+1)*(2+z*zb)*z**2/(2*z*zb+1)**4

End Subroutine Model_Sphmetric

Subroutine Model_Twist (nn, time, J, omega)
!-----------------------------------------------------------------------
! Set \omega = -\bar\eth J/2
!-----------------------------------------------------------------------

   Implicit None

   Integer, Intent (in) :: nn
   Double Precision, Intent (in) :: time
   Double Complex, Dimension (nn,nn,2), Intent (in) :: J
   Double Complex, Dimension (nn,nn,2), Intent (out) :: omega

!   t_hat = time + 5. * p_model &
!         * log((u_0 + time - t_0) / (u_0 + t_ - R_a))

!   omega(:,:,1) = - eps * R_a / (t_hat(:,:,1) - t_0 - R_a) &
!   *0.5*(1+z*zb)*z*(-3*z**2*zb**2-12*z*zb-4+4*z**3*zb**3)/(2+z*zb)**5

!   omega(:,:,2) = - eps * R_a / (t_hat(:,:,2) - t_0 - R_a) &
!   *0.5*(1+z*zb)*z*(12*z**2*zb**2+3*z*zb-4+4*z**3*zb**3)/(2*z*zb+1)**5

   omega(:,:,1) = - eps * R_a / (t_hat(:,:,1) - t_0 - R_a) &
                   * 4.0*z*(z*zb - 1.0)/(z*zb + 1.0)**2

   omega(:,:,2) = - eps * R_a / (t_hat(:,:,2) - t_0 - R_a) &
                   * 4.0*z*(z*zb - 1.0)/(z*zb + 1.0)**2

End Subroutine Model_Twist

Subroutine Model_rho_lambda (nn, time, rhol)
!-----------------------------------------------------------------------
! \rho_{,\lambda} on the horizon
!-----------------------------------------------------------------------

   Implicit None

   Integer, Intent (in) :: nn
   Double Precision, Intent (in) :: time
   Double Precision, Dimension (nn,nn,2), Intent (out) :: rhol

   Double Complex, Dimension (nn,nn,2) :: tmp

   Double Precision :: r_M

   r_M = 2. * Mass

!   t_hat = time + 5. * p_model &
!         * log((u_0 + time - t_0) / (u_0 + t_ - R_a))

   rhol(:,:,1) = 0.0
   rhol(:,:,2) = 0.0
      
!   rhol(:,:,1) = eps * R_a / r_M ** 2 &
!   * log((time - t_0 - R_a) / (t_ - t_0 - R_a)) &
!   *(1+z*zb)**3*(2*z**3*zb**3-21*z**2*zb**2+12*z*zb+4)/(2+z*zb)**6

!   rhol(:,:,2) = eps * R_a / r_M ** 2 &
!   * log((time - t_0 - R_a) / (t_ - t_0 - R_a)) &
!   *(z*zb+1)**3*(4*z**3*zb**3+12*z**2*zb**2-21*z*zb+2)/(2*z*zb+1)**6
      
End Subroutine Model_rho_lambda

Subroutine Model_J_lambda (nn, time, Jl)
!-----------------------------------------------------------------------
! \J_{,\lambda} on the horizon
!-----------------------------------------------------------------------

   Implicit None

   Integer, Intent (in) :: nn
   Double Precision, Intent (in) :: time
   Double Complex, Dimension (nn,nn,2), Intent (out) :: Jl

   Double Precision :: r_M

   r_M = 2. * Mass

   Jl(:,:,1) = cmplx(0.0, 0.0)
   Jl(:,:,2) = cmplx(0.0, 0.0)

!   Jl(:,:,1) = -1.5d0 * eps * R_a / r_M ** 2 &
!   * log((time - t_0 - R_a) / (t_ - t_0 - R_a)) &
!   *z**2*(1+z*zb)*(11*z**3*zb**3+8*z**2*zb**2-12*z*zb-8)/(2+z*zb)**6
 
!   Jl(:,:,2) = 1.5d0 * eps * R_a / r_M ** 2 &
!   * log ((time - t_0 - R_a) / (t_ - t_0 - R_a)) &
!   *z**2*(z*zb+1)*(8*z**3*zb**3+12*z**2*zb**2-8*z*zb-11)/(2*z*zb+1)**6

End Subroutine Model_J_lambda

Subroutine Model_J_r (nn, time, Jr)
!-----------------------------------------------------------------------
! \J_{,r} on the horizon
!-----------------------------------------------------------------------

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

   ! Initial values of J (in two levels)

   Call Model_Sphmetric (nn, time, J)

   !-----------------------------------------------------------------------
   ! Initial values of \rho and \dot \rho
   !-----------------------------------------------------------------------

   R = 1.0
   Rdot = 0.0

   !-----------------------------------------------------------------------
   ! Initial values of \rho_{,\lambda}
   !-----------------------------------------------------------------------

   Call Model_rho_lambda (nn, time, Rl)

   !-----------------------------------------------------------------------
   ! Set \omega = -\bar\eth J/2
   !-----------------------------------------------------------------------

   Call Model_Twist (nn, time, J, omega)

   !-----------------------------------------------------------------------
   ! Initial values of J_{,\lambda})
   !-----------------------------------------------------------------------

   Call Model_J_lambda (nn, time, Jl)

End Subroutine Model_Initial

End Module Model
