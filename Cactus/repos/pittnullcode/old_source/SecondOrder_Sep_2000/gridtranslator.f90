module gridtranslator
!====================
! purpose: communicate nullcode grid structure to ODE solver and
!          cummunicate back results

use prec ! set precision-parameter wp
implicit none

integer,       allocatable, dimension(:), Private :: iperm, rank
real(kind=wp), allocatable, dimension(:), Public  :: theta, cos_theta
real(kind=wp), allocatable, dimension(:), Public  :: uhatofu, omvec
integer, Public       :: n_x       ! how many ODEs to solve 
integer, Public       :: n_u       ! ... at how many time steps 
real(kind=wp), Public :: t_inc     ! time increment 
integer               :: flipper = 1
contains

subroutine stereogrid_to_thetavect(n_x, debug)
!============================================
! purpose: generate a 1D array of theta-values from the p-q arrays
!          and compute a ranking array for it.
 
use prec
use Model
use affine, only: sint, cost
implicit none

integer, intent(out) :: n_x        ! the number of theta values
integer, intent(in)  :: debug      ! debug level
integer, parameter   :: kflag = 2  ! return permutation vector sorted in
                                   ! increasing order and also sort input
integer              :: ier = 0    ! error indicator
integer              :: i, i1, ns  ! counters
integer              :: nn
real(kind=wp)        :: halfpi, pi, h, x
logical              :: equidist_theta 
!! == executable statements == !!

nn   = size(q,1)
if (model_switch < 2) then
   n_x  = 2*nn*nn
else
   n_x =  2*nn + 1
endif

allocate(theta(n_x),iperm(n_x),rank(n_x),cos_theta(n_x),cost(n_x),sint(n_x))

if (model_switch < 2) then
   Do ns = 1, 2            ! first North - then South
      Do i = 1, nn
         Do i1 = 1, nn
            if(ns == 1) then ! North
               cos_theta(nn*(i-1)+i1)         = (2.0_wp-pp(i,i1))/pp(i,i1)
               theta    (nn*(i-1)+i1)         = acos(cos_theta(nn*(i-1)+i1))
            else             ! South
               cos_theta(nn*(i-1)+i1 + nn*nn) = (pp(i,i1)-2.0_wp)/pp(i,i1)
               theta(nn*(i-1)+i1 + nn*nn) =acos(cos_theta(nn*(i-1)+i1 + nn*nn))
            endif
         End Do
      End Do
   End Do

   call dpsort(theta, n_x, iperm, kflag, ier) ! generate ranking array

   Do i = 1, n_x   
      Do i1 = 1, n_x
         if ( iperm(i1) == i ) rank(i) = i1
      End Do
   End Do

else
   equidist_theta = .false. 
   if(equidist_theta) then
      ! initialize spatial grid points for error checking run
      halfpi = asin(1.0_wp)
      pi     = 2.0_wp*halfpi
      h      = pi/dble(n_x-1)
      do i = 1, n_x
         theta(i) = h * dble(i-1)
         cos_theta(i) = cos(theta(i))
      enddo
   else
      ! initialize spatial grid points for full axisymmetric run
      h = 1.0/dble(nn)
      cos_theta(nn+1) = 0.0_wp
      do i = 1, nn
         cos_theta(i)      = -h * dble(nn - i + 1)
         cos_theta(nn+1+i) =  h * dble(i)
         theta(i) = acos(cos_theta(i))
      enddo
      cos_theta(1)   = -1.0_wp
      cos_theta(n_x) =  1.0_wp	 
   endif
endif

cost = cos_theta
sint = sqrt(1.0_wp - cos_theta**2)

open(86, file = 'cos_theta_input')
do i = 1, n_x
   write(86,*) cos_theta(i)
enddo
close(86)

return
999  format(3I3,2E13.5)
end subroutine stereogrid_to_thetavect

subroutine Spheroid_RJ(time, nn, Rout, Jout)
!===========================================
! purpose: calculate the area radius R an J for timestep it

use affine
use Model

implicit none

integer,          intent(in)                       :: nn
real(kind=wp),    intent(in)                       :: time
real(kind=wp),    intent(out), dimension (nn,nn,2) :: Rout
complex(kind=wp), intent(out), dimension (nn,nn,2) :: Jout
complex(kind=wp)                                   :: phase
real(kind=wp)                                      :: thetaJ_N, thetaJ_S
real(kind=wp)                                      :: rhatN, rhatS, d_th_max
real(kind=wp)                                      :: uh_N, uh_S, &
                                                    & sig_half_N, sig_half_S,&
                                                    & one_over_2m
integer                                            :: indexN, indexS
integer                                            :: i, i1   ! counters
integer                                            :: middle_nn
logical                                            :: odd_nn
!! === executable statements === !!

if (modulo(nn,2) == 1) then
   odd_nn = .true.
   middle_nn = (nn + 1)/2 
else
   middle_nn = -42 ! never reached!
endif

!write(88,99) time, uhatofu(i_equator)

one_over_2m = 1.0_wp / (2.0_wp * mass)

Do i = 1, nn
   Do i1 = 1, nn

      indexN = rank( nn*(i-1) + i1 )
      indexS = rank( nn*(i-1) + i1 + nn*nn) 
  
      uh_N =   uhatofu(indexN)
      uh_S =   uhatofu(indexS)
      sig_half_N =  0.5_wp*sigma(indexN)
      sig_half_S =  0.5_wp*sigma(indexS)

      t_hat(i1,i,1) = uh_N - u_zero(indexN) + t_0
      t_hat(i1,i,2) = uh_S - u_zero(indexS) + t_0

      rhatN = sqrt(abs((uh_N + sig_half_N)*  &
            &          (uh_N - sig_half_N)))
            
      rhatS = sqrt(abs((uh_S + sig_half_S)*  &
           &           (uh_S - sig_half_S)))
      
      Rout(i,i1,1) = rhatN * omvec(indexN) * one_over_2m
      Rout(i,i1,2) = rhatS * omvec(indexS) * one_over_2m

      thetaJ_N =  0.5_wp*(                                            &
           &          (uh_N - sig_half_N)/   &
           &          (uh_N + sig_half_N)    &
           &      -                                                   &
           &          (uh_N + sig_half_N)/   &
           &          (uh_N - sig_half_N)    &
           &            )

      thetaJ_S =  0.5_wp*(                                            &
           &          (uh_S - sig_half_S)/   &
           &          (uh_S + sig_half_S)    &
           &      -                                                   &
           &          (uh_S + sig_half_S)/   &
           &          (uh_S - sig_half_S)    &
           &            )


      if ((i == middle_nn).and.(i1 == middle_nn)) then
         phase = 1.0_wp
      else 
         phase = cmplx(q(i,i1), p(i,i1))/cmplx(q(i,i1), -p(i,i1))
         ! note: the N/S phase difference comes from the different semantics of
         !       p/q w.r.t. theta 
      endif

      Jout(i,i1,1) = cmplx(thetaJ_N, 0, kind(phase)) * phase * flipper
      Jout(i,i1,2) = cmplx(thetaJ_S, 0, kind(phase)) * phase * flipper
   End Do
End Do

return
99 format(2E13.5)
end subroutine Spheroid_RJ

subroutine Spheroid_Rdot(nn, Rdot)
!=================================
! purpose: calculate Rdot = \dot\rho for timestep it

use affine
use Model

implicit none

integer,          intent(in)                        :: nn
real(kind=wp),    intent(out), dimension (nn,nn,2)  :: Rdot

real(kind=wp),    dimension (2*nn*nn)               :: duhat_du
real(kind=wp)                                       :: rhatN,  rhatS, &
                   &    uhN, uhS, p_omN, p_omS, sigN, sigS, sig2N, sig2S, &   
                   &    rhat_dotN, rhat_dotS, omdotN, omdotS
real(kind=wp)                                       :: dummy_t = 0.0_wp
integer                                             :: indexN, indexS
integer                                             :: i, i1   ! counters
!! === executable statements === !!
 
duhat_du =  udot(dummy_t, uhatofu)

Do i = 1, nn
   Do i1 = 1, nn

      indexN = rank( nn*(i-1) + i1 )
      indexS = rank( nn*(i-1) + i1 + nn*nn) 
      
      uhN = uhatofu(indexN)
      uhS = uhatofu(indexS)

      p_omN = p_omega(indexN) 
      p_omS = p_omega(indexS)

      sigN = sigma(indexN) 
      sigS = sigma(indexS)

      sig2N = sigN**2
      sig2S = sigS**2

      rhatN = sqrt(abs((uhN + 0.5_wp*sigN)*  &
            &          (uhN - 0.5_wp*sigN)))
            
      rhatS = sqrt(abs((uhS + 0.5_wp*sigS)*  &
            &          (uhS - 0.5_wp*sigS)))

      rhat_dotN = uhN / rhatN
      rhat_dotS = uhS / rhatS

      omdotN = omega_dot(sig2N,p_omN,uhN)
      omdotS = omega_dot(sig2S,p_omS,uhS)

      Rdot(i,i1,1)      = (rhatN * omdotN + rhat_dotN * omvec(indexN)) &
                          *  duhat_du(indexN) / (2.0_wp * mass)
      Rdot(i,i1,2)      = (rhatN * omdotN + rhat_dotS * omvec(indexS)) &
                          *  duhat_du(indexS) / (2.0_wp * mass)
   End Do
End Do

return
end subroutine Spheroid_Rdot


subroutine Spheroid_RRdotCgamma(mass, Rout, Rdotout, Cgammaout, &
                             &  R2Cgammaout, lngammadot, Cgammadot, gugg)
!========================================================================
! purpose: calculate the area radius R, Rdot and quantities related to
! gamma for the current timestep
! note that gamma := Cgamma * sin(theta)^2;
! here we compute Cgamma and derived quantities in order not to spoil
! dividing by gamma later on 

! modules used:
use affine,         only: omega_dot, udot, sigma, sigma2, p_omega
use model,          only: eps

implicit none
! input parameters
!     mass        = horizon mass parameter (asymptotic Schwarzschild value)
real(kind=wp), intent(in)                 :: mass 

! output parameters
!     Rout        = the radius function
!     Rdotout     = d/du Rout
!     Cgammaout   = Gamma := gamma / sin(theta)^2
!     R2Cgammaout = R^2/Cgamma
!     lngammadot  = d/du ln(gamma) =  d/du ln(Cgamma)  
!     Cgammadot   = d/du Cgamma
!     gugg        = [d/du ln(gamma)]/gamma  
real(kind=wp), intent(out), dimension (:) ::   &
  &  Rout, Rdotout, Cgammaout, R2Cgammaout, lngammadot, Cgammadot, gugg

! local variables
real(kind=wp), dimension (:), allocatable ::   &
              & omdot, rhat, sighalf, rhatdot, duhat_du
real(kind=wp)                             :: dummy_t
integer                                   :: i   ! counter
integer                                   :: mx  ! # of grid points

!< === executable statements === >!

mx = size(Rout)
allocate(omdot(mx), rhat(mx), sighalf(mx), rhatdot(mx), duhat_du(mx))

Do i = 1, mx
   omdot(i) = omega_dot(sigma2(i), p_omega(i), uhatofu(i))
end do

duhat_du   = udot(dummy_t, uhatofu)
sighalf    = 0.5_wp*sigma
sighalf(1) = 0.0_wp
sighalf(mx)= 0.0_wp
 
rhat       = sqrt(abs( (uhatofu - sighalf) * (uhatofu + sighalf) ))
rhat(1)    = -uhatofu(1)
rhat(mx)   = -uhatofu(mx)
Rout       = rhat * omvec

rhatdot    = uhatofu / rhat  ! < 0
rhatdot(1) = -1.0_wp
rhatdot(mx)= -1.0_wp
Rdotout    = (rhat * omdot + rhatdot * omvec) * duhat_du

sighalf    = sighalf*flipper
Cgammaout  = (uhatofu + sighalf) / (uhatofu - sighalf)     ! Capital Gamma
Cgammadot  = duhat_du*(-2.0 * sighalf)/(uhatofu - sighalf)**2 ! d/d u (-"-)

R2Cgammaout= Cgammaout / Rout**2 ! Cgamma / r^2

lngammadot =  duhat_du*(-2.0_wp * sighalf) / (uhatofu**2 - sighalf**2)

gugg = flipper*eps/sqrt(1.0_wp + eps * cos_theta**2)**3 ! -sigma/sin2theta
gugg = duhat_du * gugg / (uhatofu + sighalf)**2         ! gamma,u/gamma^2

deallocate(omdot, rhat, sighalf, rhatdot, duhat_du)
return
end subroutine Spheroid_RRdotCgamma


end module gridtranslator
