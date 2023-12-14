module affine
!M===========
use prec
use Model, only: eps

! global variables, that we have to distribute more cleverly
double precision, dimension (:), allocatable  :: theta, sint, cost, sigma2, &
                                              &  sigma, p_omega, t_zero,    &
                                              &  u_zero, L_inf, th_minus
real(kind=wp) :: t_end           ! maximum value of affine parameter
real(kind=wp) :: t_minus         ! when do we initialize
real(kind=wp) :: b               ! value of b in spheroidal geometry
real(kind=wp) :: b2              ! assigned b*b in init
real(kind=wp) :: rinf            ! parameter R_inf in conformal factor Omega
integer       :: i_equator       ! index of the point which is next to pi/2
contains
! subroutines:  init
! functions:    omega, uzero_t, uzero_i, udot


double precision function omega(sig2,p_om,uh)
!F===============================================
implicit none

double precision, intent(in) :: sig2, p_om, uh

omega = - rinf/( uh + sig2/(12.0_wp*(p_om - uh)) )
end function omega

double precision function omega_dot(sig2,p_om,uh)
!F==================================================
implicit none

double precision, intent(in) :: sig2, p_om, uh

omega_dot = (12.0_wp*rinf*(sig2 + 12.0_wp*(-p_om + uh)**2)) /  &
&                         (sig2 + 12.0_wp*( p_om - uh)*uh)**2

end function omega_dot


subroutine init(parvec)
!S=====================
!purpose: initialization of model dependent parameters and arrays

use Model
implicit none

double precision, intent(in), dimension(:) :: parvec ! holds theta vals
integer        :: i, n_x
real(kind=wp)  :: pi_half

!executable statements

! shorthand
b2 = b*b
pi_half = asin(1.0_wp)

! transfer data from Roberto's module Model
rinf = R_inf
t_minus = t_

!allocate arrays
n_x = size(parvec)
write(10, *) 'n_x =', n_x

allocate(theta (n_x))
allocate(sigma2(n_x))
allocate(sigma (n_x))
allocate(p_omega   (n_x))
allocate(t_zero(n_x))
allocate(u_zero(n_x))
allocate(L_inf (n_x))
allocate(th_minus(n_x))

theta = acos(cost)

do i = 1, n_x 
  sigma2(i) = ( eps*(1.0_wp - cost(i)**2) )**2       / &
             !--------------------------------
                (1.0_wp + eps*cost(i)**2 )**3

  sigma(i) = eps *  (1.0_wp - cost(i)**2)            / &
                  !--------------------------------
                    sqrt((1.0_wp + eps*cost(i)**2 )**3)

if(L_inf_switch == 1) then 
  L_inf(i) = a_inf/sqrt( (1.0_wp + eps*cost(i)**2 )**3 )
else
 L_inf(i) = 1.0_wp
endif

 u_zero(i) = -(1.0_wp + eps*(1.0_wp - 0.5_wp*(1.0_wp - cost(i)**2))) &
          &                        /sqrt( (1.0_wp + eps*cost(i)**2 )**3 )
end do

i_equator = sum(minloc(abs(theta - pi_half))) ! minloc() has only 1 entry
write(10, *) 'theta - pi/2 @ numerical equator =', theta(i_equator) - pi_half


th_minus = t_
!p_omega = maxval(sigma)/sqrt(13.0_wp) + 0.00001_wp
!p_omega = sigma/3.0_wp  + 0.00001_wp
!p_omega = maxval(sigma)/3.5_wp
p_omega = p_model
t_zero = t_0
end subroutine init


function uzero_t(th)
implicit none

real(kind=wp), intent(in) :: th
!local:
real(kind=wp) :: cost_, cos2t_, temp, nenner, uzero_t

! abbreviations

cost_  = cos(th)
cos2t_ = cos(2.0_wp*th)

temp   = (1.0_wp +        b2 - cos2t_ + b2*cos2t_)/b2
nenner = - sqrt(2.0d0)*b*b2*temp**1.5

temp   = (1.0_wp + 3.0_wp*b2 - cos2t_ + b2*cos2t_)
uzero_t = temp/nenner

return 
end function uzero_t


function uzero_i(itheta)
implicit none

integer, intent(in) :: itheta
!local:
real(kind=wp) :: cos2t, first, sec, nenner, uzero_i

! abbreviations

cos2t  = cos(2.0_wp*theta(itheta))
first  = (1.0_wp +        b2 - cos2t + b2*cos2t)/b2
sec    = (1.0_wp + 3.0_wp*b2 - cos2t + b2*cos2t)
nenner = - sqrt(2.0d0)*b*b2*first**1.5

uzero_i = sec/nenner

end function uzero_i

function udot(t,y)
! udot := d u_hat/du

implicit none

real(kind=wp),               intent(in) :: t
real(kind=wp), dimension(:), intent(in) :: y
real(kind=wp), dimension(size(y))       :: udot
integer :: i
real(kind=wp)  :: uh, mu, mu2, xminus, xplus, ex

write(10,*) 'udot called @ time = ', t

do i=1, size(y)
   uh = y(i)
   mu2 = 13.0_wp*p_omega(i)**2 - sigma2(i)
   if(mu2 < 0.0_wp) then
      write(*, *)  'in integration for uhat(u) :  mu^2 ='
      write(*, 99) mu2
      write(10,*)  'in integration for uhat(u) :  mu^2 ='
      write(10,99) mu2
      STOP 'mu^2 < 0'
   endif
   mu = sqrt(mu2)
   xminus = (2.0_wp*uh - 5.0_wp*p_omega(i) - mu)**2 
   xplus  = (2.0_wp*uh - 5.0_wp*p_omega(i) + mu)**2
   ex     = 2.0_wp*p_omega(i)/mu
   
   udot(i) = (sigma2(i) + 12.0_wp*p_omega(i)*uh - 12.0_wp*uh**2)**2  *     &
        (xminus**(-ex-1.0_wp)*xplus**(ex-1.0_wp))/9.0_wp

   udot(i) = udot(i)*L_inf(i)
end do

return
99 format(1E13.5)
end function udot


function uprime_equ(uh)
! uprime_equ := d u/du_hat at equator = theta = pi/2

implicit none

real(kind=wp), intent(in) :: uh
real(kind=wp)             :: udot,  uprime_equ
integer                   :: i
real(kind=wp)  :: mu, mu2, xminus, xplus, ex, sigma2_e, p_omega_e 

sigma2_e = sigma2(i_equator)
p_omega_e    = p_omega   (i_equator)

mu2 = 13.0_wp*p_omega_e**2 - sigma2_e
if(mu2 < 0.0_wp) then
   write(*, *)  'in integration for u(uhat) :  mu^2 ='
   write(*, 99)  mu2
   write(10,*)  'in integration for u(uhat) :  mu^2 ='
   write(10,99)  mu2
   STOP 'mu^2 < 0'
endif
mu = sqrt(mu2)
xminus = (2.0_wp*uh - 5.0_wp*p_omega_e - mu)**2 
xplus  = (2.0_wp*uh - 5.0_wp*p_omega_e + mu)**2
ex     = 2.0_wp*p_omega_e/mu
   
udot = (sigma2_e + 12.0_wp*p_omega_e*uh - 12.0_wp*uh**2)**2  *     &
     (xminus**(-ex-1.0_wp)*xplus**(ex-1.0_wp))/9.0_wp

udot = udot * L_inf(i_equator)

uprime_equ = 1.0_wp/udot ! du /du_hat

return
99 format(1E13.5)
end function uprime_equ


end module affine

