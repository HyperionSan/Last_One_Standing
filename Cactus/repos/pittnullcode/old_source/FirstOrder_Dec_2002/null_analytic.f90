module null_analytic
contains
subroutine null_data

!!! Initialize J.

   use null_grid
   use null_vars
   use null_params
   use particle, only : p_mass, p_zn
   use model, only: Mass
   
   implicit none

   integer i,j,k
   double complex b1,b2,c1,c2,t2,t3,znorth,zbnorth,z0,zb0
   double precision P0,ppnorth
   print *
   print *, 'null_analytic: set initial data for J on outgoing null surface'

   print *, 'ID_CASE = ', ID_CASE
   select case (ID_CASE)

!!!  set an incoming pulse
     case(2)
       print *, 'null_analytic: ID_CASE = 2, ID_AMP = ', ID_AMP
       jnn = (0.d0, 0.d0)
       do i=1, nx-1
         if(rb(i).le.R_right.and.rb(i).ge.R_left) then
               jnn(:,:,i) = z**2/(1.+z*zb)**2 * &
                     (1.-R_left/rb(i))**4*(1.-R_right/rb(i))**4*ID_AMP
          end if
       end do

!!!  set J according to Newtonian limit
     case(1)
       print *, 'null_analytic: ID_CASE = 1, set J to Newtonian limit'

!!!  set J = 0
     case(0)
       print *, 'null_analytic: ID_CASE = 0, set J = 0'
       jnn = (0.d0, 0.d0)

     case DEFAULT
       print*,'You must set a value for ID_CASE of 0, 1 or 2 in null.in,'
       print*,'so that null_analytic.f90 can function.'
       print*,'Program aborting.'
       STOP

     END select
    print* 
end subroutine null_data

end module null_analytic
