module null_analytic
contains
subroutine null_data

   use null_grid
   use null_vars
   use null_params
   
   implicit none

   integer i

   jns = (0., 0.)
   jnn = (0., 0.)
   bnn = 0.
   bns = 0.
   uns = (0., 0.)
   unn = (0., 0.)
   do i = 1, nx-1
      wns(:,:,i) = -2. * mass/rb(i)**2
      wnn(:,:,i) = -2. * mass/rb(i)**2
   end do
   wns(:,:,nx) = 0.
   wns(:,:,nx) = 0.
 
end subroutine null_data

end module null_analytic
