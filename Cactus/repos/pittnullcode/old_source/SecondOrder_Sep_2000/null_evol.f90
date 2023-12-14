module null_evol

   implicit none

   double precision, dimension (:,:), allocatable, save, private :: &
       & xc, beta, dx_beta, w, k, dx_k, dx_w, f1, f2

   double complex, dimension (:,:), allocatable, save, private :: &
       & u, eth_u, dx_u, eth_dx_u, eth_beta, eth2_beta, j, n_j, h_c, &
       & jave, dx_jave, d2x_jave, test, eth_j, ethb_j, eth_k, eth_bu, &
       & eth_ethb_beta,fac, eth_dx_bu, dx_j, eth_dx_j, ethb_dx_j, jnew, &
       & n1, n2, n3, n4, n5, n6, n7, p1, p2, p3, p4

   integer, dimension (:,:), allocatable, save, private :: mask

contains

subroutine null_evol_allocate

   use null_grid, only : nn

   allocate ( xc(nn,nn), beta(nn,nn), w(nn,nn), dx_w(nn,nn), u(nn&
         & ,nn), eth_u(nn,nn), dx_u(nn,nn), eth_dx_u(nn,nn),&
         & eth_beta(nn,nn), eth2_beta(nn,nn), j(nn,nn), n_j(nn,nn),&
         & h_c(nn,nn), jave(nn,nn), dx_jave(nn,nn), d2x_jave(nn,nn),&
         & eth_ethb_beta(nn,nn), mask(nn,nn),dx_beta(nn,nn), &
         & fac(nn,nn), dx_j(nn,nn), dx_k(nn,nn), eth_k(nn,nn)&
         & ,eth_j(nn,nn), ethb_j(nn,nn), eth_bu(nn,nn),eth_dx_bu(nn&
         & ,nn),k(nn,nn), eth_dx_j(nn,nn), ethb_dx_j(nn,nn), test(nn&
         & -2,nn-2), jnew(nn,nn), f1(nn,nn), f2(nn,nn)  )

   allocate (n1(nn,nn), n2(nn,nn), n3(nn,nn), n4(nn,nn))
   allocate (n5(nn,nn), n6(nn,nn), n7(nn,nn))
   allocate (p1(nn,nn), p2(nn,nn), p3(nn,nn), p4(nn,nn))

end subroutine null_evol_allocate

subroutine null_j (i, jns, bns, uns, wns, jos, bos, uos, wos, masks, e2bm)

   use null_grid, only : nn, nx, dd, dx, dt, rwt, x, xh, pi
   use null_mask
   use null_eth
   use null_params, only : disip

   implicit none

   integer,                             intent (in)    :: i
   double complex,   dimension (:,:,:), intent (inout) :: jns
   double precision, dimension (:,:,:), intent (in)    :: bns
   double complex,   dimension (:,:,:), intent (in)    :: uns
   double precision, dimension (:,:,:), intent (in)    :: wns
   double complex,   dimension (:,:,:), intent (in)    :: jos
   double precision, dimension (:,:,:), intent (in)    :: bos
   double complex,   dimension (:,:,:), intent (in)    :: uos
   double precision, dimension (:,:,:), intent (in)    :: wos
   integer,          dimension (:,:),   intent (in)    :: masks
   double precision, dimension (:,:),   intent (in)    :: e2bm

   integer cran

   call null_remask (i, masks+1, mask)

   xc = xh(i-1)

   u = 0.25 * (uns(:,:,i-1) + uos(:,:,i-1) + uns(:,:,i-2) + uos(:,:,i))

   dx_u = 0.5 * (uns(:,:,i-1) - uns(:,:,i-2) + uos(:,:,i) - uos(:,:,i-1)) / dx

   beta = 0.25*(bns(:,:,i-1) + bos(:,:,i) )

   w = 0.5 * (wns(:,:,i-1) + wos(:,:,i) )

   if(i.ne.nx) then

      dx_beta = 0.5 * (bns(:,:,i-1) - bns(:,:,i-2) + bos(:,:,i+1)&
           & - bos(:,:,i)) / dx

      dx_w = 0.5 * (wns(:,:,i-1) - wns(:,:,i-2) + wos(:,:,i+1)   -&
           & wos(:,:,i)) / dx

      dx_j = 0.5 * (jns(:,:,i-1) - jns(:,:,i-2) + jos(:,:,i+1)   -&
           & jos(:,:,i)) / dx

   else

      dx_beta = 0.5 * (bns(:,:,i-1) - bns(:,:,i-2) + 0.5 * (3. *&
            & bos(:,:,i) - 4. * bos(:,:,i-1) + bos(:,:,i-2))) / dx


      dx_w = 0.5 * (wns(:,:,i-1) - wns(:,:,i-2) + 0.5 * (3. * wos(:&
            & ,:,i) - 4. * wos(:,:,i-1) + wos(:,:,i-2))) / dx

      dx_j = 0.5 * (jns(:,:,i-1) - jns(:,:,i-2) + 0.5 * (3. * jos(:&
            & ,:,i) - 4. * jos(:,:,i-1) + jos(:,:,i-2))) / dx

   end if
    
!!!!!mpi 22/11/99. need to hit by e^(2b_m) to w and dx_w

   w = e2bm * w
   dx_w = e2bm * dx_w
	    
!!!!!mpi    

   call null_d1 (eth_u, u, 1, 1)
   call null_d1 (eth_dx_u, dx_u, 1, 1)
   call null_d1 (eth_beta, dcmplx(beta), 0, 1)
   call null_d2 (eth2_beta, dcmplx(beta), 0, 1, 1)
   call null_d2 (eth_ethb_beta, dcmplx(beta), 0, 1, -1)
   call null_d1 (eth_bu, conjg(u), -1, 1)
   call null_d1 (eth_dx_bu, conjg(dx_u), -1, 1)

   j =  0.5*jns(:,:,i-1)+0.5*jos(:,:,i) 

   k = sqrt (1. + dble(j) ** 2 + dimag(j) ** 2)

   call null_d1 (eth_k, dcmplx(k), 0, 1)
   call null_d1 (eth_j, j, 2, 1)
   call null_d1 (ethb_j, j, 2, -1)
   call null_d1 (eth_dx_j, dx_j, 2, 1)
   call null_d1 (ethb_dx_j, dx_j, 2, -1)

   dx_k = dble( dx_j * conjg(j) ) / k 
   fac = ( - dx_k * j/k + dx_j )  

   ! non-linear terms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   n1 = - exp(2*beta) * (1.-xc)/(xc*rwt) * (    &
         &       k * ( eth_j*conjg(eth_beta) - ethb_j*eth_beta )    &
         &     + j * ( eth_beta * conjg(ethb_j) + ethb_j *conjg(eth_beta)    &
         &                   - 2.* eth_k * conjg(eth_beta) ) )

   n2 = - ( ethb_j * (xc * (1. - xc) * 0.5 * dx_u + u) +   &
         &        eth_j * (xc * (1. - xc) * conjg(dx_u) * 0.5 + conjg(u) ) )  

   n3 =  (1. - k) * (xc * (1. - xc) * eth_dx_u  + 2. * eth_u)     &
         & - j * (xc * (1. - xc) * eth_dx_bu + 2. * eth_bu)    

   n4 = 0.5 * exp(-2*beta) * xc**3 * (1.-xc) * rwt * (    &
         &       k**2 * dx_u**2 + j**2 * conjg(dx_u)**2      &
         &       + 2. * j * k * (dble(dx_u)**2 + dimag(dx_u)**2 ) )    

   n5 = - xc * (1. - xc) * dx_j * dble(eth_bu)


   n6 = +  xc * (1. - xc) * ( conjg(u) * eth_j * dble( j * conjg(dx_j) ) &
         & + u * ethb_j * (dimag(j) * dble(dx_j) - dble(j) * dimag(dx_j)) &
         & + ( conjg(u) * ethb_j - conjg(ethb_j) * u )    &
         & *  ( j * dx_k - k * dx_j ) &
         & - ethb_dx_j * u - eth_dx_j * conjg(u) )   

   n7 = + xc * (1. - xc) * (dx_j * k - dx_k * j) *  (      &
         &                2. * dimag( conjg(u) * ( ethb_j - eth_k ) )    &
         &                - k * 2. * dimag ( eth_bu )      &
         &                + 2. * dimag ( j *conjg(eth_u) ) )         

   ! p2*j/r

   p2 = + j * (1. - xc) / (xc * rwt) * exp(2. * beta) *  (    &
         &       - 2. * k * ( eth_ethb_beta + dble(eth_beta)**2    &
         &                   + dimag(eth_beta)**2  )       &
         &       + 2. * dble( ethb_j * conjg(eth_beta) )    &
         &       + 2. * dble( conjg(j) * ( eth2_beta + eth_beta**2 ) )     &
         &       - 2. * dble( eth_k * conjg(eth_beta) ) )    

   ! p3*j/r

   p3 =  + j * dble( xc * (1. - xc) * eth_dx_bu + 2. * eth_bu )   

   ! p4*j/r

   p4 =  - j * xc**3 * (1.-xc) * rwt * exp(-2. * beta) * 0.5 * (     &
         &        k * ( dble(dx_u)**2 + dimag(dx_u)**2 ) +      &
         &        dble( j * conjg(dx_u)**2 ) )  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   n_j = n1 + n2 + n3 + n4 + n5 + n6 + n7 + p2 + p3 + p4

   ! the rhs ...

   h_c = - xc * (1. - xc) * eth_dx_u - 2. * eth_u + 2 * (1. - xc) /  &
         & (rwt *xc) * exp (2. * beta) * (eth_beta * eth_beta +  &
         & eth2_beta) - xc * (1. - xc) * dx_w * j - w * j + n_j 

!!!!mpi 11/22/99 need to hit by e2bm to get things right

   h_c = e2bm * h_c

!!!!mpi

   ! interpolate the $f = x j$ values
   !crank initial

   f1 =dt*( xc * (1.-xc) * dx_w + w)
   f2 =dt* (1.-xc)**2 * (( 1.-xc) + xc * w )

   jns(:,:,i) = (1-mask)* jns(:,:,i) + jos(:,:,i) * mask

   if (i.ne.nx) then

!       jave= .5 * ( x(i)*jos(:,:,i) + x(i-1)*jos(:,:,i-1) ) +   ( 1.&
!           & -x(i) ) * ( x(i+1)*jos(:,:,i+1) - x(i-1)*jos(:,:,i-1) )&
!           & *.5 / dx + ( 1.-x(i-1) ) * ( x(i)*jos(:,:,i) - x(i-2)&
!           & *jos(:,:,i-2) ) *.25 / dx

      jave= .5 * ( x(i)*jos(:,:,i) + x(i-1)*jos(:,:,i-1) ) +  ( 1.&
 & - xc ) * ( (1.-0.75*disip)* (x(i)*jos(:,:,i) - x(i-1)*jos(:,:,i-1)) &
 & + 0.25*disip * (x(i+1)*jos(:,:,i+1) - x(i-2) *jos(:,:,i-2) ) ) / dx 


!          jave =0.5*(x(i)*jos(:,:,i)+x(i-1)*jos(:,:,i-1) ) + (1.-xc)&
!               & * ( x(i)*jos(:,:,i) - x(i-1)*jos(:,:,i-1) ) / dx


      d2x_jave= ( x(i+1)*jos(:,:,i+1)-2.*x(i)*jos(:,:,i)+x(i-1)&
            & *jos(:,:,i-1) )

      do cran = 1, 60

         jnew = jns(:,:,i) 


         p1 = j * (1. - xc)*(xc*2.*dble((jns(:,:,i) + jns(:,:,i-1)&
               & -jos(:,:,i) - jos(:,:,i-1) )/(2.*dt) * conjg(fac) )&
               & - 8. * dx_beta * (1.-xc+xc*w) )

         !	dx_jave = dx_j

         jns(:,:,i) = jns(:,:,i)*(1-mask) + mask * ( - x(i-1) *&
               & jns(:,:,i-1) * (1. - f1 * .25) * ( 0.5 - (1.-xc)/dx&
               & + 0.5 * f2/dx**2 ) + x(i-2) * jns(:,:,i-2) * f2 *&
               & .25/dx**2 + jave * (1. + f1*.25) + f2 * d2x_jave +&
               & 0.5 * (h_c + p1) * dt ) / ( x(i) * (  (1. - f1 *&
               & .25) * ( 0.5 + (1.-xc)/dx) - 0.25 * f2/dx**2 ) )

         ! new $j$ grid value

         test(:,:)= jnew(2:nn-1,2:nn-1)-jns(2:nn-1,2:nn-1,i) 

         if(maxval(abs(test) ).lt..01*dt*dd**2) exit

      end do

   else

      do cran = 1, 60 

         jnew = jns(:,:,i) 

         p1 = j * (1. - xc)*(xc*2.*dble((jns(:,:,i) + jns(:,:,i-1)&
               & -jos(:,:,i) - jos(:,:,i-1) )/(2.*dt) * conjg(fac) )&
               & - 8. * dx_beta * (1.-xc+xc*w) )


         jave =0.5*(x(i)*jos(:,:,i)+x(i-1)*jos(:,:,i-1) ) + (1.-xc)&
               & * ( x(i)*jos(:,:,i) - x(i-1)*jos(:,:,i-1) ) / dx

         !  dx_jave = dx_j

         d2x_jave= ( 2.*x(i)*jos(:,:,i)-5.*x(i-1)*jos(:,:,i-1)+4.&
               & *x(i-2)*jos(:,:,i-2) -x(i-3)*jos(:,:,i-3) )

         jns(:,:,i) = jns(:,:,i)*(1-mask) + mask * ( - x(i-1) *&
               & jns(:,:,i-1) * (1. - f1 * .25) * ( 0.5 - (1.-xc)/dx&
               & + 0.5 * f2/dx**2 ) + x(i-2) * jns(:,:,i-2) * f2 *&
               & .25/dx**2 + jave * (1. + f1*.25) + f2 * d2x_jave +&
               & 0.5 * (h_c + p1) * dt ) / ( x(i) * (   (1. - f1 *&
               & .25) * ( 0.5 + (1.-xc)/dx) - 0.25 * f2/dx**2 ) )     

         test(:,:)=jnew(2:nn-1,2:nn-1) - jns(2:nn-1,2:nn-1,i) 

         if(maxval(abs(test) ).lt..01*dt*dd**2) exit

      end do

   end if

end subroutine null_j

end module null_evol
