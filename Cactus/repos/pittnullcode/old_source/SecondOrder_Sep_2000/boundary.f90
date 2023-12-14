module boundary

contains

subroutine boundarexp

   use hdrive
   use null_vars
   use null_grid
   use horizon
   use null_params

   implicit none
	
   integer i, minr, maxr, pr, ret, pos_n, pos_s, jj, l, ck, cran, posit

   double precision, dimension (:,:), allocatable::x0s,x0n,drs,drn,drhs,drhn
   double precision,dimension(:,:), allocatable::dxs,dxn,dxhn,dxhs
   double complex, dimension (:,:), allocatable :: q
   double precision func, r_m
   double precision u, rinv, rhinv, d1, d2, fac_n, fac_s
   double precision normj, normu, normq, normw, normb
  
   double precision, dimension (:,:,:), allocatable:: tbn,tbs,twn,tws
   double complex, dimension (:,:,:), allocatable:: ww, ww_r
 
 
   allocate ( x0s(nn,nn), x0n(nn,nn), drs(nn,nn), drn(nn,nn), & 
     &        drhs(nn,nn), drhn(nn,nn), q(nn,nn),dxhs(nn,nn), &
     &        dxhn(nn,nn), dxn(nn,nn),dxs(nn,nn), ww(nn,nn,2), ww_r(nn,nn,2) )

   print*, 'call hdriver'
   call hdriver (time, dt, dd, nn, it, nt)

   r_m = 2.*mass

   e2bmn = e2beta(:,:,1)
   e2bms = e2beta(:,:,2)

   ! to get the first point out
   maskn = int(( r_m*rhonew(:,:,1)/(1.+r_m*rhonew(:,:,1)) - 0.5)/dx ) + 1  +  1  
   masks = int(( r_m*rhonew(:,:,2)/(1.+r_m*rhonew(:,:,2)) - 0.5)/dx ) + 1  +  1    

   minr = min(minval(masks),minval(maskn))
   maxr = max(maxval(masks),maxval(maskn))

   write(*,'(a,2I4)')  ' minr, maxr in boundarexp:', minr, maxr

   ww = (vnew - r_m*rhonew)/ (r_m*rhonew)**2
   ww_r = (vrnew + 1. - 2. * r_m*rhonew* ww )/ (r_m*rhonew)**2

   do jj=1, nn
      do l=1, nn
          
         pos_s = masks(jj,l)-6
         pos_n = maskn(jj,l)-6
             
         fac_s = 1.
         fac_n = 1.
         do i = nx, 1, -1
             
            if (i.lt.pos_s) fac_s = 0.
            if (i.lt.pos_n) fac_n = 0.
             
            drs(jj,l)  = rb(i) - r_m*rhonew(jj,l,2)
            drn(jj,l)  = rb(i) - r_m*rhonew(jj,l,1)

            drhs(jj,l) = rbh(i) - r_m*rhonew(jj,l,2)
            drhn(jj,l) = rbh(i) - r_m*rhonew(jj,l,1)

            dxs(jj,l) = x(i)- r_m*rhonew(jj,l,2)/(1.+r_m*rhonew(jj,l,2))
            dxn(jj,l) = x(i)- r_m*rhonew(jj,l,1)/(1.+r_m*rhonew(jj,l,1))
            dxhs(jj,l) = xh(i)- r_m*rhonew(jj,l,2)/(1.+r_m*rhonew(jj,l,2))
            dxhn(jj,l) = xh(i)- r_m*rhonew(jj,l,1)/(1.+r_m*rhonew(jj,l,1))
                 

            jns(jj,l,i) = (jnew(jj,l,2) + jrnew(jj,l,2) *drs(jj,l) )*fac_s &
                              + (1.-fac_s)*jns(jj,l,pos_s) 
            jnn(jj,l,i) = (jnew(jj,l,1) + jrnew(jj,l,1) *drn(jj,l) )*fac_n &
                              + (1.-fac_n)*jnn(jj,l,pos_n) 

            bns(jj,l,i) = (betanew(jj,l,2) + betarnew(jj,l,2) *drs(jj,l) )*fac_s &
                              + (1.-fac_s)*bns(jj,l,pos_s) 
            bnn(jj,l,i) = (betanew(jj,l,1) + betarnew(jj,l,1) *drn(jj,l) )*fac_n &
                              + (1.-fac_n)*bnn(jj,l,pos_n)   
         
            uns(jj,l,i) = (unew(jj,l,2) + urnew(jj,l,2) *drhs(jj,l) )*&
                            fac_s/e2bms(jj,l) &
                              + (1.-fac_s)*uns(jj,l,pos_s)  
            unn(jj,l,i) = (unew(jj,l,1) + urnew(jj,l,1) *drhs(jj,l) )*&
                           fac_n/e2bmn(jj,l) &
                              + (1.-fac_n)*unn(jj,l,pos_n) 
                              
            wns(jj,l,i) = (ww(jj,l,2) + ww_r(jj,l,2) *drs(jj,l) )*&
                           fac_s/e2bms(jj,l) &
                              + (1.-fac_s)*wns(jj,l,pos_s)  
            wnn(jj,l,i) = (ww(jj,l,2) + ww_r(jj,l,2) *drn(jj,l) )*&
                           fac_n/e2bmn(jj,l) &
                              + (1.-fac_n)*wnn(jj,l,pos_n)   
                
	  end do	
		                      
      end do
   end do

!note we took care of renormalizations!!
! (in j will work things from the inside)
	
   deallocate ( x0s, x0n, drs, drn, drhs, drhn,dxs,dxn,dxhs,dxhn, ww, ww_r )
   deallocate ( q )

end subroutine

end module
