module boundary

   implicit none

   double precision, dimension (:,:,:), allocatable :: rbdry, xbdry

   double precision, dimension (:,:), allocatable, private :: &
      drn, drhn, dxn, dxhn

   double precision, dimension (:,:,:), allocatable, private :: ww, ww_r
 
contains

subroutine boundary_allocate

   use null_grid, only : nn
   implicit none

   allocate ( rbdry(nn,nn,2), xbdry(nn,nn,2) )

   allocate ( drn(nn,nn), drhn(nn,nn), &
              dxn(nn,nn), dxhn(nn,nn), &
              ww(nn,nn,2), ww_r(nn,nn,2) )

end subroutine boundary_allocate

subroutine boundary_expansion ( mass )

   use null_vars
   use null_grid
   use horizon
   use null_params
   use null_eth
   use ascii_io
   use point_dump

   implicit none

   double precision, intent (in) :: mass

   integer i, minr, maxr, ret, pos_n, pos_s, jj, l
   double precision r_m, fac_n, fac_s
   double precision rmin, rmax
  
   ! "rbdry" is r at the boundary, which for the fission problem is
   ! r = r_M * rho

   ! the value of r_M on the horizon (lambda = 0)

   r_m = 2. * mass

   rbdry = r_m * rhonew
   rmin = minval(rbdry)
   rmax = maxval(rbdry)

   write(*,*) 'minval(r) @ boundary =', rmin
   write(*,*) 'maxval(r) @ boundary =', rmax
!   if (minval(rbdry) .le. rwt) then
!      write(*,*) 'rwt = ', rwt
!      write(*,*) 'mrb = ', minval(rbdry)
!      write (*,*) 'minval(r) @ boundary is less than rwt, hence:'
!      write (*,*) 'RADIAL GRID DOES NOT COVER EXTERIOR SPACETIME'
!      stop 'in BOUNDARY module'
!   end if

   ! xbdry is x at the horizon, x = r / (rwt + r)

   xbdry = rbdry / (rwt + rbdry)

   ! to get the first point out [watch out for the +2 business]

   if (horizon_on_grid) then   ! note the nint versus the int
! consider removing the +1. The reason -- the horizon is falling in
! however slowly in the close limit. Right now the two filled points 
! are just outside the horizon and just outside that. It would be
! better to bracket the horizon. This is not a big problem now becuase
! the horizon should be moving in *very* slowly.
     maskn = nint((xbdry(:,:,1) - 0.5)/dx) + 1 + NIntPts
   else
     maskn = int((xbdry(:,:,1) - 0.5)/dx) + 2 + NIntPts
   end if
   minr = minval(maskn)
   maxr = maxval(maskn)

   write(*,'(a,I4)') ' min(mask) @ boundary: ', minr
   write(*,'(a,I4)') ' max(mask) @ boundary: ', maxr
   write(*,'(a,e16.5)') ' r[min(mask)]: ', rb(minr)
   write(*,'(a,e16.5)') ' r[max(mask)]: ', rb(maxr)

   ! We will come back here later: exactly how many points is enough? 
   ! What is this? "Remask" adds one, hum...
   ! For now, assume shit will happen if a point with I<=4 is used. 
   ! Better safe than sorry.
! consider worst case scenario of horion postion versus neighboring ang
! and then neighboring time step. Worst case angular 2 radial grispoints
! ( that is we can take it as a sing of lack of resolution if the radius
! of the horizon varies by more than dx for a neighboring angular gridpoints
! however the angular profile is probably something like +- 2 dp and dq)
! now we also will not want the radius to move mor than dx in one time step
! but worst case scenario is 2 because of int math. Hence the net worst case
! scenario (looking back to the previos slice) is 4 and since the startup needs
! an extra 2 we are left with 6  ... the original. In the close limit this
! should reduce to 2 since we will demand that the horizon does not change its
! radial gripoint (approximately).
! NOTE there is no check that the above requirments are satisfied, but with
! refinement they will eventully be. 

!   if (minr .le. 5) then 
   if (minr < 5) then 
      write (*,*) 'min(mask) @ boundary does not leave enough interior buffer radial points'
      stop 'in BOUNDARY module'
   end if

   ww = (vnew - rbdry) / rbdry**2
   ww_r = (vrnew - 1. - 2. * rbdry * ww )/ rbdry**2

   do jj=1, nn
      do l=1, nn
          
         pos_n = maskn(jj,l)-4

         ! write(*,*) jj, l, pos_s, minr, maxr, rwt

         fac_n = 1.
         do i = nx, 1, -1
             
            if (i.lt.pos_n) fac_n = 0.
             
            drn(jj,l)  = rb(i) - rbdry(jj,l,1)

            drhn(jj,l) = rbh(i) - rbdry(jj,l,1)

            dxn(jj,l) = x(i) - xbdry(jj,l,1)

            dxhn(jj,l) = xh(i) - xbdry(jj,l,1)

            jnn(jj,l,i) = (jnew(jj,l,1) + jrnew(jj,l,1) *drn(jj,l) )*fac_n &
                        + (1.-fac_n)*jnn(jj,l,pos_n) 

!! my 'nu' stuff
            nunn(jj,l,i) = (cnunew(jj,l,1) + cnurnew(jj,l,1) *drn(jj,l) )*fac_n &
                              + (1.-fac_n)*nunn(jj,l,pos_n)   
!! end my 'nu' stuff 
            
            bnn(jj,l,i) = (betanew(jj,l,1) + betarnew(jj,l,1) *drn(jj,l) )*fac_n &
                              + (1.-fac_n)*bnn(jj,l,pos_n)   
         
            cbnn(jj,l,i) = (cBnew(jj,l,1) + cBrnew(jj,l,1) *drn(jj,l) )*fac_n &
                              + (1.-fac_n)*cbnn(jj,l,pos_n)   
         
            cknn(jj,l,i) = (cKnew(jj,l,1) + cKrnew(jj,l,1) *drn(jj,l) )*fac_n &
                              + (1.-fac_n)*cknn(jj,l,pos_n)   
         
            unn(jj,l,i) = (unew(jj,l,1) + urnew(jj,l,1) *drhn(jj,l) )*&
                           fac_n + (1.-fac_n)*unn(jj,l,pos_n) 
                              
            wnn(jj,l,i) = (ww(jj,l,1) + ww_r(jj,l,1) *drn(jj,l) )*&
                           fac_n + (1.-fac_n)*wnn(jj,l,pos_n)   
                
	  end do	
      end do
   end do

!   call pointdump(time, abs(jnew(:,:,1))    , 'JHor')
!   call pointdump(time, abs(jrnew(:,:,1))    , 'JrHor')
!   call pointdump(time, dble(betanew(:,:,1)) , 'BHor')
!   call pointdump(time, dble(betarnew(:,:,1)), 'BrHor')
!   call pointdump(time, abs(unew(:,:,1))    , 'UHor')
!   call pointdump(time, abs(urnew(:,:,1))   , 'UrHor')
!   call pointdump(time, dble(ww(:,:,1))      , 'WHor')
!   call pointdump(time, dble(ww_r(:,:,1))    , 'WrHor')
!   call pointdump(time, dble(rhonew(:,:,1))  , 'RHOHor')
!   call pointdump(time, dble(rholnew(:,:,1)) , 'RHOlHor')
!   call pointdump(time, dble(omeganew(:,:,1)) , 'OmegaHor')
!   call pointdump(time, abs(qnew(:,:,1)) , 'QHor')
!   call pointdump(time, abs(cBnew(:,:,1)) , 'CBHor')
!   call pointdump(time, abs(cBrnew(:,:,1)) , 'CBrHor')
!   call pointdump(time, abs(ethrhol(:,:,1)) , 'ethRHOlHor')
   ! dump fields at the boundary

   if (mod(it-it_start, iot2) ==0) then
      ret = gft_write ('jjb_B',  time,   dble(jnew(:,:,1)*conjg(jnew(:,:,1))))
      ret = gft_write ('uub_B',  time,   dble(unew(:,:,1)*conjg(unew(:,:,1))))
      ret = gft_write ('b_B',    time,     betanew(:,:,1) )
      ret = gft_write ('oob_B',    time, dble(omeganew(:,:,1)*conjg(omeganew(:,:,1))))
      ret = gft_write ('b_B',    time,     betanew(:,:,1) )
      ret = gft_write ('w_B',    time,          ww(:,:,1))
      ret = gft_write ('r_B',    time,      rhonew(:,:,1)-1.)
   end if
   if (boundary_output) then
      ret = gft_write ('rej_B',  time,   dble(jnew(:,:,1)))
      ret = gft_write ('imj_B',  time,  dimag(jnew(:,:,1)))
      ret = gft_write ('reu_B',  time,   dble(unew(:,:,1)))
      ret = gft_write ('imu_B',  time,  dimag(unew(:,:,1)))
      ret = gft_write ('b_B',    time,     betanew(:,:,1) )

      ret = gft_write ('w_B',    time,          ww(:,:,1))

      ret = gft_write ('rejr_B', time,  dble(jrnew(:,:,1)))
      ret = gft_write ('imjr_B', time, dimag(jrnew(:,:,1)))
      ret = gft_write ('reur_B', time,  dble(urnew(:,:,1)))
      ret = gft_write ('imur_B', time, dimag(urnew(:,:,1)))

      ret = gft_write ('br_B',   time,    betarnew(:,:,1) )
      ret = gft_write ('wr_B',   time,        ww_r(:,:,1))

      ret = gft_write ('r_B',    time,      rhonew(:,:,1)-1.)

      ret = gft_write ('k_B',    time,        knew(:,:,1)-1.)
      ret = gft_write ('kr_B',   time,       krnew(:,:,1))

      ret = gft_write ('reK_B', time,  dble(cKnew(:,:,1)))
      ret = gft_write ('imK_B', time, dimag(cKnew(:,:,1)))
   end if

end subroutine boundary_expansion

end module
