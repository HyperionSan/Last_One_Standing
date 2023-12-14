program null

   use null_grid
   use null_analytic
   use null_code
   use null_io
   use null_vars
   use null_news
   use coordtrans
   use ascii_io
   use boundary

!for the news
   use null_newsvarnorth
   use null_newsvarsouth
   use null_coortranvars
!end news modules

   implicit none

!news output....
   logical output
   integer pis, coun,ko,  kk,  poutx, pouty, px2, py2, ret
   double precision dtprime, outu, outb, outw, outj

!------ for bondi time tracking
   double precision centraltime,  dtnew, dtold

!------ for the news and transfomation to bondi time and angle frame.

   double precision, allocatable :: newsn(:,:), newss(:,:)
   double precision, allocatable, dimension(:,:):: uscrin, uscris
   double precision, allocatable, dimension(:,:)::zeta,zetab,jbon,nbon,jbos,nbos
   double precision, allocatable, dimension(:,:):: bonditime, dubonew, dubold

!end news output

   call null_read_params
   call null_allocate
   call null_setup_grid

!correct news output after transformation
   allocate (newsn(nn,nn), newss(nn,nn), uscrin(nn,nn), uscris(nn,nn) )
   allocate(zeta(nn,nn), zetab(nn,nn), nbon(nn,nn), bonditime(nn,nn), &
     &     jbon(nn,nn), dubonew(nn,nn), dubold(nn,nn),nbos(nn,nn),jbos(nn,nn)  )
   open(unit=91, file="jb1.dat", status="unknown")
   open(unit=92, file="jb2.dat", status="unknown") 

!***************for output
   coun = 20 
   ko = int( (0.8 - 0.5)/dx ) +1
   poutx = int( 1.5/dd ) + 3
   pouty = poutx 
   px2 = int(0.25*7./dd) + 3

   centraltime = time
   dtnew = 0. 
   dtold = 0.
   dubonew = 0.
   dubold = 0.
   bonditime = time
   zeta  = qs + (0.,1.) * ps
   zetab = qs - (0.,1.) * ps

!***************end for output
!for the transformation of the news

   print*, 'start time', time
   print*, 'delta t', dt

   call null_data
   call null_copy
   call null_io_ascii

   do it = 1, nt
     print*, '*** time step it =', it, '@ time =', time, '***'

     call boundarexp
     call null_evolve

!now calculate the news....
!	    call news(newsn, .5*(unn(:,:,nx)+unn(:,:,nx-1)),.5*(uon(:,:,nx)+ &
!     &        uon(:,:,nx-1)), bnn(:,:,nx), bon(:,:,nx),  &
!     &      newss,.5*(uns(:,:,nx)+uns(:,:,nx-1)),.5*(uos(:,:,nx)+ &
!     &        uos(:,:,nx-1)),bns(:,:,nx), bos(:,:,nx-1) )
!
!	    uscrin = 0.25*(unn(:,:,nx)+unn(:,:,nx-1)+uon(:,:,nx)+uon(:,:,nx-1))
!            uscris = 0.25*(uns(:,:,nx)+uns(:,:,nx-1)+uos(:,:,nx)+uos(:,:,nx-1))
!	
!	      call null_coord(newsn, newss,uscrin, uscris, &
!     &                    .5* ( bnn(:,:,nx)+bon(:,:,nx) ) , + &
!     &                       .5* ( bns(:,:,nx)+bos(:,:,nx) )  )
!
!	      call null_coord_s(newsn, newss,uscrin, uscris, &
!     &                    .5* ( bnn(:,:,nx)+bon(:,:,nx) ) , + &
!     &                       .5* ( bns(:,:,nx)+bos(:,:,nx) )  )

!****************************
!for output the right time!!!!
!     if (it.ne.0) then
!        dubonew = utimebon
!        dtnew = dt
!	 bonditime = bonditime + 0.5 * (dubonew + dubold)
!	 centraltime = centraltime + 0.5 * (dtnew + dtold)
!        dubold = dubonew
!        dtold = dtnew
!     end if
! ***************************

!end now calculate the news

      call null_copy
      call null_io_ascii

      time = time + dt

!finally output the news.....
!     if (mod(it,it_skip).eq.0) then
!     nbon = newsout * zetab/(zeta + 0. *dd**3) 
!     nbos = newsout_s * zetab/(zeta + 0. *dd**3) 
!     jbon = zetab/(zeta + 0.2 *dd**2) * jns(:,:,nx) 
!     write(91,*) bonditime(poutx,pouty),dble(nbon(poutx,pouty)), &
!     &     dimag(nbon(poutx,pouty))
!     write(92,*) bonditime(px2, pouty),maxval(abs(dble(nbon*maskscri))), &
!     &     maxval(abs(dimag(nbon*maskscri)))
!     call flush(91, kk)
!     call flush(92,kk)
!     ret = gft_write ('renews_n',time,(/nn,nn/),2,dble(abs(maskscri)*nbon) )
!     ret = gft_write ('imnews_n',time,(/nn,nn/),2,dimag(abs(maskscri)*nbon))
!     ret = gft_write ('renews_s',time,(/nn,nn/),2,dble(abs(maskscri_s)*nbos) )
!     ret = gft_write ('imnews_s',time,(/nn,nn/),2,dimag(abs(maskscri_s)*nbos))
!     end if
!end output the right news>.....

   end do

end program null
