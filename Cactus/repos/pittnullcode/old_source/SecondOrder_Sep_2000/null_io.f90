module null_io

contains

subroutine null_io_ascii

   use null_grid
   use null_vars
!  use null_newsvarnorth

   use null_params
   use ascii_io

   implicit none

   character*10 gfile, rfile, ufile, sfile, jfile, kfile,vfile
   integer :: pis, ret, i, posit, itout
   integer, parameter :: iunit = 82

RETURN ! die, you shit, die!

   posit = int((nn-5)/2) + 3
   pis = maxval(masks) 
 
   if (time.le.-.9437) then
      itout = it_skip 
   else
      itout = 1
   end if
       
   if (mod(it,itout).eq.0) then

      open (unit = iunit, file = 'dat/sizes.dat', position = 'append')
      write (unit = iunit, fmt = '(i5,6(1x,e14.8))') it, time, &
            max(maxval(abs(jnn)), maxval(abs(jns))), &
            max(maxval(abs(bnn)), maxval(abs(bns))), &
            max(maxval(abs(unn)), maxval(abs(uns))), &
            max(maxval(abs(wnn)), maxval(abs(wns)))
      close (unit = iunit)
	    
!     open (unit = 83, file = 'dat/sizesb.dat', position = 'append')
!     write (unit = 83, fmt = '(i5,6(1x,e14.8))') it, time, &
!           max(maxval(abs(jexpns)), maxval(abs(jexpnn))), &
!           max(maxval(abs(bexpnn)), maxval(abs(bexpns))), &
!           max(maxval(abs(uexpnn)), maxval(abs(uexpns))), &
!           max(maxval(abs(wexpnn)), maxval(abs(wexpns)))	    
!     close (unit = 83)

      write (*,'(i5,6(1x,e16.8))') it, time, &
            max(maxval(abs(jnn)), maxval(abs(jns))), &
            max(maxval(abs(bnn)), maxval(abs(bns))), &
            max(maxval(abs(unn)), maxval(abs(uns))), &
            max(maxval(abs(wnn)), maxval(abs(wns)))

   end if

   20 format(2(3x,e16.8))

   if (mod(it,itout) .eq. 0) then

      write (gfile,'(i10)') 1000000000 + it
      write (ufile,'(i10)') 1000000000 + it
      write (vfile,'(i10)') 1000000000 + it
      write (jfile,'(i10)') 1000000000 + it
      write (kfile,'(i10)') 1000000000 + it
      
      gfile(1:5) = 'dat/g'
      jfile(1:5) = 'dat/j'
      kfile(1:5) = 'dat/k'
      ufile(1:5) = 'dat/u'
      vfile(1:5) = 'dat/v'

      open (30, file = gfile, status = 'unknown') 
      open (33, file = jfile, status = 'unknown') 
      open (34, file = kfile, status = 'unknown') 
      open (35, file = ufile, status = 'unknown') 
      open (36, file = vfile, status = 'unknown') 

!     do i = masks(3,posit)-2, nx 
!           write (30,*) x(i)-rns(3,posit)/(1.+rns(3,posit)), wns(3,posit,i)
!           write (33,*) x(i)-rns(3,posit)/(1.+rns(3,posit)), dble(jns(3,posit,i))
!           write (34,*) x(i)-rns(3,posit)/(1.+rns(3,posit)), dimag(jns(3,posit,i))
!           write (35,*) x(i)-rns(3,posit)/(1.+rns(3,posit)), dble(uns(3,posit,i))
!           write (36,*) x(i)-rns(3,posit)/(1.+rns(3,posit)), dimag(uns(3,posit,i))
!     end do

      close (30)
      close (33)
      close (34)
      close (35)
      close (36)

      write (60, '(a)') 'splot "' // gfile(5:10) // '"'

      write (rfile,'(i10)') 1000000000 + it
      rfile(1:5) = 'dat/r'
      open (32, file = rfile, status = 'unknown') 

      do i = 2, nx-1
         write (32,*) x(i), bns(3,posit,i)
      end do
      close (32)

      write (62, '(a)') 'splot "' // rfile(5:10) // '"'

   end if

   if (null_output .and. mod(it,itout) .eq. 0) then

      ret = gft_write ('rj_100',time, (/nn, nx/), 2, dble(jns(3,:,:)) )
      ret = gft_write ('b_100',time,  (/nn, nx/), 2, (bns(3,:,:)) )

      ret = gft_write ('dat/rej',time, (/nn, nn/), 2, dble(jns(:,:,nx)) )
      ret = gft_write ('dat/imj',time, (/nn, nn/), 2, dimag(jns(:,:,nx)) )
      ret = gft_write ('dat/reu',time, (/nn, nn/), 2, dble(uns(:,:,nx)) )
      ret = gft_write ('dat/imu',time, (/nn, nn/), 2, dimag(uns(:,:,nx)) )
      ret = gft_write ('dat/bet',time, (/nn, nn/), 2, (bns(:,:,nx)) )
      ret = gft_write ('dat/w'  ,time, (/nn, nn/), 2, (wns(:,:,nx)) )
      ret = gft_write ('dat/rej_m',time, (/nn, nn/), 2, dble(jns(:,:,pis)) )
      ret = gft_write ('dat/imj_m',time, (/nn, nn/), 2, dimag(jns(:,:,pis)) )
      ret = gft_write ('dat/reu_m',time, (/nn, nn/), 2, dble(uns(:,:,pis)) )
      ret = gft_write ('dat/imu_m',time, (/nn, nn/), 2, dimag(uns(:,:,pis)) )
      ret = gft_write ('dat/bet_m',time, (/nn, nn/), 2, (bns(:,:,pis)) )
      ret = gft_write ('dat/w_m'  ,time, (/nn, nn/), 2, (wns(:,:,pis)) )
       
      ret = gft_write ('dat/rej_n',time, (/nn, nn/), 2, dble(jnn(:,:,nx)) )
      ret = gft_write ('dat/imj_n',time, (/nn, nn/), 2, dimag(jnn(:,:,nx)) )
      ret = gft_write ('dat/reu_n',time, (/nn, nn/), 2, dble(unn(:,:,nx)) )
      ret = gft_write ('dat/imu_n',time, (/nn, nn/), 2, dimag(unn(:,:,nx)) )
      ret = gft_write ('dat/bet_n',time, (/nn, nn/), 2, (bnn(:,:,nx)) )
      ret = gft_write ('dat/w_n'  ,time, (/nn, nn/), 2, (wnn(:,:,nx)) )
      ret = gft_write ('dat/rej_m_n',time, (/nn, nn/), 2, dble(jnn(:,:,pis)) )
      ret = gft_write ('dat/imj_m_n',time, (/nn, nn/), 2, dimag(jnn(:,:,pis)) )
      ret = gft_write ('dat/reu_m_n',time, (/nn, nn/), 2, dble(unn(:,:,pis)) )
      ret = gft_write ('dat/imu_m_n',time, (/nn, nn/), 2, dimag(unn(:,:,pis)) )
      ret = gft_write ('dat/bet_m_n',time, (/nn, nn/), 2, (bnn(:,:,pis)) )
      ret = gft_write ('dat/w_m_n'  ,time, (/nn, nn/), 2, (wnn(:,:,pis)) )

   end if

end subroutine null_io_ascii

end module null_io
