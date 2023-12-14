module null_io
contains
subroutine null_io_ascii

   use null_grid
   use null_vars
   use null_params
   use ascii_io

   implicit none

   integer            :: posit1, posit2, ret
   integer, parameter :: iunit = 82
   integer, save      :: index
   logical, parameter :: radial_output = .false.

   posit1 = 3
   posit2 = int((nn-5)/2) + 3

   if (mod(it - (it_start + 1), iot0) .eq. 0) then

      open (unit = iunit, file = 'sizes.dat', position = 'append')
      write (unit = iunit, fmt = '(i5,6(1x,e15.8))') it, time, &
            maxval(abs(jnn)), &
            maxval(abs(bnn)), &
            maxval(abs(unn)), &
            maxval(abs(wnn))
      close (unit = iunit)
	    
      write (*,'(i5,6(1x,e16.8))') it, time, &
            maxval(abs(jnn)), &
            maxval(abs(bnn)), &
            maxval(abs(unn)), &
            maxval(abs(wnn))

   end if

!   if (radial_output .and. mod(it - (it_start + 1), iot0) .eq. 0) then
! 
!      call  line_io_ascii(nx, index, 'bnsL',        bns (posit1, posit2, : ))
!      call  line_io_ascii(nx, index, 'wnsL',        wns (posit1, posit2, : ))
!      call  line_io_ascii(nx, index, 'jnsLr',  dble(jns (posit1, posit2, : )))
!      call  line_io_ascii(nx, index, 'nunsLr', dble(nuns(posit1, posit2, : )))
!      call  line_io_ascii(nx, index, 'unsLr',  dble(uns (posit1, posit2, : )))
!
!      index = index + 1
!   end if

   if (mod(it  - it_start, iot2) == 0 ) then
      print *, 'dumping fields at Scri'      
      ret = gft_write ('rJ_nI',  time,   dble(jnn(:,:,nx)))
      ret = gft_write ('iJ_nI',  time,  dimag(jnn(:,:,nx)))
      ret = gft_write ('JJB_nI',  time,  dble(jnn(:,:,nx)*conjg(jnn(:,:,nx))))
      
!      ret = gft_write ('rnu_nI', time,  dble(nunn(:,:,nx)))
!      ret = gft_write ('inu_nI', time, dimag(nunn(:,:,nx)))
      ret = gft_write ('nnb_nI', time, dble(nunn(:,:,nx)*conjg(nunn(:,:,nx))))
      
      !ret = gft_write ('rK_nxI',  time,  dble(cknn(:,:,nx)))
      !ret = gft_write ('iK_nxI',  time, dimag(cknn(:,:,nx)))
      ret = gft_write ('KKB_nxI',  time, dble(cknn(:,:,nx)*conjg(cknn(:,:,nx))))
      
      ret = gft_write ('b_nI',   time,        bnn(:,:,nx) )
      
     ! ret = gft_write ('rB_nI',  time,  dble(cbnn(:,:,nx)))
     ! ret = gft_write ('iB_nI',  time, dimag(cbnn(:,:,nx)))
      ret = gft_write ('BBb_nI',  time,  dble(cbnn(:,:,nx)*conjg(cbnn(:,:,nx))))
      
      ret = gft_write ('rU_nI',  time,   dble(unn(:,:,nx)))
      ret = gft_write ('iU_nI',  time,  dimag(unn(:,:,nx)))
      ret = gft_write ('UUb_nI',  time,  dble(unn(:,:,nx)*conjg(unn(:,:,nx))))
      ret = gft_write ('W_nI',   time,        wnn(:,:,nx) )
   end if

end subroutine null_io_ascii


subroutine line_io_ascii(n, index, filnam, line)

   use null_grid
   use null_vars
   use null_params
   use ascii_io

   implicit none

   integer, intent (in) :: n, index
   character*(*) filnam
   double precision, dimension (nx), intent (in) :: line

   integer :: i, iunit = 11

   open (unit = iunit, file = filnam // '.gnu', status = 'unknown', &
         position = 'append')
   write (unit = iunit, fmt = '(a,i4)')  &
         'plot "' // filnam // '.dat" index ', index
   close (unit = iunit)

   open (unit = iunit, file = filnam // '.dat', status = 'unknown', &
         position = 'append')
   write (unit = iunit, fmt = '("# ",e21.13e3)') time

   do i = 2, n-1
      write (unit = iunit, fmt = '(f10.6,1x,e17.8)') x(i), line(i)
   end do
   write (unit = iunit, fmt = '("")')  !\__ tell gnuplot that slice ends here
   write (unit = iunit, fmt = '("")')  !/ 
   close (unit = iunit)

end subroutine line_io_ascii
end module null_io
