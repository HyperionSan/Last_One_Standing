module ascii_io

contains

integer function gft_write (name, time, data2d)

   implicit none

   character*(*),                     intent (in) :: name
   double precision,                  intent (in) :: time
   double precision, dimension (:,:), intent (in) :: data2d

   integer, parameter :: iunit = 82
   integer :: i, j, n1, n2

   n1 = size(data2d,1)
   n2 = size(data2d,2)

   open (unit = iunit, file = name // '.dat', position = 'append')

   write (unit = iunit, fmt = '(a2,e21.13e3)') '# ', time
   do j = 1, n1
      do i = 1, n2
         write (unit = iunit, fmt = '(e21.13e3)') data2d(j,i)
      end do
      write (unit = iunit, fmt = '("")')
   end do
   write (unit = iunit, fmt = '("")')

   close (unit = iunit)

   gft_write = 0

end function gft_write

integer function gft_close_all()
   gft_close_all = 0
end function gft_close_all

integer function gft_write_mask (name, time, data2d, mask)

   implicit none

   character*(*),                     intent (in) :: name
   double precision,                  intent (in) :: time
   double precision, dimension (:,:), intent (in) :: data2d
   double precision, dimension (:,:), intent (in) :: mask

   integer, parameter :: iunit = 82
   integer :: i, j, nq, np, num
   double precision :: q, p, dq, dp

   nq = size(data2d,1)
   np = size(data2d,2)
   
   dq = 2. / dble(nq-5)
   dp = 2. / dble(np-5)

   open (unit = iunit, file = name // '.dat', position = 'append')

   write (unit = iunit, fmt = '(a2,e21.13e3)') '# ', time

   do j = 1, nq
      p = -1. + (j - 3) * dp
      num = 0
      do i = 1, np
         q = -1. + (i - 3) * dq
         if (mask(j,i) .ne. 0.0d0) then
            write (unit = iunit, fmt = '(2f10.6,1x,e21.13e3)') q, p, data2d(j,i)
            num = num + 1
         end if
      end do
      if (num .ne. 0) write (unit = iunit, fmt = '("")')
   end do

   do i = 1, np
      q = -1. + (i - 3) * dq
      num = 0
      do j = 1, nq
         p = -1. + (j - 3) * dp
         if (mask(j,i) .ne. 0) then
            write (unit = iunit, fmt = '(2f10.6,1x,e21.13e3)') q, p, data2d(i,j)
            num = num + 1
         end if
      end do
      if (num .ne. 0) write (unit = iunit, fmt = '("")')
   end do

   write (unit = iunit, fmt = '("")')

   close (unit = iunit)

   gft_write_mask = 0

end function gft_write_mask

end module ascii_io
