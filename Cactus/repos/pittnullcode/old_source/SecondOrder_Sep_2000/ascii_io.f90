module ascii_io

contains

integer function gft_write (name, time, array_shape, rank, data2d)

   implicit none

   character*(*),                     intent (in) :: name
   double precision,                  intent (in) :: time
   integer, dimension (:),            intent (in) :: array_shape
   integer,                           intent (in) :: rank
   double precision, dimension (:,:), intent (in) :: data2d

   integer, parameter :: iunit = 82
   integer :: i, j, n1, n2
   double precision :: localtime

   if (rank .ne. 2) then
      return
   end if

   localtime = time

   n1 = array_shape(1)
   n2 = array_shape(2)

   open (unit = iunit, file = name // '.dat', position = 'append')

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

end module ascii_io
