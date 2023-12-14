module axisymmetricio

   implicit none

contains


subroutine axi_outr(mx, index, filnam, gridfunc)

   implicit none

   integer, intent (in) :: mx, index
   character*(*) filnam
   double precision, dimension (mx), intent (in) :: gridfunc

   integer :: i, iunit = 11
   double precision :: dx, x

   open (unit = iunit, file = filnam // '.gnu', status = 'unknown', &
         position = 'append')
   write (unit = iunit, fmt = '(a,i4)')  &
         'plot "' // filnam // '.axidat" index ', index
   close (unit = iunit)

   open (unit = iunit, file = filnam // '.axidat', status = 'unknown', &
         position = 'append')

!  write (unit = iunit, fmt = '("# ",e21.13e3)') time

   dx = 2.0 / dble(mx - 1)

   do i = 1, mx
         x = -1.0 + (i-1)*dx
         write (unit = iunit, fmt = '(f10.6,1x,e17.8)') x,  gridfunc(i)
   end do
   write (unit = iunit, fmt = '("")')
   write (unit = iunit, fmt = '("")')

   close (unit = iunit)
   
end subroutine axi_outr

end module axisymmetricio
