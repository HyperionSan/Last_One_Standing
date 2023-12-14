module geomview
contains

subroutine zmeshr(nn, it, scale, offset, filnam, func)

   implicit none

   integer, intent (in) :: nn, it
   character*(*), intent (in) :: filnam
   double precision, intent (in) :: scale, offset, func(nn,nn)

   integer :: i, j, iunit = 11, namlen
   character*10 :: zmeshf

   namlen = len(filnam)
   do i = 1, namlen
      write (zmeshf(i:i),'(a)') filnam(i:i)
   end do
   write (zmeshf(namlen+1:namlen+6),'(i5)') 10000 + it
   zmeshf(namlen+1:namlen+1) = '_'
   open (unit = iunit, file = zmeshf(1:namlen+6), status = 'unknown')
   write (unit = iunit, fmt = '("ZMESH ",i4,1x,i4)') nn, nn
   do j = 1, nn
      do i = 1, nn
         write (unit = iunit, fmt = '(e13.3)') offset + scale * func(i,j)
      end do
   end do
   close (unit = iunit)
   
end subroutine zmeshr

subroutine zmeshc(nn, it, scale, offset, filnam, func)

   implicit none

   integer, intent (in) :: nn, it
   character*(*), intent (in) :: filnam
   double precision, intent (in) :: scale, offset
   double complex, intent (in) :: func(nn,nn)

   call zmeshr(nn, it, scale, offset, filnam // 'r',  dble(func))
   call zmeshr(nn, it, scale, offset, filnam // 'i', dimag(func))

end subroutine zmeshc

end module geomview
