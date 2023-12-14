module null_mask

contains

subroutine null_remask (k, masks, mask)

   implicit none

   integer,                  intent (in)  :: k
   integer, dimension (:,:), intent (in)  :: masks
   integer, dimension (:,:), intent (out) :: mask

   integer n1, n2, i, j

   n1 = size(masks,1)
   n2 = size(masks,2)

   ! we could use a where, but the performance of the where sucks...

   do j = 1, n2
      do i = 1, n1
         if (k .ge. masks(i,j)) then
            mask(i,j) = 1
         else
            mask(i,j) = 0
         end if
      end do
   end do

end subroutine null_remask

end module null_mask
