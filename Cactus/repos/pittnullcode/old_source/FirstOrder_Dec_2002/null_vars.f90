module null_vars

   implicit none

   double complex,   dimension (:,:,:), allocatable, save :: &
       jnn,  jon,  &
      nunn, nuon, &
      cbnn, cbon, &
      cknn, ckon, &
       unn,  uon

   double precision, dimension (:,:,:), allocatable, save :: &
      bnn, bon, &
      wnn, won

   integer,          dimension (:,:),   allocatable, save :: maskn

contains

subroutine null_vars_allocate

   use null_grid, only : nn, nx

   allocate ( jnn(nn,nn,nx),  jon(nn,nn,nx))
   allocate ( nunn(nn,nn,nx), nuon(nn,nn,nx))
   allocate ( cknn(nn,nn,nx), ckon(nn,nn,nx))
   allocate ( bnn(nn,nn,nx),  bon(nn,nn,nx))
   allocate ( cbnn(nn,nn,nx), cbon(nn,nn,nx))
   allocate ( unn(nn,nn,nx),  uon(nn,nn,nx))
   allocate ( wnn(nn,nn,nx),  won(nn,nn,nx))

   allocate ( maskn(nn,nn))

end subroutine null_vars_allocate

end module null_vars
