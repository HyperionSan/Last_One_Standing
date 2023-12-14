module null_vars

   implicit none

   double complex,   dimension (:,:,:), allocatable, save :: jnn, jon, jns, jos
   double precision, dimension (:,:,:), allocatable, save :: bnn, bon, bns, bos
   double complex,   dimension (:,:,:), allocatable, save :: unn, uon, uns, uos
   double precision, dimension (:,:,:), allocatable, save :: wnn, won, wns, wos
   double precision, dimension (:,:),   allocatable, save :: e2bmn, e2bms
   integer,          dimension (:,:),   allocatable, save :: maskn, masks


   double complex,   dimension (:,:), allocatable, save :: &
          cj0s, cj0n, cj1n, cj1s, cj2n, cj2s
   double precision, dimension (:,:), allocatable, save :: &
          beta0s, beta0n, beta1n, beta1s,  beta2s, beta2n
   double complex,   dimension (:,:), allocatable, save :: &
          cu0s, cu0n, cu1n, cu1s, cu2n, cu2s
   double precision, dimension (:,:), allocatable, save :: &
          cw0s, cw0n, cw1n, cw1s, cw2n, cw2s

contains

subroutine null_vars_allocate

   use null_grid, only : nn, nx

   allocate (jns(nn,nn,nx), jos(nn,nn,nx), jnn(nn,nn,nx), jon(nn,nn,nx), &
         bns(nn,nn,nx), bos(nn,nn,nx), bnn(nn,nn,nx), bon(nn,nn,nx), &
         uns(nn,nn,nx), uos(nn,nn,nx), unn(nn,nn,nx), uon(nn,nn,nx), &
         wns(nn,nn,nx), wos(nn,nn,nx), wnn(nn,nn,nx), won(nn,nn,nx), &
         cj0s(nn,nn), cj0n(nn,nn), cj1n(nn,nn), cj1s(nn,nn), &
         cj2n(nn,nn), cj2s(nn,nn), &
         beta0s(nn,nn), beta0n(nn,nn), beta1n(nn,nn), beta1s(nn,nn), &
         beta2s(nn,nn), beta2n(nn,nn),&
         cu0s(nn,nn), cu0n(nn,nn), cu1n(nn,nn), cu1s(nn,nn), &
         cu2n(nn,nn), cu2s(nn,nn),  &
         cw0s(nn,nn), cw0n(nn,nn), cw1n(nn,nn), cw1s(nn,nn), &
         cw2n(nn,nn),cw2s(nn,nn), &
         masks(nn,nn), maskn(nn,nn), e2bmn(nn,nn), e2bms(nn,nn) )

end subroutine null_vars_allocate

end module null_vars
