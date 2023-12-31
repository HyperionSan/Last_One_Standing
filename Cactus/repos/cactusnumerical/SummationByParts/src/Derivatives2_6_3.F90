#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine SBP_deriv2_gf_6_3 ( var, ni, nj, nk, dir, bb, gsize, delta, dvar2 )

  use All_Coeffs_mod

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(IN) :: ni, nj, nk
  CCTK_REAL, dimension(ni,nj,nk), intent(IN) :: var
  CCTK_INT, intent(IN) :: dir
  CCTK_INT, intent(IN) :: bb(2)
  CCTK_INT, intent(IN) :: gsize
  CCTK_REAL, intent(IN) :: delta
  CCTK_REAL, dimension(ni,nj,nk), intent(OUT) :: dvar2

  CCTK_REAL, dimension(4), save :: a
  CCTK_REAL, dimension(9,6), save :: q 
  CCTK_REAL :: idel2

  CCTK_INT :: il, ir, jl, jr, kl, kr

  logical, save :: first = .true.

  if ( first ) then
    call coeffs_2_6_3 ( a, q )
    first = .false.
  end if

  idel2 = 1.0_wp / delta**2

  direction: select case (dir)
  case (0) direction
    if ( bb(1) == 0 ) then
      il = 1 + gsize
    else
      dvar2(1,:,:) = ( q(1,1) * var(1,:,:) + q(2,1) * var(2,:,:)  + &
                       q(3,1) * var(3,:,:) + q(4,1) * var(4,:,:) + &
                       q(5,1) * var(5,:,:) + q(6,1) * var(6,:,:) ) * idel2
      dvar2(2,:,:) = ( q(1,2) * var(1,:,:) + q(2,2) * var(2,:,:) + &
                       q(3,2) * var(3,:,:) + q(4,2) * var(4,:,:) + &
                       q(5,2) * var(5,:,:) + q(6,2) * var(6,:,:) ) * idel2
      dvar2(3,:,:) = ( q(1,3) * var(1,:,:) + q(2,3) * var(2,:,:) + &
                       q(3,3) * var(3,:,:) + q(4,3) * var(4,:,:) + &
                       q(5,3) * var(5,:,:) + q(6,3) * var(6,:,:) ) * idel2
      dvar2(4,:,:) = ( q(1,4) * var(1,:,:) + q(2,4) * var(2,:,:) + &
                       q(3,4) * var(3,:,:) + q(4,4) * var(4,:,:) + &
                       q(5,4) * var(5,:,:) + q(6,4) * var(6,:,:) + &
                       q(7,4) * var(7,:,:) ) * idel2
      dvar2(5,:,:) = ( q(1,5) * var(1,:,:) + q(2,5) * var(2,:,:) + &
                       q(3,5) * var(3,:,:) + q(4,5) * var(4,:,:) + &
                       q(5,5) * var(5,:,:) + q(6,5) * var(6,:,:) + &
                       q(7,5) * var(7,:,:) + q(8,5) * var(8,:,:) ) * idel2
      dvar2(6,:,:) = ( q(1,6) * var(1,:,:) + q(2,6) * var(2,:,:) + &
                       q(3,6) * var(3,:,:) + q(4,6) * var(4,:,:) + &
                       q(5,6) * var(5,:,:) + q(6,6) * var(6,:,:) + &
                       q(7,6) * var(7,:,:) + q(8,6) * var(8,:,:) + &
                       q(9,6) * var(9,:,:) ) * idel2
      il = 7
    end if
    if ( bb(2) == 0 ) then
      ir = ni - gsize
    else
      dvar2(ni-5,:,:) = ( q(1,6) * var(ni,:,:) + &
                          q(2,6) * var(ni-1,:,:) + &
                          q(3,6) * var(ni-2,:,:) + &
                          q(4,6) * var(ni-3,:,:) + &
                          q(5,6) * var(ni-4,:,:) + &
                          q(6,6) * var(ni-5,:,:) + &
                          q(7,6) * var(ni-6,:,:) + &
                          q(8,6) * var(ni-7,:,:) + &
                          q(9,6) * var(ni-8,:,:) ) * idel2
      dvar2(ni-4,:,:) = ( q(1,5) * var(ni,:,:) + &
                          q(2,5) * var(ni-1,:,:) + &
                          q(3,5) * var(ni-2,:,:) + &
                          q(4,5) * var(ni-3,:,:) + &
                          q(5,5) * var(ni-4,:,:) + &
                          q(6,5) * var(ni-5,:,:) + &
                          q(7,5) * var(ni-6,:,:) + &
                          q(8,5) * var(ni-7,:,:) ) * idel2
      dvar2(ni-3,:,:) = ( q(1,4) * var(ni,:,:) + &
                          q(2,4) * var(ni-1,:,:) + &
                          q(3,4) * var(ni-2,:,:) + &
                          q(4,4) * var(ni-3,:,:) + &
                          q(5,4) * var(ni-4,:,:) + &
                          q(6,4) * var(ni-5,:,:) + &
                          q(7,4) * var(ni-6,:,:) ) * idel2
      dvar2(ni-2,:,:) = ( q(1,3) * var(ni,:,:) + &
                          q(2,3) * var(ni-1,:,:) + &
                          q(3,3) * var(ni-2,:,:) + &
                          q(4,3) * var(ni-3,:,:) + &
                          q(5,3) * var(ni-4,:,:) + &
                          q(6,3) * var(ni-5,:,:) ) * idel2
      dvar2(ni-1,:,:) = ( q(1,2) * var(ni,:,:) + &
                          q(2,2) * var(ni-1,:,:) + &
                          q(3,2) * var(ni-2,:,:) + &
                          q(4,2) * var(ni-3,:,:) + &
                          q(5,2) * var(ni-4,:,:) + &
                          q(6,2) * var(ni-5,:,:) ) * idel2
      dvar2(ni,:,:)   = ( q(1,1) * var(ni,:,:) + &
                          q(2,1) * var(ni-1,:,:) + &
                          q(3,1) * var(ni-2,:,:) + &
                          q(4,1) * var(ni-3,:,:) + &
                          q(5,1) * var(ni-4,:,:) + &
                          q(6,1) * var(ni-5,:,:) ) * idel2
      ir = ni - 6
    end if
    dvar2(il:ir,:,:) = ( a(1) * var(il:ir,:,:) + &
                         a(2) * ( var(il+1:ir+1,:,:) + &
                                  var(il-1:ir-1,:,:) ) + &
                         a(3) * ( var(il+2:ir+2,:,:) + &
                                  var(il-2:ir-2,:,:) ) + &
                         a(4) * ( var(il+3:ir+3,:,:) + &
                                  var(il-3:ir-3,:,:) ) ) * idel2
  case (1) direction
    if ( zero_derivs_y /= 0 ) then
      dvar2 = zero
    else
      if ( bb(1) == 0 ) then
        jl = 1 + gsize
      else
        dvar2(:,1,:) = ( q(1,1) * var(:,1,:) + q(2,1) * var(:,2,:)  + &
                         q(3,1) * var(:,3,:) + q(4,1) * var(:,4,:) + &
                         q(5,1) * var(:,5,:) + q(6,1) * var(:,6,:) ) * idel2
        dvar2(:,2,:) = ( q(1,2) * var(:,1,:) + q(2,2) * var(:,2,:) + &
                         q(3,2) * var(:,3,:) + q(4,2) * var(:,4,:) + &
                         q(5,2) * var(:,5,:) + q(6,2) * var(:,6,:) ) * idel2
        dvar2(:,3,:) = ( q(1,3) * var(:,1,:) + q(2,3) * var(:,2,:) + &
                         q(3,3) * var(:,3,:) + q(4,3) * var(:,4,:) + &
                         q(5,3) * var(:,5,:) + q(6,3) * var(:,6,:) ) * idel2
        dvar2(:,4,:) = ( q(1,4) * var(:,1,:) + q(2,4) * var(:,2,:) + &
                         q(3,4) * var(:,3,:) + q(4,4) * var(:,4,:) + &
                         q(5,4) * var(:,5,:) + q(6,4) * var(:,6,:) + &
                         q(7,4) * var(:,7,:) ) * idel2
        dvar2(:,5,:) = ( q(1,5) * var(:,1,:) + q(2,5) * var(:,2,:) + &
                         q(3,5) * var(:,3,:) + q(4,5) * var(:,4,:) + &
                         q(5,5) * var(:,5,:) + q(6,5) * var(:,6,:) + &
                         q(7,5) * var(:,7,:) + q(8,5) * var(:,8,:) ) * idel2
        dvar2(:,6,:) = ( q(1,6) * var(:,1,:) + q(2,6) * var(:,2,:) + &
                         q(3,6) * var(:,3,:) + q(4,6) * var(:,4,:) + &
                         q(5,6) * var(:,5,:) + q(6,6) * var(:,6,:) + &
                         q(7,6) * var(:,7,:) + q(8,6) * var(:,8,:) + &
                         q(9,6) * var(:,9,:) ) * idel2
        jl = 7
      end if
      if ( bb(2) == 0 ) then
        jr = nj - gsize
      else
        dvar2(:,nj-5,:) = ( q(1,6) * var(:,nj,:) + &
                            q(2,6) * var(:,nj-1,:) + &
                            q(3,6) * var(:,nj-2,:) + &
                            q(4,6) * var(:,nj-3,:) + &
                            q(5,6) * var(:,nj-4,:) + &
                            q(6,6) * var(:,nj-5,:) + &
                            q(7,6) * var(:,nj-6,:) + &
                            q(8,6) * var(:,nj-7,:) + &
                            q(9,6) * var(:,nj-8,:) ) * idel2
        dvar2(:,nj-4,:) = ( q(1,5) * var(:,nj,:) + &
                            q(2,5) * var(:,nj-1,:) + &
                            q(3,5) * var(:,nj-2,:) + &
                            q(4,5) * var(:,nj-3,:) + &
                            q(5,5) * var(:,nj-4,:) + &
                            q(6,5) * var(:,nj-5,:) + &
                            q(7,5) * var(:,nj-6,:) + &
                            q(8,5) * var(:,nj-7,:) ) * idel2
        dvar2(:,nj-3,:) = ( q(1,4) * var(:,nj,:) + &
                            q(2,4) * var(:,nj-1,:) + &
                            q(3,4) * var(:,nj-2,:) + &
                            q(4,4) * var(:,nj-3,:) + &
                            q(5,4) * var(:,nj-4,:) + &
                            q(6,4) * var(:,nj-5,:) + &
                            q(7,4) * var(:,nj-6,:) ) * idel2
        dvar2(:,nj-2,:) = ( q(1,3) * var(:,nj,:) + &
                            q(2,3) * var(:,nj-1,:) + &
                            q(3,3) * var(:,nj-2,:) + &
                            q(4,3) * var(:,nj-3,:) + &
                            q(5,3) * var(:,nj-4,:) + &
                            q(6,3) * var(:,nj-5,:) ) * idel2
        dvar2(:,nj-1,:) = ( q(1,2) * var(:,nj,:) + &
                            q(2,2) * var(:,nj-1,:) + &
                            q(3,2) * var(:,nj-2,:) + &
                            q(4,2) * var(:,nj-3,:) + &
                            q(5,2) * var(:,nj-4,:) + &
                            q(6,2) * var(:,nj-5,:) ) * idel2
        dvar2(:,nj,:)   = ( q(1,1) * var(:,nj,:) + &
                            q(2,1) * var(:,nj-1,:) + &
                            q(3,1) * var(:,nj-2,:) + &
                            q(4,1) * var(:,nj-3,:) + &
                            q(5,1) * var(:,nj-4,:) + &
                            q(6,1) * var(:,nj-5,:) ) * idel2
        jr = nj - 6
      end if
      dvar2(:,jl:jr,:) = ( a(1) * var(:,jl:jr,:) + &
                           a(2) * ( var(:,jl+1:jr+1,:) + &
                                    var(:,jl-1:jr-1,:) ) + &
                           a(3) * ( var(:,jl+2:jr+2,:) + &
                                    var(:,jl-2:jr-2,:) ) + &
                           a(4) * ( var(:,jl+3:jr+3,:) + &
                                    var(:,jl-3:jr-3,:) ) ) * idel2

    end if
  case (2) direction
    if ( zero_derivs_z /= 0 ) then
      dvar2 = zero
    else
      if ( bb(1) == 0 ) then
        kl = 1 + gsize
      else
        dvar2(:,:,1) = ( q(1,1) * var(:,:,1) + q(2,1) * var(:,:,2)  + &
                         q(3,1) * var(:,:,3) + q(4,1) * var(:,:,4) + &
                         q(5,1) * var(:,:,5) + q(6,1) * var(:,:,6) ) * idel2
        dvar2(:,:,2) = ( q(1,2) * var(:,:,1) + q(2,2) * var(:,:,2) + &
                         q(3,2) * var(:,:,3) + q(4,2) * var(:,:,4) + &
                         q(5,2) * var(:,:,5) + q(6,2) * var(:,:,6) ) * idel2
        dvar2(:,:,3) = ( q(1,3) * var(:,:,1) + q(2,3) * var(:,:,2) + &
                         q(3,3) * var(:,:,3) + q(4,3) * var(:,:,4) + &
                         q(5,3) * var(:,:,5) + q(6,3) * var(:,:,6) ) * idel2
        dvar2(:,:,4) = ( q(1,4) * var(:,:,1) + q(2,4) * var(:,:,2) + &
                         q(3,4) * var(:,:,3) + q(4,4) * var(:,:,4) + &
                         q(5,4) * var(:,:,5) + q(6,4) * var(:,:,6) + &
                         q(7,4) * var(:,:,7) ) * idel2
        dvar2(:,:,5) = ( q(1,5) * var(:,:,1) + q(2,5) * var(:,:,2) + &
                         q(3,5) * var(:,:,3) + q(4,5) * var(:,:,4) + &
                         q(5,5) * var(:,:,5) + q(6,5) * var(:,:,6) + &
                         q(7,5) * var(:,:,7) + q(8,5) * var(:,:,8) ) * idel2
        dvar2(:,:,6) = ( q(1,6) * var(:,:,1) + q(2,6) * var(:,:,2) + &
                         q(3,6) * var(:,:,3) + q(4,6) * var(:,:,4) + &
                         q(5,6) * var(:,:,5) + q(6,6) * var(:,:,6) + &
                         q(7,6) * var(:,:,7) + q(8,6) * var(:,:,8) + &
                         q(9,6) * var(:,:,9) ) * idel2
        kl = 7
      end if
      if ( bb(2) == 0 ) then
        kr = nk - gsize
      else 
        dvar2(:,:,nk-5) = ( q(1,6) * var(:,:,nk) + &
                            q(2,6) * var(:,:,nk-1) + &
                            q(3,6) * var(:,:,nk-2) + &
                            q(4,6) * var(:,:,nk-3) + &
                            q(5,6) * var(:,:,nk-4) + &
                            q(6,6) * var(:,:,nk-5) + &
                            q(7,6) * var(:,:,nk-6) + &
                            q(8,6) * var(:,:,nk-7) + &
                            q(9,6) * var(:,:,nk-8) ) * idel2
        dvar2(:,:,nk-4) = ( q(1,5) * var(:,:,nk) + &
                            q(2,5) * var(:,:,nk-1) + &
                            q(3,5) * var(:,:,nk-2) + &
                            q(4,5) * var(:,:,nk-3) + &
                            q(5,5) * var(:,:,nk-4) + &
                            q(6,5) * var(:,:,nk-5) + &
                            q(7,5) * var(:,:,nk-6) + &
                            q(8,5) * var(:,:,nk-7) ) * idel2
        dvar2(:,:,nk-3) = ( q(1,4) * var(:,:,nk) + &
                            q(2,4) * var(:,:,nk-1) + &
                            q(3,4) * var(:,:,nk-2) + &
                            q(4,4) * var(:,:,nk-3) + &
                            q(5,4) * var(:,:,nk-4) + &
                            q(6,4) * var(:,:,nk-5) + &
                            q(7,4) * var(:,:,nk-6) ) * idel2
        dvar2(:,:,nk-2) = ( q(1,3) * var(:,:,nk) + &
                            q(2,3) * var(:,:,nk-1) + &
                            q(3,3) * var(:,:,nk-2) + &
                            q(4,3) * var(:,:,nk-3) + &
                            q(5,3) * var(:,:,nk-4) + &
                            q(6,3) * var(:,:,nk-5) ) * idel2
        dvar2(:,:,nk-1) = ( q(1,2) * var(:,:,nk) + &
                            q(2,2) * var(:,:,nk-1) + &
                            q(3,2) * var(:,:,nk-2) + &
                            q(4,2) * var(:,:,nk-3) + &
                            q(5,2) * var(:,:,nk-4) + &
                            q(6,2) * var(:,:,nk-5) ) * idel2
        dvar2(:,:,nk)   = ( q(1,1) * var(:,:,nk) + &
                            q(2,1) * var(:,:,nk-1) + &
                            q(3,1) * var(:,:,nk-2) + &
                            q(4,1) * var(:,:,nk-3) + &
                            q(5,1) * var(:,:,nk-4) + &
                            q(6,1) * var(:,:,nk-5) ) * idel2
        kr = nk - 6
      end if
      dvar2(:,:,kl:kr) = ( a(1) * var(:,:,kl:kr) + &
                           a(2) * ( var(:,:,kl+1:kr+1) + &
                                    var(:,:,kl-1:kr-1) ) + &
                           a(3) * ( var(:,:,kl+2:kr+2) + &
                                    var(:,:,kl-2:kr-2) ) + &
                           a(4) * ( var(:,:,kl+3:kr+3) + &
                                    var(:,:,kl-3:kr-3) ) ) * idel2
    end if
  end select direction
end subroutine SBP_deriv2_gf_6_3
