#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine TestLoopFortran_all(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_LOOP3_ALL_DECLARE(all3)
  integer   :: i,j,k
  CCTK_REAL :: fsum
  
  call CCTK_INFO("TestLoopFortran_all")
  fsum = 0
  !$OMP PARALLEL &
  !$OMP CCTK_LOOP3_ALL_OMP_PRIVATE(all3, i,j,k) &
  !$OMP reduction(+: fsum) &
  !$OMP default(none) shared(cctkgh, cctk_ash, cctk_lsh, cctk_tile_min, cctk_tile_max, r)
  CCTK_LOOP3_ALL(all3, i,j,k)
     fsum = fsum + r(i,j,k)
  CCTK_ENDLOOP3_ALL(all3)
  !$OMP END PARALLEL
  fsum_all = fsum
end subroutine TestLoopFortran_all



subroutine TestLoopFortran_int(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_LOOP3_INT_DECLARE(int3)
  integer   :: i,j,k
  CCTK_REAL :: fsum
  
  call CCTK_INFO("TestLoopFortran_int")
  fsum = 0
  !$OMP PARALLEL reduction(+: fsum) CCTK_LOOP3_INT_OMP_PRIVATE(int3, i,j,k) &
  !$OMP default(none) shared(cctkgh, cctk_ash, cctk_lsh, cctk_tile_min, cctk_tile_max, r)
  CCTK_LOOP3_INT(int3, i,j,k)
     fsum = fsum + r(i,j,k)
  CCTK_ENDLOOP3_INT(int3)
  !$OMP END PARALLEL
  fsum_int = fsum
end subroutine TestLoopFortran_int



subroutine TestLoopFortran_bnd(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_LOOP3_BND_DECLARE(bnd3)
  integer   :: i,j,k
  integer   :: n_i,n_j,n_k
  CCTK_REAL :: fsum
  
  call CCTK_INFO("TestLoopFortran_bnd")
  fsum = 0
  !$omp parallel CCTK_LOOP3_BND_OMP_PRIVATE(bnd3, i,j,k, n_i,n_j,n_k) &
  !$OMP reduction(+: fsum) &
  !$OMP default(none) shared(cctkgh, cctk_ash, cctk_lsh, cctk_tile_min, cctk_tile_max, r)
  CCTK_LOOP3_BND(bnd3, i,j,k, n_i,n_j,n_k)
     fsum = fsum + r(i,j,k)
  CCTK_ENDLOOP3_BND(bnd3)
  !$omp end parallel
  fsum_bnd = fsum
end subroutine TestLoopFortran_bnd



subroutine TestLoopFortran_intbnd(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_LOOP3_INTBND_DECLARE(intbnd3)
  integer   :: i,j,k
  integer   :: n_i,n_j,n_k
  CCTK_REAL :: fsum
  
  call CCTK_INFO("TestLoopFortran_intbnd")
  fsum = 0
  !$OMP PARALLEL &
  !$OMP default(none) shared(cctkgh, cctk_ash, cctk_lsh, cctk_tile_min, cctk_tile_max, r) &
  !$OMP reduction(+: fsum) &
  !$omp CCTK_LOOP3_INTBND_OMP_PRIVATE(intbnd3, i,j,k, n_i,n_j,n_k)
  CCTK_LOOP3_INTBND(intbnd3, i,j,k, n_i,n_j,n_k)
     fsum = fsum + r(i,j,k)
  CCTK_ENDLOOP3_INTBND(intbnd3)
  !$OMP END PARALLEL
  fsum_intbnd = fsum
end subroutine TestLoopFortran_intbnd



subroutine TestLoopFortran(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  call TestLoopFortran_all(CCTK_PASS_FTOF)
  call TestLoopFortran_int(CCTK_PASS_FTOF)
  call TestLoopFortran_bnd(CCTK_PASS_FTOF)
  call TestLoopFortran_intbnd(CCTK_PASS_FTOF)
end subroutine TestLoopFortran
