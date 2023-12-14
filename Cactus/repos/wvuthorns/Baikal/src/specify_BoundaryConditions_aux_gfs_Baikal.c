#include "stdio.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Faces.h"
#include "util_Table.h"
/*
 * Set `flat` boundary conditions on BSSN constraints, similar to what Lean does.
 */
void specify_BoundaryConditions_aux_gfs_Baikal(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS_specify_BoundaryConditions_aux_gfs_Baikal;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;

  const CCTK_INT bndsize = FD_order / 2 + 1;  // <- bndsize = number of ghostzones

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "Baikal::HGF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::HGF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "Baikal::MU0GF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::MU0GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "Baikal::MU1GF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::MU1GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "Baikal::MU2GF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::MU2GF!");
}
