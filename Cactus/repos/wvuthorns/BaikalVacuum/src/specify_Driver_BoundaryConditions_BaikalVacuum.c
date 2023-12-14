#include "stdio.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Faces.h"
#include "util_Table.h"
/*
 * Set all boundary conditions for PreSync.
 * BSSN constraints are set to use `flat` boundary conditions, similar to what Lean does.
 * 
 * BSSN evolved variables are set to use `none` boundary conditions, as these are set via NewRad.
 * 
 * Since we choose NewRad boundary conditions, we must register all
 * evolved gridfunctions to have boundary type "none". This is because
 * NewRad directly modifies the RHSs.
 * 
 * This code is based on Kranc's McLachlan/ML_BSSN/src/Boundaries.cc code.
 */
void specify_Driver_BoundaryConditions_BaikalVacuum(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS_specify_Driver_BoundaryConditions_BaikalVacuum;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  const CCTK_INT bndsize = FD_order / 2 + 1;  // <- bndsize = number of ghostzones

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "BaikalVacuum::HGF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::HGF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "BaikalVacuum::MU0GF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::MU0GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "BaikalVacuum::MU1GF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::MU1GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "BaikalVacuum::MU2GF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::MU2GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::aDD00GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::aDD00GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::aDD01GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::aDD01GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::aDD02GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::aDD02GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::aDD11GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::aDD11GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::aDD12GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::aDD12GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::aDD22GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::aDD22GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::alphaGF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::alphaGF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::betU0GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::betU0GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::betU1GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::betU1GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::betU2GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::betU2GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::cfGF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::cfGF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::hDD00GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::hDD00GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::hDD01GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::hDD01GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::hDD02GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::hDD02GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::hDD11GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::hDD11GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::hDD12GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::hDD12GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::hDD22GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::hDD22GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::lambdaU0GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::lambdaU0GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::lambdaU1GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::lambdaU1GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::lambdaU2GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::lambdaU2GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::trKGF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::trKGF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::vetU0GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::vetU0GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::vetU1GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::vetU1GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::vetU2GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for BaikalVacuum::vetU2GF!");
}
