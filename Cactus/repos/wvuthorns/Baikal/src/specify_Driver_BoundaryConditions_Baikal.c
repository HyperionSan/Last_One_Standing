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
void specify_Driver_BoundaryConditions_Baikal(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS_specify_Driver_BoundaryConditions_Baikal;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  const CCTK_INT bndsize = FD_order / 2 + 1;  // <- bndsize = number of ghostzones

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "Baikal::HGF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::HGF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "Baikal::MU0GF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::MU0GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "Baikal::MU1GF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::MU1GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "Baikal::MU2GF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::MU2GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD00GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::aDD00GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD01GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::aDD01GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD02GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::aDD02GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD11GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::aDD11GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD12GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::aDD12GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD22GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::aDD22GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::alphaGF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::alphaGF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::betU0GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::betU0GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::betU1GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::betU1GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::betU2GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::betU2GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::cfGF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::cfGF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD00GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::hDD00GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD01GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::hDD01GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD02GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::hDD02GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD11GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::hDD11GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD12GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::hDD12GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD22GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::hDD22GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::lambdaU0GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::lambdaU0GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::lambdaU1GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::lambdaU1GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::lambdaU2GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::lambdaU2GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::trKGF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::trKGF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::vetU0GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::vetU0GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::vetU1GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::vetU1GF!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::vetU2GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for Baikal::vetU2GF!");
}
