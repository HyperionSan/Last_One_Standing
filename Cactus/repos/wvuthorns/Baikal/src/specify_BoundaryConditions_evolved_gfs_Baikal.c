#include "stdio.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
/*
 * Set `none` boundary conditions on BSSN evolved variables, as these are set via NewRad.
 * 
 * Since we choose NewRad boundary conditions, we must register all
 *   evolved gridfunctions to have boundary type "none". This is because
 *   NewRad directly modifies the RHSs.
 * 
 * This code is based on Kranc's McLachlan/ML_BSSN/src/Boundaries.cc code.
 */
void specify_BoundaryConditions_evolved_gfs_Baikal(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS_specify_BoundaryConditions_evolved_gfs_Baikal;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD00GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::aDD00GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD01GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::aDD01GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD02GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::aDD02GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD11GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::aDD11GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD12GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::aDD12GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD22GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::aDD22GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::alphaGF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::alphaGF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::betU0GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::betU0GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::betU1GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::betU1GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::betU2GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::betU2GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::cfGF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::cfGF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD00GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::hDD00GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD01GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::hDD01GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD02GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::hDD02GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD11GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::hDD11GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD12GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::hDD12GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD22GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::hDD22GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::lambdaU0GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::lambdaU0GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::lambdaU1GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::lambdaU1GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::lambdaU2GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::lambdaU2GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::trKGF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::trKGF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::vetU0GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::vetU0GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::vetU1GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::vetU1GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::vetU2GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::vetU2GF!");
}
