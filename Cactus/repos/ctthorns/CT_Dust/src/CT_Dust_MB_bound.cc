/*  File produced by Kranc */

#define KRANC_C

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Kranc.hh"
#include "Differencing.h"
#include "loopcontrol.h"
#include "gsl/gsl_sf_hyperg.h"

namespace CT_Dust {

extern "C" void CT_Dust_MB_bound_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % CT_Dust_MB_bound_calc_every != CT_Dust_MB_bound_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ADMBase::curv","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ADMBase::curv.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ADMBase::metric","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ADMBase::metric.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CT_Dust::CT_eps","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CT_Dust::CT_eps.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CT_Dust::CT_prs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CT_Dust::CT_prs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CT_Dust::CT_rho","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CT_Dust::CT_rho.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CT_Dust::CT_rhodp","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CT_Dust::CT_rhodp.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CT_Dust::CT_u","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CT_Dust::CT_u.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CT_Dust::CT_W","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CT_Dust::CT_W.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_curv","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN::ML_curv.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_gamma","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN::ML_gamma.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_log_confac","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN::ML_log_confac.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_metric","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN::ML_metric.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_trace_curv","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN::ML_trace_curv.");
  return;
}

static void CT_Dust_MB_bound_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Include user-supplied include files */
  /* Initialise finite differencing variables */
  const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;
  const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;
  const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;
  const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;
  const ptrdiff_t cctkLbnd1 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[0];
  const ptrdiff_t cctkLbnd2 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[1];
  const ptrdiff_t cctkLbnd3 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[2];
  const CCTK_REAL t CCTK_ATTRIBUTE_UNUSED = cctk_time;
  const CCTK_REAL cctkOriginSpace1 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(0);
  const CCTK_REAL cctkOriginSpace2 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(1);
  const CCTK_REAL cctkOriginSpace3 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(2);
  const CCTK_REAL dt CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_TIME;
  const CCTK_REAL dx CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(2);
  const CCTK_REAL dxi CCTK_ATTRIBUTE_UNUSED = pow(dx,-1);
  const CCTK_REAL dyi CCTK_ATTRIBUTE_UNUSED = pow(dy,-1);
  const CCTK_REAL dzi CCTK_ATTRIBUTE_UNUSED = pow(dz,-1);
  const CCTK_REAL khalf CCTK_ATTRIBUTE_UNUSED = 0.5;
  const CCTK_REAL kthird CCTK_ATTRIBUTE_UNUSED = 
    0.333333333333333333333333333333;
  const CCTK_REAL ktwothird CCTK_ATTRIBUTE_UNUSED = 
    0.666666666666666666666666666667;
  const CCTK_REAL kfourthird CCTK_ATTRIBUTE_UNUSED = 
    1.33333333333333333333333333333;
  const CCTK_REAL hdxi CCTK_ATTRIBUTE_UNUSED = 0.5*dxi;
  const CCTK_REAL hdyi CCTK_ATTRIBUTE_UNUSED = 0.5*dyi;
  const CCTK_REAL hdzi CCTK_ATTRIBUTE_UNUSED = 0.5*dzi;
  /* Initialize predefined quantities */
  const CCTK_REAL p1o1024dx CCTK_ATTRIBUTE_UNUSED = 0.0009765625*pow(dx,-1);
  const CCTK_REAL p1o1024dy CCTK_ATTRIBUTE_UNUSED = 0.0009765625*pow(dy,-1);
  const CCTK_REAL p1o1024dz CCTK_ATTRIBUTE_UNUSED = 0.0009765625*pow(dz,-1);
  const CCTK_REAL p1o120dx CCTK_ATTRIBUTE_UNUSED = 0.00833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL p1o120dy CCTK_ATTRIBUTE_UNUSED = 0.00833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL p1o120dz CCTK_ATTRIBUTE_UNUSED = 0.00833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL p1o12dx CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL p1o12dy CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL p1o12dz CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL p1o144dxdy CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o144dxdz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o144dydz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o1680dx CCTK_ATTRIBUTE_UNUSED = 0.000595238095238095238095238095238*pow(dx,-1);
  const CCTK_REAL p1o1680dy CCTK_ATTRIBUTE_UNUSED = 0.000595238095238095238095238095238*pow(dy,-1);
  const CCTK_REAL p1o1680dz CCTK_ATTRIBUTE_UNUSED = 0.000595238095238095238095238095238*pow(dz,-1);
  const CCTK_REAL p1o16dx CCTK_ATTRIBUTE_UNUSED = 0.0625*pow(dx,-1);
  const CCTK_REAL p1o16dy CCTK_ATTRIBUTE_UNUSED = 0.0625*pow(dy,-1);
  const CCTK_REAL p1o16dz CCTK_ATTRIBUTE_UNUSED = 0.0625*pow(dz,-1);
  const CCTK_REAL p1o180dx2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dx,-2);
  const CCTK_REAL p1o180dy2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dy,-2);
  const CCTK_REAL p1o180dz2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dz,-2);
  const CCTK_REAL p1o24dx CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL p1o24dy CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL p1o24dz CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*pow(dz,-1);
  const CCTK_REAL p1o256dx CCTK_ATTRIBUTE_UNUSED = 0.00390625*pow(dx,-1);
  const CCTK_REAL p1o256dy CCTK_ATTRIBUTE_UNUSED = 0.00390625*pow(dy,-1);
  const CCTK_REAL p1o256dz CCTK_ATTRIBUTE_UNUSED = 0.00390625*pow(dz,-1);
  const CCTK_REAL p1o2dx CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dx,-1);
  const CCTK_REAL p1o2dy CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dy,-1);
  const CCTK_REAL p1o2dz CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dz,-1);
  const CCTK_REAL p1o3600dxdy CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o3600dxdz CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o3600dydz CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dx CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1);
  const CCTK_REAL p1o4dxdy CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o4dxdz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dy CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dy,-1);
  const CCTK_REAL p1o4dydz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dz,-1);
  const CCTK_REAL p1o5040dx2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*pow(dx,-2);
  const CCTK_REAL p1o5040dy2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*pow(dy,-2);
  const CCTK_REAL p1o5040dz2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*pow(dz,-2);
  const CCTK_REAL p1o560dx CCTK_ATTRIBUTE_UNUSED = 0.00178571428571428571428571428571*pow(dx,-1);
  const CCTK_REAL p1o560dy CCTK_ATTRIBUTE_UNUSED = 0.00178571428571428571428571428571*pow(dy,-1);
  const CCTK_REAL p1o560dz CCTK_ATTRIBUTE_UNUSED = 0.00178571428571428571428571428571*pow(dz,-1);
  const CCTK_REAL p1o60dx CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL p1o60dy CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL p1o60dz CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dz,-1);
  const CCTK_REAL p1o64dx CCTK_ATTRIBUTE_UNUSED = 0.015625*pow(dx,-1);
  const CCTK_REAL p1o64dy CCTK_ATTRIBUTE_UNUSED = 0.015625*pow(dy,-1);
  const CCTK_REAL p1o64dz CCTK_ATTRIBUTE_UNUSED = 0.015625*pow(dz,-1);
  const CCTK_REAL p1o705600dxdy CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o705600dxdz CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o705600dydz CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o840dx CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*pow(dx,-1);
  const CCTK_REAL p1o840dy CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*pow(dy,-1);
  const CCTK_REAL p1o840dz CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*pow(dz,-1);
  const CCTK_REAL p1odx CCTK_ATTRIBUTE_UNUSED = pow(dx,-1);
  const CCTK_REAL p1odx2 CCTK_ATTRIBUTE_UNUSED = pow(dx,-2);
  const CCTK_REAL p1ody CCTK_ATTRIBUTE_UNUSED = pow(dy,-1);
  const CCTK_REAL p1ody2 CCTK_ATTRIBUTE_UNUSED = pow(dy,-2);
  const CCTK_REAL p1odz CCTK_ATTRIBUTE_UNUSED = pow(dz,-1);
  const CCTK_REAL p1odz2 CCTK_ATTRIBUTE_UNUSED = pow(dz,-2);
  const CCTK_REAL pm1o120dx CCTK_ATTRIBUTE_UNUSED = -0.00833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL pm1o120dy CCTK_ATTRIBUTE_UNUSED = -0.00833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL pm1o120dz CCTK_ATTRIBUTE_UNUSED = -0.00833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dx,-2);
  const CCTK_REAL pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dy,-2);
  const CCTK_REAL pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dz,-2);
  const CCTK_REAL pm1o2dx CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dx,-1);
  const CCTK_REAL pm1o2dy CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dy,-1);
  const CCTK_REAL pm1o2dz CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dz,-1);
  const CCTK_REAL pm1o4dx CCTK_ATTRIBUTE_UNUSED = -0.25*pow(dx,-1);
  const CCTK_REAL pm1o4dy CCTK_ATTRIBUTE_UNUSED = -0.25*pow(dy,-1);
  const CCTK_REAL pm1o4dz CCTK_ATTRIBUTE_UNUSED = -0.25*pow(dz,-1);
  const CCTK_REAL pm1o60dx CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL pm1o60dy CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL pm1o60dz CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dz,-1);
  const CCTK_REAL pm1o840dx CCTK_ATTRIBUTE_UNUSED = -0.00119047619047619047619047619048*pow(dx,-1);
  const CCTK_REAL pm1o840dy CCTK_ATTRIBUTE_UNUSED = -0.00119047619047619047619047619048*pow(dy,-1);
  const CCTK_REAL pm1o840dz CCTK_ATTRIBUTE_UNUSED = -0.00119047619047619047619047619048*pow(dz,-1);
  /* Jacobian variable pointers */
  const bool usejacobian1 = (!CCTK_IsFunctionAliased("MultiPatch_GetMap") || MultiPatch_GetMap(cctkGH) != jacobian_identity_map)
                        && strlen(jacobian_group) > 0;
  const bool usejacobian = assume_use_jacobian>=0 ? assume_use_jacobian : usejacobian1;
  if (usejacobian && (strlen(jacobian_derivative_group) == 0))
  {
    CCTK_WARN(CCTK_WARN_ALERT, "GenericFD::jacobian_group and GenericFD::jacobian_derivative_group must both be set to valid group names");
  }
  
  const CCTK_REAL* restrict jacobian_ptrs[9];
  if (usejacobian) GroupDataPointers(cctkGH, jacobian_group,
                                                9, jacobian_ptrs);
  
  const CCTK_REAL* restrict const J11 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[0] : 0;
  const CCTK_REAL* restrict const J12 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[1] : 0;
  const CCTK_REAL* restrict const J13 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[2] : 0;
  const CCTK_REAL* restrict const J21 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[3] : 0;
  const CCTK_REAL* restrict const J22 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[4] : 0;
  const CCTK_REAL* restrict const J23 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[5] : 0;
  const CCTK_REAL* restrict const J31 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[6] : 0;
  const CCTK_REAL* restrict const J32 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[7] : 0;
  const CCTK_REAL* restrict const J33 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[8] : 0;
  
  const CCTK_REAL* restrict jacobian_determinant_ptrs[1] CCTK_ATTRIBUTE_UNUSED;
  if (usejacobian && strlen(jacobian_determinant_group) > 0) GroupDataPointers(cctkGH, jacobian_determinant_group,
                                                1, jacobian_determinant_ptrs);
  
  const CCTK_REAL* restrict const detJ CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_determinant_ptrs[0] : 0;
  
  const CCTK_REAL* restrict jacobian_inverse_ptrs[9] CCTK_ATTRIBUTE_UNUSED;
  if (usejacobian && strlen(jacobian_inverse_group) > 0) GroupDataPointers(cctkGH, jacobian_inverse_group,
                                                9, jacobian_inverse_ptrs);
  
  const CCTK_REAL* restrict const iJ11 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[0] : 0;
  const CCTK_REAL* restrict const iJ12 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[1] : 0;
  const CCTK_REAL* restrict const iJ13 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[2] : 0;
  const CCTK_REAL* restrict const iJ21 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[3] : 0;
  const CCTK_REAL* restrict const iJ22 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[4] : 0;
  const CCTK_REAL* restrict const iJ23 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[5] : 0;
  const CCTK_REAL* restrict const iJ31 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[6] : 0;
  const CCTK_REAL* restrict const iJ32 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[7] : 0;
  const CCTK_REAL* restrict const iJ33 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[8] : 0;
  
  const CCTK_REAL* restrict jacobian_derivative_ptrs[18] CCTK_ATTRIBUTE_UNUSED;
  if (usejacobian) GroupDataPointers(cctkGH, jacobian_derivative_group,
                                      18, jacobian_derivative_ptrs);
  
  const CCTK_REAL* restrict const dJ111 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[0] : 0;
  const CCTK_REAL* restrict const dJ112 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[1] : 0;
  const CCTK_REAL* restrict const dJ113 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[2] : 0;
  const CCTK_REAL* restrict const dJ122 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[3] : 0;
  const CCTK_REAL* restrict const dJ123 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[4] : 0;
  const CCTK_REAL* restrict const dJ133 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[5] : 0;
  const CCTK_REAL* restrict const dJ211 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[6] : 0;
  const CCTK_REAL* restrict const dJ212 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[7] : 0;
  const CCTK_REAL* restrict const dJ213 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[8] : 0;
  const CCTK_REAL* restrict const dJ222 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[9] : 0;
  const CCTK_REAL* restrict const dJ223 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[10] : 0;
  const CCTK_REAL* restrict const dJ233 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[11] : 0;
  const CCTK_REAL* restrict const dJ311 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[12] : 0;
  const CCTK_REAL* restrict const dJ312 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[13] : 0;
  const CCTK_REAL* restrict const dJ313 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[14] : 0;
  const CCTK_REAL* restrict const dJ322 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[15] : 0;
  const CCTK_REAL* restrict const dJ323 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[16] : 0;
  const CCTK_REAL* restrict const dJ333 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[17] : 0;
  /* Assign local copies of arrays functions */
  
  
  /* Calculate temporaries and arrays functions */
  /* Copy local copies back to grid functions */
  /* Loop over the grid points */
  const int imin0=imin[0];
  const int imin1=imin[1];
  const int imin2=imin[2];
  const int imax0=imax[0];
  const int imax1=imax[1];
  const int imax2=imax[2];
  #pragma omp parallel
  CCTK_LOOP3(CT_Dust_MB_bound,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL alpL CCTK_ATTRIBUTE_UNUSED = alp[index];
    CCTK_REAL betaxL CCTK_ATTRIBUTE_UNUSED = betax[index];
    CCTK_REAL betayL CCTK_ATTRIBUTE_UNUSED = betay[index];
    CCTK_REAL betazL CCTK_ATTRIBUTE_UNUSED = betaz[index];
    CCTK_REAL gxxL CCTK_ATTRIBUTE_UNUSED = gxx[index];
    CCTK_REAL gxyL CCTK_ATTRIBUTE_UNUSED = gxy[index];
    CCTK_REAL gxzL CCTK_ATTRIBUTE_UNUSED = gxz[index];
    CCTK_REAL gyyL CCTK_ATTRIBUTE_UNUSED = gyy[index];
    CCTK_REAL gyzL CCTK_ATTRIBUTE_UNUSED = gyz[index];
    CCTK_REAL gzzL CCTK_ATTRIBUTE_UNUSED = gzz[index];
    CCTK_REAL kxxL CCTK_ATTRIBUTE_UNUSED = kxx[index];
    CCTK_REAL kxyL CCTK_ATTRIBUTE_UNUSED = kxy[index];
    CCTK_REAL kxzL CCTK_ATTRIBUTE_UNUSED = kxz[index];
    CCTK_REAL kyyL CCTK_ATTRIBUTE_UNUSED = kyy[index];
    CCTK_REAL kyzL CCTK_ATTRIBUTE_UNUSED = kyz[index];
    CCTK_REAL kzzL CCTK_ATTRIBUTE_UNUSED = kzz[index];
    CCTK_REAL phiL CCTK_ATTRIBUTE_UNUSED = phi[index];
    CCTK_REAL rhoL CCTK_ATTRIBUTE_UNUSED = rho[index];
    CCTK_REAL trKL CCTK_ATTRIBUTE_UNUSED = trK[index];
    CCTK_REAL u1L CCTK_ATTRIBUTE_UNUSED = u1[index];
    CCTK_REAL u2L CCTK_ATTRIBUTE_UNUSED = u2[index];
    CCTK_REAL u3L CCTK_ATTRIBUTE_UNUSED = u3[index];
    CCTK_REAL xL CCTK_ATTRIBUTE_UNUSED = x[index];
    CCTK_REAL yL CCTK_ATTRIBUTE_UNUSED = y[index];
    CCTK_REAL zL CCTK_ATTRIBUTE_UNUSED = z[index];
    
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    
    switch (fdOrder)
    {
      case 2:
      {
        break;
      }
      
      case 4:
      {
        break;
      }
      
      case 6:
      {
        break;
      }
      
      case 8:
      {
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    gxxL = pow(-1 + 
      3*pow(pflrw_H0,2)*pow(Lambda,-1),0.666666666666666666666666666667)*pow(sinh(0.86602540378443864676372317075*Lambda*t*pow(Lambda*pow(pflrw_H0,-2),-0.5)),1.33333333333333333333333333333);
    
    gxyL = 0;
    
    gxzL = 0;
    
    gyyL = pow(-1 + 
      3*pow(pflrw_H0,2)*pow(Lambda,-1),0.666666666666666666666666666667)*pow(sinh(0.86602540378443864676372317075*Lambda*t*pow(Lambda*pow(pflrw_H0,-2),-0.5)),1.33333333333333333333333333333);
    
    gyzL = 0;
    
    gzzL = 0.000277777777777777777777777777778*pow(-1 + 
      3*pow(pflrw_H0,2)*pow(Lambda,-1),0.666666666666666666666666666667)*pow(-60 
      + (15*pow(pflrw_H0,2)*(pow(xL,2) + 
      pow(yL,2))*pow(Lambda*pow(pflrw_H0,-6)*pow(Lambda - 
      3*pow(pflrw_H0,2),2),0.333333333333333333333333333333) + 
      36*gsl_sf_hyperg_2F1(0.333333333333333333333333333333,0.833333333333333333333333333333,1.83333333333333333333333333333,pow(tanh(0.86602540378443864676372317075*t*pow(Lambda,0.5)),2))*pow(tanh(0.86602540378443864676372317075*t*pow(Lambda,0.5)),0.666666666666666666666666666667))*(-1 
      + 
      sin(zL)),2)*pow(sinh(0.86602540378443864676372317075*Lambda*t*pow(Lambda*pow(pflrw_H0,-2),-0.5)),1.33333333333333333333333333333);
    
    CCTK_REAL detg CCTK_ATTRIBUTE_UNUSED = 2*gxyL*gxzL*gyzL - 
      gzzL*pow(gxyL,2) + gyyL*(gxxL*gzzL - pow(gxzL,2)) - gxxL*pow(gyzL,2);
    
    kxxL = 
      -0.577350269189625764509148780502*Lambda*cosh(0.86602540378443864676372317075*Lambda*t*pow(Lambda*pow(pflrw_H0,-2),-0.5))*pow(Lambda*pow(pflrw_H0,-2),-0.5)*pow(-1 
      + 
      3*pow(pflrw_H0,2)*pow(Lambda,-1),0.666666666666666666666666666667)*pow(sinh(0.86602540378443864676372317075*Lambda*t*pow(Lambda*pow(pflrw_H0,-2),-0.5)),0.333333333333333333333333333333);
    
    kxyL = 0;
    
    kxzL = 0;
    
    kyyL = 
      -0.577350269189625764509148780502*Lambda*cosh(0.86602540378443864676372317075*Lambda*t*pow(Lambda*pow(pflrw_H0,-2),-0.5))*pow(Lambda*pow(pflrw_H0,-2),-0.5)*pow(-1 
      + 
      3*pow(pflrw_H0,2)*pow(Lambda,-1),0.666666666666666666666666666667)*pow(sinh(0.86602540378443864676372317075*Lambda*t*pow(Lambda*pow(pflrw_H0,-2),-0.5)),0.333333333333333333333333333333);
    
    kyzL = 0;
    
    kzzL = -0.00555555555555555555555555555556*pow(-1 + 
      3*pow(pflrw_H0,2)*pow(Lambda,-1),0.666666666666666666666666666667)*pow(sinh(0.86602540378443864676372317075*Lambda*t*pow(Lambda*pow(pflrw_H0,-2),-0.5)),0.333333333333333333333333333333)*(60 
      + (-15*pow(pflrw_H0,2)*(pow(xL,2) + 
      pow(yL,2))*pow(Lambda*pow(pflrw_H0,-6)*pow(Lambda - 
      3*pow(pflrw_H0,2),2),0.333333333333333333333333333333) - 
      36*gsl_sf_hyperg_2F1(0.333333333333333333333333333333,0.833333333333333333333333333333,1.83333333333333333333333333333,pow(tanh(0.86602540378443864676372317075*t*pow(Lambda,0.5)),2))*pow(tanh(0.86602540378443864676372317075*t*pow(Lambda,0.5)),0.666666666666666666666666666667))*(-1 
      + 
      sin(zL)))*(1.73205080756887729352744634151*Lambda*cosh(0.86602540378443864676372317075*Lambda*t*pow(Lambda*pow(pflrw_H0,-2),-0.5))*pow(Lambda*pow(pflrw_H0,-2),-0.5)*(1 
      + (-0.25*pow(pflrw_H0,2)*(pow(xL,2) + 
      pow(yL,2))*pow(Lambda*pow(pflrw_H0,-6)*pow(Lambda - 
      3*pow(pflrw_H0,2),2),0.333333333333333333333333333333) - 
      0.6*gsl_sf_hyperg_2F1(0.333333333333333333333333333333,0.833333333333333333333333333333,1.83333333333333333333333333333,pow(tanh(0.86602540378443864676372317075*t*pow(Lambda,0.5)),2))*pow(tanh(0.86602540378443864676372317075*t*pow(Lambda,0.5)),0.666666666666666666666666666667))*(-1 
      + sin(zL))) + 0.519615242270663188058233902452*pow(Lambda,0.5)*(-5 + 
      3*gsl_sf_hyperg_2F1(0.333333333333333333333333333333,0.833333333333333333333333333333,1.83333333333333333333333333333,pow(tanh(0.86602540378443864676372317075*t*pow(Lambda,0.5)),2))*pow(pow(pow(cosh(0.86602540378443864676372317075*t*pow(Lambda,0.5)),-1),2),0.333333333333333333333333333333))*pow(pow(pow(cosh(0.86602540378443864676372317075*t*pow(Lambda,0.5)),-1),2),0.666666666666666666666666666667)*pow(tanh(0.86602540378443864676372317075*t*pow(Lambda,0.5)),-0.333333333333333333333333333333)*(-1 
      + 
      sin(zL))*sinh(0.86602540378443864676372317075*Lambda*t*pow(Lambda*pow(pflrw_H0,-2),-0.5)));
    
    CCTK_REAL det4g CCTK_ATTRIBUTE_UNUSED = -(alpL*detg);
    
    CCTK_REAL gu11 CCTK_ATTRIBUTE_UNUSED = (gyyL*gzzL - 
      pow(gyzL,2))*pow(detg,-1);
    
    CCTK_REAL gu12 CCTK_ATTRIBUTE_UNUSED = (gxzL*gyzL - 
      gxyL*gzzL)*pow(detg,-1);
    
    CCTK_REAL gu13 CCTK_ATTRIBUTE_UNUSED = (-(gxzL*gyyL) + 
      gxyL*gyzL)*pow(detg,-1);
    
    CCTK_REAL gu22 CCTK_ATTRIBUTE_UNUSED = (gxxL*gzzL - 
      pow(gxzL,2))*pow(detg,-1);
    
    CCTK_REAL gu23 CCTK_ATTRIBUTE_UNUSED = (gxyL*gxzL - 
      gxxL*gyzL)*pow(detg,-1);
    
    CCTK_REAL gu33 CCTK_ATTRIBUTE_UNUSED = (gxxL*gyyL - 
      pow(gxyL,2))*pow(detg,-1);
    
    CCTK_REAL f1 CCTK_ATTRIBUTE_UNUSED = (betaxL*u1L + betayL*u2L + 
      betazL*u3L)*pow(alpL,-2);
    
    CCTK_REAL f2 CCTK_ATTRIBUTE_UNUSED = pow(alpL,-2)*(1 + 
      2*(u1L*(u2L*gu12 + u3L*gu13) + u2L*u3L*gu23) + gu11*pow(u1L,2) + 
      gu22*pow(u2L,2) + gu33*pow(u3L,2));
    
    CCTK_REAL u0 CCTK_ATTRIBUTE_UNUSED = 0.5*(f1 + pow(4*f2 + 
      pow(f1,2),0.5));
    
    CCTK_REAL WL CCTK_ATTRIBUTE_UNUSED = u0*pow(-det4g,0.5);
    
    rhoL = 
      0.198943678864869169711104704216*Lambda*pow(pow(sinh(0.86602540378443864676372317075*Lambda*t*pow(Lambda*pow(pflrw_H0,-2),-0.5)),-1),2)*pow(60 
      + (-15*pow(pflrw_H0,2)*(pow(xL,2) + 
      pow(yL,2))*pow(Lambda*pow(pflrw_H0,-6)*pow(Lambda - 
      3*pow(pflrw_H0,2),2),0.333333333333333333333333333333) - 
      36*gsl_sf_hyperg_2F1(0.333333333333333333333333333333,0.833333333333333333333333333333,1.83333333333333333333333333333,pow(tanh(0.86602540378443864676372317075*t*pow(Lambda,0.5)),2))*pow(tanh(0.86602540378443864676372317075*t*pow(Lambda,0.5)),0.666666666666666666666666666667))*(-1 
      + sin(zL)),-1)*(12 - 3*pow(pflrw_H0,2)*(pow(xL,2) + 
      pow(yL,2))*pow(Lambda*pow(pflrw_H0,-6)*pow(Lambda - 
      3*pow(pflrw_H0,2),2),0.333333333333333333333333333333)*(-1 + sin(zL)));
    
    CCTK_REAL epsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL prsL CCTK_ATTRIBUTE_UNUSED = rhoL*w;
    
    u1L = 0;
    
    u2L = 0;
    
    u3L = 0;
    
    CCTK_REAL rhodpL CCTK_ATTRIBUTE_UNUSED = rhoL*pow(detg,0.5);
    
    phiL = pow(detg,-0.166666666666666666666666666667);
    
    CCTK_REAL gt11L CCTK_ATTRIBUTE_UNUSED = gxxL*pow(phiL,2);
    
    CCTK_REAL gt12L CCTK_ATTRIBUTE_UNUSED = gxyL*pow(phiL,2);
    
    CCTK_REAL gt13L CCTK_ATTRIBUTE_UNUSED = gxzL*pow(phiL,2);
    
    CCTK_REAL gt22L CCTK_ATTRIBUTE_UNUSED = gyyL*pow(phiL,2);
    
    CCTK_REAL gt23L CCTK_ATTRIBUTE_UNUSED = gyzL*pow(phiL,2);
    
    CCTK_REAL gt33L CCTK_ATTRIBUTE_UNUSED = gzzL*pow(phiL,2);
    
    trKL = kxxL*gu11 + kyyL*gu22 + 2*(kxyL*gu12 + kxzL*gu13 + kyzL*gu23) + 
      kzzL*gu33;
    
    CCTK_REAL At11L CCTK_ATTRIBUTE_UNUSED = (kxxL - 
      0.333333333333333333333333333333*gxxL*trKL)*pow(phiL,2);
    
    CCTK_REAL At12L CCTK_ATTRIBUTE_UNUSED = (kxyL - 
      0.333333333333333333333333333333*gxyL*trKL)*pow(phiL,2);
    
    CCTK_REAL At13L CCTK_ATTRIBUTE_UNUSED = (kxzL - 
      0.333333333333333333333333333333*gxzL*trKL)*pow(phiL,2);
    
    CCTK_REAL At22L CCTK_ATTRIBUTE_UNUSED = (kyyL - 
      0.333333333333333333333333333333*gyyL*trKL)*pow(phiL,2);
    
    CCTK_REAL At23L CCTK_ATTRIBUTE_UNUSED = (kyzL - 
      0.333333333333333333333333333333*gyzL*trKL)*pow(phiL,2);
    
    CCTK_REAL At33L CCTK_ATTRIBUTE_UNUSED = (kzzL - 
      0.333333333333333333333333333333*gzzL*trKL)*pow(phiL,2);
    
    CCTK_REAL Xt1L CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Xt2L CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Xt3L CCTK_ATTRIBUTE_UNUSED = 0;
    /* Copy local copies back to grid functions */
    At11[index] = At11L;
    At12[index] = At12L;
    At13[index] = At13L;
    At22[index] = At22L;
    At23[index] = At23L;
    At33[index] = At33L;
    eps[index] = epsL;
    gt11[index] = gt11L;
    gt12[index] = gt12L;
    gt13[index] = gt13L;
    gt22[index] = gt22L;
    gt23[index] = gt23L;
    gt33[index] = gt33L;
    gxx[index] = gxxL;
    gxy[index] = gxyL;
    gxz[index] = gxzL;
    gyy[index] = gyyL;
    gyz[index] = gyzL;
    gzz[index] = gzzL;
    kxx[index] = kxxL;
    kxy[index] = kxyL;
    kxz[index] = kxzL;
    kyy[index] = kyyL;
    kyz[index] = kyzL;
    kzz[index] = kzzL;
    phi[index] = phiL;
    prs[index] = prsL;
    rho[index] = rhoL;
    rhodp[index] = rhodpL;
    trK[index] = trKL;
    u1[index] = u1L;
    u2[index] = u2L;
    u3[index] = u3L;
    W[index] = WL;
    Xt1[index] = Xt1L;
    Xt2[index] = Xt2L;
    Xt3[index] = Xt3L;
  }
  CCTK_ENDLOOP3(CT_Dust_MB_bound);
}
extern "C" void CT_Dust_MB_bound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering CT_Dust_MB_bound_Body");
  }
  if (cctk_iteration % CT_Dust_MB_bound_calc_every != CT_Dust_MB_bound_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ADMBase::curv",
    "ADMBase::lapse",
    "ADMBase::metric",
    "ADMBase::shift",
    "CT_Dust::CT_eps",
    "CT_Dust::CT_prs",
    "CT_Dust::CT_rho",
    "CT_Dust::CT_rhodp",
    "CT_Dust::CT_u",
    "CT_Dust::CT_W",
    "grid::coordinates",
    "ML_BSSN::ML_curv",
    "ML_BSSN::ML_gamma",
    "ML_BSSN::ML_log_confac",
    "ML_BSSN::ML_metric",
    "ML_BSSN::ML_trace_curv"};
  AssertGroupStorage(cctkGH, "CT_Dust_MB_bound", 16, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      break;
    }
    
    case 4:
    {
      break;
    }
    
    case 6:
    {
      break;
    }
    
    case 8:
    {
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverBoundaryWithGhosts(cctkGH, CT_Dust_MB_bound_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving CT_Dust_MB_bound_Body");
  }
}

} // namespace CT_Dust
