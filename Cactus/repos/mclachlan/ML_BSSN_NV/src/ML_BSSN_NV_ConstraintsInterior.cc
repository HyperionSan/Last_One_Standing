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

namespace ML_BSSN_NV {

extern "C" void ML_BSSN_NV_ConstraintsInterior_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_BSSN_NV_ConstraintsInterior_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_NV_ConstraintsInterior_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_NV_ConstraintsInterior_calc_every != ML_BSSN_NV_ConstraintsInterior_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NV::ML_cons_Gamma","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN_NV::ML_cons_Gamma.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NV::ML_Ham","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN_NV::ML_Ham.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NV::ML_mom","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN_NV::ML_mom.");
  return;
}

static void ML_BSSN_NV_ConstraintsInterior_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  const CCTK_REAL p1o180dx2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dx,-2);
  const CCTK_REAL p1o180dy2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dy,-2);
  const CCTK_REAL p1o180dz2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dz,-2);
  const CCTK_REAL p1o24dx CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL p1o24dy CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL p1o24dz CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*pow(dz,-1);
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
  const CCTK_REAL p1odx2 CCTK_ATTRIBUTE_UNUSED = pow(dx,-2);
  const CCTK_REAL p1ody2 CCTK_ATTRIBUTE_UNUSED = pow(dy,-2);
  const CCTK_REAL p1odz2 CCTK_ATTRIBUTE_UNUSED = pow(dz,-2);
  const CCTK_REAL pm1o120dx CCTK_ATTRIBUTE_UNUSED = -0.00833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL pm1o120dy CCTK_ATTRIBUTE_UNUSED = -0.00833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL pm1o120dz CCTK_ATTRIBUTE_UNUSED = -0.00833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL pm1o12dx CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dx,-2);
  const CCTK_REAL pm1o12dy CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dy,-2);
  const CCTK_REAL pm1o12dz CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dz,-2);
  const CCTK_REAL pm1o16dx CCTK_ATTRIBUTE_UNUSED = -0.0625*pow(dx,-1);
  const CCTK_REAL pm1o16dy CCTK_ATTRIBUTE_UNUSED = -0.0625*pow(dy,-1);
  const CCTK_REAL pm1o16dz CCTK_ATTRIBUTE_UNUSED = -0.0625*pow(dz,-1);
  const CCTK_REAL pm1o256dx CCTK_ATTRIBUTE_UNUSED = -0.00390625*pow(dx,-1);
  const CCTK_REAL pm1o256dy CCTK_ATTRIBUTE_UNUSED = -0.00390625*pow(dy,-1);
  const CCTK_REAL pm1o256dz CCTK_ATTRIBUTE_UNUSED = -0.00390625*pow(dz,-1);
  const CCTK_REAL pm1o2dx CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dx,-1);
  const CCTK_REAL pm1o2dy CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dy,-1);
  const CCTK_REAL pm1o2dz CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dz,-1);
  const CCTK_REAL pm1o4dx CCTK_ATTRIBUTE_UNUSED = -0.25*pow(dx,-1);
  const CCTK_REAL pm1o4dy CCTK_ATTRIBUTE_UNUSED = -0.25*pow(dy,-1);
  const CCTK_REAL pm1o4dz CCTK_ATTRIBUTE_UNUSED = -0.25*pow(dz,-1);
  const CCTK_REAL pm1o60dx CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL pm1o60dy CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL pm1o60dz CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dz,-1);
  const CCTK_REAL pm1o6dx CCTK_ATTRIBUTE_UNUSED = -0.166666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL pm1o6dy CCTK_ATTRIBUTE_UNUSED = -0.166666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL pm1o6dz CCTK_ATTRIBUTE_UNUSED = -0.166666666666666666666666666667*pow(dz,-1);
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
  CCTK_LOOP3(ML_BSSN_NV_ConstraintsInterior,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL alphaL CCTK_ATTRIBUTE_UNUSED = alpha[index];
    CCTK_REAL At11L CCTK_ATTRIBUTE_UNUSED = At11[index];
    CCTK_REAL At12L CCTK_ATTRIBUTE_UNUSED = At12[index];
    CCTK_REAL At13L CCTK_ATTRIBUTE_UNUSED = At13[index];
    CCTK_REAL At22L CCTK_ATTRIBUTE_UNUSED = At22[index];
    CCTK_REAL At23L CCTK_ATTRIBUTE_UNUSED = At23[index];
    CCTK_REAL At33L CCTK_ATTRIBUTE_UNUSED = At33[index];
    CCTK_REAL beta1L CCTK_ATTRIBUTE_UNUSED = beta1[index];
    CCTK_REAL beta2L CCTK_ATTRIBUTE_UNUSED = beta2[index];
    CCTK_REAL beta3L CCTK_ATTRIBUTE_UNUSED = beta3[index];
    CCTK_REAL gt11L CCTK_ATTRIBUTE_UNUSED = gt11[index];
    CCTK_REAL gt12L CCTK_ATTRIBUTE_UNUSED = gt12[index];
    CCTK_REAL gt13L CCTK_ATTRIBUTE_UNUSED = gt13[index];
    CCTK_REAL gt22L CCTK_ATTRIBUTE_UNUSED = gt22[index];
    CCTK_REAL gt23L CCTK_ATTRIBUTE_UNUSED = gt23[index];
    CCTK_REAL gt33L CCTK_ATTRIBUTE_UNUSED = gt33[index];
    CCTK_REAL phiL CCTK_ATTRIBUTE_UNUSED = phi[index];
    CCTK_REAL trKL CCTK_ATTRIBUTE_UNUSED = trK[index];
    CCTK_REAL Xt1L CCTK_ATTRIBUTE_UNUSED = Xt1[index];
    CCTK_REAL Xt2L CCTK_ATTRIBUTE_UNUSED = Xt2[index];
    CCTK_REAL Xt3L CCTK_ATTRIBUTE_UNUSED = Xt3[index];
    
    CCTK_REAL eTttL, eTtxL, eTtyL, eTtzL, eTxxL, eTxyL, eTxzL, eTyyL, eTyzL, eTzzL CCTK_ATTRIBUTE_UNUSED ;
    
    if (assume_stress_energy_state>=0 ? assume_stress_energy_state : *stress_energy_state)
    {
      eTttL = eTtt[index];
      eTtxL = eTtx[index];
      eTtyL = eTty[index];
      eTtzL = eTtz[index];
      eTxxL = eTxx[index];
      eTxyL = eTxy[index];
      eTxzL = eTxz[index];
      eTyyL = eTyy[index];
      eTyzL = eTyz[index];
      eTzzL = eTzz[index];
    }
    else
    {
      eTttL = 0.;
      eTtxL = 0.;
      eTtyL = 0.;
      eTtzL = 0.;
      eTxxL = 0.;
      eTxyL = 0.;
      eTxzL = 0.;
      eTyyL = 0.;
      eTyzL = 0.;
      eTzzL = 0.;
    }
    
    CCTK_REAL dJ111L, dJ112L, dJ113L, dJ122L, dJ123L, dJ133L, dJ211L, dJ212L, dJ213L, dJ222L, dJ223L, dJ233L, dJ311L, dJ312L, dJ313L, dJ322L, dJ323L, dJ333L, J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L CCTK_ATTRIBUTE_UNUSED ;
    
    if (usejacobian)
    {
      dJ111L = dJ111[index];
      dJ112L = dJ112[index];
      dJ113L = dJ113[index];
      dJ122L = dJ122[index];
      dJ123L = dJ123[index];
      dJ133L = dJ133[index];
      dJ211L = dJ211[index];
      dJ212L = dJ212[index];
      dJ213L = dJ213[index];
      dJ222L = dJ222[index];
      dJ223L = dJ223[index];
      dJ233L = dJ233[index];
      dJ311L = dJ311[index];
      dJ312L = dJ312[index];
      dJ313L = dJ313[index];
      dJ322L = dJ322[index];
      dJ323L = dJ323[index];
      dJ333L = dJ333[index];
      J11L = J11[index];
      J12L = J12[index];
      J13L = J13[index];
      J21L = J21[index];
      J22L = J22[index];
      J23L = J23[index];
      J31L = J31[index];
      J32L = J32[index];
      J33L = J33[index];
    }
    /* Include user supplied include files */
    /* Precompute derivatives */
    CCTK_REAL PDstandardNth1At11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2At11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3At11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1At12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2At12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3At12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1At13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2At13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3At13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1At22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2At22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3At22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1At23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2At23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3At23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1At33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2At33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3At33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3Xt3 CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDstandardNth1At11 = PDstandardNthfdOrder21(&At11[index]);
        PDstandardNth2At11 = PDstandardNthfdOrder22(&At11[index]);
        PDstandardNth3At11 = PDstandardNthfdOrder23(&At11[index]);
        PDstandardNth1At12 = PDstandardNthfdOrder21(&At12[index]);
        PDstandardNth2At12 = PDstandardNthfdOrder22(&At12[index]);
        PDstandardNth3At12 = PDstandardNthfdOrder23(&At12[index]);
        PDstandardNth1At13 = PDstandardNthfdOrder21(&At13[index]);
        PDstandardNth2At13 = PDstandardNthfdOrder22(&At13[index]);
        PDstandardNth3At13 = PDstandardNthfdOrder23(&At13[index]);
        PDstandardNth1At22 = PDstandardNthfdOrder21(&At22[index]);
        PDstandardNth2At22 = PDstandardNthfdOrder22(&At22[index]);
        PDstandardNth3At22 = PDstandardNthfdOrder23(&At22[index]);
        PDstandardNth1At23 = PDstandardNthfdOrder21(&At23[index]);
        PDstandardNth2At23 = PDstandardNthfdOrder22(&At23[index]);
        PDstandardNth3At23 = PDstandardNthfdOrder23(&At23[index]);
        PDstandardNth1At33 = PDstandardNthfdOrder21(&At33[index]);
        PDstandardNth2At33 = PDstandardNthfdOrder22(&At33[index]);
        PDstandardNth3At33 = PDstandardNthfdOrder23(&At33[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder21(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder22(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder23(&gt11[index]);
        PDstandardNth11gt11 = PDstandardNthfdOrder211(&gt11[index]);
        PDstandardNth22gt11 = PDstandardNthfdOrder222(&gt11[index]);
        PDstandardNth33gt11 = PDstandardNthfdOrder233(&gt11[index]);
        PDstandardNth12gt11 = PDstandardNthfdOrder212(&gt11[index]);
        PDstandardNth13gt11 = PDstandardNthfdOrder213(&gt11[index]);
        PDstandardNth23gt11 = PDstandardNthfdOrder223(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder21(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder22(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder23(&gt12[index]);
        PDstandardNth11gt12 = PDstandardNthfdOrder211(&gt12[index]);
        PDstandardNth22gt12 = PDstandardNthfdOrder222(&gt12[index]);
        PDstandardNth33gt12 = PDstandardNthfdOrder233(&gt12[index]);
        PDstandardNth12gt12 = PDstandardNthfdOrder212(&gt12[index]);
        PDstandardNth13gt12 = PDstandardNthfdOrder213(&gt12[index]);
        PDstandardNth23gt12 = PDstandardNthfdOrder223(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder21(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder22(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder23(&gt13[index]);
        PDstandardNth11gt13 = PDstandardNthfdOrder211(&gt13[index]);
        PDstandardNth22gt13 = PDstandardNthfdOrder222(&gt13[index]);
        PDstandardNth33gt13 = PDstandardNthfdOrder233(&gt13[index]);
        PDstandardNth12gt13 = PDstandardNthfdOrder212(&gt13[index]);
        PDstandardNth13gt13 = PDstandardNthfdOrder213(&gt13[index]);
        PDstandardNth23gt13 = PDstandardNthfdOrder223(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder21(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder22(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder23(&gt22[index]);
        PDstandardNth11gt22 = PDstandardNthfdOrder211(&gt22[index]);
        PDstandardNth22gt22 = PDstandardNthfdOrder222(&gt22[index]);
        PDstandardNth33gt22 = PDstandardNthfdOrder233(&gt22[index]);
        PDstandardNth12gt22 = PDstandardNthfdOrder212(&gt22[index]);
        PDstandardNth13gt22 = PDstandardNthfdOrder213(&gt22[index]);
        PDstandardNth23gt22 = PDstandardNthfdOrder223(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder21(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder22(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder23(&gt23[index]);
        PDstandardNth11gt23 = PDstandardNthfdOrder211(&gt23[index]);
        PDstandardNth22gt23 = PDstandardNthfdOrder222(&gt23[index]);
        PDstandardNth33gt23 = PDstandardNthfdOrder233(&gt23[index]);
        PDstandardNth12gt23 = PDstandardNthfdOrder212(&gt23[index]);
        PDstandardNth13gt23 = PDstandardNthfdOrder213(&gt23[index]);
        PDstandardNth23gt23 = PDstandardNthfdOrder223(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder21(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder22(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder23(&gt33[index]);
        PDstandardNth11gt33 = PDstandardNthfdOrder211(&gt33[index]);
        PDstandardNth22gt33 = PDstandardNthfdOrder222(&gt33[index]);
        PDstandardNth33gt33 = PDstandardNthfdOrder233(&gt33[index]);
        PDstandardNth12gt33 = PDstandardNthfdOrder212(&gt33[index]);
        PDstandardNth13gt33 = PDstandardNthfdOrder213(&gt33[index]);
        PDstandardNth23gt33 = PDstandardNthfdOrder223(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder21(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder22(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder23(&phi[index]);
        PDstandardNth11phi = PDstandardNthfdOrder211(&phi[index]);
        PDstandardNth22phi = PDstandardNthfdOrder222(&phi[index]);
        PDstandardNth33phi = PDstandardNthfdOrder233(&phi[index]);
        PDstandardNth12phi = PDstandardNthfdOrder212(&phi[index]);
        PDstandardNth13phi = PDstandardNthfdOrder213(&phi[index]);
        PDstandardNth23phi = PDstandardNthfdOrder223(&phi[index]);
        PDstandardNth1trK = PDstandardNthfdOrder21(&trK[index]);
        PDstandardNth2trK = PDstandardNthfdOrder22(&trK[index]);
        PDstandardNth3trK = PDstandardNthfdOrder23(&trK[index]);
        PDstandardNth1Xt1 = PDstandardNthfdOrder21(&Xt1[index]);
        PDstandardNth2Xt1 = PDstandardNthfdOrder22(&Xt1[index]);
        PDstandardNth3Xt1 = PDstandardNthfdOrder23(&Xt1[index]);
        PDstandardNth1Xt2 = PDstandardNthfdOrder21(&Xt2[index]);
        PDstandardNth2Xt2 = PDstandardNthfdOrder22(&Xt2[index]);
        PDstandardNth3Xt2 = PDstandardNthfdOrder23(&Xt2[index]);
        PDstandardNth1Xt3 = PDstandardNthfdOrder21(&Xt3[index]);
        PDstandardNth2Xt3 = PDstandardNthfdOrder22(&Xt3[index]);
        PDstandardNth3Xt3 = PDstandardNthfdOrder23(&Xt3[index]);
        break;
      }
      
      case 4:
      {
        PDstandardNth1At11 = PDstandardNthfdOrder41(&At11[index]);
        PDstandardNth2At11 = PDstandardNthfdOrder42(&At11[index]);
        PDstandardNth3At11 = PDstandardNthfdOrder43(&At11[index]);
        PDstandardNth1At12 = PDstandardNthfdOrder41(&At12[index]);
        PDstandardNth2At12 = PDstandardNthfdOrder42(&At12[index]);
        PDstandardNth3At12 = PDstandardNthfdOrder43(&At12[index]);
        PDstandardNth1At13 = PDstandardNthfdOrder41(&At13[index]);
        PDstandardNth2At13 = PDstandardNthfdOrder42(&At13[index]);
        PDstandardNth3At13 = PDstandardNthfdOrder43(&At13[index]);
        PDstandardNth1At22 = PDstandardNthfdOrder41(&At22[index]);
        PDstandardNth2At22 = PDstandardNthfdOrder42(&At22[index]);
        PDstandardNth3At22 = PDstandardNthfdOrder43(&At22[index]);
        PDstandardNth1At23 = PDstandardNthfdOrder41(&At23[index]);
        PDstandardNth2At23 = PDstandardNthfdOrder42(&At23[index]);
        PDstandardNth3At23 = PDstandardNthfdOrder43(&At23[index]);
        PDstandardNth1At33 = PDstandardNthfdOrder41(&At33[index]);
        PDstandardNth2At33 = PDstandardNthfdOrder42(&At33[index]);
        PDstandardNth3At33 = PDstandardNthfdOrder43(&At33[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder41(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder42(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder43(&gt11[index]);
        PDstandardNth11gt11 = PDstandardNthfdOrder411(&gt11[index]);
        PDstandardNth22gt11 = PDstandardNthfdOrder422(&gt11[index]);
        PDstandardNth33gt11 = PDstandardNthfdOrder433(&gt11[index]);
        PDstandardNth12gt11 = PDstandardNthfdOrder412(&gt11[index]);
        PDstandardNth13gt11 = PDstandardNthfdOrder413(&gt11[index]);
        PDstandardNth23gt11 = PDstandardNthfdOrder423(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder41(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder42(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder43(&gt12[index]);
        PDstandardNth11gt12 = PDstandardNthfdOrder411(&gt12[index]);
        PDstandardNth22gt12 = PDstandardNthfdOrder422(&gt12[index]);
        PDstandardNth33gt12 = PDstandardNthfdOrder433(&gt12[index]);
        PDstandardNth12gt12 = PDstandardNthfdOrder412(&gt12[index]);
        PDstandardNth13gt12 = PDstandardNthfdOrder413(&gt12[index]);
        PDstandardNth23gt12 = PDstandardNthfdOrder423(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder41(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder42(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder43(&gt13[index]);
        PDstandardNth11gt13 = PDstandardNthfdOrder411(&gt13[index]);
        PDstandardNth22gt13 = PDstandardNthfdOrder422(&gt13[index]);
        PDstandardNth33gt13 = PDstandardNthfdOrder433(&gt13[index]);
        PDstandardNth12gt13 = PDstandardNthfdOrder412(&gt13[index]);
        PDstandardNth13gt13 = PDstandardNthfdOrder413(&gt13[index]);
        PDstandardNth23gt13 = PDstandardNthfdOrder423(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder41(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder42(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder43(&gt22[index]);
        PDstandardNth11gt22 = PDstandardNthfdOrder411(&gt22[index]);
        PDstandardNth22gt22 = PDstandardNthfdOrder422(&gt22[index]);
        PDstandardNth33gt22 = PDstandardNthfdOrder433(&gt22[index]);
        PDstandardNth12gt22 = PDstandardNthfdOrder412(&gt22[index]);
        PDstandardNth13gt22 = PDstandardNthfdOrder413(&gt22[index]);
        PDstandardNth23gt22 = PDstandardNthfdOrder423(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder41(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder42(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder43(&gt23[index]);
        PDstandardNth11gt23 = PDstandardNthfdOrder411(&gt23[index]);
        PDstandardNth22gt23 = PDstandardNthfdOrder422(&gt23[index]);
        PDstandardNth33gt23 = PDstandardNthfdOrder433(&gt23[index]);
        PDstandardNth12gt23 = PDstandardNthfdOrder412(&gt23[index]);
        PDstandardNth13gt23 = PDstandardNthfdOrder413(&gt23[index]);
        PDstandardNth23gt23 = PDstandardNthfdOrder423(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder41(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder42(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder43(&gt33[index]);
        PDstandardNth11gt33 = PDstandardNthfdOrder411(&gt33[index]);
        PDstandardNth22gt33 = PDstandardNthfdOrder422(&gt33[index]);
        PDstandardNth33gt33 = PDstandardNthfdOrder433(&gt33[index]);
        PDstandardNth12gt33 = PDstandardNthfdOrder412(&gt33[index]);
        PDstandardNth13gt33 = PDstandardNthfdOrder413(&gt33[index]);
        PDstandardNth23gt33 = PDstandardNthfdOrder423(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder41(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder42(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder43(&phi[index]);
        PDstandardNth11phi = PDstandardNthfdOrder411(&phi[index]);
        PDstandardNth22phi = PDstandardNthfdOrder422(&phi[index]);
        PDstandardNth33phi = PDstandardNthfdOrder433(&phi[index]);
        PDstandardNth12phi = PDstandardNthfdOrder412(&phi[index]);
        PDstandardNth13phi = PDstandardNthfdOrder413(&phi[index]);
        PDstandardNth23phi = PDstandardNthfdOrder423(&phi[index]);
        PDstandardNth1trK = PDstandardNthfdOrder41(&trK[index]);
        PDstandardNth2trK = PDstandardNthfdOrder42(&trK[index]);
        PDstandardNth3trK = PDstandardNthfdOrder43(&trK[index]);
        PDstandardNth1Xt1 = PDstandardNthfdOrder41(&Xt1[index]);
        PDstandardNth2Xt1 = PDstandardNthfdOrder42(&Xt1[index]);
        PDstandardNth3Xt1 = PDstandardNthfdOrder43(&Xt1[index]);
        PDstandardNth1Xt2 = PDstandardNthfdOrder41(&Xt2[index]);
        PDstandardNth2Xt2 = PDstandardNthfdOrder42(&Xt2[index]);
        PDstandardNth3Xt2 = PDstandardNthfdOrder43(&Xt2[index]);
        PDstandardNth1Xt3 = PDstandardNthfdOrder41(&Xt3[index]);
        PDstandardNth2Xt3 = PDstandardNthfdOrder42(&Xt3[index]);
        PDstandardNth3Xt3 = PDstandardNthfdOrder43(&Xt3[index]);
        break;
      }
      
      case 6:
      {
        PDstandardNth1At11 = PDstandardNthfdOrder61(&At11[index]);
        PDstandardNth2At11 = PDstandardNthfdOrder62(&At11[index]);
        PDstandardNth3At11 = PDstandardNthfdOrder63(&At11[index]);
        PDstandardNth1At12 = PDstandardNthfdOrder61(&At12[index]);
        PDstandardNth2At12 = PDstandardNthfdOrder62(&At12[index]);
        PDstandardNth3At12 = PDstandardNthfdOrder63(&At12[index]);
        PDstandardNth1At13 = PDstandardNthfdOrder61(&At13[index]);
        PDstandardNth2At13 = PDstandardNthfdOrder62(&At13[index]);
        PDstandardNth3At13 = PDstandardNthfdOrder63(&At13[index]);
        PDstandardNth1At22 = PDstandardNthfdOrder61(&At22[index]);
        PDstandardNth2At22 = PDstandardNthfdOrder62(&At22[index]);
        PDstandardNth3At22 = PDstandardNthfdOrder63(&At22[index]);
        PDstandardNth1At23 = PDstandardNthfdOrder61(&At23[index]);
        PDstandardNth2At23 = PDstandardNthfdOrder62(&At23[index]);
        PDstandardNth3At23 = PDstandardNthfdOrder63(&At23[index]);
        PDstandardNth1At33 = PDstandardNthfdOrder61(&At33[index]);
        PDstandardNth2At33 = PDstandardNthfdOrder62(&At33[index]);
        PDstandardNth3At33 = PDstandardNthfdOrder63(&At33[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder61(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder62(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder63(&gt11[index]);
        PDstandardNth11gt11 = PDstandardNthfdOrder611(&gt11[index]);
        PDstandardNth22gt11 = PDstandardNthfdOrder622(&gt11[index]);
        PDstandardNth33gt11 = PDstandardNthfdOrder633(&gt11[index]);
        PDstandardNth12gt11 = PDstandardNthfdOrder612(&gt11[index]);
        PDstandardNth13gt11 = PDstandardNthfdOrder613(&gt11[index]);
        PDstandardNth23gt11 = PDstandardNthfdOrder623(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder61(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder62(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder63(&gt12[index]);
        PDstandardNth11gt12 = PDstandardNthfdOrder611(&gt12[index]);
        PDstandardNth22gt12 = PDstandardNthfdOrder622(&gt12[index]);
        PDstandardNth33gt12 = PDstandardNthfdOrder633(&gt12[index]);
        PDstandardNth12gt12 = PDstandardNthfdOrder612(&gt12[index]);
        PDstandardNth13gt12 = PDstandardNthfdOrder613(&gt12[index]);
        PDstandardNth23gt12 = PDstandardNthfdOrder623(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder61(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder62(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder63(&gt13[index]);
        PDstandardNth11gt13 = PDstandardNthfdOrder611(&gt13[index]);
        PDstandardNth22gt13 = PDstandardNthfdOrder622(&gt13[index]);
        PDstandardNth33gt13 = PDstandardNthfdOrder633(&gt13[index]);
        PDstandardNth12gt13 = PDstandardNthfdOrder612(&gt13[index]);
        PDstandardNth13gt13 = PDstandardNthfdOrder613(&gt13[index]);
        PDstandardNth23gt13 = PDstandardNthfdOrder623(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder61(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder62(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder63(&gt22[index]);
        PDstandardNth11gt22 = PDstandardNthfdOrder611(&gt22[index]);
        PDstandardNth22gt22 = PDstandardNthfdOrder622(&gt22[index]);
        PDstandardNth33gt22 = PDstandardNthfdOrder633(&gt22[index]);
        PDstandardNth12gt22 = PDstandardNthfdOrder612(&gt22[index]);
        PDstandardNth13gt22 = PDstandardNthfdOrder613(&gt22[index]);
        PDstandardNth23gt22 = PDstandardNthfdOrder623(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder61(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder62(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder63(&gt23[index]);
        PDstandardNth11gt23 = PDstandardNthfdOrder611(&gt23[index]);
        PDstandardNth22gt23 = PDstandardNthfdOrder622(&gt23[index]);
        PDstandardNth33gt23 = PDstandardNthfdOrder633(&gt23[index]);
        PDstandardNth12gt23 = PDstandardNthfdOrder612(&gt23[index]);
        PDstandardNth13gt23 = PDstandardNthfdOrder613(&gt23[index]);
        PDstandardNth23gt23 = PDstandardNthfdOrder623(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder61(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder62(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder63(&gt33[index]);
        PDstandardNth11gt33 = PDstandardNthfdOrder611(&gt33[index]);
        PDstandardNth22gt33 = PDstandardNthfdOrder622(&gt33[index]);
        PDstandardNth33gt33 = PDstandardNthfdOrder633(&gt33[index]);
        PDstandardNth12gt33 = PDstandardNthfdOrder612(&gt33[index]);
        PDstandardNth13gt33 = PDstandardNthfdOrder613(&gt33[index]);
        PDstandardNth23gt33 = PDstandardNthfdOrder623(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder61(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder62(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder63(&phi[index]);
        PDstandardNth11phi = PDstandardNthfdOrder611(&phi[index]);
        PDstandardNth22phi = PDstandardNthfdOrder622(&phi[index]);
        PDstandardNth33phi = PDstandardNthfdOrder633(&phi[index]);
        PDstandardNth12phi = PDstandardNthfdOrder612(&phi[index]);
        PDstandardNth13phi = PDstandardNthfdOrder613(&phi[index]);
        PDstandardNth23phi = PDstandardNthfdOrder623(&phi[index]);
        PDstandardNth1trK = PDstandardNthfdOrder61(&trK[index]);
        PDstandardNth2trK = PDstandardNthfdOrder62(&trK[index]);
        PDstandardNth3trK = PDstandardNthfdOrder63(&trK[index]);
        PDstandardNth1Xt1 = PDstandardNthfdOrder61(&Xt1[index]);
        PDstandardNth2Xt1 = PDstandardNthfdOrder62(&Xt1[index]);
        PDstandardNth3Xt1 = PDstandardNthfdOrder63(&Xt1[index]);
        PDstandardNth1Xt2 = PDstandardNthfdOrder61(&Xt2[index]);
        PDstandardNth2Xt2 = PDstandardNthfdOrder62(&Xt2[index]);
        PDstandardNth3Xt2 = PDstandardNthfdOrder63(&Xt2[index]);
        PDstandardNth1Xt3 = PDstandardNthfdOrder61(&Xt3[index]);
        PDstandardNth2Xt3 = PDstandardNthfdOrder62(&Xt3[index]);
        PDstandardNth3Xt3 = PDstandardNthfdOrder63(&Xt3[index]);
        break;
      }
      
      case 8:
      {
        PDstandardNth1At11 = PDstandardNthfdOrder81(&At11[index]);
        PDstandardNth2At11 = PDstandardNthfdOrder82(&At11[index]);
        PDstandardNth3At11 = PDstandardNthfdOrder83(&At11[index]);
        PDstandardNth1At12 = PDstandardNthfdOrder81(&At12[index]);
        PDstandardNth2At12 = PDstandardNthfdOrder82(&At12[index]);
        PDstandardNth3At12 = PDstandardNthfdOrder83(&At12[index]);
        PDstandardNth1At13 = PDstandardNthfdOrder81(&At13[index]);
        PDstandardNth2At13 = PDstandardNthfdOrder82(&At13[index]);
        PDstandardNth3At13 = PDstandardNthfdOrder83(&At13[index]);
        PDstandardNth1At22 = PDstandardNthfdOrder81(&At22[index]);
        PDstandardNth2At22 = PDstandardNthfdOrder82(&At22[index]);
        PDstandardNth3At22 = PDstandardNthfdOrder83(&At22[index]);
        PDstandardNth1At23 = PDstandardNthfdOrder81(&At23[index]);
        PDstandardNth2At23 = PDstandardNthfdOrder82(&At23[index]);
        PDstandardNth3At23 = PDstandardNthfdOrder83(&At23[index]);
        PDstandardNth1At33 = PDstandardNthfdOrder81(&At33[index]);
        PDstandardNth2At33 = PDstandardNthfdOrder82(&At33[index]);
        PDstandardNth3At33 = PDstandardNthfdOrder83(&At33[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder81(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder82(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder83(&gt11[index]);
        PDstandardNth11gt11 = PDstandardNthfdOrder811(&gt11[index]);
        PDstandardNth22gt11 = PDstandardNthfdOrder822(&gt11[index]);
        PDstandardNth33gt11 = PDstandardNthfdOrder833(&gt11[index]);
        PDstandardNth12gt11 = PDstandardNthfdOrder812(&gt11[index]);
        PDstandardNth13gt11 = PDstandardNthfdOrder813(&gt11[index]);
        PDstandardNth23gt11 = PDstandardNthfdOrder823(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder81(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder82(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder83(&gt12[index]);
        PDstandardNth11gt12 = PDstandardNthfdOrder811(&gt12[index]);
        PDstandardNth22gt12 = PDstandardNthfdOrder822(&gt12[index]);
        PDstandardNth33gt12 = PDstandardNthfdOrder833(&gt12[index]);
        PDstandardNth12gt12 = PDstandardNthfdOrder812(&gt12[index]);
        PDstandardNth13gt12 = PDstandardNthfdOrder813(&gt12[index]);
        PDstandardNth23gt12 = PDstandardNthfdOrder823(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder81(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder82(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder83(&gt13[index]);
        PDstandardNth11gt13 = PDstandardNthfdOrder811(&gt13[index]);
        PDstandardNth22gt13 = PDstandardNthfdOrder822(&gt13[index]);
        PDstandardNth33gt13 = PDstandardNthfdOrder833(&gt13[index]);
        PDstandardNth12gt13 = PDstandardNthfdOrder812(&gt13[index]);
        PDstandardNth13gt13 = PDstandardNthfdOrder813(&gt13[index]);
        PDstandardNth23gt13 = PDstandardNthfdOrder823(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder81(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder82(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder83(&gt22[index]);
        PDstandardNth11gt22 = PDstandardNthfdOrder811(&gt22[index]);
        PDstandardNth22gt22 = PDstandardNthfdOrder822(&gt22[index]);
        PDstandardNth33gt22 = PDstandardNthfdOrder833(&gt22[index]);
        PDstandardNth12gt22 = PDstandardNthfdOrder812(&gt22[index]);
        PDstandardNth13gt22 = PDstandardNthfdOrder813(&gt22[index]);
        PDstandardNth23gt22 = PDstandardNthfdOrder823(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder81(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder82(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder83(&gt23[index]);
        PDstandardNth11gt23 = PDstandardNthfdOrder811(&gt23[index]);
        PDstandardNth22gt23 = PDstandardNthfdOrder822(&gt23[index]);
        PDstandardNth33gt23 = PDstandardNthfdOrder833(&gt23[index]);
        PDstandardNth12gt23 = PDstandardNthfdOrder812(&gt23[index]);
        PDstandardNth13gt23 = PDstandardNthfdOrder813(&gt23[index]);
        PDstandardNth23gt23 = PDstandardNthfdOrder823(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder81(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder82(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder83(&gt33[index]);
        PDstandardNth11gt33 = PDstandardNthfdOrder811(&gt33[index]);
        PDstandardNth22gt33 = PDstandardNthfdOrder822(&gt33[index]);
        PDstandardNth33gt33 = PDstandardNthfdOrder833(&gt33[index]);
        PDstandardNth12gt33 = PDstandardNthfdOrder812(&gt33[index]);
        PDstandardNth13gt33 = PDstandardNthfdOrder813(&gt33[index]);
        PDstandardNth23gt33 = PDstandardNthfdOrder823(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder81(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder82(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder83(&phi[index]);
        PDstandardNth11phi = PDstandardNthfdOrder811(&phi[index]);
        PDstandardNth22phi = PDstandardNthfdOrder822(&phi[index]);
        PDstandardNth33phi = PDstandardNthfdOrder833(&phi[index]);
        PDstandardNth12phi = PDstandardNthfdOrder812(&phi[index]);
        PDstandardNth13phi = PDstandardNthfdOrder813(&phi[index]);
        PDstandardNth23phi = PDstandardNthfdOrder823(&phi[index]);
        PDstandardNth1trK = PDstandardNthfdOrder81(&trK[index]);
        PDstandardNth2trK = PDstandardNthfdOrder82(&trK[index]);
        PDstandardNth3trK = PDstandardNthfdOrder83(&trK[index]);
        PDstandardNth1Xt1 = PDstandardNthfdOrder81(&Xt1[index]);
        PDstandardNth2Xt1 = PDstandardNthfdOrder82(&Xt1[index]);
        PDstandardNth3Xt1 = PDstandardNthfdOrder83(&Xt1[index]);
        PDstandardNth1Xt2 = PDstandardNthfdOrder81(&Xt2[index]);
        PDstandardNth2Xt2 = PDstandardNthfdOrder82(&Xt2[index]);
        PDstandardNth3Xt2 = PDstandardNthfdOrder83(&Xt2[index]);
        PDstandardNth1Xt3 = PDstandardNthfdOrder81(&Xt3[index]);
        PDstandardNth2Xt3 = PDstandardNthfdOrder82(&Xt3[index]);
        PDstandardNth3Xt3 = PDstandardNthfdOrder83(&Xt3[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL JacPDstandardNth11gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1At11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1At12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1At13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1At22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1At23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1At33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2At11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2At12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2At13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2At22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2At23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2At33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3At11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3At12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3At13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3At22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3At23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3At33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3Xt3 CCTK_ATTRIBUTE_UNUSED;
    
    if (usejacobian)
    {
      JacPDstandardNth1At11 = J11L*PDstandardNth1At11 + 
        J21L*PDstandardNth2At11 + J31L*PDstandardNth3At11;
      
      JacPDstandardNth1At12 = J11L*PDstandardNth1At12 + 
        J21L*PDstandardNth2At12 + J31L*PDstandardNth3At12;
      
      JacPDstandardNth1At13 = J11L*PDstandardNth1At13 + 
        J21L*PDstandardNth2At13 + J31L*PDstandardNth3At13;
      
      JacPDstandardNth1At22 = J11L*PDstandardNth1At22 + 
        J21L*PDstandardNth2At22 + J31L*PDstandardNth3At22;
      
      JacPDstandardNth1At23 = J11L*PDstandardNth1At23 + 
        J21L*PDstandardNth2At23 + J31L*PDstandardNth3At23;
      
      JacPDstandardNth1At33 = J11L*PDstandardNth1At33 + 
        J21L*PDstandardNth2At33 + J31L*PDstandardNth3At33;
      
      JacPDstandardNth1gt11 = J11L*PDstandardNth1gt11 + 
        J21L*PDstandardNth2gt11 + J31L*PDstandardNth3gt11;
      
      JacPDstandardNth1gt12 = J11L*PDstandardNth1gt12 + 
        J21L*PDstandardNth2gt12 + J31L*PDstandardNth3gt12;
      
      JacPDstandardNth1gt13 = J11L*PDstandardNth1gt13 + 
        J21L*PDstandardNth2gt13 + J31L*PDstandardNth3gt13;
      
      JacPDstandardNth1gt22 = J11L*PDstandardNth1gt22 + 
        J21L*PDstandardNth2gt22 + J31L*PDstandardNth3gt22;
      
      JacPDstandardNth1gt23 = J11L*PDstandardNth1gt23 + 
        J21L*PDstandardNth2gt23 + J31L*PDstandardNth3gt23;
      
      JacPDstandardNth1gt33 = J11L*PDstandardNth1gt33 + 
        J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33;
      
      JacPDstandardNth1phi = J11L*PDstandardNth1phi + J21L*PDstandardNth2phi 
        + J31L*PDstandardNth3phi;
      
      JacPDstandardNth1trK = J11L*PDstandardNth1trK + J21L*PDstandardNth2trK 
        + J31L*PDstandardNth3trK;
      
      JacPDstandardNth1Xt1 = J11L*PDstandardNth1Xt1 + J21L*PDstandardNth2Xt1 
        + J31L*PDstandardNth3Xt1;
      
      JacPDstandardNth1Xt2 = J11L*PDstandardNth1Xt2 + J21L*PDstandardNth2Xt2 
        + J31L*PDstandardNth3Xt2;
      
      JacPDstandardNth1Xt3 = J11L*PDstandardNth1Xt3 + J21L*PDstandardNth2Xt3 
        + J31L*PDstandardNth3Xt3;
      
      JacPDstandardNth2At11 = J12L*PDstandardNth1At11 + 
        J22L*PDstandardNth2At11 + J32L*PDstandardNth3At11;
      
      JacPDstandardNth2At12 = J12L*PDstandardNth1At12 + 
        J22L*PDstandardNth2At12 + J32L*PDstandardNth3At12;
      
      JacPDstandardNth2At13 = J12L*PDstandardNth1At13 + 
        J22L*PDstandardNth2At13 + J32L*PDstandardNth3At13;
      
      JacPDstandardNth2At22 = J12L*PDstandardNth1At22 + 
        J22L*PDstandardNth2At22 + J32L*PDstandardNth3At22;
      
      JacPDstandardNth2At23 = J12L*PDstandardNth1At23 + 
        J22L*PDstandardNth2At23 + J32L*PDstandardNth3At23;
      
      JacPDstandardNth2At33 = J12L*PDstandardNth1At33 + 
        J22L*PDstandardNth2At33 + J32L*PDstandardNth3At33;
      
      JacPDstandardNth2gt11 = J12L*PDstandardNth1gt11 + 
        J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11;
      
      JacPDstandardNth2gt12 = J12L*PDstandardNth1gt12 + 
        J22L*PDstandardNth2gt12 + J32L*PDstandardNth3gt12;
      
      JacPDstandardNth2gt13 = J12L*PDstandardNth1gt13 + 
        J22L*PDstandardNth2gt13 + J32L*PDstandardNth3gt13;
      
      JacPDstandardNth2gt22 = J12L*PDstandardNth1gt22 + 
        J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22;
      
      JacPDstandardNth2gt23 = J12L*PDstandardNth1gt23 + 
        J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23;
      
      JacPDstandardNth2gt33 = J12L*PDstandardNth1gt33 + 
        J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33;
      
      JacPDstandardNth2phi = J12L*PDstandardNth1phi + J22L*PDstandardNth2phi 
        + J32L*PDstandardNth3phi;
      
      JacPDstandardNth2trK = J12L*PDstandardNth1trK + J22L*PDstandardNth2trK 
        + J32L*PDstandardNth3trK;
      
      JacPDstandardNth2Xt1 = J12L*PDstandardNth1Xt1 + J22L*PDstandardNth2Xt1 
        + J32L*PDstandardNth3Xt1;
      
      JacPDstandardNth2Xt2 = J12L*PDstandardNth1Xt2 + J22L*PDstandardNth2Xt2 
        + J32L*PDstandardNth3Xt2;
      
      JacPDstandardNth2Xt3 = J12L*PDstandardNth1Xt3 + J22L*PDstandardNth2Xt3 
        + J32L*PDstandardNth3Xt3;
      
      JacPDstandardNth3At11 = J13L*PDstandardNth1At11 + 
        J23L*PDstandardNth2At11 + J33L*PDstandardNth3At11;
      
      JacPDstandardNth3At12 = J13L*PDstandardNth1At12 + 
        J23L*PDstandardNth2At12 + J33L*PDstandardNth3At12;
      
      JacPDstandardNth3At13 = J13L*PDstandardNth1At13 + 
        J23L*PDstandardNth2At13 + J33L*PDstandardNth3At13;
      
      JacPDstandardNth3At22 = J13L*PDstandardNth1At22 + 
        J23L*PDstandardNth2At22 + J33L*PDstandardNth3At22;
      
      JacPDstandardNth3At23 = J13L*PDstandardNth1At23 + 
        J23L*PDstandardNth2At23 + J33L*PDstandardNth3At23;
      
      JacPDstandardNth3At33 = J13L*PDstandardNth1At33 + 
        J23L*PDstandardNth2At33 + J33L*PDstandardNth3At33;
      
      JacPDstandardNth3gt11 = J13L*PDstandardNth1gt11 + 
        J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = J13L*PDstandardNth1gt12 + 
        J23L*PDstandardNth2gt12 + J33L*PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = J13L*PDstandardNth1gt13 + 
        J23L*PDstandardNth2gt13 + J33L*PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = J13L*PDstandardNth1gt22 + 
        J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = J13L*PDstandardNth1gt23 + 
        J23L*PDstandardNth2gt23 + J33L*PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = J13L*PDstandardNth1gt33 + 
        J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33;
      
      JacPDstandardNth3phi = J13L*PDstandardNth1phi + J23L*PDstandardNth2phi 
        + J33L*PDstandardNth3phi;
      
      JacPDstandardNth3trK = J13L*PDstandardNth1trK + J23L*PDstandardNth2trK 
        + J33L*PDstandardNth3trK;
      
      JacPDstandardNth3Xt1 = J13L*PDstandardNth1Xt1 + J23L*PDstandardNth2Xt1 
        + J33L*PDstandardNth3Xt1;
      
      JacPDstandardNth3Xt2 = J13L*PDstandardNth1Xt2 + J23L*PDstandardNth2Xt2 
        + J33L*PDstandardNth3Xt2;
      
      JacPDstandardNth3Xt3 = J13L*PDstandardNth1Xt3 + J23L*PDstandardNth2Xt3 
        + J33L*PDstandardNth3Xt3;
      
      JacPDstandardNth11gt11 = dJ111L*PDstandardNth1gt11 + 
        2*(J11L*(J21L*PDstandardNth12gt11 + J31L*PDstandardNth13gt11) + 
        J21L*J31L*PDstandardNth23gt11) + dJ211L*PDstandardNth2gt11 + 
        dJ311L*PDstandardNth3gt11 + PDstandardNth11gt11*pow(J11L,2) + 
        PDstandardNth22gt11*pow(J21L,2) + PDstandardNth33gt11*pow(J31L,2);
      
      JacPDstandardNth11gt12 = dJ111L*PDstandardNth1gt12 + 
        2*(J11L*(J21L*PDstandardNth12gt12 + J31L*PDstandardNth13gt12) + 
        J21L*J31L*PDstandardNth23gt12) + dJ211L*PDstandardNth2gt12 + 
        dJ311L*PDstandardNth3gt12 + PDstandardNth11gt12*pow(J11L,2) + 
        PDstandardNth22gt12*pow(J21L,2) + PDstandardNth33gt12*pow(J31L,2);
      
      JacPDstandardNth11gt13 = dJ111L*PDstandardNth1gt13 + 
        2*(J11L*(J21L*PDstandardNth12gt13 + J31L*PDstandardNth13gt13) + 
        J21L*J31L*PDstandardNth23gt13) + dJ211L*PDstandardNth2gt13 + 
        dJ311L*PDstandardNth3gt13 + PDstandardNth11gt13*pow(J11L,2) + 
        PDstandardNth22gt13*pow(J21L,2) + PDstandardNth33gt13*pow(J31L,2);
      
      JacPDstandardNth11gt22 = dJ111L*PDstandardNth1gt22 + 
        2*(J11L*(J21L*PDstandardNth12gt22 + J31L*PDstandardNth13gt22) + 
        J21L*J31L*PDstandardNth23gt22) + dJ211L*PDstandardNth2gt22 + 
        dJ311L*PDstandardNth3gt22 + PDstandardNth11gt22*pow(J11L,2) + 
        PDstandardNth22gt22*pow(J21L,2) + PDstandardNth33gt22*pow(J31L,2);
      
      JacPDstandardNth11gt23 = dJ111L*PDstandardNth1gt23 + 
        2*(J11L*(J21L*PDstandardNth12gt23 + J31L*PDstandardNth13gt23) + 
        J21L*J31L*PDstandardNth23gt23) + dJ211L*PDstandardNth2gt23 + 
        dJ311L*PDstandardNth3gt23 + PDstandardNth11gt23*pow(J11L,2) + 
        PDstandardNth22gt23*pow(J21L,2) + PDstandardNth33gt23*pow(J31L,2);
      
      JacPDstandardNth11gt33 = dJ111L*PDstandardNth1gt33 + 
        2*(J11L*(J21L*PDstandardNth12gt33 + J31L*PDstandardNth13gt33) + 
        J21L*J31L*PDstandardNth23gt33) + dJ211L*PDstandardNth2gt33 + 
        dJ311L*PDstandardNth3gt33 + PDstandardNth11gt33*pow(J11L,2) + 
        PDstandardNth22gt33*pow(J21L,2) + PDstandardNth33gt33*pow(J31L,2);
      
      JacPDstandardNth11phi = dJ111L*PDstandardNth1phi + 
        2*(J11L*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi) + 
        J21L*J31L*PDstandardNth23phi) + dJ211L*PDstandardNth2phi + 
        dJ311L*PDstandardNth3phi + PDstandardNth11phi*pow(J11L,2) + 
        PDstandardNth22phi*pow(J21L,2) + PDstandardNth33phi*pow(J31L,2);
      
      JacPDstandardNth22gt11 = dJ122L*PDstandardNth1gt11 + 
        2*(J12L*(J22L*PDstandardNth12gt11 + J32L*PDstandardNth13gt11) + 
        J22L*J32L*PDstandardNth23gt11) + dJ222L*PDstandardNth2gt11 + 
        dJ322L*PDstandardNth3gt11 + PDstandardNth11gt11*pow(J12L,2) + 
        PDstandardNth22gt11*pow(J22L,2) + PDstandardNth33gt11*pow(J32L,2);
      
      JacPDstandardNth22gt12 = dJ122L*PDstandardNth1gt12 + 
        2*(J12L*(J22L*PDstandardNth12gt12 + J32L*PDstandardNth13gt12) + 
        J22L*J32L*PDstandardNth23gt12) + dJ222L*PDstandardNth2gt12 + 
        dJ322L*PDstandardNth3gt12 + PDstandardNth11gt12*pow(J12L,2) + 
        PDstandardNth22gt12*pow(J22L,2) + PDstandardNth33gt12*pow(J32L,2);
      
      JacPDstandardNth22gt13 = dJ122L*PDstandardNth1gt13 + 
        2*(J12L*(J22L*PDstandardNth12gt13 + J32L*PDstandardNth13gt13) + 
        J22L*J32L*PDstandardNth23gt13) + dJ222L*PDstandardNth2gt13 + 
        dJ322L*PDstandardNth3gt13 + PDstandardNth11gt13*pow(J12L,2) + 
        PDstandardNth22gt13*pow(J22L,2) + PDstandardNth33gt13*pow(J32L,2);
      
      JacPDstandardNth22gt22 = dJ122L*PDstandardNth1gt22 + 
        2*(J12L*(J22L*PDstandardNth12gt22 + J32L*PDstandardNth13gt22) + 
        J22L*J32L*PDstandardNth23gt22) + dJ222L*PDstandardNth2gt22 + 
        dJ322L*PDstandardNth3gt22 + PDstandardNth11gt22*pow(J12L,2) + 
        PDstandardNth22gt22*pow(J22L,2) + PDstandardNth33gt22*pow(J32L,2);
      
      JacPDstandardNth22gt23 = dJ122L*PDstandardNth1gt23 + 
        2*(J12L*(J22L*PDstandardNth12gt23 + J32L*PDstandardNth13gt23) + 
        J22L*J32L*PDstandardNth23gt23) + dJ222L*PDstandardNth2gt23 + 
        dJ322L*PDstandardNth3gt23 + PDstandardNth11gt23*pow(J12L,2) + 
        PDstandardNth22gt23*pow(J22L,2) + PDstandardNth33gt23*pow(J32L,2);
      
      JacPDstandardNth22gt33 = dJ122L*PDstandardNth1gt33 + 
        2*(J12L*(J22L*PDstandardNth12gt33 + J32L*PDstandardNth13gt33) + 
        J22L*J32L*PDstandardNth23gt33) + dJ222L*PDstandardNth2gt33 + 
        dJ322L*PDstandardNth3gt33 + PDstandardNth11gt33*pow(J12L,2) + 
        PDstandardNth22gt33*pow(J22L,2) + PDstandardNth33gt33*pow(J32L,2);
      
      JacPDstandardNth22phi = dJ122L*PDstandardNth1phi + 
        2*(J12L*(J22L*PDstandardNth12phi + J32L*PDstandardNth13phi) + 
        J22L*J32L*PDstandardNth23phi) + dJ222L*PDstandardNth2phi + 
        dJ322L*PDstandardNth3phi + PDstandardNth11phi*pow(J12L,2) + 
        PDstandardNth22phi*pow(J22L,2) + PDstandardNth33phi*pow(J32L,2);
      
      JacPDstandardNth33gt11 = dJ133L*PDstandardNth1gt11 + 
        2*(J13L*(J23L*PDstandardNth12gt11 + J33L*PDstandardNth13gt11) + 
        J23L*J33L*PDstandardNth23gt11) + dJ233L*PDstandardNth2gt11 + 
        dJ333L*PDstandardNth3gt11 + PDstandardNth11gt11*pow(J13L,2) + 
        PDstandardNth22gt11*pow(J23L,2) + PDstandardNth33gt11*pow(J33L,2);
      
      JacPDstandardNth33gt12 = dJ133L*PDstandardNth1gt12 + 
        2*(J13L*(J23L*PDstandardNth12gt12 + J33L*PDstandardNth13gt12) + 
        J23L*J33L*PDstandardNth23gt12) + dJ233L*PDstandardNth2gt12 + 
        dJ333L*PDstandardNth3gt12 + PDstandardNth11gt12*pow(J13L,2) + 
        PDstandardNth22gt12*pow(J23L,2) + PDstandardNth33gt12*pow(J33L,2);
      
      JacPDstandardNth33gt13 = dJ133L*PDstandardNth1gt13 + 
        2*(J13L*(J23L*PDstandardNth12gt13 + J33L*PDstandardNth13gt13) + 
        J23L*J33L*PDstandardNth23gt13) + dJ233L*PDstandardNth2gt13 + 
        dJ333L*PDstandardNth3gt13 + PDstandardNth11gt13*pow(J13L,2) + 
        PDstandardNth22gt13*pow(J23L,2) + PDstandardNth33gt13*pow(J33L,2);
      
      JacPDstandardNth33gt22 = dJ133L*PDstandardNth1gt22 + 
        2*(J13L*(J23L*PDstandardNth12gt22 + J33L*PDstandardNth13gt22) + 
        J23L*J33L*PDstandardNth23gt22) + dJ233L*PDstandardNth2gt22 + 
        dJ333L*PDstandardNth3gt22 + PDstandardNth11gt22*pow(J13L,2) + 
        PDstandardNth22gt22*pow(J23L,2) + PDstandardNth33gt22*pow(J33L,2);
      
      JacPDstandardNth33gt23 = dJ133L*PDstandardNth1gt23 + 
        2*(J13L*(J23L*PDstandardNth12gt23 + J33L*PDstandardNth13gt23) + 
        J23L*J33L*PDstandardNth23gt23) + dJ233L*PDstandardNth2gt23 + 
        dJ333L*PDstandardNth3gt23 + PDstandardNth11gt23*pow(J13L,2) + 
        PDstandardNth22gt23*pow(J23L,2) + PDstandardNth33gt23*pow(J33L,2);
      
      JacPDstandardNth33gt33 = dJ133L*PDstandardNth1gt33 + 
        2*(J13L*(J23L*PDstandardNth12gt33 + J33L*PDstandardNth13gt33) + 
        J23L*J33L*PDstandardNth23gt33) + dJ233L*PDstandardNth2gt33 + 
        dJ333L*PDstandardNth3gt33 + PDstandardNth11gt33*pow(J13L,2) + 
        PDstandardNth22gt33*pow(J23L,2) + PDstandardNth33gt33*pow(J33L,2);
      
      JacPDstandardNth33phi = dJ133L*PDstandardNth1phi + 
        2*(J13L*(J23L*PDstandardNth12phi + J33L*PDstandardNth13phi) + 
        J23L*J33L*PDstandardNth23phi) + dJ233L*PDstandardNth2phi + 
        dJ333L*PDstandardNth3phi + PDstandardNth11phi*pow(J13L,2) + 
        PDstandardNth22phi*pow(J23L,2) + PDstandardNth33phi*pow(J33L,2);
      
      JacPDstandardNth12gt11 = J12L*(J11L*PDstandardNth11gt11 + 
        J21L*PDstandardNth12gt11 + J31L*PDstandardNth13gt11) + 
        J11L*(J22L*PDstandardNth12gt11 + J32L*PDstandardNth13gt11) + 
        dJ112L*PDstandardNth1gt11 + J22L*(J21L*PDstandardNth22gt11 + 
        J31L*PDstandardNth23gt11) + dJ212L*PDstandardNth2gt11 + 
        J32L*(J21L*PDstandardNth23gt11 + J31L*PDstandardNth33gt11) + 
        dJ312L*PDstandardNth3gt11;
      
      JacPDstandardNth12gt12 = J12L*(J11L*PDstandardNth11gt12 + 
        J21L*PDstandardNth12gt12 + J31L*PDstandardNth13gt12) + 
        J11L*(J22L*PDstandardNth12gt12 + J32L*PDstandardNth13gt12) + 
        dJ112L*PDstandardNth1gt12 + J22L*(J21L*PDstandardNth22gt12 + 
        J31L*PDstandardNth23gt12) + dJ212L*PDstandardNth2gt12 + 
        J32L*(J21L*PDstandardNth23gt12 + J31L*PDstandardNth33gt12) + 
        dJ312L*PDstandardNth3gt12;
      
      JacPDstandardNth12gt13 = J12L*(J11L*PDstandardNth11gt13 + 
        J21L*PDstandardNth12gt13 + J31L*PDstandardNth13gt13) + 
        J11L*(J22L*PDstandardNth12gt13 + J32L*PDstandardNth13gt13) + 
        dJ112L*PDstandardNth1gt13 + J22L*(J21L*PDstandardNth22gt13 + 
        J31L*PDstandardNth23gt13) + dJ212L*PDstandardNth2gt13 + 
        J32L*(J21L*PDstandardNth23gt13 + J31L*PDstandardNth33gt13) + 
        dJ312L*PDstandardNth3gt13;
      
      JacPDstandardNth12gt22 = J12L*(J11L*PDstandardNth11gt22 + 
        J21L*PDstandardNth12gt22 + J31L*PDstandardNth13gt22) + 
        J11L*(J22L*PDstandardNth12gt22 + J32L*PDstandardNth13gt22) + 
        dJ112L*PDstandardNth1gt22 + J22L*(J21L*PDstandardNth22gt22 + 
        J31L*PDstandardNth23gt22) + dJ212L*PDstandardNth2gt22 + 
        J32L*(J21L*PDstandardNth23gt22 + J31L*PDstandardNth33gt22) + 
        dJ312L*PDstandardNth3gt22;
      
      JacPDstandardNth12gt23 = J12L*(J11L*PDstandardNth11gt23 + 
        J21L*PDstandardNth12gt23 + J31L*PDstandardNth13gt23) + 
        J11L*(J22L*PDstandardNth12gt23 + J32L*PDstandardNth13gt23) + 
        dJ112L*PDstandardNth1gt23 + J22L*(J21L*PDstandardNth22gt23 + 
        J31L*PDstandardNth23gt23) + dJ212L*PDstandardNth2gt23 + 
        J32L*(J21L*PDstandardNth23gt23 + J31L*PDstandardNth33gt23) + 
        dJ312L*PDstandardNth3gt23;
      
      JacPDstandardNth12gt33 = J12L*(J11L*PDstandardNth11gt33 + 
        J21L*PDstandardNth12gt33 + J31L*PDstandardNth13gt33) + 
        J11L*(J22L*PDstandardNth12gt33 + J32L*PDstandardNth13gt33) + 
        dJ112L*PDstandardNth1gt33 + J22L*(J21L*PDstandardNth22gt33 + 
        J31L*PDstandardNth23gt33) + dJ212L*PDstandardNth2gt33 + 
        J32L*(J21L*PDstandardNth23gt33 + J31L*PDstandardNth33gt33) + 
        dJ312L*PDstandardNth3gt33;
      
      JacPDstandardNth12phi = J12L*(J11L*PDstandardNth11phi + 
        J21L*PDstandardNth12phi + J31L*PDstandardNth13phi) + 
        J11L*(J22L*PDstandardNth12phi + J32L*PDstandardNth13phi) + 
        dJ112L*PDstandardNth1phi + J22L*(J21L*PDstandardNth22phi + 
        J31L*PDstandardNth23phi) + dJ212L*PDstandardNth2phi + 
        J32L*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi) + 
        dJ312L*PDstandardNth3phi;
      
      JacPDstandardNth13gt11 = J13L*(J11L*PDstandardNth11gt11 + 
        J21L*PDstandardNth12gt11 + J31L*PDstandardNth13gt11) + 
        J11L*(J23L*PDstandardNth12gt11 + J33L*PDstandardNth13gt11) + 
        dJ113L*PDstandardNth1gt11 + J23L*(J21L*PDstandardNth22gt11 + 
        J31L*PDstandardNth23gt11) + dJ213L*PDstandardNth2gt11 + 
        J33L*(J21L*PDstandardNth23gt11 + J31L*PDstandardNth33gt11) + 
        dJ313L*PDstandardNth3gt11;
      
      JacPDstandardNth13gt12 = J13L*(J11L*PDstandardNth11gt12 + 
        J21L*PDstandardNth12gt12 + J31L*PDstandardNth13gt12) + 
        J11L*(J23L*PDstandardNth12gt12 + J33L*PDstandardNth13gt12) + 
        dJ113L*PDstandardNth1gt12 + J23L*(J21L*PDstandardNth22gt12 + 
        J31L*PDstandardNth23gt12) + dJ213L*PDstandardNth2gt12 + 
        J33L*(J21L*PDstandardNth23gt12 + J31L*PDstandardNth33gt12) + 
        dJ313L*PDstandardNth3gt12;
      
      JacPDstandardNth13gt13 = J13L*(J11L*PDstandardNth11gt13 + 
        J21L*PDstandardNth12gt13 + J31L*PDstandardNth13gt13) + 
        J11L*(J23L*PDstandardNth12gt13 + J33L*PDstandardNth13gt13) + 
        dJ113L*PDstandardNth1gt13 + J23L*(J21L*PDstandardNth22gt13 + 
        J31L*PDstandardNth23gt13) + dJ213L*PDstandardNth2gt13 + 
        J33L*(J21L*PDstandardNth23gt13 + J31L*PDstandardNth33gt13) + 
        dJ313L*PDstandardNth3gt13;
      
      JacPDstandardNth13gt22 = J13L*(J11L*PDstandardNth11gt22 + 
        J21L*PDstandardNth12gt22 + J31L*PDstandardNth13gt22) + 
        J11L*(J23L*PDstandardNth12gt22 + J33L*PDstandardNth13gt22) + 
        dJ113L*PDstandardNth1gt22 + J23L*(J21L*PDstandardNth22gt22 + 
        J31L*PDstandardNth23gt22) + dJ213L*PDstandardNth2gt22 + 
        J33L*(J21L*PDstandardNth23gt22 + J31L*PDstandardNth33gt22) + 
        dJ313L*PDstandardNth3gt22;
      
      JacPDstandardNth13gt23 = J13L*(J11L*PDstandardNth11gt23 + 
        J21L*PDstandardNth12gt23 + J31L*PDstandardNth13gt23) + 
        J11L*(J23L*PDstandardNth12gt23 + J33L*PDstandardNth13gt23) + 
        dJ113L*PDstandardNth1gt23 + J23L*(J21L*PDstandardNth22gt23 + 
        J31L*PDstandardNth23gt23) + dJ213L*PDstandardNth2gt23 + 
        J33L*(J21L*PDstandardNth23gt23 + J31L*PDstandardNth33gt23) + 
        dJ313L*PDstandardNth3gt23;
      
      JacPDstandardNth13gt33 = J13L*(J11L*PDstandardNth11gt33 + 
        J21L*PDstandardNth12gt33 + J31L*PDstandardNth13gt33) + 
        J11L*(J23L*PDstandardNth12gt33 + J33L*PDstandardNth13gt33) + 
        dJ113L*PDstandardNth1gt33 + J23L*(J21L*PDstandardNth22gt33 + 
        J31L*PDstandardNth23gt33) + dJ213L*PDstandardNth2gt33 + 
        J33L*(J21L*PDstandardNth23gt33 + J31L*PDstandardNth33gt33) + 
        dJ313L*PDstandardNth3gt33;
      
      JacPDstandardNth13phi = J13L*(J11L*PDstandardNth11phi + 
        J21L*PDstandardNth12phi + J31L*PDstandardNth13phi) + 
        J11L*(J23L*PDstandardNth12phi + J33L*PDstandardNth13phi) + 
        dJ113L*PDstandardNth1phi + J23L*(J21L*PDstandardNth22phi + 
        J31L*PDstandardNth23phi) + dJ213L*PDstandardNth2phi + 
        J33L*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi) + 
        dJ313L*PDstandardNth3phi;
      
      JacPDstandardNth21gt11 = J12L*(J11L*PDstandardNth11gt11 + 
        J21L*PDstandardNth12gt11 + J31L*PDstandardNth13gt11) + 
        J11L*(J22L*PDstandardNth12gt11 + J32L*PDstandardNth13gt11) + 
        dJ112L*PDstandardNth1gt11 + J22L*(J21L*PDstandardNth22gt11 + 
        J31L*PDstandardNth23gt11) + dJ212L*PDstandardNth2gt11 + 
        J32L*(J21L*PDstandardNth23gt11 + J31L*PDstandardNth33gt11) + 
        dJ312L*PDstandardNth3gt11;
      
      JacPDstandardNth21gt12 = J12L*(J11L*PDstandardNth11gt12 + 
        J21L*PDstandardNth12gt12 + J31L*PDstandardNth13gt12) + 
        J11L*(J22L*PDstandardNth12gt12 + J32L*PDstandardNth13gt12) + 
        dJ112L*PDstandardNth1gt12 + J22L*(J21L*PDstandardNth22gt12 + 
        J31L*PDstandardNth23gt12) + dJ212L*PDstandardNth2gt12 + 
        J32L*(J21L*PDstandardNth23gt12 + J31L*PDstandardNth33gt12) + 
        dJ312L*PDstandardNth3gt12;
      
      JacPDstandardNth21gt13 = J12L*(J11L*PDstandardNth11gt13 + 
        J21L*PDstandardNth12gt13 + J31L*PDstandardNth13gt13) + 
        J11L*(J22L*PDstandardNth12gt13 + J32L*PDstandardNth13gt13) + 
        dJ112L*PDstandardNth1gt13 + J22L*(J21L*PDstandardNth22gt13 + 
        J31L*PDstandardNth23gt13) + dJ212L*PDstandardNth2gt13 + 
        J32L*(J21L*PDstandardNth23gt13 + J31L*PDstandardNth33gt13) + 
        dJ312L*PDstandardNth3gt13;
      
      JacPDstandardNth21gt22 = J12L*(J11L*PDstandardNth11gt22 + 
        J21L*PDstandardNth12gt22 + J31L*PDstandardNth13gt22) + 
        J11L*(J22L*PDstandardNth12gt22 + J32L*PDstandardNth13gt22) + 
        dJ112L*PDstandardNth1gt22 + J22L*(J21L*PDstandardNth22gt22 + 
        J31L*PDstandardNth23gt22) + dJ212L*PDstandardNth2gt22 + 
        J32L*(J21L*PDstandardNth23gt22 + J31L*PDstandardNth33gt22) + 
        dJ312L*PDstandardNth3gt22;
      
      JacPDstandardNth21gt23 = J12L*(J11L*PDstandardNth11gt23 + 
        J21L*PDstandardNth12gt23 + J31L*PDstandardNth13gt23) + 
        J11L*(J22L*PDstandardNth12gt23 + J32L*PDstandardNth13gt23) + 
        dJ112L*PDstandardNth1gt23 + J22L*(J21L*PDstandardNth22gt23 + 
        J31L*PDstandardNth23gt23) + dJ212L*PDstandardNth2gt23 + 
        J32L*(J21L*PDstandardNth23gt23 + J31L*PDstandardNth33gt23) + 
        dJ312L*PDstandardNth3gt23;
      
      JacPDstandardNth21gt33 = J12L*(J11L*PDstandardNth11gt33 + 
        J21L*PDstandardNth12gt33 + J31L*PDstandardNth13gt33) + 
        J11L*(J22L*PDstandardNth12gt33 + J32L*PDstandardNth13gt33) + 
        dJ112L*PDstandardNth1gt33 + J22L*(J21L*PDstandardNth22gt33 + 
        J31L*PDstandardNth23gt33) + dJ212L*PDstandardNth2gt33 + 
        J32L*(J21L*PDstandardNth23gt33 + J31L*PDstandardNth33gt33) + 
        dJ312L*PDstandardNth3gt33;
      
      JacPDstandardNth21phi = J12L*(J11L*PDstandardNth11phi + 
        J21L*PDstandardNth12phi + J31L*PDstandardNth13phi) + 
        J11L*(J22L*PDstandardNth12phi + J32L*PDstandardNth13phi) + 
        dJ112L*PDstandardNth1phi + J22L*(J21L*PDstandardNth22phi + 
        J31L*PDstandardNth23phi) + dJ212L*PDstandardNth2phi + 
        J32L*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi) + 
        dJ312L*PDstandardNth3phi;
      
      JacPDstandardNth23gt11 = J13L*(J12L*PDstandardNth11gt11 + 
        J22L*PDstandardNth12gt11 + J32L*PDstandardNth13gt11) + 
        J12L*(J23L*PDstandardNth12gt11 + J33L*PDstandardNth13gt11) + 
        dJ123L*PDstandardNth1gt11 + J23L*(J22L*PDstandardNth22gt11 + 
        J32L*PDstandardNth23gt11) + dJ223L*PDstandardNth2gt11 + 
        J33L*(J22L*PDstandardNth23gt11 + J32L*PDstandardNth33gt11) + 
        dJ323L*PDstandardNth3gt11;
      
      JacPDstandardNth23gt12 = J13L*(J12L*PDstandardNth11gt12 + 
        J22L*PDstandardNth12gt12 + J32L*PDstandardNth13gt12) + 
        J12L*(J23L*PDstandardNth12gt12 + J33L*PDstandardNth13gt12) + 
        dJ123L*PDstandardNth1gt12 + J23L*(J22L*PDstandardNth22gt12 + 
        J32L*PDstandardNth23gt12) + dJ223L*PDstandardNth2gt12 + 
        J33L*(J22L*PDstandardNth23gt12 + J32L*PDstandardNth33gt12) + 
        dJ323L*PDstandardNth3gt12;
      
      JacPDstandardNth23gt13 = J13L*(J12L*PDstandardNth11gt13 + 
        J22L*PDstandardNth12gt13 + J32L*PDstandardNth13gt13) + 
        J12L*(J23L*PDstandardNth12gt13 + J33L*PDstandardNth13gt13) + 
        dJ123L*PDstandardNth1gt13 + J23L*(J22L*PDstandardNth22gt13 + 
        J32L*PDstandardNth23gt13) + dJ223L*PDstandardNth2gt13 + 
        J33L*(J22L*PDstandardNth23gt13 + J32L*PDstandardNth33gt13) + 
        dJ323L*PDstandardNth3gt13;
      
      JacPDstandardNth23gt22 = J13L*(J12L*PDstandardNth11gt22 + 
        J22L*PDstandardNth12gt22 + J32L*PDstandardNth13gt22) + 
        J12L*(J23L*PDstandardNth12gt22 + J33L*PDstandardNth13gt22) + 
        dJ123L*PDstandardNth1gt22 + J23L*(J22L*PDstandardNth22gt22 + 
        J32L*PDstandardNth23gt22) + dJ223L*PDstandardNth2gt22 + 
        J33L*(J22L*PDstandardNth23gt22 + J32L*PDstandardNth33gt22) + 
        dJ323L*PDstandardNth3gt22;
      
      JacPDstandardNth23gt23 = J13L*(J12L*PDstandardNth11gt23 + 
        J22L*PDstandardNth12gt23 + J32L*PDstandardNth13gt23) + 
        J12L*(J23L*PDstandardNth12gt23 + J33L*PDstandardNth13gt23) + 
        dJ123L*PDstandardNth1gt23 + J23L*(J22L*PDstandardNth22gt23 + 
        J32L*PDstandardNth23gt23) + dJ223L*PDstandardNth2gt23 + 
        J33L*(J22L*PDstandardNth23gt23 + J32L*PDstandardNth33gt23) + 
        dJ323L*PDstandardNth3gt23;
      
      JacPDstandardNth23gt33 = J13L*(J12L*PDstandardNth11gt33 + 
        J22L*PDstandardNth12gt33 + J32L*PDstandardNth13gt33) + 
        J12L*(J23L*PDstandardNth12gt33 + J33L*PDstandardNth13gt33) + 
        dJ123L*PDstandardNth1gt33 + J23L*(J22L*PDstandardNth22gt33 + 
        J32L*PDstandardNth23gt33) + dJ223L*PDstandardNth2gt33 + 
        J33L*(J22L*PDstandardNth23gt33 + J32L*PDstandardNth33gt33) + 
        dJ323L*PDstandardNth3gt33;
      
      JacPDstandardNth23phi = J13L*(J12L*PDstandardNth11phi + 
        J22L*PDstandardNth12phi + J32L*PDstandardNth13phi) + 
        J12L*(J23L*PDstandardNth12phi + J33L*PDstandardNth13phi) + 
        dJ123L*PDstandardNth1phi + J23L*(J22L*PDstandardNth22phi + 
        J32L*PDstandardNth23phi) + dJ223L*PDstandardNth2phi + 
        J33L*(J22L*PDstandardNth23phi + J32L*PDstandardNth33phi) + 
        dJ323L*PDstandardNth3phi;
      
      JacPDstandardNth31gt11 = J13L*(J11L*PDstandardNth11gt11 + 
        J21L*PDstandardNth12gt11 + J31L*PDstandardNth13gt11) + 
        J11L*(J23L*PDstandardNth12gt11 + J33L*PDstandardNth13gt11) + 
        dJ113L*PDstandardNth1gt11 + J23L*(J21L*PDstandardNth22gt11 + 
        J31L*PDstandardNth23gt11) + dJ213L*PDstandardNth2gt11 + 
        J33L*(J21L*PDstandardNth23gt11 + J31L*PDstandardNth33gt11) + 
        dJ313L*PDstandardNth3gt11;
      
      JacPDstandardNth31gt12 = J13L*(J11L*PDstandardNth11gt12 + 
        J21L*PDstandardNth12gt12 + J31L*PDstandardNth13gt12) + 
        J11L*(J23L*PDstandardNth12gt12 + J33L*PDstandardNth13gt12) + 
        dJ113L*PDstandardNth1gt12 + J23L*(J21L*PDstandardNth22gt12 + 
        J31L*PDstandardNth23gt12) + dJ213L*PDstandardNth2gt12 + 
        J33L*(J21L*PDstandardNth23gt12 + J31L*PDstandardNth33gt12) + 
        dJ313L*PDstandardNth3gt12;
      
      JacPDstandardNth31gt13 = J13L*(J11L*PDstandardNth11gt13 + 
        J21L*PDstandardNth12gt13 + J31L*PDstandardNth13gt13) + 
        J11L*(J23L*PDstandardNth12gt13 + J33L*PDstandardNth13gt13) + 
        dJ113L*PDstandardNth1gt13 + J23L*(J21L*PDstandardNth22gt13 + 
        J31L*PDstandardNth23gt13) + dJ213L*PDstandardNth2gt13 + 
        J33L*(J21L*PDstandardNth23gt13 + J31L*PDstandardNth33gt13) + 
        dJ313L*PDstandardNth3gt13;
      
      JacPDstandardNth31gt22 = J13L*(J11L*PDstandardNth11gt22 + 
        J21L*PDstandardNth12gt22 + J31L*PDstandardNth13gt22) + 
        J11L*(J23L*PDstandardNth12gt22 + J33L*PDstandardNth13gt22) + 
        dJ113L*PDstandardNth1gt22 + J23L*(J21L*PDstandardNth22gt22 + 
        J31L*PDstandardNth23gt22) + dJ213L*PDstandardNth2gt22 + 
        J33L*(J21L*PDstandardNth23gt22 + J31L*PDstandardNth33gt22) + 
        dJ313L*PDstandardNth3gt22;
      
      JacPDstandardNth31gt23 = J13L*(J11L*PDstandardNth11gt23 + 
        J21L*PDstandardNth12gt23 + J31L*PDstandardNth13gt23) + 
        J11L*(J23L*PDstandardNth12gt23 + J33L*PDstandardNth13gt23) + 
        dJ113L*PDstandardNth1gt23 + J23L*(J21L*PDstandardNth22gt23 + 
        J31L*PDstandardNth23gt23) + dJ213L*PDstandardNth2gt23 + 
        J33L*(J21L*PDstandardNth23gt23 + J31L*PDstandardNth33gt23) + 
        dJ313L*PDstandardNth3gt23;
      
      JacPDstandardNth31gt33 = J13L*(J11L*PDstandardNth11gt33 + 
        J21L*PDstandardNth12gt33 + J31L*PDstandardNth13gt33) + 
        J11L*(J23L*PDstandardNth12gt33 + J33L*PDstandardNth13gt33) + 
        dJ113L*PDstandardNth1gt33 + J23L*(J21L*PDstandardNth22gt33 + 
        J31L*PDstandardNth23gt33) + dJ213L*PDstandardNth2gt33 + 
        J33L*(J21L*PDstandardNth23gt33 + J31L*PDstandardNth33gt33) + 
        dJ313L*PDstandardNth3gt33;
      
      JacPDstandardNth31phi = J13L*(J11L*PDstandardNth11phi + 
        J21L*PDstandardNth12phi + J31L*PDstandardNth13phi) + 
        J11L*(J23L*PDstandardNth12phi + J33L*PDstandardNth13phi) + 
        dJ113L*PDstandardNth1phi + J23L*(J21L*PDstandardNth22phi + 
        J31L*PDstandardNth23phi) + dJ213L*PDstandardNth2phi + 
        J33L*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi) + 
        dJ313L*PDstandardNth3phi;
      
      JacPDstandardNth32gt11 = J13L*(J12L*PDstandardNth11gt11 + 
        J22L*PDstandardNth12gt11 + J32L*PDstandardNth13gt11) + 
        J12L*(J23L*PDstandardNth12gt11 + J33L*PDstandardNth13gt11) + 
        dJ123L*PDstandardNth1gt11 + J23L*(J22L*PDstandardNth22gt11 + 
        J32L*PDstandardNth23gt11) + dJ223L*PDstandardNth2gt11 + 
        J33L*(J22L*PDstandardNth23gt11 + J32L*PDstandardNth33gt11) + 
        dJ323L*PDstandardNth3gt11;
      
      JacPDstandardNth32gt12 = J13L*(J12L*PDstandardNth11gt12 + 
        J22L*PDstandardNth12gt12 + J32L*PDstandardNth13gt12) + 
        J12L*(J23L*PDstandardNth12gt12 + J33L*PDstandardNth13gt12) + 
        dJ123L*PDstandardNth1gt12 + J23L*(J22L*PDstandardNth22gt12 + 
        J32L*PDstandardNth23gt12) + dJ223L*PDstandardNth2gt12 + 
        J33L*(J22L*PDstandardNth23gt12 + J32L*PDstandardNth33gt12) + 
        dJ323L*PDstandardNth3gt12;
      
      JacPDstandardNth32gt13 = J13L*(J12L*PDstandardNth11gt13 + 
        J22L*PDstandardNth12gt13 + J32L*PDstandardNth13gt13) + 
        J12L*(J23L*PDstandardNth12gt13 + J33L*PDstandardNth13gt13) + 
        dJ123L*PDstandardNth1gt13 + J23L*(J22L*PDstandardNth22gt13 + 
        J32L*PDstandardNth23gt13) + dJ223L*PDstandardNth2gt13 + 
        J33L*(J22L*PDstandardNth23gt13 + J32L*PDstandardNth33gt13) + 
        dJ323L*PDstandardNth3gt13;
      
      JacPDstandardNth32gt22 = J13L*(J12L*PDstandardNth11gt22 + 
        J22L*PDstandardNth12gt22 + J32L*PDstandardNth13gt22) + 
        J12L*(J23L*PDstandardNth12gt22 + J33L*PDstandardNth13gt22) + 
        dJ123L*PDstandardNth1gt22 + J23L*(J22L*PDstandardNth22gt22 + 
        J32L*PDstandardNth23gt22) + dJ223L*PDstandardNth2gt22 + 
        J33L*(J22L*PDstandardNth23gt22 + J32L*PDstandardNth33gt22) + 
        dJ323L*PDstandardNth3gt22;
      
      JacPDstandardNth32gt23 = J13L*(J12L*PDstandardNth11gt23 + 
        J22L*PDstandardNth12gt23 + J32L*PDstandardNth13gt23) + 
        J12L*(J23L*PDstandardNth12gt23 + J33L*PDstandardNth13gt23) + 
        dJ123L*PDstandardNth1gt23 + J23L*(J22L*PDstandardNth22gt23 + 
        J32L*PDstandardNth23gt23) + dJ223L*PDstandardNth2gt23 + 
        J33L*(J22L*PDstandardNth23gt23 + J32L*PDstandardNth33gt23) + 
        dJ323L*PDstandardNth3gt23;
      
      JacPDstandardNth32gt33 = J13L*(J12L*PDstandardNth11gt33 + 
        J22L*PDstandardNth12gt33 + J32L*PDstandardNth13gt33) + 
        J12L*(J23L*PDstandardNth12gt33 + J33L*PDstandardNth13gt33) + 
        dJ123L*PDstandardNth1gt33 + J23L*(J22L*PDstandardNth22gt33 + 
        J32L*PDstandardNth23gt33) + dJ223L*PDstandardNth2gt33 + 
        J33L*(J22L*PDstandardNth23gt33 + J32L*PDstandardNth33gt33) + 
        dJ323L*PDstandardNth3gt33;
      
      JacPDstandardNth32phi = J13L*(J12L*PDstandardNth11phi + 
        J22L*PDstandardNth12phi + J32L*PDstandardNth13phi) + 
        J12L*(J23L*PDstandardNth12phi + J33L*PDstandardNth13phi) + 
        dJ123L*PDstandardNth1phi + J23L*(J22L*PDstandardNth22phi + 
        J32L*PDstandardNth23phi) + dJ223L*PDstandardNth2phi + 
        J33L*(J22L*PDstandardNth23phi + J32L*PDstandardNth33phi) + 
        dJ323L*PDstandardNth3phi;
    }
    else
    {
      JacPDstandardNth1At11 = PDstandardNth1At11;
      
      JacPDstandardNth1At12 = PDstandardNth1At12;
      
      JacPDstandardNth1At13 = PDstandardNth1At13;
      
      JacPDstandardNth1At22 = PDstandardNth1At22;
      
      JacPDstandardNth1At23 = PDstandardNth1At23;
      
      JacPDstandardNth1At33 = PDstandardNth1At33;
      
      JacPDstandardNth1gt11 = PDstandardNth1gt11;
      
      JacPDstandardNth1gt12 = PDstandardNth1gt12;
      
      JacPDstandardNth1gt13 = PDstandardNth1gt13;
      
      JacPDstandardNth1gt22 = PDstandardNth1gt22;
      
      JacPDstandardNth1gt23 = PDstandardNth1gt23;
      
      JacPDstandardNth1gt33 = PDstandardNth1gt33;
      
      JacPDstandardNth1phi = PDstandardNth1phi;
      
      JacPDstandardNth1trK = PDstandardNth1trK;
      
      JacPDstandardNth1Xt1 = PDstandardNth1Xt1;
      
      JacPDstandardNth1Xt2 = PDstandardNth1Xt2;
      
      JacPDstandardNth1Xt3 = PDstandardNth1Xt3;
      
      JacPDstandardNth2At11 = PDstandardNth2At11;
      
      JacPDstandardNth2At12 = PDstandardNth2At12;
      
      JacPDstandardNth2At13 = PDstandardNth2At13;
      
      JacPDstandardNth2At22 = PDstandardNth2At22;
      
      JacPDstandardNth2At23 = PDstandardNth2At23;
      
      JacPDstandardNth2At33 = PDstandardNth2At33;
      
      JacPDstandardNth2gt11 = PDstandardNth2gt11;
      
      JacPDstandardNth2gt12 = PDstandardNth2gt12;
      
      JacPDstandardNth2gt13 = PDstandardNth2gt13;
      
      JacPDstandardNth2gt22 = PDstandardNth2gt22;
      
      JacPDstandardNth2gt23 = PDstandardNth2gt23;
      
      JacPDstandardNth2gt33 = PDstandardNth2gt33;
      
      JacPDstandardNth2phi = PDstandardNth2phi;
      
      JacPDstandardNth2trK = PDstandardNth2trK;
      
      JacPDstandardNth2Xt1 = PDstandardNth2Xt1;
      
      JacPDstandardNth2Xt2 = PDstandardNth2Xt2;
      
      JacPDstandardNth2Xt3 = PDstandardNth2Xt3;
      
      JacPDstandardNth3At11 = PDstandardNth3At11;
      
      JacPDstandardNth3At12 = PDstandardNth3At12;
      
      JacPDstandardNth3At13 = PDstandardNth3At13;
      
      JacPDstandardNth3At22 = PDstandardNth3At22;
      
      JacPDstandardNth3At23 = PDstandardNth3At23;
      
      JacPDstandardNth3At33 = PDstandardNth3At33;
      
      JacPDstandardNth3gt11 = PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = PDstandardNth3gt33;
      
      JacPDstandardNth3phi = PDstandardNth3phi;
      
      JacPDstandardNth3trK = PDstandardNth3trK;
      
      JacPDstandardNth3Xt1 = PDstandardNth3Xt1;
      
      JacPDstandardNth3Xt2 = PDstandardNth3Xt2;
      
      JacPDstandardNth3Xt3 = PDstandardNth3Xt3;
      
      JacPDstandardNth11gt11 = PDstandardNth11gt11;
      
      JacPDstandardNth11gt12 = PDstandardNth11gt12;
      
      JacPDstandardNth11gt13 = PDstandardNth11gt13;
      
      JacPDstandardNth11gt22 = PDstandardNth11gt22;
      
      JacPDstandardNth11gt23 = PDstandardNth11gt23;
      
      JacPDstandardNth11gt33 = PDstandardNth11gt33;
      
      JacPDstandardNth11phi = PDstandardNth11phi;
      
      JacPDstandardNth22gt11 = PDstandardNth22gt11;
      
      JacPDstandardNth22gt12 = PDstandardNth22gt12;
      
      JacPDstandardNth22gt13 = PDstandardNth22gt13;
      
      JacPDstandardNth22gt22 = PDstandardNth22gt22;
      
      JacPDstandardNth22gt23 = PDstandardNth22gt23;
      
      JacPDstandardNth22gt33 = PDstandardNth22gt33;
      
      JacPDstandardNth22phi = PDstandardNth22phi;
      
      JacPDstandardNth33gt11 = PDstandardNth33gt11;
      
      JacPDstandardNth33gt12 = PDstandardNth33gt12;
      
      JacPDstandardNth33gt13 = PDstandardNth33gt13;
      
      JacPDstandardNth33gt22 = PDstandardNth33gt22;
      
      JacPDstandardNth33gt23 = PDstandardNth33gt23;
      
      JacPDstandardNth33gt33 = PDstandardNth33gt33;
      
      JacPDstandardNth33phi = PDstandardNth33phi;
      
      JacPDstandardNth12gt11 = PDstandardNth12gt11;
      
      JacPDstandardNth12gt12 = PDstandardNth12gt12;
      
      JacPDstandardNth12gt13 = PDstandardNth12gt13;
      
      JacPDstandardNth12gt22 = PDstandardNth12gt22;
      
      JacPDstandardNth12gt23 = PDstandardNth12gt23;
      
      JacPDstandardNth12gt33 = PDstandardNth12gt33;
      
      JacPDstandardNth12phi = PDstandardNth12phi;
      
      JacPDstandardNth13gt11 = PDstandardNth13gt11;
      
      JacPDstandardNth13gt12 = PDstandardNth13gt12;
      
      JacPDstandardNth13gt13 = PDstandardNth13gt13;
      
      JacPDstandardNth13gt22 = PDstandardNth13gt22;
      
      JacPDstandardNth13gt23 = PDstandardNth13gt23;
      
      JacPDstandardNth13gt33 = PDstandardNth13gt33;
      
      JacPDstandardNth13phi = PDstandardNth13phi;
      
      JacPDstandardNth21gt11 = PDstandardNth12gt11;
      
      JacPDstandardNth21gt12 = PDstandardNth12gt12;
      
      JacPDstandardNth21gt13 = PDstandardNth12gt13;
      
      JacPDstandardNth21gt22 = PDstandardNth12gt22;
      
      JacPDstandardNth21gt23 = PDstandardNth12gt23;
      
      JacPDstandardNth21gt33 = PDstandardNth12gt33;
      
      JacPDstandardNth21phi = PDstandardNth12phi;
      
      JacPDstandardNth23gt11 = PDstandardNth23gt11;
      
      JacPDstandardNth23gt12 = PDstandardNth23gt12;
      
      JacPDstandardNth23gt13 = PDstandardNth23gt13;
      
      JacPDstandardNth23gt22 = PDstandardNth23gt22;
      
      JacPDstandardNth23gt23 = PDstandardNth23gt23;
      
      JacPDstandardNth23gt33 = PDstandardNth23gt33;
      
      JacPDstandardNth23phi = PDstandardNth23phi;
      
      JacPDstandardNth31gt11 = PDstandardNth13gt11;
      
      JacPDstandardNth31gt12 = PDstandardNth13gt12;
      
      JacPDstandardNth31gt13 = PDstandardNth13gt13;
      
      JacPDstandardNth31gt22 = PDstandardNth13gt22;
      
      JacPDstandardNth31gt23 = PDstandardNth13gt23;
      
      JacPDstandardNth31gt33 = PDstandardNth13gt33;
      
      JacPDstandardNth31phi = PDstandardNth13phi;
      
      JacPDstandardNth32gt11 = PDstandardNth23gt11;
      
      JacPDstandardNth32gt12 = PDstandardNth23gt12;
      
      JacPDstandardNth32gt13 = PDstandardNth23gt13;
      
      JacPDstandardNth32gt22 = PDstandardNth23gt22;
      
      JacPDstandardNth32gt23 = PDstandardNth23gt23;
      
      JacPDstandardNth32gt33 = PDstandardNth23gt33;
      
      JacPDstandardNth32phi = PDstandardNth23phi;
    }
    
    CCTK_REAL detgt CCTK_ATTRIBUTE_UNUSED = 1;
    
    CCTK_REAL gtu11 CCTK_ATTRIBUTE_UNUSED = (gt22L*gt33L - 
      pow(gt23L,2))*pow(detgt,-1);
    
    CCTK_REAL gtu12 CCTK_ATTRIBUTE_UNUSED = (gt13L*gt23L - 
      gt12L*gt33L)*pow(detgt,-1);
    
    CCTK_REAL gtu13 CCTK_ATTRIBUTE_UNUSED = (-(gt13L*gt22L) + 
      gt12L*gt23L)*pow(detgt,-1);
    
    CCTK_REAL gtu22 CCTK_ATTRIBUTE_UNUSED = (gt11L*gt33L - 
      pow(gt13L,2))*pow(detgt,-1);
    
    CCTK_REAL gtu23 CCTK_ATTRIBUTE_UNUSED = (gt12L*gt13L - 
      gt11L*gt23L)*pow(detgt,-1);
    
    CCTK_REAL gtu33 CCTK_ATTRIBUTE_UNUSED = (gt11L*gt22L - 
      pow(gt12L,2))*pow(detgt,-1);
    
    CCTK_REAL Gtl111 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth1gt11;
    
    CCTK_REAL Gtl112 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth2gt11;
    
    CCTK_REAL Gtl113 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth3gt11;
    
    CCTK_REAL Gtl122 CCTK_ATTRIBUTE_UNUSED = -0.5*JacPDstandardNth1gt22 + 
      JacPDstandardNth2gt12;
    
    CCTK_REAL Gtl123 CCTK_ATTRIBUTE_UNUSED = 0.5*(-JacPDstandardNth1gt23 + 
      JacPDstandardNth2gt13 + JacPDstandardNth3gt12);
    
    CCTK_REAL Gtl133 CCTK_ATTRIBUTE_UNUSED = -0.5*JacPDstandardNth1gt33 + 
      JacPDstandardNth3gt13;
    
    CCTK_REAL Gtl211 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1gt12 - 
      0.5*JacPDstandardNth2gt11;
    
    CCTK_REAL Gtl212 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth1gt22;
    
    CCTK_REAL Gtl213 CCTK_ATTRIBUTE_UNUSED = 0.5*(JacPDstandardNth1gt23 - 
      JacPDstandardNth2gt13 + JacPDstandardNth3gt12);
    
    CCTK_REAL Gtl222 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth2gt22;
    
    CCTK_REAL Gtl223 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth3gt22;
    
    CCTK_REAL Gtl233 CCTK_ATTRIBUTE_UNUSED = -0.5*JacPDstandardNth2gt33 + 
      JacPDstandardNth3gt23;
    
    CCTK_REAL Gtl311 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1gt13 - 
      0.5*JacPDstandardNth3gt11;
    
    CCTK_REAL Gtl312 CCTK_ATTRIBUTE_UNUSED = 0.5*(JacPDstandardNth1gt23 + 
      JacPDstandardNth2gt13 - JacPDstandardNth3gt12);
    
    CCTK_REAL Gtl313 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth1gt33;
    
    CCTK_REAL Gtl322 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth2gt23 - 
      0.5*JacPDstandardNth3gt22;
    
    CCTK_REAL Gtl323 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth2gt33;
    
    CCTK_REAL Gtl333 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth3gt33;
    
    CCTK_REAL Gtlu111 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu11 + Gtl112*gtu12 
      + Gtl113*gtu13;
    
    CCTK_REAL Gtlu112 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu12 + Gtl112*gtu22 
      + Gtl113*gtu23;
    
    CCTK_REAL Gtlu113 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu13 + Gtl112*gtu23 
      + Gtl113*gtu33;
    
    CCTK_REAL Gtlu121 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu11 + Gtl122*gtu12 
      + Gtl123*gtu13;
    
    CCTK_REAL Gtlu122 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu12 + Gtl122*gtu22 
      + Gtl123*gtu23;
    
    CCTK_REAL Gtlu123 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu13 + Gtl122*gtu23 
      + Gtl123*gtu33;
    
    CCTK_REAL Gtlu131 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu11 + Gtl123*gtu12 
      + Gtl133*gtu13;
    
    CCTK_REAL Gtlu132 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu12 + Gtl123*gtu22 
      + Gtl133*gtu23;
    
    CCTK_REAL Gtlu133 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu13 + Gtl123*gtu23 
      + Gtl133*gtu33;
    
    CCTK_REAL Gtlu211 CCTK_ATTRIBUTE_UNUSED = Gtl211*gtu11 + Gtl212*gtu12 
      + Gtl213*gtu13;
    
    CCTK_REAL Gtlu212 CCTK_ATTRIBUTE_UNUSED = Gtl211*gtu12 + Gtl212*gtu22 
      + Gtl213*gtu23;
    
    CCTK_REAL Gtlu213 CCTK_ATTRIBUTE_UNUSED = Gtl211*gtu13 + Gtl212*gtu23 
      + Gtl213*gtu33;
    
    CCTK_REAL Gtlu221 CCTK_ATTRIBUTE_UNUSED = Gtl212*gtu11 + Gtl222*gtu12 
      + Gtl223*gtu13;
    
    CCTK_REAL Gtlu222 CCTK_ATTRIBUTE_UNUSED = Gtl212*gtu12 + Gtl222*gtu22 
      + Gtl223*gtu23;
    
    CCTK_REAL Gtlu223 CCTK_ATTRIBUTE_UNUSED = Gtl212*gtu13 + Gtl222*gtu23 
      + Gtl223*gtu33;
    
    CCTK_REAL Gtlu231 CCTK_ATTRIBUTE_UNUSED = Gtl213*gtu11 + Gtl223*gtu12 
      + Gtl233*gtu13;
    
    CCTK_REAL Gtlu232 CCTK_ATTRIBUTE_UNUSED = Gtl213*gtu12 + Gtl223*gtu22 
      + Gtl233*gtu23;
    
    CCTK_REAL Gtlu233 CCTK_ATTRIBUTE_UNUSED = Gtl213*gtu13 + Gtl223*gtu23 
      + Gtl233*gtu33;
    
    CCTK_REAL Gtlu311 CCTK_ATTRIBUTE_UNUSED = Gtl311*gtu11 + Gtl312*gtu12 
      + Gtl313*gtu13;
    
    CCTK_REAL Gtlu312 CCTK_ATTRIBUTE_UNUSED = Gtl311*gtu12 + Gtl312*gtu22 
      + Gtl313*gtu23;
    
    CCTK_REAL Gtlu313 CCTK_ATTRIBUTE_UNUSED = Gtl311*gtu13 + Gtl312*gtu23 
      + Gtl313*gtu33;
    
    CCTK_REAL Gtlu321 CCTK_ATTRIBUTE_UNUSED = Gtl312*gtu11 + Gtl322*gtu12 
      + Gtl323*gtu13;
    
    CCTK_REAL Gtlu322 CCTK_ATTRIBUTE_UNUSED = Gtl312*gtu12 + Gtl322*gtu22 
      + Gtl323*gtu23;
    
    CCTK_REAL Gtlu323 CCTK_ATTRIBUTE_UNUSED = Gtl312*gtu13 + Gtl322*gtu23 
      + Gtl323*gtu33;
    
    CCTK_REAL Gtlu331 CCTK_ATTRIBUTE_UNUSED = Gtl313*gtu11 + Gtl323*gtu12 
      + Gtl333*gtu13;
    
    CCTK_REAL Gtlu332 CCTK_ATTRIBUTE_UNUSED = Gtl313*gtu12 + Gtl323*gtu22 
      + Gtl333*gtu23;
    
    CCTK_REAL Gtlu333 CCTK_ATTRIBUTE_UNUSED = Gtl313*gtu13 + Gtl323*gtu23 
      + Gtl333*gtu33;
    
    CCTK_REAL Gt111 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu11 + Gtl211*gtu12 + 
      Gtl311*gtu13;
    
    CCTK_REAL Gt211 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu12 + Gtl211*gtu22 + 
      Gtl311*gtu23;
    
    CCTK_REAL Gt311 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu13 + Gtl211*gtu23 + 
      Gtl311*gtu33;
    
    CCTK_REAL Gt112 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu11 + Gtl212*gtu12 + 
      Gtl312*gtu13;
    
    CCTK_REAL Gt212 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu12 + Gtl212*gtu22 + 
      Gtl312*gtu23;
    
    CCTK_REAL Gt312 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu13 + Gtl212*gtu23 + 
      Gtl312*gtu33;
    
    CCTK_REAL Gt113 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu11 + Gtl213*gtu12 + 
      Gtl313*gtu13;
    
    CCTK_REAL Gt213 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu12 + Gtl213*gtu22 + 
      Gtl313*gtu23;
    
    CCTK_REAL Gt313 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu13 + Gtl213*gtu23 + 
      Gtl313*gtu33;
    
    CCTK_REAL Gt122 CCTK_ATTRIBUTE_UNUSED = Gtl122*gtu11 + Gtl222*gtu12 + 
      Gtl322*gtu13;
    
    CCTK_REAL Gt222 CCTK_ATTRIBUTE_UNUSED = Gtl122*gtu12 + Gtl222*gtu22 + 
      Gtl322*gtu23;
    
    CCTK_REAL Gt322 CCTK_ATTRIBUTE_UNUSED = Gtl122*gtu13 + Gtl222*gtu23 + 
      Gtl322*gtu33;
    
    CCTK_REAL Gt123 CCTK_ATTRIBUTE_UNUSED = Gtl123*gtu11 + Gtl223*gtu12 + 
      Gtl323*gtu13;
    
    CCTK_REAL Gt223 CCTK_ATTRIBUTE_UNUSED = Gtl123*gtu12 + Gtl223*gtu22 + 
      Gtl323*gtu23;
    
    CCTK_REAL Gt323 CCTK_ATTRIBUTE_UNUSED = Gtl123*gtu13 + Gtl223*gtu23 + 
      Gtl323*gtu33;
    
    CCTK_REAL Gt133 CCTK_ATTRIBUTE_UNUSED = Gtl133*gtu11 + Gtl233*gtu12 + 
      Gtl333*gtu13;
    
    CCTK_REAL Gt233 CCTK_ATTRIBUTE_UNUSED = Gtl133*gtu12 + Gtl233*gtu22 + 
      Gtl333*gtu23;
    
    CCTK_REAL Gt333 CCTK_ATTRIBUTE_UNUSED = Gtl133*gtu13 + Gtl233*gtu23 + 
      Gtl333*gtu33;
    
    CCTK_REAL em4phi CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod != 
      0,pow(phiL,2),exp(-4*phiL));
    
    CCTK_REAL e4phi CCTK_ATTRIBUTE_UNUSED = pow(em4phi,-1);
    
    CCTK_REAL g11 CCTK_ATTRIBUTE_UNUSED = gt11L*e4phi;
    
    CCTK_REAL g12 CCTK_ATTRIBUTE_UNUSED = gt12L*e4phi;
    
    CCTK_REAL g13 CCTK_ATTRIBUTE_UNUSED = gt13L*e4phi;
    
    CCTK_REAL g22 CCTK_ATTRIBUTE_UNUSED = gt22L*e4phi;
    
    CCTK_REAL g23 CCTK_ATTRIBUTE_UNUSED = gt23L*e4phi;
    
    CCTK_REAL g33 CCTK_ATTRIBUTE_UNUSED = gt33L*e4phi;
    
    CCTK_REAL gu11 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu11;
    
    CCTK_REAL gu12 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu12;
    
    CCTK_REAL gu13 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu13;
    
    CCTK_REAL gu22 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu22;
    
    CCTK_REAL gu23 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu23;
    
    CCTK_REAL gu33 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu33;
    
    CCTK_REAL Xtn1 CCTK_ATTRIBUTE_UNUSED = Gt111*gtu11 + Gt122*gtu22 + 
      2*(Gt112*gtu12 + Gt113*gtu13 + Gt123*gtu23) + Gt133*gtu33;
    
    CCTK_REAL Xtn2 CCTK_ATTRIBUTE_UNUSED = Gt211*gtu11 + Gt222*gtu22 + 
      2*(Gt212*gtu12 + Gt213*gtu13 + Gt223*gtu23) + Gt233*gtu33;
    
    CCTK_REAL Xtn3 CCTK_ATTRIBUTE_UNUSED = Gt311*gtu11 + Gt322*gtu22 + 
      2*(Gt312*gtu12 + Gt313*gtu13 + Gt323*gtu23) + Gt333*gtu33;
    
    CCTK_REAL Zl1 CCTK_ATTRIBUTE_UNUSED = 0.5*(gt11L*Xt1L + gt12L*Xt2L + 
      gt13L*Xt3L - gtu11*JacPDstandardNth1gt11 - gtu12*(JacPDstandardNth1gt12 
      + JacPDstandardNth2gt11) - gtu22*JacPDstandardNth2gt12 + 
      gtu13*(-JacPDstandardNth1gt13 - JacPDstandardNth3gt11) + 
      gtu23*(-JacPDstandardNth2gt13 - JacPDstandardNth3gt12) - 
      gtu33*JacPDstandardNth3gt13);
    
    CCTK_REAL Zl2 CCTK_ATTRIBUTE_UNUSED = 0.5*(gt12L*Xt1L + gt22L*Xt2L + 
      gt23L*Xt3L - gtu11*JacPDstandardNth1gt12 - gtu12*(JacPDstandardNth1gt22 
      + JacPDstandardNth2gt12) - gtu22*JacPDstandardNth2gt22 + 
      gtu13*(-JacPDstandardNth1gt23 - JacPDstandardNth3gt12) + 
      gtu23*(-JacPDstandardNth2gt23 - JacPDstandardNth3gt22) - 
      gtu33*JacPDstandardNth3gt23);
    
    CCTK_REAL Zl3 CCTK_ATTRIBUTE_UNUSED = 0.5*(gt13L*Xt1L + gt23L*Xt2L + 
      gt33L*Xt3L - gtu11*JacPDstandardNth1gt13 - gtu12*(JacPDstandardNth1gt23 
      + JacPDstandardNth2gt13) - gtu22*JacPDstandardNth2gt23 + 
      gtu13*(-JacPDstandardNth1gt33 - JacPDstandardNth3gt13) + 
      gtu23*(-JacPDstandardNth2gt33 - JacPDstandardNth3gt23) - 
      gtu33*JacPDstandardNth3gt33);
    
    CCTK_REAL Z1 CCTK_ATTRIBUTE_UNUSED = gu11*Zl1 + gu12*Zl2 + gu13*Zl3;
    
    CCTK_REAL Z2 CCTK_ATTRIBUTE_UNUSED = gu12*Zl1 + gu22*Zl2 + gu23*Zl3;
    
    CCTK_REAL Z3 CCTK_ATTRIBUTE_UNUSED = gu13*Zl1 + gu23*Zl2 + gu33*Zl3;
    
    CCTK_REAL Rt11 CCTK_ATTRIBUTE_UNUSED = 3*(Gt111*Gtlu111 + 
      Gt112*Gtlu112 + Gt113*Gtlu113) + 2*(Gt211*Gtlu121 + Gt212*Gtlu122 + 
      Gt213*Gtlu123 + Gt311*Gtlu131 + Gt312*Gtlu132 + Gt313*Gtlu133) + 
      Gt211*Gtlu211 + Gt212*Gtlu212 + Gt213*Gtlu213 + Gt311*Gtlu311 + 
      Gt312*Gtlu312 + Gt313*Gtlu313 + gt11L*JacPDstandardNth1Xt1 + 
      gt12L*JacPDstandardNth1Xt2 + gt13L*JacPDstandardNth1Xt3 + 
      0.5*(-(gtu11*JacPDstandardNth11gt11) - gtu12*(JacPDstandardNth12gt11 + 
      JacPDstandardNth21gt11) - gtu22*JacPDstandardNth22gt11 + 
      gtu13*(-JacPDstandardNth13gt11 - JacPDstandardNth31gt11) + 
      gtu23*(-JacPDstandardNth23gt11 - JacPDstandardNth32gt11) - 
      gtu33*JacPDstandardNth33gt11) + Gtl111*Xtn1 + Gtl112*Xtn2 + 
      Gtl113*Xtn3;
    
    CCTK_REAL Rt12 CCTK_ATTRIBUTE_UNUSED = 0.5*(4*(Gt211*Gtlu221 + 
      Gt212*Gtlu222 + Gt213*Gtlu223) + 2*(Gt122*Gtlu112 + Gt123*Gtlu113 + 
      Gt111*Gtlu121 + Gt212*Gtlu121 + Gt222*Gtlu122 + Gt113*Gtlu123 + 
      Gt223*Gtlu123 + Gt312*Gtlu131 + Gt322*Gtlu132 + Gt323*Gtlu133 + 
      Gt111*Gtlu211 + Gt112*(Gtlu111 + Gtlu122 + Gtlu212) + Gt113*Gtlu213 + 
      Gt311*Gtlu231 + Gt312*Gtlu232 + Gt313*Gtlu233 + Gt311*Gtlu321 + 
      Gt312*Gtlu322 + Gt313*Gtlu323) - gtu11*JacPDstandardNth11gt12 + 
      gt12L*JacPDstandardNth1Xt1 + gt22L*JacPDstandardNth1Xt2 + 
      gt23L*JacPDstandardNth1Xt3 + gtu12*(-JacPDstandardNth12gt12 - 
      JacPDstandardNth21gt12) - gtu22*JacPDstandardNth22gt12 + 
      gt11L*JacPDstandardNth2Xt1 + gt12L*JacPDstandardNth2Xt2 + 
      gt13L*JacPDstandardNth2Xt3 + gtu13*(-JacPDstandardNth13gt12 - 
      JacPDstandardNth31gt12) + gtu23*(-JacPDstandardNth23gt12 - 
      JacPDstandardNth32gt12) - gtu33*JacPDstandardNth33gt12 + Gtl112*Xtn1 + 
      Gtl211*Xtn1 + Gtl122*Xtn2 + Gtl212*Xtn2 + Gtl123*Xtn3 + Gtl213*Xtn3);
    
    CCTK_REAL Rt13 CCTK_ATTRIBUTE_UNUSED = 0.5*(2*(Gt123*Gtlu112 + 
      Gt133*Gtlu113 + Gt213*Gtlu121 + Gt223*Gtlu122 + Gt233*Gtlu123 + 
      Gt111*Gtlu131 + Gt313*Gtlu131 + Gt112*Gtlu132 + Gt323*Gtlu132 + 
      Gt333*Gtlu133 + Gt211*Gtlu231 + Gt212*Gtlu232 + Gt213*Gtlu233 + 
      Gt111*Gtlu311 + Gt112*Gtlu312 + Gt113*(Gtlu111 + Gtlu133 + Gtlu313) + 
      Gt211*Gtlu321 + Gt212*Gtlu322 + Gt213*Gtlu323) + 4*(Gt311*Gtlu331 + 
      Gt312*Gtlu332 + Gt313*Gtlu333) - gtu11*JacPDstandardNth11gt13 + 
      gt13L*JacPDstandardNth1Xt1 + gt23L*JacPDstandardNth1Xt2 + 
      gt33L*JacPDstandardNth1Xt3 + gtu12*(-JacPDstandardNth12gt13 - 
      JacPDstandardNth21gt13) - gtu22*JacPDstandardNth22gt13 + 
      gtu13*(-JacPDstandardNth13gt13 - JacPDstandardNth31gt13) + 
      gtu23*(-JacPDstandardNth23gt13 - JacPDstandardNth32gt13) - 
      gtu33*JacPDstandardNth33gt13 + gt11L*JacPDstandardNth3Xt1 + 
      gt12L*JacPDstandardNth3Xt2 + gt13L*JacPDstandardNth3Xt3 + Gtl113*Xtn1 + 
      Gtl311*Xtn1 + Gtl123*Xtn2 + Gtl312*Xtn2 + Gtl133*Xtn3 + Gtl313*Xtn3);
    
    CCTK_REAL Rt22 CCTK_ATTRIBUTE_UNUSED = Gt112*(Gtlu121 + 2*Gtlu211) + 
      Gt122*(Gtlu122 + 2*Gtlu212) + Gt123*(Gtlu123 + 2*Gtlu213) + 
      3*(Gt212*Gtlu221 + Gt222*Gtlu222 + Gt223*Gtlu223) + 2*(Gt312*Gtlu231 + 
      Gt322*Gtlu232 + Gt323*Gtlu233) + Gt312*Gtlu321 + Gt322*Gtlu322 + 
      Gt323*Gtlu323 + gt12L*JacPDstandardNth2Xt1 + gt22L*JacPDstandardNth2Xt2 
      + gt23L*JacPDstandardNth2Xt3 + 0.5*(-(gtu11*JacPDstandardNth11gt22) - 
      gtu12*(JacPDstandardNth12gt22 + JacPDstandardNth21gt22) - 
      gtu22*JacPDstandardNth22gt22 + gtu13*(-JacPDstandardNth13gt22 - 
      JacPDstandardNth31gt22) + gtu23*(-JacPDstandardNth23gt22 - 
      JacPDstandardNth32gt22) - gtu33*JacPDstandardNth33gt22) + Gtl212*Xtn1 + 
      Gtl222*Xtn2 + Gtl223*Xtn3;
    
    CCTK_REAL Rt23 CCTK_ATTRIBUTE_UNUSED = 0.5*(2*(Gt123*Gtlu133 + 
      Gt113*Gtlu211 + Gt123*Gtlu212 + Gt133*Gtlu213 + Gt213*Gtlu221 + 
      Gt223*Gtlu222 + Gt233*Gtlu223 + Gt212*Gtlu231 + Gt313*Gtlu231 + 
      Gt222*Gtlu232 + Gt323*Gtlu232 + Gt223*Gtlu233 + Gt333*Gtlu233 + 
      Gt112*(Gtlu131 + Gtlu311) + Gt122*(Gtlu132 + Gtlu312) + Gt123*Gtlu313 + 
      Gt212*Gtlu321 + Gt222*Gtlu322 + Gt223*Gtlu323) + 4*(Gt312*Gtlu331 + 
      Gt322*Gtlu332 + Gt323*Gtlu333) - gtu11*JacPDstandardNth11gt23 + 
      gtu12*(-JacPDstandardNth12gt23 - JacPDstandardNth21gt23) - 
      gtu22*JacPDstandardNth22gt23 + gt13L*JacPDstandardNth2Xt1 + 
      gt23L*JacPDstandardNth2Xt2 + gt33L*JacPDstandardNth2Xt3 + 
      gtu13*(-JacPDstandardNth13gt23 - JacPDstandardNth31gt23) + 
      gtu23*(-JacPDstandardNth23gt23 - JacPDstandardNth32gt23) - 
      gtu33*JacPDstandardNth33gt23 + gt12L*JacPDstandardNth3Xt1 + 
      gt22L*JacPDstandardNth3Xt2 + gt23L*JacPDstandardNth3Xt3 + Gtl213*Xtn1 + 
      Gtl312*Xtn1 + Gtl223*Xtn2 + Gtl322*Xtn2 + Gtl233*Xtn3 + Gtl323*Xtn3);
    
    CCTK_REAL Rt33 CCTK_ATTRIBUTE_UNUSED = Gt113*(Gtlu131 + 2*Gtlu311) + 
      Gt123*(Gtlu132 + 2*Gtlu312) + Gt133*(Gtlu133 + 2*Gtlu313) + 
      Gt213*(Gtlu231 + 2*Gtlu321) + Gt223*(Gtlu232 + 2*Gtlu322) + 
      Gt233*(Gtlu233 + 2*Gtlu323) + 3*(Gt313*Gtlu331 + Gt323*Gtlu332 + 
      Gt333*Gtlu333) + 0.5*(-(gtu11*JacPDstandardNth11gt33) - 
      gtu12*(JacPDstandardNth12gt33 + JacPDstandardNth21gt33) - 
      gtu22*JacPDstandardNth22gt33 + gtu13*(-JacPDstandardNth13gt33 - 
      JacPDstandardNth31gt33) + gtu23*(-JacPDstandardNth23gt33 - 
      JacPDstandardNth32gt33) - gtu33*JacPDstandardNth33gt33) + 
      gt13L*JacPDstandardNth3Xt1 + gt23L*JacPDstandardNth3Xt2 + 
      gt33L*JacPDstandardNth3Xt3 + Gtl313*Xtn1 + Gtl323*Xtn2 + Gtl333*Xtn3;
    
    CCTK_REAL fac1 CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod != 
      0,-0.5*pow(phiL,-1),1);
    
    CCTK_REAL cdphi1 CCTK_ATTRIBUTE_UNUSED = fac1*JacPDstandardNth1phi;
    
    CCTK_REAL cdphi2 CCTK_ATTRIBUTE_UNUSED = fac1*JacPDstandardNth2phi;
    
    CCTK_REAL cdphi3 CCTK_ATTRIBUTE_UNUSED = fac1*JacPDstandardNth3phi;
    
    CCTK_REAL fac2 CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod != 
      0,0.5*pow(phiL,-2),0);
    
    CCTK_REAL cdphi211 CCTK_ATTRIBUTE_UNUSED = fac1*(JacPDstandardNth11phi 
      - Gt111*JacPDstandardNth1phi - Gt211*JacPDstandardNth2phi - 
      Gt311*JacPDstandardNth3phi) + fac2*pow(JacPDstandardNth1phi,2);
    
    CCTK_REAL cdphi212 CCTK_ATTRIBUTE_UNUSED = 
      fac2*JacPDstandardNth1phi*JacPDstandardNth2phi + 
      fac1*(JacPDstandardNth12phi - Gt112*JacPDstandardNth1phi - 
      Gt212*JacPDstandardNth2phi - Gt312*JacPDstandardNth3phi);
    
    CCTK_REAL cdphi213 CCTK_ATTRIBUTE_UNUSED = 
      fac2*JacPDstandardNth1phi*JacPDstandardNth3phi + 
      fac1*(JacPDstandardNth13phi - Gt113*JacPDstandardNth1phi - 
      Gt213*JacPDstandardNth2phi - Gt313*JacPDstandardNth3phi);
    
    CCTK_REAL cdphi221 CCTK_ATTRIBUTE_UNUSED = 
      fac2*JacPDstandardNth1phi*JacPDstandardNth2phi - 
      fac1*(Gt112*JacPDstandardNth1phi - JacPDstandardNth21phi + 
      Gt212*JacPDstandardNth2phi + Gt312*JacPDstandardNth3phi);
    
    CCTK_REAL cdphi222 CCTK_ATTRIBUTE_UNUSED = 
      -(fac1*(Gt122*JacPDstandardNth1phi - JacPDstandardNth22phi + 
      Gt222*JacPDstandardNth2phi + Gt322*JacPDstandardNth3phi)) + 
      fac2*pow(JacPDstandardNth2phi,2);
    
    CCTK_REAL cdphi223 CCTK_ATTRIBUTE_UNUSED = 
      fac2*JacPDstandardNth2phi*JacPDstandardNth3phi - 
      fac1*(Gt123*JacPDstandardNth1phi - JacPDstandardNth23phi + 
      Gt223*JacPDstandardNth2phi + Gt323*JacPDstandardNth3phi);
    
    CCTK_REAL cdphi231 CCTK_ATTRIBUTE_UNUSED = 
      fac2*JacPDstandardNth1phi*JacPDstandardNth3phi - 
      fac1*(Gt113*JacPDstandardNth1phi + Gt213*JacPDstandardNth2phi - 
      JacPDstandardNth31phi + Gt313*JacPDstandardNth3phi);
    
    CCTK_REAL cdphi232 CCTK_ATTRIBUTE_UNUSED = 
      fac2*JacPDstandardNth2phi*JacPDstandardNth3phi - 
      fac1*(Gt123*JacPDstandardNth1phi + Gt223*JacPDstandardNth2phi - 
      JacPDstandardNth32phi + Gt323*JacPDstandardNth3phi);
    
    CCTK_REAL cdphi233 CCTK_ATTRIBUTE_UNUSED = 
      -(fac1*(Gt133*JacPDstandardNth1phi + Gt233*JacPDstandardNth2phi - 
      JacPDstandardNth33phi + Gt333*JacPDstandardNth3phi)) + 
      fac2*pow(JacPDstandardNth3phi,2);
    
    CCTK_REAL Rphi11 CCTK_ATTRIBUTE_UNUSED = -8*gt11L*cdphi1*(cdphi2*gtu12 
      + cdphi3*gtu13) + (4 - 4*gt11L*gtu11)*pow(cdphi1,2) - 2*(cdphi211 + 
      gt11L*(cdphi211*gtu11 + (cdphi212 + cdphi221)*gtu12 + (cdphi213 + 
      cdphi231)*gtu13 + (cdphi223 + cdphi232 + 4*cdphi2*cdphi3)*gtu23 + 
      gtu22*(cdphi222 + 2*pow(cdphi2,2)) + gtu33*(cdphi233 + 
      2*pow(cdphi3,2))));
    
    CCTK_REAL Rphi12 CCTK_ATTRIBUTE_UNUSED = 4*cdphi1*cdphi2 - 2*(cdphi221 
      + gt12L*(cdphi211*gtu11 + (cdphi212 + cdphi221)*gtu12 + (cdphi213 + 
      cdphi231)*gtu13 + cdphi222*gtu22 + (cdphi223 + cdphi232)*gtu23 + 
      cdphi233*gtu33)) - 4*gt12L*(2*(cdphi1*(cdphi2*gtu12 + cdphi3*gtu13) + 
      cdphi2*cdphi3*gtu23) + gtu11*pow(cdphi1,2) + gtu22*pow(cdphi2,2) + 
      gtu33*pow(cdphi3,2));
    
    CCTK_REAL Rphi13 CCTK_ATTRIBUTE_UNUSED = 4*cdphi1*cdphi3 - 2*(cdphi231 
      + gt13L*(cdphi211*gtu11 + (cdphi212 + cdphi221)*gtu12 + (cdphi213 + 
      cdphi231)*gtu13 + cdphi222*gtu22 + (cdphi223 + cdphi232)*gtu23 + 
      cdphi233*gtu33)) - 4*gt13L*(2*(cdphi1*(cdphi2*gtu12 + cdphi3*gtu13) + 
      cdphi2*cdphi3*gtu23) + gtu11*pow(cdphi1,2) + gtu22*pow(cdphi2,2) + 
      gtu33*pow(cdphi3,2));
    
    CCTK_REAL Rphi22 CCTK_ATTRIBUTE_UNUSED = -8*gt22L*cdphi2*(cdphi1*gtu12 
      + cdphi3*gtu23) + (4 - 4*gt22L*gtu22)*pow(cdphi2,2) - 2*(cdphi222 + 
      gt22L*((cdphi212 + cdphi221)*gtu12 + (cdphi213 + cdphi231 + 
      4*cdphi1*cdphi3)*gtu13 + cdphi222*gtu22 + (cdphi223 + cdphi232)*gtu23 + 
      gtu11*(cdphi211 + 2*pow(cdphi1,2)) + gtu33*(cdphi233 + 
      2*pow(cdphi3,2))));
    
    CCTK_REAL Rphi23 CCTK_ATTRIBUTE_UNUSED = 4*cdphi2*cdphi3 - 2*(cdphi232 
      + gt23L*(cdphi211*gtu11 + (cdphi212 + cdphi221)*gtu12 + (cdphi213 + 
      cdphi231)*gtu13 + cdphi222*gtu22 + (cdphi223 + cdphi232)*gtu23 + 
      cdphi233*gtu33)) - 4*gt23L*(2*(cdphi1*(cdphi2*gtu12 + cdphi3*gtu13) + 
      cdphi2*cdphi3*gtu23) + gtu11*pow(cdphi1,2) + gtu22*pow(cdphi2,2) + 
      gtu33*pow(cdphi3,2));
    
    CCTK_REAL Rphi33 CCTK_ATTRIBUTE_UNUSED = -2*(cdphi233 + 
      gt33L*(cdphi211*gtu11 + (cdphi212 + cdphi221)*gtu12 + (cdphi213 + 
      cdphi231)*gtu13 + cdphi222*gtu22 + (cdphi223 + cdphi232)*gtu23 + 
      cdphi233*gtu33)) + 4*pow(cdphi3,2) - 4*gt33L*(2*(cdphi1*(cdphi2*gtu12 + 
      cdphi3*gtu13) + cdphi2*cdphi3*gtu23) + gtu11*pow(cdphi1,2) + 
      gtu22*pow(cdphi2,2) + gtu33*pow(cdphi3,2));
    
    CCTK_REAL Atm11 CCTK_ATTRIBUTE_UNUSED = At11L*gtu11 + At12L*gtu12 + 
      At13L*gtu13;
    
    CCTK_REAL Atm21 CCTK_ATTRIBUTE_UNUSED = At11L*gtu12 + At12L*gtu22 + 
      At13L*gtu23;
    
    CCTK_REAL Atm31 CCTK_ATTRIBUTE_UNUSED = At11L*gtu13 + At12L*gtu23 + 
      At13L*gtu33;
    
    CCTK_REAL Atm12 CCTK_ATTRIBUTE_UNUSED = At12L*gtu11 + At22L*gtu12 + 
      At23L*gtu13;
    
    CCTK_REAL Atm22 CCTK_ATTRIBUTE_UNUSED = At12L*gtu12 + At22L*gtu22 + 
      At23L*gtu23;
    
    CCTK_REAL Atm32 CCTK_ATTRIBUTE_UNUSED = At12L*gtu13 + At22L*gtu23 + 
      At23L*gtu33;
    
    CCTK_REAL Atm13 CCTK_ATTRIBUTE_UNUSED = At13L*gtu11 + At23L*gtu12 + 
      At33L*gtu13;
    
    CCTK_REAL Atm23 CCTK_ATTRIBUTE_UNUSED = At13L*gtu12 + At23L*gtu22 + 
      At33L*gtu23;
    
    CCTK_REAL Atm33 CCTK_ATTRIBUTE_UNUSED = At13L*gtu13 + At23L*gtu23 + 
      At33L*gtu33;
    
    CCTK_REAL R11 CCTK_ATTRIBUTE_UNUSED = Rphi11 + Rt11 + 
      IfThen(formulation != 0,(4*JacPDstandardNth1phi*(g12*Z2 + g13*Z3) + 
      phiL*e4phi*(JacPDstandardNth1gt11*Z1 + JacPDstandardNth2gt11*Z2 + 
      JacPDstandardNth3gt11*Z3) + 2*g11*(JacPDstandardNth1phi*Z1 - 
      JacPDstandardNth2phi*Z2 - JacPDstandardNth3phi*Z3))*pow(phiL,-1),0);
    
    CCTK_REAL R12 CCTK_ATTRIBUTE_UNUSED = Rphi12 + Rt12 + 
      IfThen(formulation != 0,(phiL*e4phi*(JacPDstandardNth1gt12*Z1 + 
      JacPDstandardNth2gt12*Z2 + JacPDstandardNth3gt12*Z3) + 
      2*(g11*JacPDstandardNth2phi*Z1 + g22*JacPDstandardNth1phi*Z2 + 
      (g23*JacPDstandardNth1phi + g13*JacPDstandardNth2phi - 
      g12*JacPDstandardNth3phi)*Z3))*pow(phiL,-1),0);
    
    CCTK_REAL R13 CCTK_ATTRIBUTE_UNUSED = Rphi13 + Rt13 + 
      IfThen(formulation != 0,(2*((g23*JacPDstandardNth1phi - 
      g13*JacPDstandardNth2phi)*Z2 + JacPDstandardNth3phi*(g11*Z1 + g12*Z2) + 
      g33*JacPDstandardNth1phi*Z3) + phiL*e4phi*(JacPDstandardNth1gt13*Z1 + 
      JacPDstandardNth2gt13*Z2 + JacPDstandardNth3gt13*Z3))*pow(phiL,-1),0);
    
    CCTK_REAL R22 CCTK_ATTRIBUTE_UNUSED = Rphi22 + Rt22 + 
      IfThen(formulation != 0,e4phi*(JacPDstandardNth1gt22*Z1 + 
      JacPDstandardNth2gt22*Z2 + JacPDstandardNth3gt22*Z3) + 
      (4*JacPDstandardNth2phi*(g12*Z1 + g23*Z3) - 
      2*g22*(JacPDstandardNth1phi*Z1 - JacPDstandardNth2phi*Z2 + 
      JacPDstandardNth3phi*Z3))*pow(phiL,-1),0);
    
    CCTK_REAL R23 CCTK_ATTRIBUTE_UNUSED = Rphi23 + Rt23 + 
      IfThen(formulation != 0,((phiL*e4phi*JacPDstandardNth1gt23 - 
      2*g23*JacPDstandardNth1phi + 2*(g13*JacPDstandardNth2phi + 
      g12*JacPDstandardNth3phi))*Z1 + (phiL*e4phi*JacPDstandardNth2gt23 + 
      2*g22*JacPDstandardNth3phi)*Z2 + (2*g33*JacPDstandardNth2phi + 
      phiL*e4phi*JacPDstandardNth3gt23)*Z3)*pow(phiL,-1),0);
    
    CCTK_REAL R33 CCTK_ATTRIBUTE_UNUSED = Rphi33 + Rt33 + 
      IfThen(formulation != 0,e4phi*(JacPDstandardNth1gt33*Z1 + 
      JacPDstandardNth2gt33*Z2 + JacPDstandardNth3gt33*Z3) + 
      (4*JacPDstandardNth3phi*(g13*Z1 + g23*Z2) - 
      2*g33*(JacPDstandardNth1phi*Z1 + JacPDstandardNth2phi*Z2 - 
      JacPDstandardNth3phi*Z3))*pow(phiL,-1),0);
    
    CCTK_REAL rho CCTK_ATTRIBUTE_UNUSED = pow(alphaL,-2)*(eTttL - 
      2*(beta2L*eTtyL + beta3L*eTtzL) + 2*(beta1L*(-eTtxL + beta2L*eTxyL + 
      beta3L*eTxzL) + beta2L*beta3L*eTyzL) + eTxxL*pow(beta1L,2) + 
      eTyyL*pow(beta2L,2) + eTzzL*pow(beta3L,2));
    
    CCTK_REAL S1 CCTK_ATTRIBUTE_UNUSED = (-eTtxL + beta1L*eTxxL + 
      beta2L*eTxyL + beta3L*eTxzL)*pow(alphaL,-1);
    
    CCTK_REAL S2 CCTK_ATTRIBUTE_UNUSED = (-eTtyL + beta1L*eTxyL + 
      beta2L*eTyyL + beta3L*eTyzL)*pow(alphaL,-1);
    
    CCTK_REAL S3 CCTK_ATTRIBUTE_UNUSED = (-eTtzL + beta1L*eTxzL + 
      beta2L*eTyzL + beta3L*eTzzL)*pow(alphaL,-1);
    
    CCTK_REAL trR CCTK_ATTRIBUTE_UNUSED = gu11*R11 + gu22*R22 + 
      2*(gu12*R12 + gu13*R13 + gu23*R23) + gu33*R33;
    
    CCTK_REAL cXt1L CCTK_ATTRIBUTE_UNUSED = -Xt1L + Xtn1;
    
    CCTK_REAL cXt2L CCTK_ATTRIBUTE_UNUSED = -Xt2L + Xtn2;
    
    CCTK_REAL cXt3L CCTK_ATTRIBUTE_UNUSED = -Xt3L + Xtn3;
    
    CCTK_REAL HL CCTK_ATTRIBUTE_UNUSED = -2*Atm12*Atm21 - 2*Atm13*Atm31 - 
      2*Atm23*Atm32 - 16*Pi*rho + trR + 
      0.666666666666666666666666666667*pow(trKL,2) - pow(Atm11,2) - 
      pow(Atm22,2) - pow(Atm33,2);
    
    CCTK_REAL M1L CCTK_ATTRIBUTE_UNUSED = Gt311*(-2*At13L*gtu11 - 
      At23L*gtu12 - At33L*gtu13) + (-(At22L*Gt212) - At23L*Gt312)*gtu22 + 
      (6*At13L*cdphi2 - At22L*Gt213 - At33L*Gt312 + At23L*(-Gt212 - 
      Gt313))*gtu23 + (6*At13L*cdphi3 - At23L*Gt213 - At33L*Gt313)*gtu33 + 
      At13L*(-2*Gt323*gtu23 - Gt113*gtu33) + At11L*(6*cdphi1*gtu11 - 
      2*Gt111*gtu11 + 6*cdphi2*gtu12 - 3*Gt112*gtu12 + 6*cdphi3*gtu13 - 
      3*Gt113*gtu13 - Gt122*gtu22 - 2*Gt123*gtu23 - Gt133*gtu33) - 
      At12L*(2*Gt211*gtu11 - 6*cdphi1*gtu12 + Gt111*gtu12 + 3*Gt212*gtu12 + 
      3*Gt213*gtu13 - 6*cdphi2*gtu22 + Gt112*gtu22 + Gt222*gtu22 - 
      6*cdphi3*gtu23 + Gt113*gtu23 + 2*Gt223*gtu23 + Gt233*gtu33) + 
      At13L*((6*cdphi1 - Gt111)*gtu13 - Gt322*gtu22 - Gt112*gtu23 - 
      Gt333*gtu33) + gtu11*JacPDstandardNth1At11 + gtu12*(-(At22L*Gt211) - 
      3*At13L*Gt312 + JacPDstandardNth1At12) + gtu13*(-(At23L*Gt211) - 
      3*At13L*Gt313 + JacPDstandardNth1At13) - 
      0.666666666666666666666666666667*JacPDstandardNth1trK + 
      gtu12*JacPDstandardNth2At11 + gtu22*JacPDstandardNth2At12 + 
      gtu23*JacPDstandardNth2At13 + gtu13*JacPDstandardNth3At11 + 
      gtu23*JacPDstandardNth3At12 + gtu33*JacPDstandardNth3At13 - 8*Pi*S1;
    
    CCTK_REAL M2L CCTK_ATTRIBUTE_UNUSED = (-(At22L*Gt211) - At23L*Gt311 - 
      At13L*Gt312)*gtu11 + (6*At22L*cdphi1 - At13L*Gt322)*gtu12 + 
      (6*At23L*cdphi1 + At13L*(-Gt112 - Gt323))*gtu13 + Gt212*(-3*At22L*gtu12 
      - At23L*gtu13) + Gt312*(-3*At23L*gtu12 - At33L*gtu13) - 
      At11L*(Gt112*gtu11 + Gt122*gtu12 + Gt123*gtu13) + 6*At22L*cdphi2*gtu22 
      - 2*At22L*Gt222*gtu22 + 6*At22L*cdphi3*gtu23 + (6*At23L*cdphi2 - 
      At13L*Gt122)*gtu23 + Gt322*(-2*At23L*gtu22 - At33L*gtu23) + 
      (6*At23L*cdphi3 - At13L*Gt123)*gtu33 + Gt223*(-3*At22L*gtu23 - 
      At23L*gtu33) + Gt323*(-3*At23L*gtu23 - At33L*gtu33) + At12L*((6*cdphi1 
      - Gt111 - Gt212)*gtu11 - 3*Gt112*gtu12 + (6*cdphi2 - Gt222)*gtu12 - 
      2*Gt113*gtu13 + (6*cdphi3 - Gt223)*gtu13 - 2*Gt122*gtu22 - 
      3*Gt123*gtu23 - Gt133*gtu33) + At22L*(-2*Gt213*gtu13 - Gt233*gtu33) + 
      At23L*(-2*Gt313*gtu13 - Gt222*gtu23 - Gt333*gtu33) + 
      gtu11*JacPDstandardNth1At12 + gtu12*JacPDstandardNth1At22 + 
      gtu13*JacPDstandardNth1At23 + gtu12*JacPDstandardNth2At12 + 
      gtu22*JacPDstandardNth2At22 + gtu23*JacPDstandardNth2At23 - 
      0.666666666666666666666666666667*JacPDstandardNth2trK + 
      gtu13*JacPDstandardNth3At12 + gtu23*JacPDstandardNth3At22 + 
      gtu33*JacPDstandardNth3At23 - 8*Pi*S2;
    
    CCTK_REAL M3L CCTK_ATTRIBUTE_UNUSED = (-(At23L*Gt211) - At12L*Gt213 - 
      At33L*Gt311)*gtu11 + (-2*At23L*Gt212 - At12L*Gt223)*gtu12 + 
      (-(At12L*Gt113) - At22L*Gt213 + At23L*(6*cdphi1 - Gt313))*gtu12 + 
      (6*At33L*cdphi1 - At12L*Gt233)*gtu13 - 3*At33L*Gt313*gtu13 - 
      At11L*(Gt113*gtu11 + Gt123*gtu12 + Gt133*gtu13) + (-(At12L*Gt123) - 
      At22L*Gt223 + At23L*(6*cdphi2 - Gt323))*gtu22 + At33L*(-2*Gt312*gtu12 - 
      Gt322*gtu22) + (6*At33L*cdphi2 - At12L*Gt133)*gtu23 - 
      3*At23L*Gt223*gtu23 + (6*At23L*cdphi3 - At22L*Gt233)*gtu23 - 
      3*At33L*Gt323*gtu23 + At23L*(-3*Gt213*gtu13 - Gt222*gtu22 - 
      Gt333*gtu23) + 6*At33L*cdphi3*gtu33 - 2*At23L*Gt233*gtu33 - 
      2*At33L*Gt333*gtu33 + At13L*((6*cdphi1 - Gt111 - Gt313)*gtu11 - 
      2*Gt112*gtu12 + (6*cdphi2 - Gt323)*gtu12 - 3*Gt113*gtu13 + (6*cdphi3 - 
      Gt333)*gtu13 - Gt122*gtu22 - 3*Gt123*gtu23 - 2*Gt133*gtu33) + 
      gtu11*JacPDstandardNth1At13 + gtu12*JacPDstandardNth1At23 + 
      gtu13*JacPDstandardNth1At33 + gtu12*JacPDstandardNth2At13 + 
      gtu22*JacPDstandardNth2At23 + gtu23*JacPDstandardNth2At33 + 
      gtu13*JacPDstandardNth3At13 + gtu23*JacPDstandardNth3At23 + 
      gtu33*JacPDstandardNth3At33 - 
      0.666666666666666666666666666667*JacPDstandardNth3trK - 8*Pi*S3;
    /* Copy local copies back to grid functions */
    cXt1[index] = cXt1L;
    cXt2[index] = cXt2L;
    cXt3[index] = cXt3L;
    H[index] = HL;
    M1[index] = M1L;
    M2[index] = M2L;
    M3[index] = M3L;
  }
  CCTK_ENDLOOP3(ML_BSSN_NV_ConstraintsInterior);
}
extern "C" void ML_BSSN_NV_ConstraintsInterior(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_BSSN_NV_ConstraintsInterior
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_NV_ConstraintsInterior);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_NV_ConstraintsInterior_Body");
  }
  if (cctk_iteration % ML_BSSN_NV_ConstraintsInterior_calc_every != ML_BSSN_NV_ConstraintsInterior_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ML_BSSN_NV::ML_cons_Gamma",
    "ML_BSSN_NV::ML_curv",
    "ML_BSSN_NV::ML_Gamma",
    "ML_BSSN_NV::ML_Ham",
    "ML_BSSN_NV::ML_lapse",
    "ML_BSSN_NV::ML_log_confac",
    "ML_BSSN_NV::ML_metric",
    "ML_BSSN_NV::ML_mom",
    "ML_BSSN_NV::ML_shift",
    "ML_BSSN_NV::ML_trace_curv"};
  AssertGroupStorage(cctkGH, "ML_BSSN_NV_ConstraintsInterior", 10, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_NV_ConstraintsInterior", 1, 1, 1);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_NV_ConstraintsInterior", 2, 2, 2);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_NV_ConstraintsInterior", 3, 3, 3);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_NV_ConstraintsInterior", 4, 4, 4);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, ML_BSSN_NV_ConstraintsInterior_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_NV_ConstraintsInterior_Body");
  }
}

} // namespace ML_BSSN_NV
