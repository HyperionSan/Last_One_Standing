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

extern "C" void ML_BSSN_NV_InitialADMBase2Interior_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_BSSN_NV_InitialADMBase2Interior_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_NV_InitialADMBase2Interior_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_NV_InitialADMBase2Interior_calc_every != ML_BSSN_NV_InitialADMBase2Interior_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NV::ML_dtlapse","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN_NV::ML_dtlapse.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NV::ML_dtshift","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN_NV::ML_dtshift.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NV::ML_Gamma","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN_NV::ML_Gamma.");
  return;
}

static void ML_BSSN_NV_InitialADMBase2Interior_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(ML_BSSN_NV_InitialADMBase2Interior,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL alpL CCTK_ATTRIBUTE_UNUSED = alp[index];
    CCTK_REAL betaxL CCTK_ATTRIBUTE_UNUSED = betax[index];
    CCTK_REAL betayL CCTK_ATTRIBUTE_UNUSED = betay[index];
    CCTK_REAL betazL CCTK_ATTRIBUTE_UNUSED = betaz[index];
    CCTK_REAL dtalpL CCTK_ATTRIBUTE_UNUSED = dtalp[index];
    CCTK_REAL dtbetaxL CCTK_ATTRIBUTE_UNUSED = dtbetax[index];
    CCTK_REAL dtbetayL CCTK_ATTRIBUTE_UNUSED = dtbetay[index];
    CCTK_REAL dtbetazL CCTK_ATTRIBUTE_UNUSED = dtbetaz[index];
    CCTK_REAL gt11L CCTK_ATTRIBUTE_UNUSED = gt11[index];
    CCTK_REAL gt12L CCTK_ATTRIBUTE_UNUSED = gt12[index];
    CCTK_REAL gt13L CCTK_ATTRIBUTE_UNUSED = gt13[index];
    CCTK_REAL gt22L CCTK_ATTRIBUTE_UNUSED = gt22[index];
    CCTK_REAL gt23L CCTK_ATTRIBUTE_UNUSED = gt23[index];
    CCTK_REAL gt33L CCTK_ATTRIBUTE_UNUSED = gt33[index];
    CCTK_REAL gxxL CCTK_ATTRIBUTE_UNUSED = gxx[index];
    CCTK_REAL gxyL CCTK_ATTRIBUTE_UNUSED = gxy[index];
    CCTK_REAL gxzL CCTK_ATTRIBUTE_UNUSED = gxz[index];
    CCTK_REAL gyyL CCTK_ATTRIBUTE_UNUSED = gyy[index];
    CCTK_REAL gyzL CCTK_ATTRIBUTE_UNUSED = gyz[index];
    CCTK_REAL gzzL CCTK_ATTRIBUTE_UNUSED = gzz[index];
    CCTK_REAL phiL CCTK_ATTRIBUTE_UNUSED = phi[index];
    CCTK_REAL rL CCTK_ATTRIBUTE_UNUSED = r[index];
    
    
    CCTK_REAL J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L CCTK_ATTRIBUTE_UNUSED ;
    
    if (usejacobian)
    {
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
    CCTK_REAL PDupwindNthSymm1alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDupwindNthSymm1alp = PDupwindNthSymmfdOrder21(&alp[index]);
        PDupwindNthSymm2alp = PDupwindNthSymmfdOrder22(&alp[index]);
        PDupwindNthSymm3alp = PDupwindNthSymmfdOrder23(&alp[index]);
        PDupwindNthAnti1alp = PDupwindNthAntifdOrder21(&alp[index]);
        PDupwindNthAnti2alp = PDupwindNthAntifdOrder22(&alp[index]);
        PDupwindNthAnti3alp = PDupwindNthAntifdOrder23(&alp[index]);
        PDupwindNthSymm1betax = PDupwindNthSymmfdOrder21(&betax[index]);
        PDupwindNthSymm2betax = PDupwindNthSymmfdOrder22(&betax[index]);
        PDupwindNthSymm3betax = PDupwindNthSymmfdOrder23(&betax[index]);
        PDupwindNthAnti1betax = PDupwindNthAntifdOrder21(&betax[index]);
        PDupwindNthAnti2betax = PDupwindNthAntifdOrder22(&betax[index]);
        PDupwindNthAnti3betax = PDupwindNthAntifdOrder23(&betax[index]);
        PDupwindNthSymm1betay = PDupwindNthSymmfdOrder21(&betay[index]);
        PDupwindNthSymm2betay = PDupwindNthSymmfdOrder22(&betay[index]);
        PDupwindNthSymm3betay = PDupwindNthSymmfdOrder23(&betay[index]);
        PDupwindNthAnti1betay = PDupwindNthAntifdOrder21(&betay[index]);
        PDupwindNthAnti2betay = PDupwindNthAntifdOrder22(&betay[index]);
        PDupwindNthAnti3betay = PDupwindNthAntifdOrder23(&betay[index]);
        PDupwindNthSymm1betaz = PDupwindNthSymmfdOrder21(&betaz[index]);
        PDupwindNthSymm2betaz = PDupwindNthSymmfdOrder22(&betaz[index]);
        PDupwindNthSymm3betaz = PDupwindNthSymmfdOrder23(&betaz[index]);
        PDupwindNthAnti1betaz = PDupwindNthAntifdOrder21(&betaz[index]);
        PDupwindNthAnti2betaz = PDupwindNthAntifdOrder22(&betaz[index]);
        PDupwindNthAnti3betaz = PDupwindNthAntifdOrder23(&betaz[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder21(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder22(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder23(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder21(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder22(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder23(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder21(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder22(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder23(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder21(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder22(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder23(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder21(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder22(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder23(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder21(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder22(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder23(&gt33[index]);
        break;
      }
      
      case 4:
      {
        PDupwindNthSymm1alp = PDupwindNthSymmfdOrder41(&alp[index]);
        PDupwindNthSymm2alp = PDupwindNthSymmfdOrder42(&alp[index]);
        PDupwindNthSymm3alp = PDupwindNthSymmfdOrder43(&alp[index]);
        PDupwindNthAnti1alp = PDupwindNthAntifdOrder41(&alp[index]);
        PDupwindNthAnti2alp = PDupwindNthAntifdOrder42(&alp[index]);
        PDupwindNthAnti3alp = PDupwindNthAntifdOrder43(&alp[index]);
        PDupwindNthSymm1betax = PDupwindNthSymmfdOrder41(&betax[index]);
        PDupwindNthSymm2betax = PDupwindNthSymmfdOrder42(&betax[index]);
        PDupwindNthSymm3betax = PDupwindNthSymmfdOrder43(&betax[index]);
        PDupwindNthAnti1betax = PDupwindNthAntifdOrder41(&betax[index]);
        PDupwindNthAnti2betax = PDupwindNthAntifdOrder42(&betax[index]);
        PDupwindNthAnti3betax = PDupwindNthAntifdOrder43(&betax[index]);
        PDupwindNthSymm1betay = PDupwindNthSymmfdOrder41(&betay[index]);
        PDupwindNthSymm2betay = PDupwindNthSymmfdOrder42(&betay[index]);
        PDupwindNthSymm3betay = PDupwindNthSymmfdOrder43(&betay[index]);
        PDupwindNthAnti1betay = PDupwindNthAntifdOrder41(&betay[index]);
        PDupwindNthAnti2betay = PDupwindNthAntifdOrder42(&betay[index]);
        PDupwindNthAnti3betay = PDupwindNthAntifdOrder43(&betay[index]);
        PDupwindNthSymm1betaz = PDupwindNthSymmfdOrder41(&betaz[index]);
        PDupwindNthSymm2betaz = PDupwindNthSymmfdOrder42(&betaz[index]);
        PDupwindNthSymm3betaz = PDupwindNthSymmfdOrder43(&betaz[index]);
        PDupwindNthAnti1betaz = PDupwindNthAntifdOrder41(&betaz[index]);
        PDupwindNthAnti2betaz = PDupwindNthAntifdOrder42(&betaz[index]);
        PDupwindNthAnti3betaz = PDupwindNthAntifdOrder43(&betaz[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder41(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder42(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder43(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder41(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder42(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder43(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder41(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder42(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder43(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder41(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder42(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder43(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder41(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder42(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder43(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder41(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder42(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder43(&gt33[index]);
        break;
      }
      
      case 6:
      {
        PDupwindNthSymm1alp = PDupwindNthSymmfdOrder61(&alp[index]);
        PDupwindNthSymm2alp = PDupwindNthSymmfdOrder62(&alp[index]);
        PDupwindNthSymm3alp = PDupwindNthSymmfdOrder63(&alp[index]);
        PDupwindNthAnti1alp = PDupwindNthAntifdOrder61(&alp[index]);
        PDupwindNthAnti2alp = PDupwindNthAntifdOrder62(&alp[index]);
        PDupwindNthAnti3alp = PDupwindNthAntifdOrder63(&alp[index]);
        PDupwindNthSymm1betax = PDupwindNthSymmfdOrder61(&betax[index]);
        PDupwindNthSymm2betax = PDupwindNthSymmfdOrder62(&betax[index]);
        PDupwindNthSymm3betax = PDupwindNthSymmfdOrder63(&betax[index]);
        PDupwindNthAnti1betax = PDupwindNthAntifdOrder61(&betax[index]);
        PDupwindNthAnti2betax = PDupwindNthAntifdOrder62(&betax[index]);
        PDupwindNthAnti3betax = PDupwindNthAntifdOrder63(&betax[index]);
        PDupwindNthSymm1betay = PDupwindNthSymmfdOrder61(&betay[index]);
        PDupwindNthSymm2betay = PDupwindNthSymmfdOrder62(&betay[index]);
        PDupwindNthSymm3betay = PDupwindNthSymmfdOrder63(&betay[index]);
        PDupwindNthAnti1betay = PDupwindNthAntifdOrder61(&betay[index]);
        PDupwindNthAnti2betay = PDupwindNthAntifdOrder62(&betay[index]);
        PDupwindNthAnti3betay = PDupwindNthAntifdOrder63(&betay[index]);
        PDupwindNthSymm1betaz = PDupwindNthSymmfdOrder61(&betaz[index]);
        PDupwindNthSymm2betaz = PDupwindNthSymmfdOrder62(&betaz[index]);
        PDupwindNthSymm3betaz = PDupwindNthSymmfdOrder63(&betaz[index]);
        PDupwindNthAnti1betaz = PDupwindNthAntifdOrder61(&betaz[index]);
        PDupwindNthAnti2betaz = PDupwindNthAntifdOrder62(&betaz[index]);
        PDupwindNthAnti3betaz = PDupwindNthAntifdOrder63(&betaz[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder61(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder62(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder63(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder61(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder62(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder63(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder61(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder62(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder63(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder61(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder62(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder63(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder61(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder62(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder63(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder61(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder62(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder63(&gt33[index]);
        break;
      }
      
      case 8:
      {
        PDupwindNthSymm1alp = PDupwindNthSymmfdOrder81(&alp[index]);
        PDupwindNthSymm2alp = PDupwindNthSymmfdOrder82(&alp[index]);
        PDupwindNthSymm3alp = PDupwindNthSymmfdOrder83(&alp[index]);
        PDupwindNthAnti1alp = PDupwindNthAntifdOrder81(&alp[index]);
        PDupwindNthAnti2alp = PDupwindNthAntifdOrder82(&alp[index]);
        PDupwindNthAnti3alp = PDupwindNthAntifdOrder83(&alp[index]);
        PDupwindNthSymm1betax = PDupwindNthSymmfdOrder81(&betax[index]);
        PDupwindNthSymm2betax = PDupwindNthSymmfdOrder82(&betax[index]);
        PDupwindNthSymm3betax = PDupwindNthSymmfdOrder83(&betax[index]);
        PDupwindNthAnti1betax = PDupwindNthAntifdOrder81(&betax[index]);
        PDupwindNthAnti2betax = PDupwindNthAntifdOrder82(&betax[index]);
        PDupwindNthAnti3betax = PDupwindNthAntifdOrder83(&betax[index]);
        PDupwindNthSymm1betay = PDupwindNthSymmfdOrder81(&betay[index]);
        PDupwindNthSymm2betay = PDupwindNthSymmfdOrder82(&betay[index]);
        PDupwindNthSymm3betay = PDupwindNthSymmfdOrder83(&betay[index]);
        PDupwindNthAnti1betay = PDupwindNthAntifdOrder81(&betay[index]);
        PDupwindNthAnti2betay = PDupwindNthAntifdOrder82(&betay[index]);
        PDupwindNthAnti3betay = PDupwindNthAntifdOrder83(&betay[index]);
        PDupwindNthSymm1betaz = PDupwindNthSymmfdOrder81(&betaz[index]);
        PDupwindNthSymm2betaz = PDupwindNthSymmfdOrder82(&betaz[index]);
        PDupwindNthSymm3betaz = PDupwindNthSymmfdOrder83(&betaz[index]);
        PDupwindNthAnti1betaz = PDupwindNthAntifdOrder81(&betaz[index]);
        PDupwindNthAnti2betaz = PDupwindNthAntifdOrder82(&betaz[index]);
        PDupwindNthAnti3betaz = PDupwindNthAntifdOrder83(&betaz[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder81(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder82(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder83(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder81(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder82(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder83(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder81(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder82(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder83(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder81(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder82(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder83(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder81(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder82(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder83(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder81(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder82(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder83(&gt33[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL JacPDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3betaz CCTK_ATTRIBUTE_UNUSED;
    
    if (usejacobian)
    {
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
      
      JacPDupwindNthSymm1alp = J11L*PDupwindNthSymm1alp + 
        J21L*PDupwindNthSymm2alp + J31L*PDupwindNthSymm3alp;
      
      JacPDupwindNthSymm1betax = J11L*PDupwindNthSymm1betax + 
        J21L*PDupwindNthSymm2betax + J31L*PDupwindNthSymm3betax;
      
      JacPDupwindNthSymm1betay = J11L*PDupwindNthSymm1betay + 
        J21L*PDupwindNthSymm2betay + J31L*PDupwindNthSymm3betay;
      
      JacPDupwindNthSymm1betaz = J11L*PDupwindNthSymm1betaz + 
        J21L*PDupwindNthSymm2betaz + J31L*PDupwindNthSymm3betaz;
      
      JacPDupwindNthSymm2alp = J12L*PDupwindNthSymm1alp + 
        J22L*PDupwindNthSymm2alp + J32L*PDupwindNthSymm3alp;
      
      JacPDupwindNthSymm2betax = J12L*PDupwindNthSymm1betax + 
        J22L*PDupwindNthSymm2betax + J32L*PDupwindNthSymm3betax;
      
      JacPDupwindNthSymm2betay = J12L*PDupwindNthSymm1betay + 
        J22L*PDupwindNthSymm2betay + J32L*PDupwindNthSymm3betay;
      
      JacPDupwindNthSymm2betaz = J12L*PDupwindNthSymm1betaz + 
        J22L*PDupwindNthSymm2betaz + J32L*PDupwindNthSymm3betaz;
      
      JacPDupwindNthSymm3alp = J13L*PDupwindNthSymm1alp + 
        J23L*PDupwindNthSymm2alp + J33L*PDupwindNthSymm3alp;
      
      JacPDupwindNthSymm3betax = J13L*PDupwindNthSymm1betax + 
        J23L*PDupwindNthSymm2betax + J33L*PDupwindNthSymm3betax;
      
      JacPDupwindNthSymm3betay = J13L*PDupwindNthSymm1betay + 
        J23L*PDupwindNthSymm2betay + J33L*PDupwindNthSymm3betay;
      
      JacPDupwindNthSymm3betaz = J13L*PDupwindNthSymm1betaz + 
        J23L*PDupwindNthSymm2betaz + J33L*PDupwindNthSymm3betaz;
      
      JacPDupwindNthAnti1alp = J11L*PDupwindNthAnti1alp + 
        J21L*PDupwindNthAnti2alp + J31L*PDupwindNthAnti3alp;
      
      JacPDupwindNthAnti1betax = J11L*PDupwindNthAnti1betax + 
        J21L*PDupwindNthAnti2betax + J31L*PDupwindNthAnti3betax;
      
      JacPDupwindNthAnti1betay = J11L*PDupwindNthAnti1betay + 
        J21L*PDupwindNthAnti2betay + J31L*PDupwindNthAnti3betay;
      
      JacPDupwindNthAnti1betaz = J11L*PDupwindNthAnti1betaz + 
        J21L*PDupwindNthAnti2betaz + J31L*PDupwindNthAnti3betaz;
      
      JacPDupwindNthAnti2alp = J12L*PDupwindNthAnti1alp + 
        J22L*PDupwindNthAnti2alp + J32L*PDupwindNthAnti3alp;
      
      JacPDupwindNthAnti2betax = J12L*PDupwindNthAnti1betax + 
        J22L*PDupwindNthAnti2betax + J32L*PDupwindNthAnti3betax;
      
      JacPDupwindNthAnti2betay = J12L*PDupwindNthAnti1betay + 
        J22L*PDupwindNthAnti2betay + J32L*PDupwindNthAnti3betay;
      
      JacPDupwindNthAnti2betaz = J12L*PDupwindNthAnti1betaz + 
        J22L*PDupwindNthAnti2betaz + J32L*PDupwindNthAnti3betaz;
      
      JacPDupwindNthAnti3alp = J13L*PDupwindNthAnti1alp + 
        J23L*PDupwindNthAnti2alp + J33L*PDupwindNthAnti3alp;
      
      JacPDupwindNthAnti3betax = J13L*PDupwindNthAnti1betax + 
        J23L*PDupwindNthAnti2betax + J33L*PDupwindNthAnti3betax;
      
      JacPDupwindNthAnti3betay = J13L*PDupwindNthAnti1betay + 
        J23L*PDupwindNthAnti2betay + J33L*PDupwindNthAnti3betay;
      
      JacPDupwindNthAnti3betaz = J13L*PDupwindNthAnti1betaz + 
        J23L*PDupwindNthAnti2betaz + J33L*PDupwindNthAnti3betaz;
    }
    else
    {
      JacPDstandardNth1gt11 = PDstandardNth1gt11;
      
      JacPDstandardNth1gt12 = PDstandardNth1gt12;
      
      JacPDstandardNth1gt13 = PDstandardNth1gt13;
      
      JacPDstandardNth1gt22 = PDstandardNth1gt22;
      
      JacPDstandardNth1gt23 = PDstandardNth1gt23;
      
      JacPDstandardNth1gt33 = PDstandardNth1gt33;
      
      JacPDstandardNth2gt11 = PDstandardNth2gt11;
      
      JacPDstandardNth2gt12 = PDstandardNth2gt12;
      
      JacPDstandardNth2gt13 = PDstandardNth2gt13;
      
      JacPDstandardNth2gt22 = PDstandardNth2gt22;
      
      JacPDstandardNth2gt23 = PDstandardNth2gt23;
      
      JacPDstandardNth2gt33 = PDstandardNth2gt33;
      
      JacPDstandardNth3gt11 = PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = PDstandardNth3gt33;
      
      JacPDupwindNthSymm1alp = PDupwindNthSymm1alp;
      
      JacPDupwindNthSymm1betax = PDupwindNthSymm1betax;
      
      JacPDupwindNthSymm1betay = PDupwindNthSymm1betay;
      
      JacPDupwindNthSymm1betaz = PDupwindNthSymm1betaz;
      
      JacPDupwindNthSymm2alp = PDupwindNthSymm2alp;
      
      JacPDupwindNthSymm2betax = PDupwindNthSymm2betax;
      
      JacPDupwindNthSymm2betay = PDupwindNthSymm2betay;
      
      JacPDupwindNthSymm2betaz = PDupwindNthSymm2betaz;
      
      JacPDupwindNthSymm3alp = PDupwindNthSymm3alp;
      
      JacPDupwindNthSymm3betax = PDupwindNthSymm3betax;
      
      JacPDupwindNthSymm3betay = PDupwindNthSymm3betay;
      
      JacPDupwindNthSymm3betaz = PDupwindNthSymm3betaz;
      
      JacPDupwindNthAnti1alp = PDupwindNthAnti1alp;
      
      JacPDupwindNthAnti1betax = PDupwindNthAnti1betax;
      
      JacPDupwindNthAnti1betay = PDupwindNthAnti1betay;
      
      JacPDupwindNthAnti1betaz = PDupwindNthAnti1betaz;
      
      JacPDupwindNthAnti2alp = PDupwindNthAnti2alp;
      
      JacPDupwindNthAnti2betax = PDupwindNthAnti2betax;
      
      JacPDupwindNthAnti2betay = PDupwindNthAnti2betay;
      
      JacPDupwindNthAnti2betaz = PDupwindNthAnti2betaz;
      
      JacPDupwindNthAnti3alp = PDupwindNthAnti3alp;
      
      JacPDupwindNthAnti3betax = PDupwindNthAnti3betax;
      
      JacPDupwindNthAnti3betay = PDupwindNthAnti3betay;
      
      JacPDupwindNthAnti3betaz = PDupwindNthAnti3betaz;
    }
    
    CCTK_REAL g11 CCTK_ATTRIBUTE_UNUSED = gxxL;
    
    CCTK_REAL g12 CCTK_ATTRIBUTE_UNUSED = gxyL;
    
    CCTK_REAL g13 CCTK_ATTRIBUTE_UNUSED = gxzL;
    
    CCTK_REAL g22 CCTK_ATTRIBUTE_UNUSED = gyyL;
    
    CCTK_REAL g23 CCTK_ATTRIBUTE_UNUSED = gyzL;
    
    CCTK_REAL g33 CCTK_ATTRIBUTE_UNUSED = gzzL;
    
    CCTK_REAL detg CCTK_ATTRIBUTE_UNUSED = 2*g12*g13*g23 - g33*pow(g12,2) 
      + g22*(g11*g33 - pow(g13,2)) - g11*pow(g23,2);
    
    CCTK_REAL gu11 CCTK_ATTRIBUTE_UNUSED = pow(detg,-1)*(g22*g33 - 
      pow(g23,2));
    
    CCTK_REAL gu12 CCTK_ATTRIBUTE_UNUSED = (g13*g23 - 
      g12*g33)*pow(detg,-1);
    
    CCTK_REAL gu13 CCTK_ATTRIBUTE_UNUSED = (-(g13*g22) + 
      g12*g23)*pow(detg,-1);
    
    CCTK_REAL gu22 CCTK_ATTRIBUTE_UNUSED = pow(detg,-1)*(g11*g33 - 
      pow(g13,2));
    
    CCTK_REAL gu23 CCTK_ATTRIBUTE_UNUSED = (g12*g13 - 
      g11*g23)*pow(detg,-1);
    
    CCTK_REAL gu33 CCTK_ATTRIBUTE_UNUSED = pow(detg,-1)*(g11*g22 - 
      pow(g12,2));
    
    CCTK_REAL em4phi CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod != 
      0,pow(phiL,2),exp(-4*phiL));
    
    CCTK_REAL gtu11 CCTK_ATTRIBUTE_UNUSED = gu11*pow(em4phi,-1);
    
    CCTK_REAL gtu12 CCTK_ATTRIBUTE_UNUSED = gu12*pow(em4phi,-1);
    
    CCTK_REAL gtu13 CCTK_ATTRIBUTE_UNUSED = gu13*pow(em4phi,-1);
    
    CCTK_REAL gtu22 CCTK_ATTRIBUTE_UNUSED = gu22*pow(em4phi,-1);
    
    CCTK_REAL gtu23 CCTK_ATTRIBUTE_UNUSED = gu23*pow(em4phi,-1);
    
    CCTK_REAL gtu33 CCTK_ATTRIBUTE_UNUSED = gu33*pow(em4phi,-1);
    
    CCTK_REAL Gt111 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu11*JacPDstandardNth1gt11 + 2*(gtu12*JacPDstandardNth1gt12 + 
      gtu13*JacPDstandardNth1gt13) - gtu12*JacPDstandardNth2gt11 - 
      gtu13*JacPDstandardNth3gt11);
    
    CCTK_REAL Gt211 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu12*JacPDstandardNth1gt11 + 2*(gtu22*JacPDstandardNth1gt12 + 
      gtu23*JacPDstandardNth1gt13) - gtu22*JacPDstandardNth2gt11 - 
      gtu23*JacPDstandardNth3gt11);
    
    CCTK_REAL Gt311 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu13*JacPDstandardNth1gt11 + 2*(gtu23*JacPDstandardNth1gt12 + 
      gtu33*JacPDstandardNth1gt13) - gtu23*JacPDstandardNth2gt11 - 
      gtu33*JacPDstandardNth3gt11);
    
    CCTK_REAL Gt112 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu12*JacPDstandardNth1gt22 + gtu11*JacPDstandardNth2gt11 + 
      gtu13*(JacPDstandardNth1gt23 + JacPDstandardNth2gt13 - 
      JacPDstandardNth3gt12));
    
    CCTK_REAL Gt212 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu22*JacPDstandardNth1gt22 + gtu12*JacPDstandardNth2gt11 + 
      gtu23*(JacPDstandardNth1gt23 + JacPDstandardNth2gt13 - 
      JacPDstandardNth3gt12));
    
    CCTK_REAL Gt312 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu23*JacPDstandardNth1gt22 + gtu13*JacPDstandardNth2gt11 + 
      gtu33*(JacPDstandardNth1gt23 + JacPDstandardNth2gt13 - 
      JacPDstandardNth3gt12));
    
    CCTK_REAL Gt113 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu13*JacPDstandardNth1gt33 + gtu11*JacPDstandardNth3gt11 + 
      gtu12*(JacPDstandardNth1gt23 - JacPDstandardNth2gt13 + 
      JacPDstandardNth3gt12));
    
    CCTK_REAL Gt213 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu23*JacPDstandardNth1gt33 + gtu12*JacPDstandardNth3gt11 + 
      gtu22*(JacPDstandardNth1gt23 - JacPDstandardNth2gt13 + 
      JacPDstandardNth3gt12));
    
    CCTK_REAL Gt313 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu33*JacPDstandardNth1gt33 + gtu13*JacPDstandardNth3gt11 + 
      gtu23*(JacPDstandardNth1gt23 - JacPDstandardNth2gt13 + 
      JacPDstandardNth3gt12));
    
    CCTK_REAL Gt122 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu11*(-JacPDstandardNth1gt22 + 2*JacPDstandardNth2gt12) + 
      gtu12*JacPDstandardNth2gt22 + gtu13*(2*JacPDstandardNth2gt23 - 
      JacPDstandardNth3gt22));
    
    CCTK_REAL Gt222 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu12*(-JacPDstandardNth1gt22 + 2*JacPDstandardNth2gt12) + 
      gtu22*JacPDstandardNth2gt22 + gtu23*(2*JacPDstandardNth2gt23 - 
      JacPDstandardNth3gt22));
    
    CCTK_REAL Gt322 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu13*(-JacPDstandardNth1gt22 + 2*JacPDstandardNth2gt12) + 
      gtu23*JacPDstandardNth2gt22 + gtu33*(2*JacPDstandardNth2gt23 - 
      JacPDstandardNth3gt22));
    
    CCTK_REAL Gt123 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu13*JacPDstandardNth2gt33 + gtu11*(-JacPDstandardNth1gt23 + 
      JacPDstandardNth2gt13 + JacPDstandardNth3gt12) + 
      gtu12*JacPDstandardNth3gt22);
    
    CCTK_REAL Gt223 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu23*JacPDstandardNth2gt33 + gtu12*(-JacPDstandardNth1gt23 + 
      JacPDstandardNth2gt13 + JacPDstandardNth3gt12) + 
      gtu22*JacPDstandardNth3gt22);
    
    CCTK_REAL Gt323 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu33*JacPDstandardNth2gt33 + gtu13*(-JacPDstandardNth1gt23 + 
      JacPDstandardNth2gt13 + JacPDstandardNth3gt12) + 
      gtu23*JacPDstandardNth3gt22);
    
    CCTK_REAL Gt133 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu11*(-JacPDstandardNth1gt33 + 2*JacPDstandardNth3gt13) + 
      gtu12*(-JacPDstandardNth2gt33 + 2*JacPDstandardNth3gt23) + 
      gtu13*JacPDstandardNth3gt33);
    
    CCTK_REAL Gt233 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu12*(-JacPDstandardNth1gt33 + 2*JacPDstandardNth3gt13) + 
      gtu22*(-JacPDstandardNth2gt33 + 2*JacPDstandardNth3gt23) + 
      gtu23*JacPDstandardNth3gt33);
    
    CCTK_REAL Gt333 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu13*(-JacPDstandardNth1gt33 + 2*JacPDstandardNth3gt13) + 
      gtu23*(-JacPDstandardNth2gt33 + 2*JacPDstandardNth3gt23) + 
      gtu33*JacPDstandardNth3gt33);
    
    CCTK_REAL Xt1L CCTK_ATTRIBUTE_UNUSED = Gt111*gtu11 + Gt122*gtu22 + 
      2*(Gt112*gtu12 + Gt113*gtu13 + Gt123*gtu23) + Gt133*gtu33;
    
    CCTK_REAL Xt2L CCTK_ATTRIBUTE_UNUSED = Gt211*gtu11 + Gt222*gtu22 + 
      2*(Gt212*gtu12 + Gt213*gtu13 + Gt223*gtu23) + Gt233*gtu33;
    
    CCTK_REAL Xt3L CCTK_ATTRIBUTE_UNUSED = Gt311*gtu11 + Gt322*gtu22 + 
      2*(Gt312*gtu12 + Gt313*gtu13 + Gt323*gtu23) + Gt333*gtu33;
    
    CCTK_REAL AL CCTK_ATTRIBUTE_UNUSED = IfThen(evolveA != 0,(-dtalpL + 
      IfThen(advectLapse != 0,betaxL*JacPDupwindNthAnti1alp + 
      betayL*JacPDupwindNthAnti2alp + betazL*JacPDupwindNthAnti3alp + 
      JacPDupwindNthSymm1alp*fabs(betaxL) + 
      JacPDupwindNthSymm2alp*fabs(betayL) + 
      JacPDupwindNthSymm3alp*fabs(betazL),0))*pow(alpL,-harmonicN)*pow(harmonicF,-1),0);
    
    CCTK_REAL shiftGammaCoeffValue CCTK_ATTRIBUTE_UNUSED = 
      IfThen(useSpatialShiftGammaCoeff != 0,shiftGammaCoeff*fmin(1,exp(1 - 
      rL*pow(spatialShiftGammaCoeffRadius,-1))),shiftGammaCoeff);
    
    CCTK_REAL B1L CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL B2L CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL B3L CCTK_ATTRIBUTE_UNUSED;
    
    if (evolveB != 0)
    {
      B1L = (dtbetaxL - IfThen(advectShift != 
        0,betaxL*JacPDupwindNthAnti1betax + betayL*JacPDupwindNthAnti2betax + 
        betazL*JacPDupwindNthAnti3betax + JacPDupwindNthSymm1betax*fabs(betaxL) 
        + JacPDupwindNthSymm2betax*fabs(betayL) + 
        JacPDupwindNthSymm3betax*fabs(betazL),0))*pow(alpL,-shiftAlphaPower)*pow(shiftGammaCoeffValue,-1);
      
      B2L = (dtbetayL - IfThen(advectShift != 
        0,betaxL*JacPDupwindNthAnti1betay + betayL*JacPDupwindNthAnti2betay + 
        betazL*JacPDupwindNthAnti3betay + JacPDupwindNthSymm1betay*fabs(betaxL) 
        + JacPDupwindNthSymm2betay*fabs(betayL) + 
        JacPDupwindNthSymm3betay*fabs(betazL),0))*pow(alpL,-shiftAlphaPower)*pow(shiftGammaCoeffValue,-1);
      
      B3L = (dtbetazL - IfThen(advectShift != 
        0,betaxL*JacPDupwindNthAnti1betaz + betayL*JacPDupwindNthAnti2betaz + 
        betazL*JacPDupwindNthAnti3betaz + JacPDupwindNthSymm1betaz*fabs(betaxL) 
        + JacPDupwindNthSymm2betaz*fabs(betayL) + 
        JacPDupwindNthSymm3betaz*fabs(betazL),0))*pow(alpL,-shiftAlphaPower)*pow(shiftGammaCoeffValue,-1);
    }
    else
    {
      B1L = 0;
      
      B2L = 0;
      
      B3L = 0;
    }
    /* Copy local copies back to grid functions */
    A[index] = AL;
    B1[index] = B1L;
    B2[index] = B2L;
    B3[index] = B3L;
    Xt1[index] = Xt1L;
    Xt2[index] = Xt2L;
    Xt3[index] = Xt3L;
  }
  CCTK_ENDLOOP3(ML_BSSN_NV_InitialADMBase2Interior);
}
extern "C" void ML_BSSN_NV_InitialADMBase2Interior(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_BSSN_NV_InitialADMBase2Interior
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_NV_InitialADMBase2Interior);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_NV_InitialADMBase2Interior_Body");
  }
  if (cctk_iteration % ML_BSSN_NV_InitialADMBase2Interior_calc_every != ML_BSSN_NV_InitialADMBase2Interior_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ADMBase::dtlapse",
    "ADMBase::dtshift",
    "ADMBase::lapse",
    "ADMBase::metric",
    "ADMBase::shift",
    "grid::coordinates",
    "ML_BSSN_NV::ML_dtlapse",
    "ML_BSSN_NV::ML_dtshift",
    "ML_BSSN_NV::ML_Gamma",
    "ML_BSSN_NV::ML_log_confac",
    "ML_BSSN_NV::ML_metric"};
  AssertGroupStorage(cctkGH, "ML_BSSN_NV_InitialADMBase2Interior", 11, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_NV_InitialADMBase2Interior", 2, 2, 2);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_NV_InitialADMBase2Interior", 3, 3, 3);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_NV_InitialADMBase2Interior", 4, 4, 4);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_NV_InitialADMBase2Interior", 5, 5, 5);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, ML_BSSN_NV_InitialADMBase2Interior_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_NV_InitialADMBase2Interior_Body");
  }
}

} // namespace ML_BSSN_NV
