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

extern "C" void ML_BSSN_NV_ADMBaseInterior_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_BSSN_NV_ADMBaseInterior_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_NV_ADMBaseInterior_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_NV_ADMBaseInterior_calc_every != ML_BSSN_NV_ADMBaseInterior_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ADMBase::dtlapse","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ADMBase::dtlapse.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ADMBase::dtshift","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ADMBase::dtshift.");
  return;
}

static void ML_BSSN_NV_ADMBaseInterior_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(ML_BSSN_NV_ADMBaseInterior,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL AL CCTK_ATTRIBUTE_UNUSED = A[index];
    CCTK_REAL alphaL CCTK_ATTRIBUTE_UNUSED = alpha[index];
    CCTK_REAL B1L CCTK_ATTRIBUTE_UNUSED = B1[index];
    CCTK_REAL B2L CCTK_ATTRIBUTE_UNUSED = B2[index];
    CCTK_REAL B3L CCTK_ATTRIBUTE_UNUSED = B3[index];
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
    CCTK_REAL rL CCTK_ATTRIBUTE_UNUSED = r[index];
    CCTK_REAL ThetaL CCTK_ATTRIBUTE_UNUSED = Theta[index];
    CCTK_REAL trKL CCTK_ATTRIBUTE_UNUSED = trK[index];
    CCTK_REAL Xt1L CCTK_ATTRIBUTE_UNUSED = Xt1[index];
    CCTK_REAL Xt2L CCTK_ATTRIBUTE_UNUSED = Xt2[index];
    CCTK_REAL Xt3L CCTK_ATTRIBUTE_UNUSED = Xt3[index];
    
    
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
    CCTK_REAL PDstandardNth1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3beta3 CCTK_ATTRIBUTE_UNUSED;
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
    CCTK_REAL PDstandardNth1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3phi CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDstandardNth1alpha = PDstandardNthfdOrder21(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder22(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder23(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder21(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder22(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder23(&alpha[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder21(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder22(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder23(&alpha[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder21(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder22(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder23(&beta1[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder21(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder22(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder23(&beta1[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder21(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder22(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder23(&beta2[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder21(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder22(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder23(&beta2[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder21(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder22(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder23(&beta3[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder21(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder22(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder23(&beta3[index]);
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
        PDstandardNth1phi = PDstandardNthfdOrder21(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder22(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder23(&phi[index]);
        break;
      }
      
      case 4:
      {
        PDstandardNth1alpha = PDstandardNthfdOrder41(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder42(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder43(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder41(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder42(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder43(&alpha[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder41(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder42(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder43(&alpha[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder41(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder42(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder43(&beta1[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder41(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder42(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder43(&beta1[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder41(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder42(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder43(&beta2[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder41(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder42(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder43(&beta2[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder41(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder42(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder43(&beta3[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder41(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder42(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder43(&beta3[index]);
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
        PDstandardNth1phi = PDstandardNthfdOrder41(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder42(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder43(&phi[index]);
        break;
      }
      
      case 6:
      {
        PDstandardNth1alpha = PDstandardNthfdOrder61(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder62(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder63(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder61(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder62(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder63(&alpha[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder61(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder62(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder63(&alpha[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder61(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder62(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder63(&beta1[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder61(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder62(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder63(&beta1[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder61(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder62(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder63(&beta2[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder61(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder62(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder63(&beta2[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder61(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder62(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder63(&beta3[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder61(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder62(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder63(&beta3[index]);
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
        PDstandardNth1phi = PDstandardNthfdOrder61(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder62(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder63(&phi[index]);
        break;
      }
      
      case 8:
      {
        PDstandardNth1alpha = PDstandardNthfdOrder81(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder82(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder83(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder81(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder82(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder83(&alpha[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder81(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder82(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder83(&alpha[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder81(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder82(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder83(&beta1[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder81(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder82(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder83(&beta1[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder81(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder82(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder83(&beta2[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder81(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder82(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder83(&beta2[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder81(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder82(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder83(&beta3[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder81(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder82(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder83(&beta3[index]);
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
        PDstandardNth1phi = PDstandardNthfdOrder81(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder82(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder83(&phi[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL JacPDstandardNth1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3beta3 CCTK_ATTRIBUTE_UNUSED;
    
    if (usejacobian)
    {
      JacPDstandardNth1alpha = J11L*PDstandardNth1alpha + 
        J21L*PDstandardNth2alpha + J31L*PDstandardNth3alpha;
      
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
      
      JacPDstandardNth2alpha = J12L*PDstandardNth1alpha + 
        J22L*PDstandardNth2alpha + J32L*PDstandardNth3alpha;
      
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
      
      JacPDstandardNth3alpha = J13L*PDstandardNth1alpha + 
        J23L*PDstandardNth2alpha + J33L*PDstandardNth3alpha;
      
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
      
      JacPDupwindNthSymm1alpha = J11L*PDupwindNthSymm1alpha + 
        J21L*PDupwindNthSymm2alpha + J31L*PDupwindNthSymm3alpha;
      
      JacPDupwindNthSymm1beta1 = J11L*PDupwindNthSymm1beta1 + 
        J21L*PDupwindNthSymm2beta1 + J31L*PDupwindNthSymm3beta1;
      
      JacPDupwindNthSymm1beta2 = J11L*PDupwindNthSymm1beta2 + 
        J21L*PDupwindNthSymm2beta2 + J31L*PDupwindNthSymm3beta2;
      
      JacPDupwindNthSymm1beta3 = J11L*PDupwindNthSymm1beta3 + 
        J21L*PDupwindNthSymm2beta3 + J31L*PDupwindNthSymm3beta3;
      
      JacPDupwindNthSymm2alpha = J12L*PDupwindNthSymm1alpha + 
        J22L*PDupwindNthSymm2alpha + J32L*PDupwindNthSymm3alpha;
      
      JacPDupwindNthSymm2beta1 = J12L*PDupwindNthSymm1beta1 + 
        J22L*PDupwindNthSymm2beta1 + J32L*PDupwindNthSymm3beta1;
      
      JacPDupwindNthSymm2beta2 = J12L*PDupwindNthSymm1beta2 + 
        J22L*PDupwindNthSymm2beta2 + J32L*PDupwindNthSymm3beta2;
      
      JacPDupwindNthSymm2beta3 = J12L*PDupwindNthSymm1beta3 + 
        J22L*PDupwindNthSymm2beta3 + J32L*PDupwindNthSymm3beta3;
      
      JacPDupwindNthSymm3alpha = J13L*PDupwindNthSymm1alpha + 
        J23L*PDupwindNthSymm2alpha + J33L*PDupwindNthSymm3alpha;
      
      JacPDupwindNthSymm3beta1 = J13L*PDupwindNthSymm1beta1 + 
        J23L*PDupwindNthSymm2beta1 + J33L*PDupwindNthSymm3beta1;
      
      JacPDupwindNthSymm3beta2 = J13L*PDupwindNthSymm1beta2 + 
        J23L*PDupwindNthSymm2beta2 + J33L*PDupwindNthSymm3beta2;
      
      JacPDupwindNthSymm3beta3 = J13L*PDupwindNthSymm1beta3 + 
        J23L*PDupwindNthSymm2beta3 + J33L*PDupwindNthSymm3beta3;
      
      JacPDupwindNthAnti1alpha = J11L*PDupwindNthAnti1alpha + 
        J21L*PDupwindNthAnti2alpha + J31L*PDupwindNthAnti3alpha;
      
      JacPDupwindNthAnti1beta1 = J11L*PDupwindNthAnti1beta1 + 
        J21L*PDupwindNthAnti2beta1 + J31L*PDupwindNthAnti3beta1;
      
      JacPDupwindNthAnti1beta2 = J11L*PDupwindNthAnti1beta2 + 
        J21L*PDupwindNthAnti2beta2 + J31L*PDupwindNthAnti3beta2;
      
      JacPDupwindNthAnti1beta3 = J11L*PDupwindNthAnti1beta3 + 
        J21L*PDupwindNthAnti2beta3 + J31L*PDupwindNthAnti3beta3;
      
      JacPDupwindNthAnti2alpha = J12L*PDupwindNthAnti1alpha + 
        J22L*PDupwindNthAnti2alpha + J32L*PDupwindNthAnti3alpha;
      
      JacPDupwindNthAnti2beta1 = J12L*PDupwindNthAnti1beta1 + 
        J22L*PDupwindNthAnti2beta1 + J32L*PDupwindNthAnti3beta1;
      
      JacPDupwindNthAnti2beta2 = J12L*PDupwindNthAnti1beta2 + 
        J22L*PDupwindNthAnti2beta2 + J32L*PDupwindNthAnti3beta2;
      
      JacPDupwindNthAnti2beta3 = J12L*PDupwindNthAnti1beta3 + 
        J22L*PDupwindNthAnti2beta3 + J32L*PDupwindNthAnti3beta3;
      
      JacPDupwindNthAnti3alpha = J13L*PDupwindNthAnti1alpha + 
        J23L*PDupwindNthAnti2alpha + J33L*PDupwindNthAnti3alpha;
      
      JacPDupwindNthAnti3beta1 = J13L*PDupwindNthAnti1beta1 + 
        J23L*PDupwindNthAnti2beta1 + J33L*PDupwindNthAnti3beta1;
      
      JacPDupwindNthAnti3beta2 = J13L*PDupwindNthAnti1beta2 + 
        J23L*PDupwindNthAnti2beta2 + J33L*PDupwindNthAnti3beta2;
      
      JacPDupwindNthAnti3beta3 = J13L*PDupwindNthAnti1beta3 + 
        J23L*PDupwindNthAnti2beta3 + J33L*PDupwindNthAnti3beta3;
    }
    else
    {
      JacPDstandardNth1alpha = PDstandardNth1alpha;
      
      JacPDstandardNth1gt11 = PDstandardNth1gt11;
      
      JacPDstandardNth1gt12 = PDstandardNth1gt12;
      
      JacPDstandardNth1gt13 = PDstandardNth1gt13;
      
      JacPDstandardNth1gt22 = PDstandardNth1gt22;
      
      JacPDstandardNth1gt23 = PDstandardNth1gt23;
      
      JacPDstandardNth1gt33 = PDstandardNth1gt33;
      
      JacPDstandardNth1phi = PDstandardNth1phi;
      
      JacPDstandardNth2alpha = PDstandardNth2alpha;
      
      JacPDstandardNth2gt11 = PDstandardNth2gt11;
      
      JacPDstandardNth2gt12 = PDstandardNth2gt12;
      
      JacPDstandardNth2gt13 = PDstandardNth2gt13;
      
      JacPDstandardNth2gt22 = PDstandardNth2gt22;
      
      JacPDstandardNth2gt23 = PDstandardNth2gt23;
      
      JacPDstandardNth2gt33 = PDstandardNth2gt33;
      
      JacPDstandardNth2phi = PDstandardNth2phi;
      
      JacPDstandardNth3alpha = PDstandardNth3alpha;
      
      JacPDstandardNth3gt11 = PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = PDstandardNth3gt33;
      
      JacPDstandardNth3phi = PDstandardNth3phi;
      
      JacPDupwindNthSymm1alpha = PDupwindNthSymm1alpha;
      
      JacPDupwindNthSymm1beta1 = PDupwindNthSymm1beta1;
      
      JacPDupwindNthSymm1beta2 = PDupwindNthSymm1beta2;
      
      JacPDupwindNthSymm1beta3 = PDupwindNthSymm1beta3;
      
      JacPDupwindNthSymm2alpha = PDupwindNthSymm2alpha;
      
      JacPDupwindNthSymm2beta1 = PDupwindNthSymm2beta1;
      
      JacPDupwindNthSymm2beta2 = PDupwindNthSymm2beta2;
      
      JacPDupwindNthSymm2beta3 = PDupwindNthSymm2beta3;
      
      JacPDupwindNthSymm3alpha = PDupwindNthSymm3alpha;
      
      JacPDupwindNthSymm3beta1 = PDupwindNthSymm3beta1;
      
      JacPDupwindNthSymm3beta2 = PDupwindNthSymm3beta2;
      
      JacPDupwindNthSymm3beta3 = PDupwindNthSymm3beta3;
      
      JacPDupwindNthAnti1alpha = PDupwindNthAnti1alpha;
      
      JacPDupwindNthAnti1beta1 = PDupwindNthAnti1beta1;
      
      JacPDupwindNthAnti1beta2 = PDupwindNthAnti1beta2;
      
      JacPDupwindNthAnti1beta3 = PDupwindNthAnti1beta3;
      
      JacPDupwindNthAnti2alpha = PDupwindNthAnti2alpha;
      
      JacPDupwindNthAnti2beta1 = PDupwindNthAnti2beta1;
      
      JacPDupwindNthAnti2beta2 = PDupwindNthAnti2beta2;
      
      JacPDupwindNthAnti2beta3 = PDupwindNthAnti2beta3;
      
      JacPDupwindNthAnti3alpha = PDupwindNthAnti3alpha;
      
      JacPDupwindNthAnti3beta1 = PDupwindNthAnti3beta1;
      
      JacPDupwindNthAnti3beta2 = PDupwindNthAnti3beta2;
      
      JacPDupwindNthAnti3beta3 = PDupwindNthAnti3beta3;
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
    
    CCTK_REAL em4phi CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod != 
      0,pow(phiL,2),exp(-4*phiL));
    
    CCTK_REAL gu11 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu11;
    
    CCTK_REAL gu12 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu12;
    
    CCTK_REAL gu13 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu13;
    
    CCTK_REAL gu22 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu22;
    
    CCTK_REAL gu23 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu23;
    
    CCTK_REAL gu33 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu33;
    
    CCTK_REAL fac1 CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod != 
      0,-0.5*pow(phiL,-1),1);
    
    CCTK_REAL cdphi1 CCTK_ATTRIBUTE_UNUSED = fac1*JacPDstandardNth1phi;
    
    CCTK_REAL cdphi2 CCTK_ATTRIBUTE_UNUSED = fac1*JacPDstandardNth2phi;
    
    CCTK_REAL cdphi3 CCTK_ATTRIBUTE_UNUSED = fac1*JacPDstandardNth3phi;
    
    CCTK_REAL dotalpha CCTK_ATTRIBUTE_UNUSED = -(harmonicF*IfThen(evolveA 
      != 0,AL,trKL + (-1 + alphaL)*alphaDriver - IfThen(formulation != 
      0,2*ThetaL,0))*pow(alphaL,harmonicN));
    
    CCTK_REAL betaDriverValue CCTK_ATTRIBUTE_UNUSED = 
      IfThen(useSpatialBetaDriver != 
      0,betaDriver*spatialBetaDriverRadius*pow(fmax(rL,spatialBetaDriverRadius),-1),betaDriver);
    
    CCTK_REAL shiftGammaCoeffValue CCTK_ATTRIBUTE_UNUSED = 
      IfThen(useSpatialShiftGammaCoeff != 0,shiftGammaCoeff*fmin(1,exp(1 - 
      rL*pow(spatialShiftGammaCoeffRadius,-1))),shiftGammaCoeff);
    
    CCTK_REAL ddetgt1 CCTK_ATTRIBUTE_UNUSED = gtu11*JacPDstandardNth1gt11 
      + gtu22*JacPDstandardNth1gt22 + 2*(gtu12*JacPDstandardNth1gt12 + 
      gtu13*JacPDstandardNth1gt13 + gtu23*JacPDstandardNth1gt23) + 
      gtu33*JacPDstandardNth1gt33;
    
    CCTK_REAL ddetgt2 CCTK_ATTRIBUTE_UNUSED = gtu11*JacPDstandardNth2gt11 
      + gtu22*JacPDstandardNth2gt22 + 2*(gtu12*JacPDstandardNth2gt12 + 
      gtu13*JacPDstandardNth2gt13 + gtu23*JacPDstandardNth2gt23) + 
      gtu33*JacPDstandardNth2gt33;
    
    CCTK_REAL ddetgt3 CCTK_ATTRIBUTE_UNUSED = gtu11*JacPDstandardNth3gt11 
      + gtu22*JacPDstandardNth3gt22 + 2*(gtu12*JacPDstandardNth3gt12 + 
      gtu13*JacPDstandardNth3gt13 + gtu23*JacPDstandardNth3gt23) + 
      gtu33*JacPDstandardNth3gt33;
    
    CCTK_REAL dotbeta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL dotbeta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL dotbeta3 CCTK_ATTRIBUTE_UNUSED;
    
    if (shiftFormulation == 0)
    {
      dotbeta1 = shiftGammaCoeffValue*IfThen(evolveB != 0,B1L,Xt1L - 
        beta1L*betaDriverValue)*pow(alphaL,shiftAlphaPower);
      
      dotbeta2 = shiftGammaCoeffValue*IfThen(evolveB != 0,B2L,Xt2L - 
        beta2L*betaDriverValue)*pow(alphaL,shiftAlphaPower);
      
      dotbeta3 = shiftGammaCoeffValue*IfThen(evolveB != 0,B3L,Xt3L - 
        beta3L*betaDriverValue)*pow(alphaL,shiftAlphaPower);
    }
    else
    {
      dotbeta1 = -(alphaL*(gu11*(JacPDstandardNth1alpha + alphaL*(2*cdphi1 + 
        0.5*ddetgt1 - gtu11*JacPDstandardNth1gt11 - 
        gtu12*(JacPDstandardNth1gt12 + JacPDstandardNth2gt11) - 
        gtu22*JacPDstandardNth2gt12 + gtu13*(-JacPDstandardNth1gt13 - 
        JacPDstandardNth3gt11) + gtu23*(-JacPDstandardNth2gt13 - 
        JacPDstandardNth3gt12) - gtu33*JacPDstandardNth3gt13)) + 
        gu12*(JacPDstandardNth2alpha + alphaL*(2*cdphi2 + 0.5*ddetgt2 - 
        gtu11*JacPDstandardNth1gt12 - gtu12*(JacPDstandardNth1gt22 + 
        JacPDstandardNth2gt12) - gtu22*JacPDstandardNth2gt22 + 
        gtu13*(-JacPDstandardNth1gt23 - JacPDstandardNth3gt12) + 
        gtu23*(-JacPDstandardNth2gt23 - JacPDstandardNth3gt22) - 
        gtu33*JacPDstandardNth3gt23)) + gu13*(JacPDstandardNth3alpha + 
        alphaL*(2*cdphi3 + 0.5*ddetgt3 - gtu11*JacPDstandardNth1gt13 - 
        gtu12*(JacPDstandardNth1gt23 + JacPDstandardNth2gt13) - 
        gtu22*JacPDstandardNth2gt23 + gtu13*(-JacPDstandardNth1gt33 - 
        JacPDstandardNth3gt13) + gtu23*(-JacPDstandardNth2gt33 - 
        JacPDstandardNth3gt23) - gtu33*JacPDstandardNth3gt33))));
      
      dotbeta2 = -(alphaL*(gu12*(JacPDstandardNth1alpha + alphaL*(2*cdphi1 + 
        0.5*ddetgt1 - gtu11*JacPDstandardNth1gt11 - 
        gtu12*(JacPDstandardNth1gt12 + JacPDstandardNth2gt11) - 
        gtu22*JacPDstandardNth2gt12 + gtu13*(-JacPDstandardNth1gt13 - 
        JacPDstandardNth3gt11) + gtu23*(-JacPDstandardNth2gt13 - 
        JacPDstandardNth3gt12) - gtu33*JacPDstandardNth3gt13)) + 
        gu22*(JacPDstandardNth2alpha + alphaL*(2*cdphi2 + 0.5*ddetgt2 - 
        gtu11*JacPDstandardNth1gt12 - gtu12*(JacPDstandardNth1gt22 + 
        JacPDstandardNth2gt12) - gtu22*JacPDstandardNth2gt22 + 
        gtu13*(-JacPDstandardNth1gt23 - JacPDstandardNth3gt12) + 
        gtu23*(-JacPDstandardNth2gt23 - JacPDstandardNth3gt22) - 
        gtu33*JacPDstandardNth3gt23)) + gu23*(JacPDstandardNth3alpha + 
        alphaL*(2*cdphi3 + 0.5*ddetgt3 - gtu11*JacPDstandardNth1gt13 - 
        gtu12*(JacPDstandardNth1gt23 + JacPDstandardNth2gt13) - 
        gtu22*JacPDstandardNth2gt23 + gtu13*(-JacPDstandardNth1gt33 - 
        JacPDstandardNth3gt13) + gtu23*(-JacPDstandardNth2gt33 - 
        JacPDstandardNth3gt23) - gtu33*JacPDstandardNth3gt33))));
      
      dotbeta3 = -(alphaL*(gu13*(JacPDstandardNth1alpha + alphaL*(2*cdphi1 + 
        0.5*ddetgt1 - gtu11*JacPDstandardNth1gt11 - 
        gtu12*(JacPDstandardNth1gt12 + JacPDstandardNth2gt11) - 
        gtu22*JacPDstandardNth2gt12 + gtu13*(-JacPDstandardNth1gt13 - 
        JacPDstandardNth3gt11) + gtu23*(-JacPDstandardNth2gt13 - 
        JacPDstandardNth3gt12) - gtu33*JacPDstandardNth3gt13)) + 
        gu23*(JacPDstandardNth2alpha + alphaL*(2*cdphi2 + 0.5*ddetgt2 - 
        gtu11*JacPDstandardNth1gt12 - gtu12*(JacPDstandardNth1gt22 + 
        JacPDstandardNth2gt12) - gtu22*JacPDstandardNth2gt22 + 
        gtu13*(-JacPDstandardNth1gt23 - JacPDstandardNth3gt12) + 
        gtu23*(-JacPDstandardNth2gt23 - JacPDstandardNth3gt22) - 
        gtu33*JacPDstandardNth3gt23)) + gu33*(JacPDstandardNth3alpha + 
        alphaL*(2*cdphi3 + 0.5*ddetgt3 - gtu11*JacPDstandardNth1gt13 - 
        gtu12*(JacPDstandardNth1gt23 + JacPDstandardNth2gt13) - 
        gtu22*JacPDstandardNth2gt23 + gtu13*(-JacPDstandardNth1gt33 - 
        JacPDstandardNth3gt13) + gtu23*(-JacPDstandardNth2gt33 - 
        JacPDstandardNth3gt23) - gtu33*JacPDstandardNth3gt33))));
    }
    
    CCTK_REAL dtalpL CCTK_ATTRIBUTE_UNUSED = dotalpha + IfThen(advectLapse 
      != 0,beta1L*JacPDupwindNthAnti1alpha + beta2L*JacPDupwindNthAnti2alpha 
      + beta3L*JacPDupwindNthAnti3alpha + 
      JacPDupwindNthSymm1alpha*fabs(beta1L) + 
      JacPDupwindNthSymm2alpha*fabs(beta2L) + 
      JacPDupwindNthSymm3alpha*fabs(beta3L),0);
    
    CCTK_REAL dtbetaxL CCTK_ATTRIBUTE_UNUSED = dotbeta1 + 
      IfThen(advectShift != 0,beta1L*JacPDupwindNthAnti1beta1 + 
      beta2L*JacPDupwindNthAnti2beta1 + beta3L*JacPDupwindNthAnti3beta1 + 
      JacPDupwindNthSymm1beta1*fabs(beta1L) + 
      JacPDupwindNthSymm2beta1*fabs(beta2L) + 
      JacPDupwindNthSymm3beta1*fabs(beta3L),0);
    
    CCTK_REAL dtbetayL CCTK_ATTRIBUTE_UNUSED = dotbeta2 + 
      IfThen(advectShift != 0,beta1L*JacPDupwindNthAnti1beta2 + 
      beta2L*JacPDupwindNthAnti2beta2 + beta3L*JacPDupwindNthAnti3beta2 + 
      JacPDupwindNthSymm1beta2*fabs(beta1L) + 
      JacPDupwindNthSymm2beta2*fabs(beta2L) + 
      JacPDupwindNthSymm3beta2*fabs(beta3L),0);
    
    CCTK_REAL dtbetazL CCTK_ATTRIBUTE_UNUSED = dotbeta3 + 
      IfThen(advectShift != 0,beta1L*JacPDupwindNthAnti1beta3 + 
      beta2L*JacPDupwindNthAnti2beta3 + beta3L*JacPDupwindNthAnti3beta3 + 
      JacPDupwindNthSymm1beta3*fabs(beta1L) + 
      JacPDupwindNthSymm2beta3*fabs(beta2L) + 
      JacPDupwindNthSymm3beta3*fabs(beta3L),0);
    /* Copy local copies back to grid functions */
    dtalp[index] = dtalpL;
    dtbetax[index] = dtbetaxL;
    dtbetay[index] = dtbetayL;
    dtbetaz[index] = dtbetazL;
  }
  CCTK_ENDLOOP3(ML_BSSN_NV_ADMBaseInterior);
}
extern "C" void ML_BSSN_NV_ADMBaseInterior(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_BSSN_NV_ADMBaseInterior
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_NV_ADMBaseInterior);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_NV_ADMBaseInterior_Body");
  }
  if (cctk_iteration % ML_BSSN_NV_ADMBaseInterior_calc_every != ML_BSSN_NV_ADMBaseInterior_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ADMBase::dtlapse",
    "ADMBase::dtshift",
    "grid::coordinates",
    "ML_BSSN_NV::ML_dtlapse",
    "ML_BSSN_NV::ML_dtshift",
    "ML_BSSN_NV::ML_Gamma",
    "ML_BSSN_NV::ML_lapse",
    "ML_BSSN_NV::ML_log_confac",
    "ML_BSSN_NV::ML_metric",
    "ML_BSSN_NV::ML_shift",
    "ML_BSSN_NV::ML_Theta",
    "ML_BSSN_NV::ML_trace_curv"};
  AssertGroupStorage(cctkGH, "ML_BSSN_NV_ADMBaseInterior", 12, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_NV_ADMBaseInterior", 2, 2, 2);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_NV_ADMBaseInterior", 3, 3, 3);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_NV_ADMBaseInterior", 4, 4, 4);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_NV_ADMBaseInterior", 5, 5, 5);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, ML_BSSN_NV_ADMBaseInterior_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_NV_ADMBaseInterior_Body");
  }
}

} // namespace ML_BSSN_NV
