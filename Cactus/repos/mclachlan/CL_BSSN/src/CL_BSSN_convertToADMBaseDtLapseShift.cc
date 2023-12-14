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
#include "vectors.h"

namespace CL_BSSN {

extern "C" void CL_BSSN_convertToADMBaseDtLapseShift_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CL_BSSN_convertToADMBaseDtLapseShift_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(CL_BSSN_convertToADMBaseDtLapseShift_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % CL_BSSN_convertToADMBaseDtLapseShift_calc_every != CL_BSSN_convertToADMBaseDtLapseShift_calc_offset)
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

static void CL_BSSN_convertToADMBaseDtLapseShift_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  const CCTK_REAL_VEC t CCTK_ATTRIBUTE_UNUSED = ToReal(cctk_time);
  const CCTK_REAL_VEC cctkOriginSpace1 CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_ORIGIN_SPACE(0));
  const CCTK_REAL_VEC cctkOriginSpace2 CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_ORIGIN_SPACE(1));
  const CCTK_REAL_VEC cctkOriginSpace3 CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_ORIGIN_SPACE(2));
  const CCTK_REAL_VEC dt CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_TIME);
  const CCTK_REAL_VEC dx CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_SPACE(0));
  const CCTK_REAL_VEC dy CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_SPACE(1));
  const CCTK_REAL_VEC dz CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_SPACE(2));
  const CCTK_REAL_VEC dxi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dx);
  const CCTK_REAL_VEC dyi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dy);
  const CCTK_REAL_VEC dzi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dz);
  const CCTK_REAL_VEC khalf CCTK_ATTRIBUTE_UNUSED = ToReal(0.5);
  const CCTK_REAL_VEC kthird CCTK_ATTRIBUTE_UNUSED = 
    ToReal(0.333333333333333333333333333333);
  const CCTK_REAL_VEC ktwothird CCTK_ATTRIBUTE_UNUSED = 
    ToReal(0.666666666666666666666666666667);
  const CCTK_REAL_VEC kfourthird CCTK_ATTRIBUTE_UNUSED = 
    ToReal(1.33333333333333333333333333333);
  const CCTK_REAL_VEC hdxi CCTK_ATTRIBUTE_UNUSED = 
    kmul(dxi,ToReal(0.5));
  const CCTK_REAL_VEC hdyi CCTK_ATTRIBUTE_UNUSED = 
    kmul(dyi,ToReal(0.5));
  const CCTK_REAL_VEC hdzi CCTK_ATTRIBUTE_UNUSED = 
    kmul(dzi,ToReal(0.5));
  /* Initialize predefined quantities */
  const CCTK_REAL_VEC p1o1024dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0009765625),dx);
  const CCTK_REAL_VEC p1o1024dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0009765625),dy);
  const CCTK_REAL_VEC p1o1024dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0009765625),dz);
  const CCTK_REAL_VEC p1o120dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00833333333333333333333333333333),dx);
  const CCTK_REAL_VEC p1o120dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00833333333333333333333333333333),dy);
  const CCTK_REAL_VEC p1o120dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00833333333333333333333333333333),dz);
  const CCTK_REAL_VEC p1o12dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dx);
  const CCTK_REAL_VEC p1o12dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dy);
  const CCTK_REAL_VEC p1o12dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dz);
  const CCTK_REAL_VEC p1o144dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dx,dy));
  const CCTK_REAL_VEC p1o144dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dx,dz));
  const CCTK_REAL_VEC p1o144dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dy,dz));
  const CCTK_REAL_VEC p1o1680dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000595238095238095238095238095238),dx);
  const CCTK_REAL_VEC p1o1680dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000595238095238095238095238095238),dy);
  const CCTK_REAL_VEC p1o1680dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000595238095238095238095238095238),dz);
  const CCTK_REAL_VEC p1o16dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0625),dx);
  const CCTK_REAL_VEC p1o16dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0625),dy);
  const CCTK_REAL_VEC p1o16dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0625),dz);
  const CCTK_REAL_VEC p1o180dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dx,dx));
  const CCTK_REAL_VEC p1o180dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dy,dy));
  const CCTK_REAL_VEC p1o180dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dz,dz));
  const CCTK_REAL_VEC p1o24dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dx);
  const CCTK_REAL_VEC p1o24dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dy);
  const CCTK_REAL_VEC p1o24dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dz);
  const CCTK_REAL_VEC p1o256dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00390625),dx);
  const CCTK_REAL_VEC p1o256dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00390625),dy);
  const CCTK_REAL_VEC p1o256dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00390625),dz);
  const CCTK_REAL_VEC p1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dx);
  const CCTK_REAL_VEC p1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dy);
  const CCTK_REAL_VEC p1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dz);
  const CCTK_REAL_VEC p1o3600dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dx,dy));
  const CCTK_REAL_VEC p1o3600dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dx,dz));
  const CCTK_REAL_VEC p1o3600dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dy,dz));
  const CCTK_REAL_VEC p1o4dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dx);
  const CCTK_REAL_VEC p1o4dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dx,dy));
  const CCTK_REAL_VEC p1o4dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dx,dz));
  const CCTK_REAL_VEC p1o4dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dy);
  const CCTK_REAL_VEC p1o4dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dy,dz));
  const CCTK_REAL_VEC p1o4dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dz);
  const CCTK_REAL_VEC p1o5040dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dx,dx));
  const CCTK_REAL_VEC p1o5040dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dy,dy));
  const CCTK_REAL_VEC p1o5040dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dz,dz));
  const CCTK_REAL_VEC p1o560dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00178571428571428571428571428571),dx);
  const CCTK_REAL_VEC p1o560dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00178571428571428571428571428571),dy);
  const CCTK_REAL_VEC p1o560dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00178571428571428571428571428571),dz);
  const CCTK_REAL_VEC p1o60dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0166666666666666666666666666667),dx);
  const CCTK_REAL_VEC p1o60dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0166666666666666666666666666667),dy);
  const CCTK_REAL_VEC p1o60dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0166666666666666666666666666667),dz);
  const CCTK_REAL_VEC p1o64dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.015625),dx);
  const CCTK_REAL_VEC p1o64dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.015625),dy);
  const CCTK_REAL_VEC p1o64dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.015625),dz);
  const CCTK_REAL_VEC p1o705600dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dx,dy));
  const CCTK_REAL_VEC p1o705600dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dx,dz));
  const CCTK_REAL_VEC p1o705600dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dy,dz));
  const CCTK_REAL_VEC p1o840dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dx);
  const CCTK_REAL_VEC p1o840dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dy);
  const CCTK_REAL_VEC p1o840dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dz);
  const CCTK_REAL_VEC p1odx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dx);
  const CCTK_REAL_VEC p1odx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dx,dx));
  const CCTK_REAL_VEC p1ody CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dy);
  const CCTK_REAL_VEC p1ody2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dy,dy));
  const CCTK_REAL_VEC p1odz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dz);
  const CCTK_REAL_VEC p1odz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dz,dz));
  const CCTK_REAL_VEC pm1o120dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00833333333333333333333333333333),dx);
  const CCTK_REAL_VEC pm1o120dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00833333333333333333333333333333),dy);
  const CCTK_REAL_VEC pm1o120dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00833333333333333333333333333333),dz);
  const CCTK_REAL_VEC pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dx,dx));
  const CCTK_REAL_VEC pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dy,dy));
  const CCTK_REAL_VEC pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dz,dz));
  const CCTK_REAL_VEC pm1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dx);
  const CCTK_REAL_VEC pm1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dy);
  const CCTK_REAL_VEC pm1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dz);
  const CCTK_REAL_VEC pm1o4dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.25),dx);
  const CCTK_REAL_VEC pm1o4dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.25),dy);
  const CCTK_REAL_VEC pm1o4dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.25),dz);
  const CCTK_REAL_VEC pm1o60dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0166666666666666666666666666667),dx);
  const CCTK_REAL_VEC pm1o60dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0166666666666666666666666666667),dy);
  const CCTK_REAL_VEC pm1o60dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0166666666666666666666666666667),dz);
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
  CCTK_LOOP3STR(CL_BSSN_convertToADMBaseDtLapseShift,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alphaL CCTK_ATTRIBUTE_UNUSED = vec_load(alpha[index]);
    CCTK_REAL_VEC B1L CCTK_ATTRIBUTE_UNUSED = vec_load(B1[index]);
    CCTK_REAL_VEC B2L CCTK_ATTRIBUTE_UNUSED = vec_load(B2[index]);
    CCTK_REAL_VEC B3L CCTK_ATTRIBUTE_UNUSED = vec_load(B3[index]);
    CCTK_REAL_VEC beta1L CCTK_ATTRIBUTE_UNUSED = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L CCTK_ATTRIBUTE_UNUSED = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L CCTK_ATTRIBUTE_UNUSED = vec_load(beta3[index]);
    CCTK_REAL_VEC trKL CCTK_ATTRIBUTE_UNUSED = vec_load(trK[index]);
    
    
    CCTK_REAL_VEC J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L CCTK_ATTRIBUTE_UNUSED ;
    
    if (usejacobian)
    {
      J11L = vec_load(J11[index]);
      J12L = vec_load(J12[index]);
      J13L = vec_load(J13[index]);
      J21L = vec_load(J21[index]);
      J22L = vec_load(J22[index]);
      J23L = vec_load(J23[index]);
      J31L = vec_load(J31[index]);
      J32L = vec_load(J32[index]);
      J33L = vec_load(J33[index]);
    }
    /* Include user supplied include files */
    /* Precompute derivatives */
    CCTK_REAL_VEC PDupwindNthAnti1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3beta3 CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder21(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder21(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder22(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder22(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder23(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder23(&alpha[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder21(&beta1[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder21(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder22(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder22(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder23(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder23(&beta1[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder21(&beta2[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder21(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder22(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder22(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder23(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder23(&beta2[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder21(&beta3[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder21(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder22(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder22(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder23(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder23(&beta3[index]);
        break;
      }
      
      case 4:
      {
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder41(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder41(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder42(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder42(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder43(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder43(&alpha[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder41(&beta1[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder41(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder42(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder42(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder43(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder43(&beta1[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder41(&beta2[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder41(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder42(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder42(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder43(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder43(&beta2[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder41(&beta3[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder41(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder42(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder42(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder43(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder43(&beta3[index]);
        break;
      }
      
      case 6:
      {
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder61(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder61(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder62(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder62(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder63(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder63(&alpha[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder61(&beta1[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder61(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder62(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder62(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder63(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder63(&beta1[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder61(&beta2[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder61(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder62(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder62(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder63(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder63(&beta2[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder61(&beta3[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder61(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder62(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder62(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder63(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder63(&beta3[index]);
        break;
      }
      
      case 8:
      {
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder81(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder81(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder82(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder82(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder83(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder83(&alpha[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder81(&beta1[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder81(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder82(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder82(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder83(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder83(&beta1[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder81(&beta2[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder81(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder82(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder82(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder83(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder83(&beta2[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder81(&beta3[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder81(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder82(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder82(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder83(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder83(&beta3[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = kisgn(beta1L);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = kisgn(beta2L);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = kisgn(beta3L);
    
    CCTK_REAL_VEC JacPDupwindNthAnti1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3beta3 CCTK_ATTRIBUTE_UNUSED;
    
    if (usejacobian)
    {
      JacPDupwindNthAnti1alpha = 
        kmadd(J11L,PDupwindNthAnti1alpha,kmadd(J21L,PDupwindNthAnti2alpha,kmul(J31L,PDupwindNthAnti3alpha)));
      
      JacPDupwindNthAnti1beta1 = 
        kmadd(J11L,PDupwindNthAnti1beta1,kmadd(J21L,PDupwindNthAnti2beta1,kmul(J31L,PDupwindNthAnti3beta1)));
      
      JacPDupwindNthAnti1beta2 = 
        kmadd(J11L,PDupwindNthAnti1beta2,kmadd(J21L,PDupwindNthAnti2beta2,kmul(J31L,PDupwindNthAnti3beta2)));
      
      JacPDupwindNthAnti1beta3 = 
        kmadd(J11L,PDupwindNthAnti1beta3,kmadd(J21L,PDupwindNthAnti2beta3,kmul(J31L,PDupwindNthAnti3beta3)));
      
      JacPDupwindNthSymm1alpha = 
        kmadd(J11L,PDupwindNthSymm1alpha,kmadd(J21L,PDupwindNthSymm2alpha,kmul(J31L,PDupwindNthSymm3alpha)));
      
      JacPDupwindNthSymm1beta1 = 
        kmadd(J11L,PDupwindNthSymm1beta1,kmadd(J21L,PDupwindNthSymm2beta1,kmul(J31L,PDupwindNthSymm3beta1)));
      
      JacPDupwindNthSymm1beta2 = 
        kmadd(J11L,PDupwindNthSymm1beta2,kmadd(J21L,PDupwindNthSymm2beta2,kmul(J31L,PDupwindNthSymm3beta2)));
      
      JacPDupwindNthSymm1beta3 = 
        kmadd(J11L,PDupwindNthSymm1beta3,kmadd(J21L,PDupwindNthSymm2beta3,kmul(J31L,PDupwindNthSymm3beta3)));
      
      JacPDupwindNthAnti2alpha = 
        kmadd(J12L,PDupwindNthAnti1alpha,kmadd(J22L,PDupwindNthAnti2alpha,kmul(J32L,PDupwindNthAnti3alpha)));
      
      JacPDupwindNthAnti2beta1 = 
        kmadd(J12L,PDupwindNthAnti1beta1,kmadd(J22L,PDupwindNthAnti2beta1,kmul(J32L,PDupwindNthAnti3beta1)));
      
      JacPDupwindNthAnti2beta2 = 
        kmadd(J12L,PDupwindNthAnti1beta2,kmadd(J22L,PDupwindNthAnti2beta2,kmul(J32L,PDupwindNthAnti3beta2)));
      
      JacPDupwindNthAnti2beta3 = 
        kmadd(J12L,PDupwindNthAnti1beta3,kmadd(J22L,PDupwindNthAnti2beta3,kmul(J32L,PDupwindNthAnti3beta3)));
      
      JacPDupwindNthSymm2alpha = 
        kmadd(J12L,PDupwindNthSymm1alpha,kmadd(J22L,PDupwindNthSymm2alpha,kmul(J32L,PDupwindNthSymm3alpha)));
      
      JacPDupwindNthSymm2beta1 = 
        kmadd(J12L,PDupwindNthSymm1beta1,kmadd(J22L,PDupwindNthSymm2beta1,kmul(J32L,PDupwindNthSymm3beta1)));
      
      JacPDupwindNthSymm2beta2 = 
        kmadd(J12L,PDupwindNthSymm1beta2,kmadd(J22L,PDupwindNthSymm2beta2,kmul(J32L,PDupwindNthSymm3beta2)));
      
      JacPDupwindNthSymm2beta3 = 
        kmadd(J12L,PDupwindNthSymm1beta3,kmadd(J22L,PDupwindNthSymm2beta3,kmul(J32L,PDupwindNthSymm3beta3)));
      
      JacPDupwindNthAnti3alpha = 
        kmadd(J13L,PDupwindNthAnti1alpha,kmadd(J23L,PDupwindNthAnti2alpha,kmul(J33L,PDupwindNthAnti3alpha)));
      
      JacPDupwindNthAnti3beta1 = 
        kmadd(J13L,PDupwindNthAnti1beta1,kmadd(J23L,PDupwindNthAnti2beta1,kmul(J33L,PDupwindNthAnti3beta1)));
      
      JacPDupwindNthAnti3beta2 = 
        kmadd(J13L,PDupwindNthAnti1beta2,kmadd(J23L,PDupwindNthAnti2beta2,kmul(J33L,PDupwindNthAnti3beta2)));
      
      JacPDupwindNthAnti3beta3 = 
        kmadd(J13L,PDupwindNthAnti1beta3,kmadd(J23L,PDupwindNthAnti2beta3,kmul(J33L,PDupwindNthAnti3beta3)));
      
      JacPDupwindNthSymm3alpha = 
        kmadd(J13L,PDupwindNthSymm1alpha,kmadd(J23L,PDupwindNthSymm2alpha,kmul(J33L,PDupwindNthSymm3alpha)));
      
      JacPDupwindNthSymm3beta1 = 
        kmadd(J13L,PDupwindNthSymm1beta1,kmadd(J23L,PDupwindNthSymm2beta1,kmul(J33L,PDupwindNthSymm3beta1)));
      
      JacPDupwindNthSymm3beta2 = 
        kmadd(J13L,PDupwindNthSymm1beta2,kmadd(J23L,PDupwindNthSymm2beta2,kmul(J33L,PDupwindNthSymm3beta2)));
      
      JacPDupwindNthSymm3beta3 = 
        kmadd(J13L,PDupwindNthSymm1beta3,kmadd(J23L,PDupwindNthSymm2beta3,kmul(J33L,PDupwindNthSymm3beta3)));
    }
    else
    {
      JacPDupwindNthAnti1alpha = PDupwindNthAnti1alpha;
      
      JacPDupwindNthAnti1beta1 = PDupwindNthAnti1beta1;
      
      JacPDupwindNthAnti1beta2 = PDupwindNthAnti1beta2;
      
      JacPDupwindNthAnti1beta3 = PDupwindNthAnti1beta3;
      
      JacPDupwindNthSymm1alpha = PDupwindNthSymm1alpha;
      
      JacPDupwindNthSymm1beta1 = PDupwindNthSymm1beta1;
      
      JacPDupwindNthSymm1beta2 = PDupwindNthSymm1beta2;
      
      JacPDupwindNthSymm1beta3 = PDupwindNthSymm1beta3;
      
      JacPDupwindNthAnti2alpha = PDupwindNthAnti2alpha;
      
      JacPDupwindNthAnti2beta1 = PDupwindNthAnti2beta1;
      
      JacPDupwindNthAnti2beta2 = PDupwindNthAnti2beta2;
      
      JacPDupwindNthAnti2beta3 = PDupwindNthAnti2beta3;
      
      JacPDupwindNthSymm2alpha = PDupwindNthSymm2alpha;
      
      JacPDupwindNthSymm2beta1 = PDupwindNthSymm2beta1;
      
      JacPDupwindNthSymm2beta2 = PDupwindNthSymm2beta2;
      
      JacPDupwindNthSymm2beta3 = PDupwindNthSymm2beta3;
      
      JacPDupwindNthAnti3alpha = PDupwindNthAnti3alpha;
      
      JacPDupwindNthAnti3beta1 = PDupwindNthAnti3beta1;
      
      JacPDupwindNthAnti3beta2 = PDupwindNthAnti3beta2;
      
      JacPDupwindNthAnti3beta3 = PDupwindNthAnti3beta3;
      
      JacPDupwindNthSymm3alpha = PDupwindNthSymm3alpha;
      
      JacPDupwindNthSymm3beta1 = PDupwindNthSymm3beta1;
      
      JacPDupwindNthSymm3beta2 = PDupwindNthSymm3beta2;
      
      JacPDupwindNthSymm3beta3 = PDupwindNthSymm3beta3;
    }
    
    CCTK_REAL_VEC dtalpL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(beta1L,JacPDupwindNthAnti1alpha,kmadd(beta2L,JacPDupwindNthAnti2alpha,kmadd(beta3L,JacPDupwindNthAnti3alpha,kmadd(JacPDupwindNthSymm1alpha,kfabs(beta1L),kmadd(JacPDupwindNthSymm2alpha,kfabs(beta2L),kmsub(JacPDupwindNthSymm3alpha,kfabs(beta3L),kmul(kmul(trKL,kpow(alphaL,harmonicN)),ToReal(harmonicF))))))));
    
    CCTK_REAL_VEC dtbetaxL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(beta1L,JacPDupwindNthAnti1beta1,kmadd(beta2L,JacPDupwindNthAnti2beta1,kmadd(beta3L,JacPDupwindNthAnti3beta1,kmadd(JacPDupwindNthSymm1beta1,kfabs(beta1L),kmadd(JacPDupwindNthSymm2beta1,kfabs(beta2L),kmadd(JacPDupwindNthSymm3beta1,kfabs(beta3L),kmul(B1L,ToReal(ShiftGammaCoeff))))))));
    
    CCTK_REAL_VEC dtbetayL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(beta1L,JacPDupwindNthAnti1beta2,kmadd(beta2L,JacPDupwindNthAnti2beta2,kmadd(beta3L,JacPDupwindNthAnti3beta2,kmadd(JacPDupwindNthSymm1beta2,kfabs(beta1L),kmadd(JacPDupwindNthSymm2beta2,kfabs(beta2L),kmadd(JacPDupwindNthSymm3beta2,kfabs(beta3L),kmul(B2L,ToReal(ShiftGammaCoeff))))))));
    
    CCTK_REAL_VEC dtbetazL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(beta1L,JacPDupwindNthAnti1beta3,kmadd(beta2L,JacPDupwindNthAnti2beta3,kmadd(beta3L,JacPDupwindNthAnti3beta3,kmadd(JacPDupwindNthSymm1beta3,kfabs(beta1L),kmadd(JacPDupwindNthSymm2beta3,kfabs(beta2L),kmadd(JacPDupwindNthSymm3beta3,kfabs(beta3L),kmul(B3L,ToReal(ShiftGammaCoeff))))))));
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(dtalp[index],dtalpL);
    vec_store_nta_partial(dtbetax[index],dtbetaxL);
    vec_store_nta_partial(dtbetay[index],dtbetayL);
    vec_store_nta_partial(dtbetaz[index],dtbetazL);
  }
  CCTK_ENDLOOP3STR(CL_BSSN_convertToADMBaseDtLapseShift);
}
extern "C" void CL_BSSN_convertToADMBaseDtLapseShift(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CL_BSSN_convertToADMBaseDtLapseShift
  DECLARE_CCTK_ARGUMENTS_CHECKED(CL_BSSN_convertToADMBaseDtLapseShift);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering CL_BSSN_convertToADMBaseDtLapseShift_Body");
  }
  if (cctk_iteration % CL_BSSN_convertToADMBaseDtLapseShift_calc_every != CL_BSSN_convertToADMBaseDtLapseShift_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ADMBase::dtlapse",
    "ADMBase::dtshift",
    "CL_BSSN::CL_dtshift",
    "CL_BSSN::CL_lapse",
    "CL_BSSN::CL_shift",
    "CL_BSSN::CL_trace_curv"};
  AssertGroupStorage(cctkGH, "CL_BSSN_convertToADMBaseDtLapseShift", 6, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "CL_BSSN_convertToADMBaseDtLapseShift", 2, 2, 2);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "CL_BSSN_convertToADMBaseDtLapseShift", 3, 3, 3);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "CL_BSSN_convertToADMBaseDtLapseShift", 4, 4, 4);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "CL_BSSN_convertToADMBaseDtLapseShift", 5, 5, 5);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, CL_BSSN_convertToADMBaseDtLapseShift_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving CL_BSSN_convertToADMBaseDtLapseShift_Body");
  }
}

} // namespace CL_BSSN
