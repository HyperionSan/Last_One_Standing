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

extern "C" void CL_BSSN_RHS1_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CL_BSSN_RHS1_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(CL_BSSN_RHS1_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % CL_BSSN_RHS1_calc_every != CL_BSSN_RHS1_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CL_BSSN::CL_cons_dmetric","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CL_BSSN::CL_cons_dmetric.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CL_BSSN::CL_dmetricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CL_BSSN::CL_dmetricrhs.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CL_BSSN::CL_Gammarhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CL_BSSN::CL_Gammarhs.");
  return;
}

static void CL_BSSN_RHS1_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3STR(CL_BSSN_RHS1,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alphaL CCTK_ATTRIBUTE_UNUSED = vec_load(alpha[index]);
    CCTK_REAL_VEC At11L CCTK_ATTRIBUTE_UNUSED = vec_load(At11[index]);
    CCTK_REAL_VEC At12L CCTK_ATTRIBUTE_UNUSED = vec_load(At12[index]);
    CCTK_REAL_VEC At13L CCTK_ATTRIBUTE_UNUSED = vec_load(At13[index]);
    CCTK_REAL_VEC At22L CCTK_ATTRIBUTE_UNUSED = vec_load(At22[index]);
    CCTK_REAL_VEC At23L CCTK_ATTRIBUTE_UNUSED = vec_load(At23[index]);
    CCTK_REAL_VEC At33L CCTK_ATTRIBUTE_UNUSED = vec_load(At33[index]);
    CCTK_REAL_VEC beta1L CCTK_ATTRIBUTE_UNUSED = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L CCTK_ATTRIBUTE_UNUSED = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L CCTK_ATTRIBUTE_UNUSED = vec_load(beta3[index]);
    CCTK_REAL_VEC cdgt111L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt111[index]);
    CCTK_REAL_VEC cdgt112L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt112[index]);
    CCTK_REAL_VEC cdgt113L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt113[index]);
    CCTK_REAL_VEC cdgt122L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt122[index]);
    CCTK_REAL_VEC cdgt123L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt123[index]);
    CCTK_REAL_VEC cdgt133L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt133[index]);
    CCTK_REAL_VEC cdgt211L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt211[index]);
    CCTK_REAL_VEC cdgt212L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt212[index]);
    CCTK_REAL_VEC cdgt213L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt213[index]);
    CCTK_REAL_VEC cdgt222L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt222[index]);
    CCTK_REAL_VEC cdgt223L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt223[index]);
    CCTK_REAL_VEC cdgt233L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt233[index]);
    CCTK_REAL_VEC cdgt311L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt311[index]);
    CCTK_REAL_VEC cdgt312L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt312[index]);
    CCTK_REAL_VEC cdgt313L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt313[index]);
    CCTK_REAL_VEC cdgt322L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt322[index]);
    CCTK_REAL_VEC cdgt323L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt323[index]);
    CCTK_REAL_VEC cdgt333L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(cdgt333[index]);
    CCTK_REAL_VEC dalpha1L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(dalpha1[index]);
    CCTK_REAL_VEC dalpha2L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(dalpha2[index]);
    CCTK_REAL_VEC dalpha3L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(dalpha3[index]);
    CCTK_REAL_VEC dbeta11L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(dbeta11[index]);
    CCTK_REAL_VEC dbeta12L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(dbeta12[index]);
    CCTK_REAL_VEC dbeta13L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(dbeta13[index]);
    CCTK_REAL_VEC dbeta21L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(dbeta21[index]);
    CCTK_REAL_VEC dbeta22L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(dbeta22[index]);
    CCTK_REAL_VEC dbeta23L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(dbeta23[index]);
    CCTK_REAL_VEC dbeta31L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(dbeta31[index]);
    CCTK_REAL_VEC dbeta32L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(dbeta32[index]);
    CCTK_REAL_VEC dbeta33L CCTK_ATTRIBUTE_UNUSED = 
      vec_load(dbeta33[index]);
    CCTK_REAL_VEC dgt111L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt111[index]);
    CCTK_REAL_VEC dgt112L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt112[index]);
    CCTK_REAL_VEC dgt113L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt113[index]);
    CCTK_REAL_VEC dgt122L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt122[index]);
    CCTK_REAL_VEC dgt123L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt123[index]);
    CCTK_REAL_VEC dgt133L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt133[index]);
    CCTK_REAL_VEC dgt211L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt211[index]);
    CCTK_REAL_VEC dgt212L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt212[index]);
    CCTK_REAL_VEC dgt213L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt213[index]);
    CCTK_REAL_VEC dgt222L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt222[index]);
    CCTK_REAL_VEC dgt223L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt223[index]);
    CCTK_REAL_VEC dgt233L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt233[index]);
    CCTK_REAL_VEC dgt311L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt311[index]);
    CCTK_REAL_VEC dgt312L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt312[index]);
    CCTK_REAL_VEC dgt313L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt313[index]);
    CCTK_REAL_VEC dgt322L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt322[index]);
    CCTK_REAL_VEC dgt323L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt323[index]);
    CCTK_REAL_VEC dgt333L CCTK_ATTRIBUTE_UNUSED = vec_load(dgt333[index]);
    CCTK_REAL_VEC dphi1L CCTK_ATTRIBUTE_UNUSED = vec_load(dphi1[index]);
    CCTK_REAL_VEC dphi2L CCTK_ATTRIBUTE_UNUSED = vec_load(dphi2[index]);
    CCTK_REAL_VEC dphi3L CCTK_ATTRIBUTE_UNUSED = vec_load(dphi3[index]);
    CCTK_REAL_VEC gt11L CCTK_ATTRIBUTE_UNUSED = vec_load(gt11[index]);
    CCTK_REAL_VEC gt12L CCTK_ATTRIBUTE_UNUSED = vec_load(gt12[index]);
    CCTK_REAL_VEC gt13L CCTK_ATTRIBUTE_UNUSED = vec_load(gt13[index]);
    CCTK_REAL_VEC gt22L CCTK_ATTRIBUTE_UNUSED = vec_load(gt22[index]);
    CCTK_REAL_VEC gt23L CCTK_ATTRIBUTE_UNUSED = vec_load(gt23[index]);
    CCTK_REAL_VEC gt33L CCTK_ATTRIBUTE_UNUSED = vec_load(gt33[index]);
    CCTK_REAL_VEC trKL CCTK_ATTRIBUTE_UNUSED = vec_load(trK[index]);
    CCTK_REAL_VEC Xt1L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt1[index]);
    CCTK_REAL_VEC Xt2L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt2[index]);
    CCTK_REAL_VEC Xt3L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt3[index]);
    
    CCTK_REAL_VEC eTtxL, eTtyL, eTtzL, eTxxL, eTxyL, eTxzL, eTyyL, eTyzL, eTzzL CCTK_ATTRIBUTE_UNUSED ;
    
    if (assume_stress_energy_state>=0 ? assume_stress_energy_state : *stress_energy_state)
    {
      eTtxL = vec_load(eTtx[index]);
      eTtyL = vec_load(eTty[index]);
      eTtzL = vec_load(eTtz[index]);
      eTxxL = vec_load(eTxx[index]);
      eTxyL = vec_load(eTxy[index]);
      eTxzL = vec_load(eTxz[index]);
      eTyyL = vec_load(eTyy[index]);
      eTyzL = vec_load(eTyz[index]);
      eTzzL = vec_load(eTzz[index]);
    }
    else
    {
      eTtxL = ToReal(0.);
      eTtyL = ToReal(0.);
      eTtzL = ToReal(0.);
      eTxxL = ToReal(0.);
      eTxyL = ToReal(0.);
      eTxzL = ToReal(0.);
      eTyyL = ToReal(0.);
      eTyzL = ToReal(0.);
      eTzzL = ToReal(0.);
    }
    
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
    CCTK_REAL_VEC PDstandardNth1At11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2At11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3At11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1At12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2At12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3At12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1At13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2At13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3At13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1At22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2At22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3At22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1At23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2At23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3At23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1At33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2At33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3At33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dbeta11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dbeta11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dbeta11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dbeta12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dbeta12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dbeta12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dbeta13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dbeta13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dbeta13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dbeta21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dbeta21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dbeta21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dbeta22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dbeta22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dbeta22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dbeta23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dbeta23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dbeta23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dbeta31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dbeta31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dbeta31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dbeta32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dbeta32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dbeta32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dbeta33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dbeta33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dbeta33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3Xt3 CCTK_ATTRIBUTE_UNUSED;
    
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
        PDstandardNth1dbeta11 = PDstandardNthfdOrder21(&dbeta11[index]);
        PDstandardNth2dbeta11 = PDstandardNthfdOrder22(&dbeta11[index]);
        PDstandardNth3dbeta11 = PDstandardNthfdOrder23(&dbeta11[index]);
        PDstandardNth1dbeta12 = PDstandardNthfdOrder21(&dbeta12[index]);
        PDstandardNth2dbeta12 = PDstandardNthfdOrder22(&dbeta12[index]);
        PDstandardNth3dbeta12 = PDstandardNthfdOrder23(&dbeta12[index]);
        PDstandardNth1dbeta13 = PDstandardNthfdOrder21(&dbeta13[index]);
        PDstandardNth2dbeta13 = PDstandardNthfdOrder22(&dbeta13[index]);
        PDstandardNth3dbeta13 = PDstandardNthfdOrder23(&dbeta13[index]);
        PDstandardNth1dbeta21 = PDstandardNthfdOrder21(&dbeta21[index]);
        PDstandardNth2dbeta21 = PDstandardNthfdOrder22(&dbeta21[index]);
        PDstandardNth3dbeta21 = PDstandardNthfdOrder23(&dbeta21[index]);
        PDstandardNth1dbeta22 = PDstandardNthfdOrder21(&dbeta22[index]);
        PDstandardNth2dbeta22 = PDstandardNthfdOrder22(&dbeta22[index]);
        PDstandardNth3dbeta22 = PDstandardNthfdOrder23(&dbeta22[index]);
        PDstandardNth1dbeta23 = PDstandardNthfdOrder21(&dbeta23[index]);
        PDstandardNth2dbeta23 = PDstandardNthfdOrder22(&dbeta23[index]);
        PDstandardNth3dbeta23 = PDstandardNthfdOrder23(&dbeta23[index]);
        PDstandardNth1dbeta31 = PDstandardNthfdOrder21(&dbeta31[index]);
        PDstandardNth2dbeta31 = PDstandardNthfdOrder22(&dbeta31[index]);
        PDstandardNth3dbeta31 = PDstandardNthfdOrder23(&dbeta31[index]);
        PDstandardNth1dbeta32 = PDstandardNthfdOrder21(&dbeta32[index]);
        PDstandardNth2dbeta32 = PDstandardNthfdOrder22(&dbeta32[index]);
        PDstandardNth3dbeta32 = PDstandardNthfdOrder23(&dbeta32[index]);
        PDstandardNth1dbeta33 = PDstandardNthfdOrder21(&dbeta33[index]);
        PDstandardNth2dbeta33 = PDstandardNthfdOrder22(&dbeta33[index]);
        PDstandardNth3dbeta33 = PDstandardNthfdOrder23(&dbeta33[index]);
        PDupwindNthAnti1dgt111 = PDupwindNthAntifdOrder21(&dgt111[index]);
        PDupwindNthSymm1dgt111 = PDupwindNthSymmfdOrder21(&dgt111[index]);
        PDupwindNthAnti2dgt111 = PDupwindNthAntifdOrder22(&dgt111[index]);
        PDupwindNthSymm2dgt111 = PDupwindNthSymmfdOrder22(&dgt111[index]);
        PDupwindNthAnti3dgt111 = PDupwindNthAntifdOrder23(&dgt111[index]);
        PDupwindNthSymm3dgt111 = PDupwindNthSymmfdOrder23(&dgt111[index]);
        PDupwindNthAnti1dgt112 = PDupwindNthAntifdOrder21(&dgt112[index]);
        PDupwindNthSymm1dgt112 = PDupwindNthSymmfdOrder21(&dgt112[index]);
        PDupwindNthAnti2dgt112 = PDupwindNthAntifdOrder22(&dgt112[index]);
        PDupwindNthSymm2dgt112 = PDupwindNthSymmfdOrder22(&dgt112[index]);
        PDupwindNthAnti3dgt112 = PDupwindNthAntifdOrder23(&dgt112[index]);
        PDupwindNthSymm3dgt112 = PDupwindNthSymmfdOrder23(&dgt112[index]);
        PDupwindNthAnti1dgt113 = PDupwindNthAntifdOrder21(&dgt113[index]);
        PDupwindNthSymm1dgt113 = PDupwindNthSymmfdOrder21(&dgt113[index]);
        PDupwindNthAnti2dgt113 = PDupwindNthAntifdOrder22(&dgt113[index]);
        PDupwindNthSymm2dgt113 = PDupwindNthSymmfdOrder22(&dgt113[index]);
        PDupwindNthAnti3dgt113 = PDupwindNthAntifdOrder23(&dgt113[index]);
        PDupwindNthSymm3dgt113 = PDupwindNthSymmfdOrder23(&dgt113[index]);
        PDupwindNthAnti1dgt122 = PDupwindNthAntifdOrder21(&dgt122[index]);
        PDupwindNthSymm1dgt122 = PDupwindNthSymmfdOrder21(&dgt122[index]);
        PDupwindNthAnti2dgt122 = PDupwindNthAntifdOrder22(&dgt122[index]);
        PDupwindNthSymm2dgt122 = PDupwindNthSymmfdOrder22(&dgt122[index]);
        PDupwindNthAnti3dgt122 = PDupwindNthAntifdOrder23(&dgt122[index]);
        PDupwindNthSymm3dgt122 = PDupwindNthSymmfdOrder23(&dgt122[index]);
        PDupwindNthAnti1dgt123 = PDupwindNthAntifdOrder21(&dgt123[index]);
        PDupwindNthSymm1dgt123 = PDupwindNthSymmfdOrder21(&dgt123[index]);
        PDupwindNthAnti2dgt123 = PDupwindNthAntifdOrder22(&dgt123[index]);
        PDupwindNthSymm2dgt123 = PDupwindNthSymmfdOrder22(&dgt123[index]);
        PDupwindNthAnti3dgt123 = PDupwindNthAntifdOrder23(&dgt123[index]);
        PDupwindNthSymm3dgt123 = PDupwindNthSymmfdOrder23(&dgt123[index]);
        PDupwindNthAnti1dgt133 = PDupwindNthAntifdOrder21(&dgt133[index]);
        PDupwindNthSymm1dgt133 = PDupwindNthSymmfdOrder21(&dgt133[index]);
        PDupwindNthAnti2dgt133 = PDupwindNthAntifdOrder22(&dgt133[index]);
        PDupwindNthSymm2dgt133 = PDupwindNthSymmfdOrder22(&dgt133[index]);
        PDupwindNthAnti3dgt133 = PDupwindNthAntifdOrder23(&dgt133[index]);
        PDupwindNthSymm3dgt133 = PDupwindNthSymmfdOrder23(&dgt133[index]);
        PDupwindNthAnti1dgt211 = PDupwindNthAntifdOrder21(&dgt211[index]);
        PDupwindNthSymm1dgt211 = PDupwindNthSymmfdOrder21(&dgt211[index]);
        PDupwindNthAnti2dgt211 = PDupwindNthAntifdOrder22(&dgt211[index]);
        PDupwindNthSymm2dgt211 = PDupwindNthSymmfdOrder22(&dgt211[index]);
        PDupwindNthAnti3dgt211 = PDupwindNthAntifdOrder23(&dgt211[index]);
        PDupwindNthSymm3dgt211 = PDupwindNthSymmfdOrder23(&dgt211[index]);
        PDupwindNthAnti1dgt212 = PDupwindNthAntifdOrder21(&dgt212[index]);
        PDupwindNthSymm1dgt212 = PDupwindNthSymmfdOrder21(&dgt212[index]);
        PDupwindNthAnti2dgt212 = PDupwindNthAntifdOrder22(&dgt212[index]);
        PDupwindNthSymm2dgt212 = PDupwindNthSymmfdOrder22(&dgt212[index]);
        PDupwindNthAnti3dgt212 = PDupwindNthAntifdOrder23(&dgt212[index]);
        PDupwindNthSymm3dgt212 = PDupwindNthSymmfdOrder23(&dgt212[index]);
        PDupwindNthAnti1dgt213 = PDupwindNthAntifdOrder21(&dgt213[index]);
        PDupwindNthSymm1dgt213 = PDupwindNthSymmfdOrder21(&dgt213[index]);
        PDupwindNthAnti2dgt213 = PDupwindNthAntifdOrder22(&dgt213[index]);
        PDupwindNthSymm2dgt213 = PDupwindNthSymmfdOrder22(&dgt213[index]);
        PDupwindNthAnti3dgt213 = PDupwindNthAntifdOrder23(&dgt213[index]);
        PDupwindNthSymm3dgt213 = PDupwindNthSymmfdOrder23(&dgt213[index]);
        PDupwindNthAnti1dgt222 = PDupwindNthAntifdOrder21(&dgt222[index]);
        PDupwindNthSymm1dgt222 = PDupwindNthSymmfdOrder21(&dgt222[index]);
        PDupwindNthAnti2dgt222 = PDupwindNthAntifdOrder22(&dgt222[index]);
        PDupwindNthSymm2dgt222 = PDupwindNthSymmfdOrder22(&dgt222[index]);
        PDupwindNthAnti3dgt222 = PDupwindNthAntifdOrder23(&dgt222[index]);
        PDupwindNthSymm3dgt222 = PDupwindNthSymmfdOrder23(&dgt222[index]);
        PDupwindNthAnti1dgt223 = PDupwindNthAntifdOrder21(&dgt223[index]);
        PDupwindNthSymm1dgt223 = PDupwindNthSymmfdOrder21(&dgt223[index]);
        PDupwindNthAnti2dgt223 = PDupwindNthAntifdOrder22(&dgt223[index]);
        PDupwindNthSymm2dgt223 = PDupwindNthSymmfdOrder22(&dgt223[index]);
        PDupwindNthAnti3dgt223 = PDupwindNthAntifdOrder23(&dgt223[index]);
        PDupwindNthSymm3dgt223 = PDupwindNthSymmfdOrder23(&dgt223[index]);
        PDupwindNthAnti1dgt233 = PDupwindNthAntifdOrder21(&dgt233[index]);
        PDupwindNthSymm1dgt233 = PDupwindNthSymmfdOrder21(&dgt233[index]);
        PDupwindNthAnti2dgt233 = PDupwindNthAntifdOrder22(&dgt233[index]);
        PDupwindNthSymm2dgt233 = PDupwindNthSymmfdOrder22(&dgt233[index]);
        PDupwindNthAnti3dgt233 = PDupwindNthAntifdOrder23(&dgt233[index]);
        PDupwindNthSymm3dgt233 = PDupwindNthSymmfdOrder23(&dgt233[index]);
        PDupwindNthAnti1dgt311 = PDupwindNthAntifdOrder21(&dgt311[index]);
        PDupwindNthSymm1dgt311 = PDupwindNthSymmfdOrder21(&dgt311[index]);
        PDupwindNthAnti2dgt311 = PDupwindNthAntifdOrder22(&dgt311[index]);
        PDupwindNthSymm2dgt311 = PDupwindNthSymmfdOrder22(&dgt311[index]);
        PDupwindNthAnti3dgt311 = PDupwindNthAntifdOrder23(&dgt311[index]);
        PDupwindNthSymm3dgt311 = PDupwindNthSymmfdOrder23(&dgt311[index]);
        PDupwindNthAnti1dgt312 = PDupwindNthAntifdOrder21(&dgt312[index]);
        PDupwindNthSymm1dgt312 = PDupwindNthSymmfdOrder21(&dgt312[index]);
        PDupwindNthAnti2dgt312 = PDupwindNthAntifdOrder22(&dgt312[index]);
        PDupwindNthSymm2dgt312 = PDupwindNthSymmfdOrder22(&dgt312[index]);
        PDupwindNthAnti3dgt312 = PDupwindNthAntifdOrder23(&dgt312[index]);
        PDupwindNthSymm3dgt312 = PDupwindNthSymmfdOrder23(&dgt312[index]);
        PDupwindNthAnti1dgt313 = PDupwindNthAntifdOrder21(&dgt313[index]);
        PDupwindNthSymm1dgt313 = PDupwindNthSymmfdOrder21(&dgt313[index]);
        PDupwindNthAnti2dgt313 = PDupwindNthAntifdOrder22(&dgt313[index]);
        PDupwindNthSymm2dgt313 = PDupwindNthSymmfdOrder22(&dgt313[index]);
        PDupwindNthAnti3dgt313 = PDupwindNthAntifdOrder23(&dgt313[index]);
        PDupwindNthSymm3dgt313 = PDupwindNthSymmfdOrder23(&dgt313[index]);
        PDupwindNthAnti1dgt322 = PDupwindNthAntifdOrder21(&dgt322[index]);
        PDupwindNthSymm1dgt322 = PDupwindNthSymmfdOrder21(&dgt322[index]);
        PDupwindNthAnti2dgt322 = PDupwindNthAntifdOrder22(&dgt322[index]);
        PDupwindNthSymm2dgt322 = PDupwindNthSymmfdOrder22(&dgt322[index]);
        PDupwindNthAnti3dgt322 = PDupwindNthAntifdOrder23(&dgt322[index]);
        PDupwindNthSymm3dgt322 = PDupwindNthSymmfdOrder23(&dgt322[index]);
        PDupwindNthAnti1dgt323 = PDupwindNthAntifdOrder21(&dgt323[index]);
        PDupwindNthSymm1dgt323 = PDupwindNthSymmfdOrder21(&dgt323[index]);
        PDupwindNthAnti2dgt323 = PDupwindNthAntifdOrder22(&dgt323[index]);
        PDupwindNthSymm2dgt323 = PDupwindNthSymmfdOrder22(&dgt323[index]);
        PDupwindNthAnti3dgt323 = PDupwindNthAntifdOrder23(&dgt323[index]);
        PDupwindNthSymm3dgt323 = PDupwindNthSymmfdOrder23(&dgt323[index]);
        PDupwindNthAnti1dgt333 = PDupwindNthAntifdOrder21(&dgt333[index]);
        PDupwindNthSymm1dgt333 = PDupwindNthSymmfdOrder21(&dgt333[index]);
        PDupwindNthAnti2dgt333 = PDupwindNthAntifdOrder22(&dgt333[index]);
        PDupwindNthSymm2dgt333 = PDupwindNthSymmfdOrder22(&dgt333[index]);
        PDupwindNthAnti3dgt333 = PDupwindNthAntifdOrder23(&dgt333[index]);
        PDupwindNthSymm3dgt333 = PDupwindNthSymmfdOrder23(&dgt333[index]);
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
        PDstandardNth1trK = PDstandardNthfdOrder21(&trK[index]);
        PDstandardNth2trK = PDstandardNthfdOrder22(&trK[index]);
        PDstandardNth3trK = PDstandardNthfdOrder23(&trK[index]);
        PDupwindNthAnti1Xt1 = PDupwindNthAntifdOrder21(&Xt1[index]);
        PDupwindNthSymm1Xt1 = PDupwindNthSymmfdOrder21(&Xt1[index]);
        PDupwindNthAnti2Xt1 = PDupwindNthAntifdOrder22(&Xt1[index]);
        PDupwindNthSymm2Xt1 = PDupwindNthSymmfdOrder22(&Xt1[index]);
        PDupwindNthAnti3Xt1 = PDupwindNthAntifdOrder23(&Xt1[index]);
        PDupwindNthSymm3Xt1 = PDupwindNthSymmfdOrder23(&Xt1[index]);
        PDupwindNthAnti1Xt2 = PDupwindNthAntifdOrder21(&Xt2[index]);
        PDupwindNthSymm1Xt2 = PDupwindNthSymmfdOrder21(&Xt2[index]);
        PDupwindNthAnti2Xt2 = PDupwindNthAntifdOrder22(&Xt2[index]);
        PDupwindNthSymm2Xt2 = PDupwindNthSymmfdOrder22(&Xt2[index]);
        PDupwindNthAnti3Xt2 = PDupwindNthAntifdOrder23(&Xt2[index]);
        PDupwindNthSymm3Xt2 = PDupwindNthSymmfdOrder23(&Xt2[index]);
        PDupwindNthAnti1Xt3 = PDupwindNthAntifdOrder21(&Xt3[index]);
        PDupwindNthSymm1Xt3 = PDupwindNthSymmfdOrder21(&Xt3[index]);
        PDupwindNthAnti2Xt3 = PDupwindNthAntifdOrder22(&Xt3[index]);
        PDupwindNthSymm2Xt3 = PDupwindNthSymmfdOrder22(&Xt3[index]);
        PDupwindNthAnti3Xt3 = PDupwindNthAntifdOrder23(&Xt3[index]);
        PDupwindNthSymm3Xt3 = PDupwindNthSymmfdOrder23(&Xt3[index]);
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
        PDstandardNth1dbeta11 = PDstandardNthfdOrder41(&dbeta11[index]);
        PDstandardNth2dbeta11 = PDstandardNthfdOrder42(&dbeta11[index]);
        PDstandardNth3dbeta11 = PDstandardNthfdOrder43(&dbeta11[index]);
        PDstandardNth1dbeta12 = PDstandardNthfdOrder41(&dbeta12[index]);
        PDstandardNth2dbeta12 = PDstandardNthfdOrder42(&dbeta12[index]);
        PDstandardNth3dbeta12 = PDstandardNthfdOrder43(&dbeta12[index]);
        PDstandardNth1dbeta13 = PDstandardNthfdOrder41(&dbeta13[index]);
        PDstandardNth2dbeta13 = PDstandardNthfdOrder42(&dbeta13[index]);
        PDstandardNth3dbeta13 = PDstandardNthfdOrder43(&dbeta13[index]);
        PDstandardNth1dbeta21 = PDstandardNthfdOrder41(&dbeta21[index]);
        PDstandardNth2dbeta21 = PDstandardNthfdOrder42(&dbeta21[index]);
        PDstandardNth3dbeta21 = PDstandardNthfdOrder43(&dbeta21[index]);
        PDstandardNth1dbeta22 = PDstandardNthfdOrder41(&dbeta22[index]);
        PDstandardNth2dbeta22 = PDstandardNthfdOrder42(&dbeta22[index]);
        PDstandardNth3dbeta22 = PDstandardNthfdOrder43(&dbeta22[index]);
        PDstandardNth1dbeta23 = PDstandardNthfdOrder41(&dbeta23[index]);
        PDstandardNth2dbeta23 = PDstandardNthfdOrder42(&dbeta23[index]);
        PDstandardNth3dbeta23 = PDstandardNthfdOrder43(&dbeta23[index]);
        PDstandardNth1dbeta31 = PDstandardNthfdOrder41(&dbeta31[index]);
        PDstandardNth2dbeta31 = PDstandardNthfdOrder42(&dbeta31[index]);
        PDstandardNth3dbeta31 = PDstandardNthfdOrder43(&dbeta31[index]);
        PDstandardNth1dbeta32 = PDstandardNthfdOrder41(&dbeta32[index]);
        PDstandardNth2dbeta32 = PDstandardNthfdOrder42(&dbeta32[index]);
        PDstandardNth3dbeta32 = PDstandardNthfdOrder43(&dbeta32[index]);
        PDstandardNth1dbeta33 = PDstandardNthfdOrder41(&dbeta33[index]);
        PDstandardNth2dbeta33 = PDstandardNthfdOrder42(&dbeta33[index]);
        PDstandardNth3dbeta33 = PDstandardNthfdOrder43(&dbeta33[index]);
        PDupwindNthAnti1dgt111 = PDupwindNthAntifdOrder41(&dgt111[index]);
        PDupwindNthSymm1dgt111 = PDupwindNthSymmfdOrder41(&dgt111[index]);
        PDupwindNthAnti2dgt111 = PDupwindNthAntifdOrder42(&dgt111[index]);
        PDupwindNthSymm2dgt111 = PDupwindNthSymmfdOrder42(&dgt111[index]);
        PDupwindNthAnti3dgt111 = PDupwindNthAntifdOrder43(&dgt111[index]);
        PDupwindNthSymm3dgt111 = PDupwindNthSymmfdOrder43(&dgt111[index]);
        PDupwindNthAnti1dgt112 = PDupwindNthAntifdOrder41(&dgt112[index]);
        PDupwindNthSymm1dgt112 = PDupwindNthSymmfdOrder41(&dgt112[index]);
        PDupwindNthAnti2dgt112 = PDupwindNthAntifdOrder42(&dgt112[index]);
        PDupwindNthSymm2dgt112 = PDupwindNthSymmfdOrder42(&dgt112[index]);
        PDupwindNthAnti3dgt112 = PDupwindNthAntifdOrder43(&dgt112[index]);
        PDupwindNthSymm3dgt112 = PDupwindNthSymmfdOrder43(&dgt112[index]);
        PDupwindNthAnti1dgt113 = PDupwindNthAntifdOrder41(&dgt113[index]);
        PDupwindNthSymm1dgt113 = PDupwindNthSymmfdOrder41(&dgt113[index]);
        PDupwindNthAnti2dgt113 = PDupwindNthAntifdOrder42(&dgt113[index]);
        PDupwindNthSymm2dgt113 = PDupwindNthSymmfdOrder42(&dgt113[index]);
        PDupwindNthAnti3dgt113 = PDupwindNthAntifdOrder43(&dgt113[index]);
        PDupwindNthSymm3dgt113 = PDupwindNthSymmfdOrder43(&dgt113[index]);
        PDupwindNthAnti1dgt122 = PDupwindNthAntifdOrder41(&dgt122[index]);
        PDupwindNthSymm1dgt122 = PDupwindNthSymmfdOrder41(&dgt122[index]);
        PDupwindNthAnti2dgt122 = PDupwindNthAntifdOrder42(&dgt122[index]);
        PDupwindNthSymm2dgt122 = PDupwindNthSymmfdOrder42(&dgt122[index]);
        PDupwindNthAnti3dgt122 = PDupwindNthAntifdOrder43(&dgt122[index]);
        PDupwindNthSymm3dgt122 = PDupwindNthSymmfdOrder43(&dgt122[index]);
        PDupwindNthAnti1dgt123 = PDupwindNthAntifdOrder41(&dgt123[index]);
        PDupwindNthSymm1dgt123 = PDupwindNthSymmfdOrder41(&dgt123[index]);
        PDupwindNthAnti2dgt123 = PDupwindNthAntifdOrder42(&dgt123[index]);
        PDupwindNthSymm2dgt123 = PDupwindNthSymmfdOrder42(&dgt123[index]);
        PDupwindNthAnti3dgt123 = PDupwindNthAntifdOrder43(&dgt123[index]);
        PDupwindNthSymm3dgt123 = PDupwindNthSymmfdOrder43(&dgt123[index]);
        PDupwindNthAnti1dgt133 = PDupwindNthAntifdOrder41(&dgt133[index]);
        PDupwindNthSymm1dgt133 = PDupwindNthSymmfdOrder41(&dgt133[index]);
        PDupwindNthAnti2dgt133 = PDupwindNthAntifdOrder42(&dgt133[index]);
        PDupwindNthSymm2dgt133 = PDupwindNthSymmfdOrder42(&dgt133[index]);
        PDupwindNthAnti3dgt133 = PDupwindNthAntifdOrder43(&dgt133[index]);
        PDupwindNthSymm3dgt133 = PDupwindNthSymmfdOrder43(&dgt133[index]);
        PDupwindNthAnti1dgt211 = PDupwindNthAntifdOrder41(&dgt211[index]);
        PDupwindNthSymm1dgt211 = PDupwindNthSymmfdOrder41(&dgt211[index]);
        PDupwindNthAnti2dgt211 = PDupwindNthAntifdOrder42(&dgt211[index]);
        PDupwindNthSymm2dgt211 = PDupwindNthSymmfdOrder42(&dgt211[index]);
        PDupwindNthAnti3dgt211 = PDupwindNthAntifdOrder43(&dgt211[index]);
        PDupwindNthSymm3dgt211 = PDupwindNthSymmfdOrder43(&dgt211[index]);
        PDupwindNthAnti1dgt212 = PDupwindNthAntifdOrder41(&dgt212[index]);
        PDupwindNthSymm1dgt212 = PDupwindNthSymmfdOrder41(&dgt212[index]);
        PDupwindNthAnti2dgt212 = PDupwindNthAntifdOrder42(&dgt212[index]);
        PDupwindNthSymm2dgt212 = PDupwindNthSymmfdOrder42(&dgt212[index]);
        PDupwindNthAnti3dgt212 = PDupwindNthAntifdOrder43(&dgt212[index]);
        PDupwindNthSymm3dgt212 = PDupwindNthSymmfdOrder43(&dgt212[index]);
        PDupwindNthAnti1dgt213 = PDupwindNthAntifdOrder41(&dgt213[index]);
        PDupwindNthSymm1dgt213 = PDupwindNthSymmfdOrder41(&dgt213[index]);
        PDupwindNthAnti2dgt213 = PDupwindNthAntifdOrder42(&dgt213[index]);
        PDupwindNthSymm2dgt213 = PDupwindNthSymmfdOrder42(&dgt213[index]);
        PDupwindNthAnti3dgt213 = PDupwindNthAntifdOrder43(&dgt213[index]);
        PDupwindNthSymm3dgt213 = PDupwindNthSymmfdOrder43(&dgt213[index]);
        PDupwindNthAnti1dgt222 = PDupwindNthAntifdOrder41(&dgt222[index]);
        PDupwindNthSymm1dgt222 = PDupwindNthSymmfdOrder41(&dgt222[index]);
        PDupwindNthAnti2dgt222 = PDupwindNthAntifdOrder42(&dgt222[index]);
        PDupwindNthSymm2dgt222 = PDupwindNthSymmfdOrder42(&dgt222[index]);
        PDupwindNthAnti3dgt222 = PDupwindNthAntifdOrder43(&dgt222[index]);
        PDupwindNthSymm3dgt222 = PDupwindNthSymmfdOrder43(&dgt222[index]);
        PDupwindNthAnti1dgt223 = PDupwindNthAntifdOrder41(&dgt223[index]);
        PDupwindNthSymm1dgt223 = PDupwindNthSymmfdOrder41(&dgt223[index]);
        PDupwindNthAnti2dgt223 = PDupwindNthAntifdOrder42(&dgt223[index]);
        PDupwindNthSymm2dgt223 = PDupwindNthSymmfdOrder42(&dgt223[index]);
        PDupwindNthAnti3dgt223 = PDupwindNthAntifdOrder43(&dgt223[index]);
        PDupwindNthSymm3dgt223 = PDupwindNthSymmfdOrder43(&dgt223[index]);
        PDupwindNthAnti1dgt233 = PDupwindNthAntifdOrder41(&dgt233[index]);
        PDupwindNthSymm1dgt233 = PDupwindNthSymmfdOrder41(&dgt233[index]);
        PDupwindNthAnti2dgt233 = PDupwindNthAntifdOrder42(&dgt233[index]);
        PDupwindNthSymm2dgt233 = PDupwindNthSymmfdOrder42(&dgt233[index]);
        PDupwindNthAnti3dgt233 = PDupwindNthAntifdOrder43(&dgt233[index]);
        PDupwindNthSymm3dgt233 = PDupwindNthSymmfdOrder43(&dgt233[index]);
        PDupwindNthAnti1dgt311 = PDupwindNthAntifdOrder41(&dgt311[index]);
        PDupwindNthSymm1dgt311 = PDupwindNthSymmfdOrder41(&dgt311[index]);
        PDupwindNthAnti2dgt311 = PDupwindNthAntifdOrder42(&dgt311[index]);
        PDupwindNthSymm2dgt311 = PDupwindNthSymmfdOrder42(&dgt311[index]);
        PDupwindNthAnti3dgt311 = PDupwindNthAntifdOrder43(&dgt311[index]);
        PDupwindNthSymm3dgt311 = PDupwindNthSymmfdOrder43(&dgt311[index]);
        PDupwindNthAnti1dgt312 = PDupwindNthAntifdOrder41(&dgt312[index]);
        PDupwindNthSymm1dgt312 = PDupwindNthSymmfdOrder41(&dgt312[index]);
        PDupwindNthAnti2dgt312 = PDupwindNthAntifdOrder42(&dgt312[index]);
        PDupwindNthSymm2dgt312 = PDupwindNthSymmfdOrder42(&dgt312[index]);
        PDupwindNthAnti3dgt312 = PDupwindNthAntifdOrder43(&dgt312[index]);
        PDupwindNthSymm3dgt312 = PDupwindNthSymmfdOrder43(&dgt312[index]);
        PDupwindNthAnti1dgt313 = PDupwindNthAntifdOrder41(&dgt313[index]);
        PDupwindNthSymm1dgt313 = PDupwindNthSymmfdOrder41(&dgt313[index]);
        PDupwindNthAnti2dgt313 = PDupwindNthAntifdOrder42(&dgt313[index]);
        PDupwindNthSymm2dgt313 = PDupwindNthSymmfdOrder42(&dgt313[index]);
        PDupwindNthAnti3dgt313 = PDupwindNthAntifdOrder43(&dgt313[index]);
        PDupwindNthSymm3dgt313 = PDupwindNthSymmfdOrder43(&dgt313[index]);
        PDupwindNthAnti1dgt322 = PDupwindNthAntifdOrder41(&dgt322[index]);
        PDupwindNthSymm1dgt322 = PDupwindNthSymmfdOrder41(&dgt322[index]);
        PDupwindNthAnti2dgt322 = PDupwindNthAntifdOrder42(&dgt322[index]);
        PDupwindNthSymm2dgt322 = PDupwindNthSymmfdOrder42(&dgt322[index]);
        PDupwindNthAnti3dgt322 = PDupwindNthAntifdOrder43(&dgt322[index]);
        PDupwindNthSymm3dgt322 = PDupwindNthSymmfdOrder43(&dgt322[index]);
        PDupwindNthAnti1dgt323 = PDupwindNthAntifdOrder41(&dgt323[index]);
        PDupwindNthSymm1dgt323 = PDupwindNthSymmfdOrder41(&dgt323[index]);
        PDupwindNthAnti2dgt323 = PDupwindNthAntifdOrder42(&dgt323[index]);
        PDupwindNthSymm2dgt323 = PDupwindNthSymmfdOrder42(&dgt323[index]);
        PDupwindNthAnti3dgt323 = PDupwindNthAntifdOrder43(&dgt323[index]);
        PDupwindNthSymm3dgt323 = PDupwindNthSymmfdOrder43(&dgt323[index]);
        PDupwindNthAnti1dgt333 = PDupwindNthAntifdOrder41(&dgt333[index]);
        PDupwindNthSymm1dgt333 = PDupwindNthSymmfdOrder41(&dgt333[index]);
        PDupwindNthAnti2dgt333 = PDupwindNthAntifdOrder42(&dgt333[index]);
        PDupwindNthSymm2dgt333 = PDupwindNthSymmfdOrder42(&dgt333[index]);
        PDupwindNthAnti3dgt333 = PDupwindNthAntifdOrder43(&dgt333[index]);
        PDupwindNthSymm3dgt333 = PDupwindNthSymmfdOrder43(&dgt333[index]);
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
        PDstandardNth1trK = PDstandardNthfdOrder41(&trK[index]);
        PDstandardNth2trK = PDstandardNthfdOrder42(&trK[index]);
        PDstandardNth3trK = PDstandardNthfdOrder43(&trK[index]);
        PDupwindNthAnti1Xt1 = PDupwindNthAntifdOrder41(&Xt1[index]);
        PDupwindNthSymm1Xt1 = PDupwindNthSymmfdOrder41(&Xt1[index]);
        PDupwindNthAnti2Xt1 = PDupwindNthAntifdOrder42(&Xt1[index]);
        PDupwindNthSymm2Xt1 = PDupwindNthSymmfdOrder42(&Xt1[index]);
        PDupwindNthAnti3Xt1 = PDupwindNthAntifdOrder43(&Xt1[index]);
        PDupwindNthSymm3Xt1 = PDupwindNthSymmfdOrder43(&Xt1[index]);
        PDupwindNthAnti1Xt2 = PDupwindNthAntifdOrder41(&Xt2[index]);
        PDupwindNthSymm1Xt2 = PDupwindNthSymmfdOrder41(&Xt2[index]);
        PDupwindNthAnti2Xt2 = PDupwindNthAntifdOrder42(&Xt2[index]);
        PDupwindNthSymm2Xt2 = PDupwindNthSymmfdOrder42(&Xt2[index]);
        PDupwindNthAnti3Xt2 = PDupwindNthAntifdOrder43(&Xt2[index]);
        PDupwindNthSymm3Xt2 = PDupwindNthSymmfdOrder43(&Xt2[index]);
        PDupwindNthAnti1Xt3 = PDupwindNthAntifdOrder41(&Xt3[index]);
        PDupwindNthSymm1Xt3 = PDupwindNthSymmfdOrder41(&Xt3[index]);
        PDupwindNthAnti2Xt3 = PDupwindNthAntifdOrder42(&Xt3[index]);
        PDupwindNthSymm2Xt3 = PDupwindNthSymmfdOrder42(&Xt3[index]);
        PDupwindNthAnti3Xt3 = PDupwindNthAntifdOrder43(&Xt3[index]);
        PDupwindNthSymm3Xt3 = PDupwindNthSymmfdOrder43(&Xt3[index]);
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
        PDstandardNth1dbeta11 = PDstandardNthfdOrder61(&dbeta11[index]);
        PDstandardNth2dbeta11 = PDstandardNthfdOrder62(&dbeta11[index]);
        PDstandardNth3dbeta11 = PDstandardNthfdOrder63(&dbeta11[index]);
        PDstandardNth1dbeta12 = PDstandardNthfdOrder61(&dbeta12[index]);
        PDstandardNth2dbeta12 = PDstandardNthfdOrder62(&dbeta12[index]);
        PDstandardNth3dbeta12 = PDstandardNthfdOrder63(&dbeta12[index]);
        PDstandardNth1dbeta13 = PDstandardNthfdOrder61(&dbeta13[index]);
        PDstandardNth2dbeta13 = PDstandardNthfdOrder62(&dbeta13[index]);
        PDstandardNth3dbeta13 = PDstandardNthfdOrder63(&dbeta13[index]);
        PDstandardNth1dbeta21 = PDstandardNthfdOrder61(&dbeta21[index]);
        PDstandardNth2dbeta21 = PDstandardNthfdOrder62(&dbeta21[index]);
        PDstandardNth3dbeta21 = PDstandardNthfdOrder63(&dbeta21[index]);
        PDstandardNth1dbeta22 = PDstandardNthfdOrder61(&dbeta22[index]);
        PDstandardNth2dbeta22 = PDstandardNthfdOrder62(&dbeta22[index]);
        PDstandardNth3dbeta22 = PDstandardNthfdOrder63(&dbeta22[index]);
        PDstandardNth1dbeta23 = PDstandardNthfdOrder61(&dbeta23[index]);
        PDstandardNth2dbeta23 = PDstandardNthfdOrder62(&dbeta23[index]);
        PDstandardNth3dbeta23 = PDstandardNthfdOrder63(&dbeta23[index]);
        PDstandardNth1dbeta31 = PDstandardNthfdOrder61(&dbeta31[index]);
        PDstandardNth2dbeta31 = PDstandardNthfdOrder62(&dbeta31[index]);
        PDstandardNth3dbeta31 = PDstandardNthfdOrder63(&dbeta31[index]);
        PDstandardNth1dbeta32 = PDstandardNthfdOrder61(&dbeta32[index]);
        PDstandardNth2dbeta32 = PDstandardNthfdOrder62(&dbeta32[index]);
        PDstandardNth3dbeta32 = PDstandardNthfdOrder63(&dbeta32[index]);
        PDstandardNth1dbeta33 = PDstandardNthfdOrder61(&dbeta33[index]);
        PDstandardNth2dbeta33 = PDstandardNthfdOrder62(&dbeta33[index]);
        PDstandardNth3dbeta33 = PDstandardNthfdOrder63(&dbeta33[index]);
        PDupwindNthAnti1dgt111 = PDupwindNthAntifdOrder61(&dgt111[index]);
        PDupwindNthSymm1dgt111 = PDupwindNthSymmfdOrder61(&dgt111[index]);
        PDupwindNthAnti2dgt111 = PDupwindNthAntifdOrder62(&dgt111[index]);
        PDupwindNthSymm2dgt111 = PDupwindNthSymmfdOrder62(&dgt111[index]);
        PDupwindNthAnti3dgt111 = PDupwindNthAntifdOrder63(&dgt111[index]);
        PDupwindNthSymm3dgt111 = PDupwindNthSymmfdOrder63(&dgt111[index]);
        PDupwindNthAnti1dgt112 = PDupwindNthAntifdOrder61(&dgt112[index]);
        PDupwindNthSymm1dgt112 = PDupwindNthSymmfdOrder61(&dgt112[index]);
        PDupwindNthAnti2dgt112 = PDupwindNthAntifdOrder62(&dgt112[index]);
        PDupwindNthSymm2dgt112 = PDupwindNthSymmfdOrder62(&dgt112[index]);
        PDupwindNthAnti3dgt112 = PDupwindNthAntifdOrder63(&dgt112[index]);
        PDupwindNthSymm3dgt112 = PDupwindNthSymmfdOrder63(&dgt112[index]);
        PDupwindNthAnti1dgt113 = PDupwindNthAntifdOrder61(&dgt113[index]);
        PDupwindNthSymm1dgt113 = PDupwindNthSymmfdOrder61(&dgt113[index]);
        PDupwindNthAnti2dgt113 = PDupwindNthAntifdOrder62(&dgt113[index]);
        PDupwindNthSymm2dgt113 = PDupwindNthSymmfdOrder62(&dgt113[index]);
        PDupwindNthAnti3dgt113 = PDupwindNthAntifdOrder63(&dgt113[index]);
        PDupwindNthSymm3dgt113 = PDupwindNthSymmfdOrder63(&dgt113[index]);
        PDupwindNthAnti1dgt122 = PDupwindNthAntifdOrder61(&dgt122[index]);
        PDupwindNthSymm1dgt122 = PDupwindNthSymmfdOrder61(&dgt122[index]);
        PDupwindNthAnti2dgt122 = PDupwindNthAntifdOrder62(&dgt122[index]);
        PDupwindNthSymm2dgt122 = PDupwindNthSymmfdOrder62(&dgt122[index]);
        PDupwindNthAnti3dgt122 = PDupwindNthAntifdOrder63(&dgt122[index]);
        PDupwindNthSymm3dgt122 = PDupwindNthSymmfdOrder63(&dgt122[index]);
        PDupwindNthAnti1dgt123 = PDupwindNthAntifdOrder61(&dgt123[index]);
        PDupwindNthSymm1dgt123 = PDupwindNthSymmfdOrder61(&dgt123[index]);
        PDupwindNthAnti2dgt123 = PDupwindNthAntifdOrder62(&dgt123[index]);
        PDupwindNthSymm2dgt123 = PDupwindNthSymmfdOrder62(&dgt123[index]);
        PDupwindNthAnti3dgt123 = PDupwindNthAntifdOrder63(&dgt123[index]);
        PDupwindNthSymm3dgt123 = PDupwindNthSymmfdOrder63(&dgt123[index]);
        PDupwindNthAnti1dgt133 = PDupwindNthAntifdOrder61(&dgt133[index]);
        PDupwindNthSymm1dgt133 = PDupwindNthSymmfdOrder61(&dgt133[index]);
        PDupwindNthAnti2dgt133 = PDupwindNthAntifdOrder62(&dgt133[index]);
        PDupwindNthSymm2dgt133 = PDupwindNthSymmfdOrder62(&dgt133[index]);
        PDupwindNthAnti3dgt133 = PDupwindNthAntifdOrder63(&dgt133[index]);
        PDupwindNthSymm3dgt133 = PDupwindNthSymmfdOrder63(&dgt133[index]);
        PDupwindNthAnti1dgt211 = PDupwindNthAntifdOrder61(&dgt211[index]);
        PDupwindNthSymm1dgt211 = PDupwindNthSymmfdOrder61(&dgt211[index]);
        PDupwindNthAnti2dgt211 = PDupwindNthAntifdOrder62(&dgt211[index]);
        PDupwindNthSymm2dgt211 = PDupwindNthSymmfdOrder62(&dgt211[index]);
        PDupwindNthAnti3dgt211 = PDupwindNthAntifdOrder63(&dgt211[index]);
        PDupwindNthSymm3dgt211 = PDupwindNthSymmfdOrder63(&dgt211[index]);
        PDupwindNthAnti1dgt212 = PDupwindNthAntifdOrder61(&dgt212[index]);
        PDupwindNthSymm1dgt212 = PDupwindNthSymmfdOrder61(&dgt212[index]);
        PDupwindNthAnti2dgt212 = PDupwindNthAntifdOrder62(&dgt212[index]);
        PDupwindNthSymm2dgt212 = PDupwindNthSymmfdOrder62(&dgt212[index]);
        PDupwindNthAnti3dgt212 = PDupwindNthAntifdOrder63(&dgt212[index]);
        PDupwindNthSymm3dgt212 = PDupwindNthSymmfdOrder63(&dgt212[index]);
        PDupwindNthAnti1dgt213 = PDupwindNthAntifdOrder61(&dgt213[index]);
        PDupwindNthSymm1dgt213 = PDupwindNthSymmfdOrder61(&dgt213[index]);
        PDupwindNthAnti2dgt213 = PDupwindNthAntifdOrder62(&dgt213[index]);
        PDupwindNthSymm2dgt213 = PDupwindNthSymmfdOrder62(&dgt213[index]);
        PDupwindNthAnti3dgt213 = PDupwindNthAntifdOrder63(&dgt213[index]);
        PDupwindNthSymm3dgt213 = PDupwindNthSymmfdOrder63(&dgt213[index]);
        PDupwindNthAnti1dgt222 = PDupwindNthAntifdOrder61(&dgt222[index]);
        PDupwindNthSymm1dgt222 = PDupwindNthSymmfdOrder61(&dgt222[index]);
        PDupwindNthAnti2dgt222 = PDupwindNthAntifdOrder62(&dgt222[index]);
        PDupwindNthSymm2dgt222 = PDupwindNthSymmfdOrder62(&dgt222[index]);
        PDupwindNthAnti3dgt222 = PDupwindNthAntifdOrder63(&dgt222[index]);
        PDupwindNthSymm3dgt222 = PDupwindNthSymmfdOrder63(&dgt222[index]);
        PDupwindNthAnti1dgt223 = PDupwindNthAntifdOrder61(&dgt223[index]);
        PDupwindNthSymm1dgt223 = PDupwindNthSymmfdOrder61(&dgt223[index]);
        PDupwindNthAnti2dgt223 = PDupwindNthAntifdOrder62(&dgt223[index]);
        PDupwindNthSymm2dgt223 = PDupwindNthSymmfdOrder62(&dgt223[index]);
        PDupwindNthAnti3dgt223 = PDupwindNthAntifdOrder63(&dgt223[index]);
        PDupwindNthSymm3dgt223 = PDupwindNthSymmfdOrder63(&dgt223[index]);
        PDupwindNthAnti1dgt233 = PDupwindNthAntifdOrder61(&dgt233[index]);
        PDupwindNthSymm1dgt233 = PDupwindNthSymmfdOrder61(&dgt233[index]);
        PDupwindNthAnti2dgt233 = PDupwindNthAntifdOrder62(&dgt233[index]);
        PDupwindNthSymm2dgt233 = PDupwindNthSymmfdOrder62(&dgt233[index]);
        PDupwindNthAnti3dgt233 = PDupwindNthAntifdOrder63(&dgt233[index]);
        PDupwindNthSymm3dgt233 = PDupwindNthSymmfdOrder63(&dgt233[index]);
        PDupwindNthAnti1dgt311 = PDupwindNthAntifdOrder61(&dgt311[index]);
        PDupwindNthSymm1dgt311 = PDupwindNthSymmfdOrder61(&dgt311[index]);
        PDupwindNthAnti2dgt311 = PDupwindNthAntifdOrder62(&dgt311[index]);
        PDupwindNthSymm2dgt311 = PDupwindNthSymmfdOrder62(&dgt311[index]);
        PDupwindNthAnti3dgt311 = PDupwindNthAntifdOrder63(&dgt311[index]);
        PDupwindNthSymm3dgt311 = PDupwindNthSymmfdOrder63(&dgt311[index]);
        PDupwindNthAnti1dgt312 = PDupwindNthAntifdOrder61(&dgt312[index]);
        PDupwindNthSymm1dgt312 = PDupwindNthSymmfdOrder61(&dgt312[index]);
        PDupwindNthAnti2dgt312 = PDupwindNthAntifdOrder62(&dgt312[index]);
        PDupwindNthSymm2dgt312 = PDupwindNthSymmfdOrder62(&dgt312[index]);
        PDupwindNthAnti3dgt312 = PDupwindNthAntifdOrder63(&dgt312[index]);
        PDupwindNthSymm3dgt312 = PDupwindNthSymmfdOrder63(&dgt312[index]);
        PDupwindNthAnti1dgt313 = PDupwindNthAntifdOrder61(&dgt313[index]);
        PDupwindNthSymm1dgt313 = PDupwindNthSymmfdOrder61(&dgt313[index]);
        PDupwindNthAnti2dgt313 = PDupwindNthAntifdOrder62(&dgt313[index]);
        PDupwindNthSymm2dgt313 = PDupwindNthSymmfdOrder62(&dgt313[index]);
        PDupwindNthAnti3dgt313 = PDupwindNthAntifdOrder63(&dgt313[index]);
        PDupwindNthSymm3dgt313 = PDupwindNthSymmfdOrder63(&dgt313[index]);
        PDupwindNthAnti1dgt322 = PDupwindNthAntifdOrder61(&dgt322[index]);
        PDupwindNthSymm1dgt322 = PDupwindNthSymmfdOrder61(&dgt322[index]);
        PDupwindNthAnti2dgt322 = PDupwindNthAntifdOrder62(&dgt322[index]);
        PDupwindNthSymm2dgt322 = PDupwindNthSymmfdOrder62(&dgt322[index]);
        PDupwindNthAnti3dgt322 = PDupwindNthAntifdOrder63(&dgt322[index]);
        PDupwindNthSymm3dgt322 = PDupwindNthSymmfdOrder63(&dgt322[index]);
        PDupwindNthAnti1dgt323 = PDupwindNthAntifdOrder61(&dgt323[index]);
        PDupwindNthSymm1dgt323 = PDupwindNthSymmfdOrder61(&dgt323[index]);
        PDupwindNthAnti2dgt323 = PDupwindNthAntifdOrder62(&dgt323[index]);
        PDupwindNthSymm2dgt323 = PDupwindNthSymmfdOrder62(&dgt323[index]);
        PDupwindNthAnti3dgt323 = PDupwindNthAntifdOrder63(&dgt323[index]);
        PDupwindNthSymm3dgt323 = PDupwindNthSymmfdOrder63(&dgt323[index]);
        PDupwindNthAnti1dgt333 = PDupwindNthAntifdOrder61(&dgt333[index]);
        PDupwindNthSymm1dgt333 = PDupwindNthSymmfdOrder61(&dgt333[index]);
        PDupwindNthAnti2dgt333 = PDupwindNthAntifdOrder62(&dgt333[index]);
        PDupwindNthSymm2dgt333 = PDupwindNthSymmfdOrder62(&dgt333[index]);
        PDupwindNthAnti3dgt333 = PDupwindNthAntifdOrder63(&dgt333[index]);
        PDupwindNthSymm3dgt333 = PDupwindNthSymmfdOrder63(&dgt333[index]);
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
        PDstandardNth1trK = PDstandardNthfdOrder61(&trK[index]);
        PDstandardNth2trK = PDstandardNthfdOrder62(&trK[index]);
        PDstandardNth3trK = PDstandardNthfdOrder63(&trK[index]);
        PDupwindNthAnti1Xt1 = PDupwindNthAntifdOrder61(&Xt1[index]);
        PDupwindNthSymm1Xt1 = PDupwindNthSymmfdOrder61(&Xt1[index]);
        PDupwindNthAnti2Xt1 = PDupwindNthAntifdOrder62(&Xt1[index]);
        PDupwindNthSymm2Xt1 = PDupwindNthSymmfdOrder62(&Xt1[index]);
        PDupwindNthAnti3Xt1 = PDupwindNthAntifdOrder63(&Xt1[index]);
        PDupwindNthSymm3Xt1 = PDupwindNthSymmfdOrder63(&Xt1[index]);
        PDupwindNthAnti1Xt2 = PDupwindNthAntifdOrder61(&Xt2[index]);
        PDupwindNthSymm1Xt2 = PDupwindNthSymmfdOrder61(&Xt2[index]);
        PDupwindNthAnti2Xt2 = PDupwindNthAntifdOrder62(&Xt2[index]);
        PDupwindNthSymm2Xt2 = PDupwindNthSymmfdOrder62(&Xt2[index]);
        PDupwindNthAnti3Xt2 = PDupwindNthAntifdOrder63(&Xt2[index]);
        PDupwindNthSymm3Xt2 = PDupwindNthSymmfdOrder63(&Xt2[index]);
        PDupwindNthAnti1Xt3 = PDupwindNthAntifdOrder61(&Xt3[index]);
        PDupwindNthSymm1Xt3 = PDupwindNthSymmfdOrder61(&Xt3[index]);
        PDupwindNthAnti2Xt3 = PDupwindNthAntifdOrder62(&Xt3[index]);
        PDupwindNthSymm2Xt3 = PDupwindNthSymmfdOrder62(&Xt3[index]);
        PDupwindNthAnti3Xt3 = PDupwindNthAntifdOrder63(&Xt3[index]);
        PDupwindNthSymm3Xt3 = PDupwindNthSymmfdOrder63(&Xt3[index]);
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
        PDstandardNth1dbeta11 = PDstandardNthfdOrder81(&dbeta11[index]);
        PDstandardNth2dbeta11 = PDstandardNthfdOrder82(&dbeta11[index]);
        PDstandardNth3dbeta11 = PDstandardNthfdOrder83(&dbeta11[index]);
        PDstandardNth1dbeta12 = PDstandardNthfdOrder81(&dbeta12[index]);
        PDstandardNth2dbeta12 = PDstandardNthfdOrder82(&dbeta12[index]);
        PDstandardNth3dbeta12 = PDstandardNthfdOrder83(&dbeta12[index]);
        PDstandardNth1dbeta13 = PDstandardNthfdOrder81(&dbeta13[index]);
        PDstandardNth2dbeta13 = PDstandardNthfdOrder82(&dbeta13[index]);
        PDstandardNth3dbeta13 = PDstandardNthfdOrder83(&dbeta13[index]);
        PDstandardNth1dbeta21 = PDstandardNthfdOrder81(&dbeta21[index]);
        PDstandardNth2dbeta21 = PDstandardNthfdOrder82(&dbeta21[index]);
        PDstandardNth3dbeta21 = PDstandardNthfdOrder83(&dbeta21[index]);
        PDstandardNth1dbeta22 = PDstandardNthfdOrder81(&dbeta22[index]);
        PDstandardNth2dbeta22 = PDstandardNthfdOrder82(&dbeta22[index]);
        PDstandardNth3dbeta22 = PDstandardNthfdOrder83(&dbeta22[index]);
        PDstandardNth1dbeta23 = PDstandardNthfdOrder81(&dbeta23[index]);
        PDstandardNth2dbeta23 = PDstandardNthfdOrder82(&dbeta23[index]);
        PDstandardNth3dbeta23 = PDstandardNthfdOrder83(&dbeta23[index]);
        PDstandardNth1dbeta31 = PDstandardNthfdOrder81(&dbeta31[index]);
        PDstandardNth2dbeta31 = PDstandardNthfdOrder82(&dbeta31[index]);
        PDstandardNth3dbeta31 = PDstandardNthfdOrder83(&dbeta31[index]);
        PDstandardNth1dbeta32 = PDstandardNthfdOrder81(&dbeta32[index]);
        PDstandardNth2dbeta32 = PDstandardNthfdOrder82(&dbeta32[index]);
        PDstandardNth3dbeta32 = PDstandardNthfdOrder83(&dbeta32[index]);
        PDstandardNth1dbeta33 = PDstandardNthfdOrder81(&dbeta33[index]);
        PDstandardNth2dbeta33 = PDstandardNthfdOrder82(&dbeta33[index]);
        PDstandardNth3dbeta33 = PDstandardNthfdOrder83(&dbeta33[index]);
        PDupwindNthAnti1dgt111 = PDupwindNthAntifdOrder81(&dgt111[index]);
        PDupwindNthSymm1dgt111 = PDupwindNthSymmfdOrder81(&dgt111[index]);
        PDupwindNthAnti2dgt111 = PDupwindNthAntifdOrder82(&dgt111[index]);
        PDupwindNthSymm2dgt111 = PDupwindNthSymmfdOrder82(&dgt111[index]);
        PDupwindNthAnti3dgt111 = PDupwindNthAntifdOrder83(&dgt111[index]);
        PDupwindNthSymm3dgt111 = PDupwindNthSymmfdOrder83(&dgt111[index]);
        PDupwindNthAnti1dgt112 = PDupwindNthAntifdOrder81(&dgt112[index]);
        PDupwindNthSymm1dgt112 = PDupwindNthSymmfdOrder81(&dgt112[index]);
        PDupwindNthAnti2dgt112 = PDupwindNthAntifdOrder82(&dgt112[index]);
        PDupwindNthSymm2dgt112 = PDupwindNthSymmfdOrder82(&dgt112[index]);
        PDupwindNthAnti3dgt112 = PDupwindNthAntifdOrder83(&dgt112[index]);
        PDupwindNthSymm3dgt112 = PDupwindNthSymmfdOrder83(&dgt112[index]);
        PDupwindNthAnti1dgt113 = PDupwindNthAntifdOrder81(&dgt113[index]);
        PDupwindNthSymm1dgt113 = PDupwindNthSymmfdOrder81(&dgt113[index]);
        PDupwindNthAnti2dgt113 = PDupwindNthAntifdOrder82(&dgt113[index]);
        PDupwindNthSymm2dgt113 = PDupwindNthSymmfdOrder82(&dgt113[index]);
        PDupwindNthAnti3dgt113 = PDupwindNthAntifdOrder83(&dgt113[index]);
        PDupwindNthSymm3dgt113 = PDupwindNthSymmfdOrder83(&dgt113[index]);
        PDupwindNthAnti1dgt122 = PDupwindNthAntifdOrder81(&dgt122[index]);
        PDupwindNthSymm1dgt122 = PDupwindNthSymmfdOrder81(&dgt122[index]);
        PDupwindNthAnti2dgt122 = PDupwindNthAntifdOrder82(&dgt122[index]);
        PDupwindNthSymm2dgt122 = PDupwindNthSymmfdOrder82(&dgt122[index]);
        PDupwindNthAnti3dgt122 = PDupwindNthAntifdOrder83(&dgt122[index]);
        PDupwindNthSymm3dgt122 = PDupwindNthSymmfdOrder83(&dgt122[index]);
        PDupwindNthAnti1dgt123 = PDupwindNthAntifdOrder81(&dgt123[index]);
        PDupwindNthSymm1dgt123 = PDupwindNthSymmfdOrder81(&dgt123[index]);
        PDupwindNthAnti2dgt123 = PDupwindNthAntifdOrder82(&dgt123[index]);
        PDupwindNthSymm2dgt123 = PDupwindNthSymmfdOrder82(&dgt123[index]);
        PDupwindNthAnti3dgt123 = PDupwindNthAntifdOrder83(&dgt123[index]);
        PDupwindNthSymm3dgt123 = PDupwindNthSymmfdOrder83(&dgt123[index]);
        PDupwindNthAnti1dgt133 = PDupwindNthAntifdOrder81(&dgt133[index]);
        PDupwindNthSymm1dgt133 = PDupwindNthSymmfdOrder81(&dgt133[index]);
        PDupwindNthAnti2dgt133 = PDupwindNthAntifdOrder82(&dgt133[index]);
        PDupwindNthSymm2dgt133 = PDupwindNthSymmfdOrder82(&dgt133[index]);
        PDupwindNthAnti3dgt133 = PDupwindNthAntifdOrder83(&dgt133[index]);
        PDupwindNthSymm3dgt133 = PDupwindNthSymmfdOrder83(&dgt133[index]);
        PDupwindNthAnti1dgt211 = PDupwindNthAntifdOrder81(&dgt211[index]);
        PDupwindNthSymm1dgt211 = PDupwindNthSymmfdOrder81(&dgt211[index]);
        PDupwindNthAnti2dgt211 = PDupwindNthAntifdOrder82(&dgt211[index]);
        PDupwindNthSymm2dgt211 = PDupwindNthSymmfdOrder82(&dgt211[index]);
        PDupwindNthAnti3dgt211 = PDupwindNthAntifdOrder83(&dgt211[index]);
        PDupwindNthSymm3dgt211 = PDupwindNthSymmfdOrder83(&dgt211[index]);
        PDupwindNthAnti1dgt212 = PDupwindNthAntifdOrder81(&dgt212[index]);
        PDupwindNthSymm1dgt212 = PDupwindNthSymmfdOrder81(&dgt212[index]);
        PDupwindNthAnti2dgt212 = PDupwindNthAntifdOrder82(&dgt212[index]);
        PDupwindNthSymm2dgt212 = PDupwindNthSymmfdOrder82(&dgt212[index]);
        PDupwindNthAnti3dgt212 = PDupwindNthAntifdOrder83(&dgt212[index]);
        PDupwindNthSymm3dgt212 = PDupwindNthSymmfdOrder83(&dgt212[index]);
        PDupwindNthAnti1dgt213 = PDupwindNthAntifdOrder81(&dgt213[index]);
        PDupwindNthSymm1dgt213 = PDupwindNthSymmfdOrder81(&dgt213[index]);
        PDupwindNthAnti2dgt213 = PDupwindNthAntifdOrder82(&dgt213[index]);
        PDupwindNthSymm2dgt213 = PDupwindNthSymmfdOrder82(&dgt213[index]);
        PDupwindNthAnti3dgt213 = PDupwindNthAntifdOrder83(&dgt213[index]);
        PDupwindNthSymm3dgt213 = PDupwindNthSymmfdOrder83(&dgt213[index]);
        PDupwindNthAnti1dgt222 = PDupwindNthAntifdOrder81(&dgt222[index]);
        PDupwindNthSymm1dgt222 = PDupwindNthSymmfdOrder81(&dgt222[index]);
        PDupwindNthAnti2dgt222 = PDupwindNthAntifdOrder82(&dgt222[index]);
        PDupwindNthSymm2dgt222 = PDupwindNthSymmfdOrder82(&dgt222[index]);
        PDupwindNthAnti3dgt222 = PDupwindNthAntifdOrder83(&dgt222[index]);
        PDupwindNthSymm3dgt222 = PDupwindNthSymmfdOrder83(&dgt222[index]);
        PDupwindNthAnti1dgt223 = PDupwindNthAntifdOrder81(&dgt223[index]);
        PDupwindNthSymm1dgt223 = PDupwindNthSymmfdOrder81(&dgt223[index]);
        PDupwindNthAnti2dgt223 = PDupwindNthAntifdOrder82(&dgt223[index]);
        PDupwindNthSymm2dgt223 = PDupwindNthSymmfdOrder82(&dgt223[index]);
        PDupwindNthAnti3dgt223 = PDupwindNthAntifdOrder83(&dgt223[index]);
        PDupwindNthSymm3dgt223 = PDupwindNthSymmfdOrder83(&dgt223[index]);
        PDupwindNthAnti1dgt233 = PDupwindNthAntifdOrder81(&dgt233[index]);
        PDupwindNthSymm1dgt233 = PDupwindNthSymmfdOrder81(&dgt233[index]);
        PDupwindNthAnti2dgt233 = PDupwindNthAntifdOrder82(&dgt233[index]);
        PDupwindNthSymm2dgt233 = PDupwindNthSymmfdOrder82(&dgt233[index]);
        PDupwindNthAnti3dgt233 = PDupwindNthAntifdOrder83(&dgt233[index]);
        PDupwindNthSymm3dgt233 = PDupwindNthSymmfdOrder83(&dgt233[index]);
        PDupwindNthAnti1dgt311 = PDupwindNthAntifdOrder81(&dgt311[index]);
        PDupwindNthSymm1dgt311 = PDupwindNthSymmfdOrder81(&dgt311[index]);
        PDupwindNthAnti2dgt311 = PDupwindNthAntifdOrder82(&dgt311[index]);
        PDupwindNthSymm2dgt311 = PDupwindNthSymmfdOrder82(&dgt311[index]);
        PDupwindNthAnti3dgt311 = PDupwindNthAntifdOrder83(&dgt311[index]);
        PDupwindNthSymm3dgt311 = PDupwindNthSymmfdOrder83(&dgt311[index]);
        PDupwindNthAnti1dgt312 = PDupwindNthAntifdOrder81(&dgt312[index]);
        PDupwindNthSymm1dgt312 = PDupwindNthSymmfdOrder81(&dgt312[index]);
        PDupwindNthAnti2dgt312 = PDupwindNthAntifdOrder82(&dgt312[index]);
        PDupwindNthSymm2dgt312 = PDupwindNthSymmfdOrder82(&dgt312[index]);
        PDupwindNthAnti3dgt312 = PDupwindNthAntifdOrder83(&dgt312[index]);
        PDupwindNthSymm3dgt312 = PDupwindNthSymmfdOrder83(&dgt312[index]);
        PDupwindNthAnti1dgt313 = PDupwindNthAntifdOrder81(&dgt313[index]);
        PDupwindNthSymm1dgt313 = PDupwindNthSymmfdOrder81(&dgt313[index]);
        PDupwindNthAnti2dgt313 = PDupwindNthAntifdOrder82(&dgt313[index]);
        PDupwindNthSymm2dgt313 = PDupwindNthSymmfdOrder82(&dgt313[index]);
        PDupwindNthAnti3dgt313 = PDupwindNthAntifdOrder83(&dgt313[index]);
        PDupwindNthSymm3dgt313 = PDupwindNthSymmfdOrder83(&dgt313[index]);
        PDupwindNthAnti1dgt322 = PDupwindNthAntifdOrder81(&dgt322[index]);
        PDupwindNthSymm1dgt322 = PDupwindNthSymmfdOrder81(&dgt322[index]);
        PDupwindNthAnti2dgt322 = PDupwindNthAntifdOrder82(&dgt322[index]);
        PDupwindNthSymm2dgt322 = PDupwindNthSymmfdOrder82(&dgt322[index]);
        PDupwindNthAnti3dgt322 = PDupwindNthAntifdOrder83(&dgt322[index]);
        PDupwindNthSymm3dgt322 = PDupwindNthSymmfdOrder83(&dgt322[index]);
        PDupwindNthAnti1dgt323 = PDupwindNthAntifdOrder81(&dgt323[index]);
        PDupwindNthSymm1dgt323 = PDupwindNthSymmfdOrder81(&dgt323[index]);
        PDupwindNthAnti2dgt323 = PDupwindNthAntifdOrder82(&dgt323[index]);
        PDupwindNthSymm2dgt323 = PDupwindNthSymmfdOrder82(&dgt323[index]);
        PDupwindNthAnti3dgt323 = PDupwindNthAntifdOrder83(&dgt323[index]);
        PDupwindNthSymm3dgt323 = PDupwindNthSymmfdOrder83(&dgt323[index]);
        PDupwindNthAnti1dgt333 = PDupwindNthAntifdOrder81(&dgt333[index]);
        PDupwindNthSymm1dgt333 = PDupwindNthSymmfdOrder81(&dgt333[index]);
        PDupwindNthAnti2dgt333 = PDupwindNthAntifdOrder82(&dgt333[index]);
        PDupwindNthSymm2dgt333 = PDupwindNthSymmfdOrder82(&dgt333[index]);
        PDupwindNthAnti3dgt333 = PDupwindNthAntifdOrder83(&dgt333[index]);
        PDupwindNthSymm3dgt333 = PDupwindNthSymmfdOrder83(&dgt333[index]);
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
        PDstandardNth1trK = PDstandardNthfdOrder81(&trK[index]);
        PDstandardNth2trK = PDstandardNthfdOrder82(&trK[index]);
        PDstandardNth3trK = PDstandardNthfdOrder83(&trK[index]);
        PDupwindNthAnti1Xt1 = PDupwindNthAntifdOrder81(&Xt1[index]);
        PDupwindNthSymm1Xt1 = PDupwindNthSymmfdOrder81(&Xt1[index]);
        PDupwindNthAnti2Xt1 = PDupwindNthAntifdOrder82(&Xt1[index]);
        PDupwindNthSymm2Xt1 = PDupwindNthSymmfdOrder82(&Xt1[index]);
        PDupwindNthAnti3Xt1 = PDupwindNthAntifdOrder83(&Xt1[index]);
        PDupwindNthSymm3Xt1 = PDupwindNthSymmfdOrder83(&Xt1[index]);
        PDupwindNthAnti1Xt2 = PDupwindNthAntifdOrder81(&Xt2[index]);
        PDupwindNthSymm1Xt2 = PDupwindNthSymmfdOrder81(&Xt2[index]);
        PDupwindNthAnti2Xt2 = PDupwindNthAntifdOrder82(&Xt2[index]);
        PDupwindNthSymm2Xt2 = PDupwindNthSymmfdOrder82(&Xt2[index]);
        PDupwindNthAnti3Xt2 = PDupwindNthAntifdOrder83(&Xt2[index]);
        PDupwindNthSymm3Xt2 = PDupwindNthSymmfdOrder83(&Xt2[index]);
        PDupwindNthAnti1Xt3 = PDupwindNthAntifdOrder81(&Xt3[index]);
        PDupwindNthSymm1Xt3 = PDupwindNthSymmfdOrder81(&Xt3[index]);
        PDupwindNthAnti2Xt3 = PDupwindNthAntifdOrder82(&Xt3[index]);
        PDupwindNthSymm2Xt3 = PDupwindNthSymmfdOrder82(&Xt3[index]);
        PDupwindNthAnti3Xt3 = PDupwindNthAntifdOrder83(&Xt3[index]);
        PDupwindNthSymm3Xt3 = PDupwindNthSymmfdOrder83(&Xt3[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = kisgn(beta1L);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = kisgn(beta2L);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = kisgn(beta3L);
    
    CCTK_REAL_VEC JacPDstandardNth1At11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1At12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1At13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1At22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1At23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1At33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dbeta11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dbeta12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dbeta13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dbeta21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dbeta22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dbeta23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dbeta31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dbeta32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dbeta33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2At11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2At12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2At13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2At22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2At23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2At33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dbeta11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dbeta12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dbeta13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dbeta21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dbeta22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dbeta23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dbeta31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dbeta32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dbeta33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3At11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3At12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3At13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3At22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3At23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3At33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dbeta11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dbeta12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dbeta13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dbeta21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dbeta22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dbeta23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dbeta31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dbeta32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dbeta33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3Xt3 CCTK_ATTRIBUTE_UNUSED;
    
    if (usejacobian)
    {
      JacPDstandardNth1At11 = 
        kmadd(J11L,PDstandardNth1At11,kmadd(J21L,PDstandardNth2At11,kmul(J31L,PDstandardNth3At11)));
      
      JacPDstandardNth1At12 = 
        kmadd(J11L,PDstandardNth1At12,kmadd(J21L,PDstandardNth2At12,kmul(J31L,PDstandardNth3At12)));
      
      JacPDstandardNth1At13 = 
        kmadd(J11L,PDstandardNth1At13,kmadd(J21L,PDstandardNth2At13,kmul(J31L,PDstandardNth3At13)));
      
      JacPDstandardNth1At22 = 
        kmadd(J11L,PDstandardNth1At22,kmadd(J21L,PDstandardNth2At22,kmul(J31L,PDstandardNth3At22)));
      
      JacPDstandardNth1At23 = 
        kmadd(J11L,PDstandardNth1At23,kmadd(J21L,PDstandardNth2At23,kmul(J31L,PDstandardNth3At23)));
      
      JacPDstandardNth1At33 = 
        kmadd(J11L,PDstandardNth1At33,kmadd(J21L,PDstandardNth2At33,kmul(J31L,PDstandardNth3At33)));
      
      JacPDstandardNth1dbeta11 = 
        kmadd(J11L,PDstandardNth1dbeta11,kmadd(J21L,PDstandardNth2dbeta11,kmul(J31L,PDstandardNth3dbeta11)));
      
      JacPDstandardNth1dbeta12 = 
        kmadd(J11L,PDstandardNth1dbeta12,kmadd(J21L,PDstandardNth2dbeta12,kmul(J31L,PDstandardNth3dbeta12)));
      
      JacPDstandardNth1dbeta13 = 
        kmadd(J11L,PDstandardNth1dbeta13,kmadd(J21L,PDstandardNth2dbeta13,kmul(J31L,PDstandardNth3dbeta13)));
      
      JacPDstandardNth1dbeta21 = 
        kmadd(J11L,PDstandardNth1dbeta21,kmadd(J21L,PDstandardNth2dbeta21,kmul(J31L,PDstandardNth3dbeta21)));
      
      JacPDstandardNth1dbeta22 = 
        kmadd(J11L,PDstandardNth1dbeta22,kmadd(J21L,PDstandardNth2dbeta22,kmul(J31L,PDstandardNth3dbeta22)));
      
      JacPDstandardNth1dbeta23 = 
        kmadd(J11L,PDstandardNth1dbeta23,kmadd(J21L,PDstandardNth2dbeta23,kmul(J31L,PDstandardNth3dbeta23)));
      
      JacPDstandardNth1dbeta31 = 
        kmadd(J11L,PDstandardNth1dbeta31,kmadd(J21L,PDstandardNth2dbeta31,kmul(J31L,PDstandardNth3dbeta31)));
      
      JacPDstandardNth1dbeta32 = 
        kmadd(J11L,PDstandardNth1dbeta32,kmadd(J21L,PDstandardNth2dbeta32,kmul(J31L,PDstandardNth3dbeta32)));
      
      JacPDstandardNth1dbeta33 = 
        kmadd(J11L,PDstandardNth1dbeta33,kmadd(J21L,PDstandardNth2dbeta33,kmul(J31L,PDstandardNth3dbeta33)));
      
      JacPDstandardNth1gt11 = 
        kmadd(J11L,PDstandardNth1gt11,kmadd(J21L,PDstandardNth2gt11,kmul(J31L,PDstandardNth3gt11)));
      
      JacPDstandardNth1gt12 = 
        kmadd(J11L,PDstandardNth1gt12,kmadd(J21L,PDstandardNth2gt12,kmul(J31L,PDstandardNth3gt12)));
      
      JacPDstandardNth1gt13 = 
        kmadd(J11L,PDstandardNth1gt13,kmadd(J21L,PDstandardNth2gt13,kmul(J31L,PDstandardNth3gt13)));
      
      JacPDstandardNth1gt22 = 
        kmadd(J11L,PDstandardNth1gt22,kmadd(J21L,PDstandardNth2gt22,kmul(J31L,PDstandardNth3gt22)));
      
      JacPDstandardNth1gt23 = 
        kmadd(J11L,PDstandardNth1gt23,kmadd(J21L,PDstandardNth2gt23,kmul(J31L,PDstandardNth3gt23)));
      
      JacPDstandardNth1gt33 = 
        kmadd(J11L,PDstandardNth1gt33,kmadd(J21L,PDstandardNth2gt33,kmul(J31L,PDstandardNth3gt33)));
      
      JacPDstandardNth1trK = 
        kmadd(J11L,PDstandardNth1trK,kmadd(J21L,PDstandardNth2trK,kmul(J31L,PDstandardNth3trK)));
      
      JacPDstandardNth2At11 = 
        kmadd(J12L,PDstandardNth1At11,kmadd(J22L,PDstandardNth2At11,kmul(J32L,PDstandardNth3At11)));
      
      JacPDstandardNth2At12 = 
        kmadd(J12L,PDstandardNth1At12,kmadd(J22L,PDstandardNth2At12,kmul(J32L,PDstandardNth3At12)));
      
      JacPDstandardNth2At13 = 
        kmadd(J12L,PDstandardNth1At13,kmadd(J22L,PDstandardNth2At13,kmul(J32L,PDstandardNth3At13)));
      
      JacPDstandardNth2At22 = 
        kmadd(J12L,PDstandardNth1At22,kmadd(J22L,PDstandardNth2At22,kmul(J32L,PDstandardNth3At22)));
      
      JacPDstandardNth2At23 = 
        kmadd(J12L,PDstandardNth1At23,kmadd(J22L,PDstandardNth2At23,kmul(J32L,PDstandardNth3At23)));
      
      JacPDstandardNth2At33 = 
        kmadd(J12L,PDstandardNth1At33,kmadd(J22L,PDstandardNth2At33,kmul(J32L,PDstandardNth3At33)));
      
      JacPDstandardNth2dbeta11 = 
        kmadd(J12L,PDstandardNth1dbeta11,kmadd(J22L,PDstandardNth2dbeta11,kmul(J32L,PDstandardNth3dbeta11)));
      
      JacPDstandardNth2dbeta12 = 
        kmadd(J12L,PDstandardNth1dbeta12,kmadd(J22L,PDstandardNth2dbeta12,kmul(J32L,PDstandardNth3dbeta12)));
      
      JacPDstandardNth2dbeta13 = 
        kmadd(J12L,PDstandardNth1dbeta13,kmadd(J22L,PDstandardNth2dbeta13,kmul(J32L,PDstandardNth3dbeta13)));
      
      JacPDstandardNth2dbeta21 = 
        kmadd(J12L,PDstandardNth1dbeta21,kmadd(J22L,PDstandardNth2dbeta21,kmul(J32L,PDstandardNth3dbeta21)));
      
      JacPDstandardNth2dbeta22 = 
        kmadd(J12L,PDstandardNth1dbeta22,kmadd(J22L,PDstandardNth2dbeta22,kmul(J32L,PDstandardNth3dbeta22)));
      
      JacPDstandardNth2dbeta23 = 
        kmadd(J12L,PDstandardNth1dbeta23,kmadd(J22L,PDstandardNth2dbeta23,kmul(J32L,PDstandardNth3dbeta23)));
      
      JacPDstandardNth2dbeta31 = 
        kmadd(J12L,PDstandardNth1dbeta31,kmadd(J22L,PDstandardNth2dbeta31,kmul(J32L,PDstandardNth3dbeta31)));
      
      JacPDstandardNth2dbeta32 = 
        kmadd(J12L,PDstandardNth1dbeta32,kmadd(J22L,PDstandardNth2dbeta32,kmul(J32L,PDstandardNth3dbeta32)));
      
      JacPDstandardNth2dbeta33 = 
        kmadd(J12L,PDstandardNth1dbeta33,kmadd(J22L,PDstandardNth2dbeta33,kmul(J32L,PDstandardNth3dbeta33)));
      
      JacPDstandardNth2gt11 = 
        kmadd(J12L,PDstandardNth1gt11,kmadd(J22L,PDstandardNth2gt11,kmul(J32L,PDstandardNth3gt11)));
      
      JacPDstandardNth2gt12 = 
        kmadd(J12L,PDstandardNth1gt12,kmadd(J22L,PDstandardNth2gt12,kmul(J32L,PDstandardNth3gt12)));
      
      JacPDstandardNth2gt13 = 
        kmadd(J12L,PDstandardNth1gt13,kmadd(J22L,PDstandardNth2gt13,kmul(J32L,PDstandardNth3gt13)));
      
      JacPDstandardNth2gt22 = 
        kmadd(J12L,PDstandardNth1gt22,kmadd(J22L,PDstandardNth2gt22,kmul(J32L,PDstandardNth3gt22)));
      
      JacPDstandardNth2gt23 = 
        kmadd(J12L,PDstandardNth1gt23,kmadd(J22L,PDstandardNth2gt23,kmul(J32L,PDstandardNth3gt23)));
      
      JacPDstandardNth2gt33 = 
        kmadd(J12L,PDstandardNth1gt33,kmadd(J22L,PDstandardNth2gt33,kmul(J32L,PDstandardNth3gt33)));
      
      JacPDstandardNth2trK = 
        kmadd(J12L,PDstandardNth1trK,kmadd(J22L,PDstandardNth2trK,kmul(J32L,PDstandardNth3trK)));
      
      JacPDstandardNth3At11 = 
        kmadd(J13L,PDstandardNth1At11,kmadd(J23L,PDstandardNth2At11,kmul(J33L,PDstandardNth3At11)));
      
      JacPDstandardNth3At12 = 
        kmadd(J13L,PDstandardNth1At12,kmadd(J23L,PDstandardNth2At12,kmul(J33L,PDstandardNth3At12)));
      
      JacPDstandardNth3At13 = 
        kmadd(J13L,PDstandardNth1At13,kmadd(J23L,PDstandardNth2At13,kmul(J33L,PDstandardNth3At13)));
      
      JacPDstandardNth3At22 = 
        kmadd(J13L,PDstandardNth1At22,kmadd(J23L,PDstandardNth2At22,kmul(J33L,PDstandardNth3At22)));
      
      JacPDstandardNth3At23 = 
        kmadd(J13L,PDstandardNth1At23,kmadd(J23L,PDstandardNth2At23,kmul(J33L,PDstandardNth3At23)));
      
      JacPDstandardNth3At33 = 
        kmadd(J13L,PDstandardNth1At33,kmadd(J23L,PDstandardNth2At33,kmul(J33L,PDstandardNth3At33)));
      
      JacPDstandardNth3dbeta11 = 
        kmadd(J13L,PDstandardNth1dbeta11,kmadd(J23L,PDstandardNth2dbeta11,kmul(J33L,PDstandardNth3dbeta11)));
      
      JacPDstandardNth3dbeta12 = 
        kmadd(J13L,PDstandardNth1dbeta12,kmadd(J23L,PDstandardNth2dbeta12,kmul(J33L,PDstandardNth3dbeta12)));
      
      JacPDstandardNth3dbeta13 = 
        kmadd(J13L,PDstandardNth1dbeta13,kmadd(J23L,PDstandardNth2dbeta13,kmul(J33L,PDstandardNth3dbeta13)));
      
      JacPDstandardNth3dbeta21 = 
        kmadd(J13L,PDstandardNth1dbeta21,kmadd(J23L,PDstandardNth2dbeta21,kmul(J33L,PDstandardNth3dbeta21)));
      
      JacPDstandardNth3dbeta22 = 
        kmadd(J13L,PDstandardNth1dbeta22,kmadd(J23L,PDstandardNth2dbeta22,kmul(J33L,PDstandardNth3dbeta22)));
      
      JacPDstandardNth3dbeta23 = 
        kmadd(J13L,PDstandardNth1dbeta23,kmadd(J23L,PDstandardNth2dbeta23,kmul(J33L,PDstandardNth3dbeta23)));
      
      JacPDstandardNth3dbeta31 = 
        kmadd(J13L,PDstandardNth1dbeta31,kmadd(J23L,PDstandardNth2dbeta31,kmul(J33L,PDstandardNth3dbeta31)));
      
      JacPDstandardNth3dbeta32 = 
        kmadd(J13L,PDstandardNth1dbeta32,kmadd(J23L,PDstandardNth2dbeta32,kmul(J33L,PDstandardNth3dbeta32)));
      
      JacPDstandardNth3dbeta33 = 
        kmadd(J13L,PDstandardNth1dbeta33,kmadd(J23L,PDstandardNth2dbeta33,kmul(J33L,PDstandardNth3dbeta33)));
      
      JacPDstandardNth3gt11 = 
        kmadd(J13L,PDstandardNth1gt11,kmadd(J23L,PDstandardNth2gt11,kmul(J33L,PDstandardNth3gt11)));
      
      JacPDstandardNth3gt12 = 
        kmadd(J13L,PDstandardNth1gt12,kmadd(J23L,PDstandardNth2gt12,kmul(J33L,PDstandardNth3gt12)));
      
      JacPDstandardNth3gt13 = 
        kmadd(J13L,PDstandardNth1gt13,kmadd(J23L,PDstandardNth2gt13,kmul(J33L,PDstandardNth3gt13)));
      
      JacPDstandardNth3gt22 = 
        kmadd(J13L,PDstandardNth1gt22,kmadd(J23L,PDstandardNth2gt22,kmul(J33L,PDstandardNth3gt22)));
      
      JacPDstandardNth3gt23 = 
        kmadd(J13L,PDstandardNth1gt23,kmadd(J23L,PDstandardNth2gt23,kmul(J33L,PDstandardNth3gt23)));
      
      JacPDstandardNth3gt33 = 
        kmadd(J13L,PDstandardNth1gt33,kmadd(J23L,PDstandardNth2gt33,kmul(J33L,PDstandardNth3gt33)));
      
      JacPDstandardNth3trK = 
        kmadd(J13L,PDstandardNth1trK,kmadd(J23L,PDstandardNth2trK,kmul(J33L,PDstandardNth3trK)));
      
      JacPDupwindNthAnti1dgt111 = 
        kmadd(J11L,PDupwindNthAnti1dgt111,kmadd(J21L,PDupwindNthAnti2dgt111,kmul(J31L,PDupwindNthAnti3dgt111)));
      
      JacPDupwindNthAnti1dgt112 = 
        kmadd(J11L,PDupwindNthAnti1dgt112,kmadd(J21L,PDupwindNthAnti2dgt112,kmul(J31L,PDupwindNthAnti3dgt112)));
      
      JacPDupwindNthAnti1dgt113 = 
        kmadd(J11L,PDupwindNthAnti1dgt113,kmadd(J21L,PDupwindNthAnti2dgt113,kmul(J31L,PDupwindNthAnti3dgt113)));
      
      JacPDupwindNthAnti1dgt122 = 
        kmadd(J11L,PDupwindNthAnti1dgt122,kmadd(J21L,PDupwindNthAnti2dgt122,kmul(J31L,PDupwindNthAnti3dgt122)));
      
      JacPDupwindNthAnti1dgt123 = 
        kmadd(J11L,PDupwindNthAnti1dgt123,kmadd(J21L,PDupwindNthAnti2dgt123,kmul(J31L,PDupwindNthAnti3dgt123)));
      
      JacPDupwindNthAnti1dgt133 = 
        kmadd(J11L,PDupwindNthAnti1dgt133,kmadd(J21L,PDupwindNthAnti2dgt133,kmul(J31L,PDupwindNthAnti3dgt133)));
      
      JacPDupwindNthAnti1dgt211 = 
        kmadd(J11L,PDupwindNthAnti1dgt211,kmadd(J21L,PDupwindNthAnti2dgt211,kmul(J31L,PDupwindNthAnti3dgt211)));
      
      JacPDupwindNthAnti1dgt212 = 
        kmadd(J11L,PDupwindNthAnti1dgt212,kmadd(J21L,PDupwindNthAnti2dgt212,kmul(J31L,PDupwindNthAnti3dgt212)));
      
      JacPDupwindNthAnti1dgt213 = 
        kmadd(J11L,PDupwindNthAnti1dgt213,kmadd(J21L,PDupwindNthAnti2dgt213,kmul(J31L,PDupwindNthAnti3dgt213)));
      
      JacPDupwindNthAnti1dgt222 = 
        kmadd(J11L,PDupwindNthAnti1dgt222,kmadd(J21L,PDupwindNthAnti2dgt222,kmul(J31L,PDupwindNthAnti3dgt222)));
      
      JacPDupwindNthAnti1dgt223 = 
        kmadd(J11L,PDupwindNthAnti1dgt223,kmadd(J21L,PDupwindNthAnti2dgt223,kmul(J31L,PDupwindNthAnti3dgt223)));
      
      JacPDupwindNthAnti1dgt233 = 
        kmadd(J11L,PDupwindNthAnti1dgt233,kmadd(J21L,PDupwindNthAnti2dgt233,kmul(J31L,PDupwindNthAnti3dgt233)));
      
      JacPDupwindNthAnti1dgt311 = 
        kmadd(J11L,PDupwindNthAnti1dgt311,kmadd(J21L,PDupwindNthAnti2dgt311,kmul(J31L,PDupwindNthAnti3dgt311)));
      
      JacPDupwindNthAnti1dgt312 = 
        kmadd(J11L,PDupwindNthAnti1dgt312,kmadd(J21L,PDupwindNthAnti2dgt312,kmul(J31L,PDupwindNthAnti3dgt312)));
      
      JacPDupwindNthAnti1dgt313 = 
        kmadd(J11L,PDupwindNthAnti1dgt313,kmadd(J21L,PDupwindNthAnti2dgt313,kmul(J31L,PDupwindNthAnti3dgt313)));
      
      JacPDupwindNthAnti1dgt322 = 
        kmadd(J11L,PDupwindNthAnti1dgt322,kmadd(J21L,PDupwindNthAnti2dgt322,kmul(J31L,PDupwindNthAnti3dgt322)));
      
      JacPDupwindNthAnti1dgt323 = 
        kmadd(J11L,PDupwindNthAnti1dgt323,kmadd(J21L,PDupwindNthAnti2dgt323,kmul(J31L,PDupwindNthAnti3dgt323)));
      
      JacPDupwindNthAnti1dgt333 = 
        kmadd(J11L,PDupwindNthAnti1dgt333,kmadd(J21L,PDupwindNthAnti2dgt333,kmul(J31L,PDupwindNthAnti3dgt333)));
      
      JacPDupwindNthAnti1Xt1 = 
        kmadd(J11L,PDupwindNthAnti1Xt1,kmadd(J21L,PDupwindNthAnti2Xt1,kmul(J31L,PDupwindNthAnti3Xt1)));
      
      JacPDupwindNthAnti1Xt2 = 
        kmadd(J11L,PDupwindNthAnti1Xt2,kmadd(J21L,PDupwindNthAnti2Xt2,kmul(J31L,PDupwindNthAnti3Xt2)));
      
      JacPDupwindNthAnti1Xt3 = 
        kmadd(J11L,PDupwindNthAnti1Xt3,kmadd(J21L,PDupwindNthAnti2Xt3,kmul(J31L,PDupwindNthAnti3Xt3)));
      
      JacPDupwindNthSymm1dgt111 = 
        kmadd(J11L,PDupwindNthSymm1dgt111,kmadd(J21L,PDupwindNthSymm2dgt111,kmul(J31L,PDupwindNthSymm3dgt111)));
      
      JacPDupwindNthSymm1dgt112 = 
        kmadd(J11L,PDupwindNthSymm1dgt112,kmadd(J21L,PDupwindNthSymm2dgt112,kmul(J31L,PDupwindNthSymm3dgt112)));
      
      JacPDupwindNthSymm1dgt113 = 
        kmadd(J11L,PDupwindNthSymm1dgt113,kmadd(J21L,PDupwindNthSymm2dgt113,kmul(J31L,PDupwindNthSymm3dgt113)));
      
      JacPDupwindNthSymm1dgt122 = 
        kmadd(J11L,PDupwindNthSymm1dgt122,kmadd(J21L,PDupwindNthSymm2dgt122,kmul(J31L,PDupwindNthSymm3dgt122)));
      
      JacPDupwindNthSymm1dgt123 = 
        kmadd(J11L,PDupwindNthSymm1dgt123,kmadd(J21L,PDupwindNthSymm2dgt123,kmul(J31L,PDupwindNthSymm3dgt123)));
      
      JacPDupwindNthSymm1dgt133 = 
        kmadd(J11L,PDupwindNthSymm1dgt133,kmadd(J21L,PDupwindNthSymm2dgt133,kmul(J31L,PDupwindNthSymm3dgt133)));
      
      JacPDupwindNthSymm1dgt211 = 
        kmadd(J11L,PDupwindNthSymm1dgt211,kmadd(J21L,PDupwindNthSymm2dgt211,kmul(J31L,PDupwindNthSymm3dgt211)));
      
      JacPDupwindNthSymm1dgt212 = 
        kmadd(J11L,PDupwindNthSymm1dgt212,kmadd(J21L,PDupwindNthSymm2dgt212,kmul(J31L,PDupwindNthSymm3dgt212)));
      
      JacPDupwindNthSymm1dgt213 = 
        kmadd(J11L,PDupwindNthSymm1dgt213,kmadd(J21L,PDupwindNthSymm2dgt213,kmul(J31L,PDupwindNthSymm3dgt213)));
      
      JacPDupwindNthSymm1dgt222 = 
        kmadd(J11L,PDupwindNthSymm1dgt222,kmadd(J21L,PDupwindNthSymm2dgt222,kmul(J31L,PDupwindNthSymm3dgt222)));
      
      JacPDupwindNthSymm1dgt223 = 
        kmadd(J11L,PDupwindNthSymm1dgt223,kmadd(J21L,PDupwindNthSymm2dgt223,kmul(J31L,PDupwindNthSymm3dgt223)));
      
      JacPDupwindNthSymm1dgt233 = 
        kmadd(J11L,PDupwindNthSymm1dgt233,kmadd(J21L,PDupwindNthSymm2dgt233,kmul(J31L,PDupwindNthSymm3dgt233)));
      
      JacPDupwindNthSymm1dgt311 = 
        kmadd(J11L,PDupwindNthSymm1dgt311,kmadd(J21L,PDupwindNthSymm2dgt311,kmul(J31L,PDupwindNthSymm3dgt311)));
      
      JacPDupwindNthSymm1dgt312 = 
        kmadd(J11L,PDupwindNthSymm1dgt312,kmadd(J21L,PDupwindNthSymm2dgt312,kmul(J31L,PDupwindNthSymm3dgt312)));
      
      JacPDupwindNthSymm1dgt313 = 
        kmadd(J11L,PDupwindNthSymm1dgt313,kmadd(J21L,PDupwindNthSymm2dgt313,kmul(J31L,PDupwindNthSymm3dgt313)));
      
      JacPDupwindNthSymm1dgt322 = 
        kmadd(J11L,PDupwindNthSymm1dgt322,kmadd(J21L,PDupwindNthSymm2dgt322,kmul(J31L,PDupwindNthSymm3dgt322)));
      
      JacPDupwindNthSymm1dgt323 = 
        kmadd(J11L,PDupwindNthSymm1dgt323,kmadd(J21L,PDupwindNthSymm2dgt323,kmul(J31L,PDupwindNthSymm3dgt323)));
      
      JacPDupwindNthSymm1dgt333 = 
        kmadd(J11L,PDupwindNthSymm1dgt333,kmadd(J21L,PDupwindNthSymm2dgt333,kmul(J31L,PDupwindNthSymm3dgt333)));
      
      JacPDupwindNthSymm1Xt1 = 
        kmadd(J11L,PDupwindNthSymm1Xt1,kmadd(J21L,PDupwindNthSymm2Xt1,kmul(J31L,PDupwindNthSymm3Xt1)));
      
      JacPDupwindNthSymm1Xt2 = 
        kmadd(J11L,PDupwindNthSymm1Xt2,kmadd(J21L,PDupwindNthSymm2Xt2,kmul(J31L,PDupwindNthSymm3Xt2)));
      
      JacPDupwindNthSymm1Xt3 = 
        kmadd(J11L,PDupwindNthSymm1Xt3,kmadd(J21L,PDupwindNthSymm2Xt3,kmul(J31L,PDupwindNthSymm3Xt3)));
      
      JacPDupwindNthAnti2dgt111 = 
        kmadd(J12L,PDupwindNthAnti1dgt111,kmadd(J22L,PDupwindNthAnti2dgt111,kmul(J32L,PDupwindNthAnti3dgt111)));
      
      JacPDupwindNthAnti2dgt112 = 
        kmadd(J12L,PDupwindNthAnti1dgt112,kmadd(J22L,PDupwindNthAnti2dgt112,kmul(J32L,PDupwindNthAnti3dgt112)));
      
      JacPDupwindNthAnti2dgt113 = 
        kmadd(J12L,PDupwindNthAnti1dgt113,kmadd(J22L,PDupwindNthAnti2dgt113,kmul(J32L,PDupwindNthAnti3dgt113)));
      
      JacPDupwindNthAnti2dgt122 = 
        kmadd(J12L,PDupwindNthAnti1dgt122,kmadd(J22L,PDupwindNthAnti2dgt122,kmul(J32L,PDupwindNthAnti3dgt122)));
      
      JacPDupwindNthAnti2dgt123 = 
        kmadd(J12L,PDupwindNthAnti1dgt123,kmadd(J22L,PDupwindNthAnti2dgt123,kmul(J32L,PDupwindNthAnti3dgt123)));
      
      JacPDupwindNthAnti2dgt133 = 
        kmadd(J12L,PDupwindNthAnti1dgt133,kmadd(J22L,PDupwindNthAnti2dgt133,kmul(J32L,PDupwindNthAnti3dgt133)));
      
      JacPDupwindNthAnti2dgt211 = 
        kmadd(J12L,PDupwindNthAnti1dgt211,kmadd(J22L,PDupwindNthAnti2dgt211,kmul(J32L,PDupwindNthAnti3dgt211)));
      
      JacPDupwindNthAnti2dgt212 = 
        kmadd(J12L,PDupwindNthAnti1dgt212,kmadd(J22L,PDupwindNthAnti2dgt212,kmul(J32L,PDupwindNthAnti3dgt212)));
      
      JacPDupwindNthAnti2dgt213 = 
        kmadd(J12L,PDupwindNthAnti1dgt213,kmadd(J22L,PDupwindNthAnti2dgt213,kmul(J32L,PDupwindNthAnti3dgt213)));
      
      JacPDupwindNthAnti2dgt222 = 
        kmadd(J12L,PDupwindNthAnti1dgt222,kmadd(J22L,PDupwindNthAnti2dgt222,kmul(J32L,PDupwindNthAnti3dgt222)));
      
      JacPDupwindNthAnti2dgt223 = 
        kmadd(J12L,PDupwindNthAnti1dgt223,kmadd(J22L,PDupwindNthAnti2dgt223,kmul(J32L,PDupwindNthAnti3dgt223)));
      
      JacPDupwindNthAnti2dgt233 = 
        kmadd(J12L,PDupwindNthAnti1dgt233,kmadd(J22L,PDupwindNthAnti2dgt233,kmul(J32L,PDupwindNthAnti3dgt233)));
      
      JacPDupwindNthAnti2dgt311 = 
        kmadd(J12L,PDupwindNthAnti1dgt311,kmadd(J22L,PDupwindNthAnti2dgt311,kmul(J32L,PDupwindNthAnti3dgt311)));
      
      JacPDupwindNthAnti2dgt312 = 
        kmadd(J12L,PDupwindNthAnti1dgt312,kmadd(J22L,PDupwindNthAnti2dgt312,kmul(J32L,PDupwindNthAnti3dgt312)));
      
      JacPDupwindNthAnti2dgt313 = 
        kmadd(J12L,PDupwindNthAnti1dgt313,kmadd(J22L,PDupwindNthAnti2dgt313,kmul(J32L,PDupwindNthAnti3dgt313)));
      
      JacPDupwindNthAnti2dgt322 = 
        kmadd(J12L,PDupwindNthAnti1dgt322,kmadd(J22L,PDupwindNthAnti2dgt322,kmul(J32L,PDupwindNthAnti3dgt322)));
      
      JacPDupwindNthAnti2dgt323 = 
        kmadd(J12L,PDupwindNthAnti1dgt323,kmadd(J22L,PDupwindNthAnti2dgt323,kmul(J32L,PDupwindNthAnti3dgt323)));
      
      JacPDupwindNthAnti2dgt333 = 
        kmadd(J12L,PDupwindNthAnti1dgt333,kmadd(J22L,PDupwindNthAnti2dgt333,kmul(J32L,PDupwindNthAnti3dgt333)));
      
      JacPDupwindNthAnti2Xt1 = 
        kmadd(J12L,PDupwindNthAnti1Xt1,kmadd(J22L,PDupwindNthAnti2Xt1,kmul(J32L,PDupwindNthAnti3Xt1)));
      
      JacPDupwindNthAnti2Xt2 = 
        kmadd(J12L,PDupwindNthAnti1Xt2,kmadd(J22L,PDupwindNthAnti2Xt2,kmul(J32L,PDupwindNthAnti3Xt2)));
      
      JacPDupwindNthAnti2Xt3 = 
        kmadd(J12L,PDupwindNthAnti1Xt3,kmadd(J22L,PDupwindNthAnti2Xt3,kmul(J32L,PDupwindNthAnti3Xt3)));
      
      JacPDupwindNthSymm2dgt111 = 
        kmadd(J12L,PDupwindNthSymm1dgt111,kmadd(J22L,PDupwindNthSymm2dgt111,kmul(J32L,PDupwindNthSymm3dgt111)));
      
      JacPDupwindNthSymm2dgt112 = 
        kmadd(J12L,PDupwindNthSymm1dgt112,kmadd(J22L,PDupwindNthSymm2dgt112,kmul(J32L,PDupwindNthSymm3dgt112)));
      
      JacPDupwindNthSymm2dgt113 = 
        kmadd(J12L,PDupwindNthSymm1dgt113,kmadd(J22L,PDupwindNthSymm2dgt113,kmul(J32L,PDupwindNthSymm3dgt113)));
      
      JacPDupwindNthSymm2dgt122 = 
        kmadd(J12L,PDupwindNthSymm1dgt122,kmadd(J22L,PDupwindNthSymm2dgt122,kmul(J32L,PDupwindNthSymm3dgt122)));
      
      JacPDupwindNthSymm2dgt123 = 
        kmadd(J12L,PDupwindNthSymm1dgt123,kmadd(J22L,PDupwindNthSymm2dgt123,kmul(J32L,PDupwindNthSymm3dgt123)));
      
      JacPDupwindNthSymm2dgt133 = 
        kmadd(J12L,PDupwindNthSymm1dgt133,kmadd(J22L,PDupwindNthSymm2dgt133,kmul(J32L,PDupwindNthSymm3dgt133)));
      
      JacPDupwindNthSymm2dgt211 = 
        kmadd(J12L,PDupwindNthSymm1dgt211,kmadd(J22L,PDupwindNthSymm2dgt211,kmul(J32L,PDupwindNthSymm3dgt211)));
      
      JacPDupwindNthSymm2dgt212 = 
        kmadd(J12L,PDupwindNthSymm1dgt212,kmadd(J22L,PDupwindNthSymm2dgt212,kmul(J32L,PDupwindNthSymm3dgt212)));
      
      JacPDupwindNthSymm2dgt213 = 
        kmadd(J12L,PDupwindNthSymm1dgt213,kmadd(J22L,PDupwindNthSymm2dgt213,kmul(J32L,PDupwindNthSymm3dgt213)));
      
      JacPDupwindNthSymm2dgt222 = 
        kmadd(J12L,PDupwindNthSymm1dgt222,kmadd(J22L,PDupwindNthSymm2dgt222,kmul(J32L,PDupwindNthSymm3dgt222)));
      
      JacPDupwindNthSymm2dgt223 = 
        kmadd(J12L,PDupwindNthSymm1dgt223,kmadd(J22L,PDupwindNthSymm2dgt223,kmul(J32L,PDupwindNthSymm3dgt223)));
      
      JacPDupwindNthSymm2dgt233 = 
        kmadd(J12L,PDupwindNthSymm1dgt233,kmadd(J22L,PDupwindNthSymm2dgt233,kmul(J32L,PDupwindNthSymm3dgt233)));
      
      JacPDupwindNthSymm2dgt311 = 
        kmadd(J12L,PDupwindNthSymm1dgt311,kmadd(J22L,PDupwindNthSymm2dgt311,kmul(J32L,PDupwindNthSymm3dgt311)));
      
      JacPDupwindNthSymm2dgt312 = 
        kmadd(J12L,PDupwindNthSymm1dgt312,kmadd(J22L,PDupwindNthSymm2dgt312,kmul(J32L,PDupwindNthSymm3dgt312)));
      
      JacPDupwindNthSymm2dgt313 = 
        kmadd(J12L,PDupwindNthSymm1dgt313,kmadd(J22L,PDupwindNthSymm2dgt313,kmul(J32L,PDupwindNthSymm3dgt313)));
      
      JacPDupwindNthSymm2dgt322 = 
        kmadd(J12L,PDupwindNthSymm1dgt322,kmadd(J22L,PDupwindNthSymm2dgt322,kmul(J32L,PDupwindNthSymm3dgt322)));
      
      JacPDupwindNthSymm2dgt323 = 
        kmadd(J12L,PDupwindNthSymm1dgt323,kmadd(J22L,PDupwindNthSymm2dgt323,kmul(J32L,PDupwindNthSymm3dgt323)));
      
      JacPDupwindNthSymm2dgt333 = 
        kmadd(J12L,PDupwindNthSymm1dgt333,kmadd(J22L,PDupwindNthSymm2dgt333,kmul(J32L,PDupwindNthSymm3dgt333)));
      
      JacPDupwindNthSymm2Xt1 = 
        kmadd(J12L,PDupwindNthSymm1Xt1,kmadd(J22L,PDupwindNthSymm2Xt1,kmul(J32L,PDupwindNthSymm3Xt1)));
      
      JacPDupwindNthSymm2Xt2 = 
        kmadd(J12L,PDupwindNthSymm1Xt2,kmadd(J22L,PDupwindNthSymm2Xt2,kmul(J32L,PDupwindNthSymm3Xt2)));
      
      JacPDupwindNthSymm2Xt3 = 
        kmadd(J12L,PDupwindNthSymm1Xt3,kmadd(J22L,PDupwindNthSymm2Xt3,kmul(J32L,PDupwindNthSymm3Xt3)));
      
      JacPDupwindNthAnti3dgt111 = 
        kmadd(J13L,PDupwindNthAnti1dgt111,kmadd(J23L,PDupwindNthAnti2dgt111,kmul(J33L,PDupwindNthAnti3dgt111)));
      
      JacPDupwindNthAnti3dgt112 = 
        kmadd(J13L,PDupwindNthAnti1dgt112,kmadd(J23L,PDupwindNthAnti2dgt112,kmul(J33L,PDupwindNthAnti3dgt112)));
      
      JacPDupwindNthAnti3dgt113 = 
        kmadd(J13L,PDupwindNthAnti1dgt113,kmadd(J23L,PDupwindNthAnti2dgt113,kmul(J33L,PDupwindNthAnti3dgt113)));
      
      JacPDupwindNthAnti3dgt122 = 
        kmadd(J13L,PDupwindNthAnti1dgt122,kmadd(J23L,PDupwindNthAnti2dgt122,kmul(J33L,PDupwindNthAnti3dgt122)));
      
      JacPDupwindNthAnti3dgt123 = 
        kmadd(J13L,PDupwindNthAnti1dgt123,kmadd(J23L,PDupwindNthAnti2dgt123,kmul(J33L,PDupwindNthAnti3dgt123)));
      
      JacPDupwindNthAnti3dgt133 = 
        kmadd(J13L,PDupwindNthAnti1dgt133,kmadd(J23L,PDupwindNthAnti2dgt133,kmul(J33L,PDupwindNthAnti3dgt133)));
      
      JacPDupwindNthAnti3dgt211 = 
        kmadd(J13L,PDupwindNthAnti1dgt211,kmadd(J23L,PDupwindNthAnti2dgt211,kmul(J33L,PDupwindNthAnti3dgt211)));
      
      JacPDupwindNthAnti3dgt212 = 
        kmadd(J13L,PDupwindNthAnti1dgt212,kmadd(J23L,PDupwindNthAnti2dgt212,kmul(J33L,PDupwindNthAnti3dgt212)));
      
      JacPDupwindNthAnti3dgt213 = 
        kmadd(J13L,PDupwindNthAnti1dgt213,kmadd(J23L,PDupwindNthAnti2dgt213,kmul(J33L,PDupwindNthAnti3dgt213)));
      
      JacPDupwindNthAnti3dgt222 = 
        kmadd(J13L,PDupwindNthAnti1dgt222,kmadd(J23L,PDupwindNthAnti2dgt222,kmul(J33L,PDupwindNthAnti3dgt222)));
      
      JacPDupwindNthAnti3dgt223 = 
        kmadd(J13L,PDupwindNthAnti1dgt223,kmadd(J23L,PDupwindNthAnti2dgt223,kmul(J33L,PDupwindNthAnti3dgt223)));
      
      JacPDupwindNthAnti3dgt233 = 
        kmadd(J13L,PDupwindNthAnti1dgt233,kmadd(J23L,PDupwindNthAnti2dgt233,kmul(J33L,PDupwindNthAnti3dgt233)));
      
      JacPDupwindNthAnti3dgt311 = 
        kmadd(J13L,PDupwindNthAnti1dgt311,kmadd(J23L,PDupwindNthAnti2dgt311,kmul(J33L,PDupwindNthAnti3dgt311)));
      
      JacPDupwindNthAnti3dgt312 = 
        kmadd(J13L,PDupwindNthAnti1dgt312,kmadd(J23L,PDupwindNthAnti2dgt312,kmul(J33L,PDupwindNthAnti3dgt312)));
      
      JacPDupwindNthAnti3dgt313 = 
        kmadd(J13L,PDupwindNthAnti1dgt313,kmadd(J23L,PDupwindNthAnti2dgt313,kmul(J33L,PDupwindNthAnti3dgt313)));
      
      JacPDupwindNthAnti3dgt322 = 
        kmadd(J13L,PDupwindNthAnti1dgt322,kmadd(J23L,PDupwindNthAnti2dgt322,kmul(J33L,PDupwindNthAnti3dgt322)));
      
      JacPDupwindNthAnti3dgt323 = 
        kmadd(J13L,PDupwindNthAnti1dgt323,kmadd(J23L,PDupwindNthAnti2dgt323,kmul(J33L,PDupwindNthAnti3dgt323)));
      
      JacPDupwindNthAnti3dgt333 = 
        kmadd(J13L,PDupwindNthAnti1dgt333,kmadd(J23L,PDupwindNthAnti2dgt333,kmul(J33L,PDupwindNthAnti3dgt333)));
      
      JacPDupwindNthAnti3Xt1 = 
        kmadd(J13L,PDupwindNthAnti1Xt1,kmadd(J23L,PDupwindNthAnti2Xt1,kmul(J33L,PDupwindNthAnti3Xt1)));
      
      JacPDupwindNthAnti3Xt2 = 
        kmadd(J13L,PDupwindNthAnti1Xt2,kmadd(J23L,PDupwindNthAnti2Xt2,kmul(J33L,PDupwindNthAnti3Xt2)));
      
      JacPDupwindNthAnti3Xt3 = 
        kmadd(J13L,PDupwindNthAnti1Xt3,kmadd(J23L,PDupwindNthAnti2Xt3,kmul(J33L,PDupwindNthAnti3Xt3)));
      
      JacPDupwindNthSymm3dgt111 = 
        kmadd(J13L,PDupwindNthSymm1dgt111,kmadd(J23L,PDupwindNthSymm2dgt111,kmul(J33L,PDupwindNthSymm3dgt111)));
      
      JacPDupwindNthSymm3dgt112 = 
        kmadd(J13L,PDupwindNthSymm1dgt112,kmadd(J23L,PDupwindNthSymm2dgt112,kmul(J33L,PDupwindNthSymm3dgt112)));
      
      JacPDupwindNthSymm3dgt113 = 
        kmadd(J13L,PDupwindNthSymm1dgt113,kmadd(J23L,PDupwindNthSymm2dgt113,kmul(J33L,PDupwindNthSymm3dgt113)));
      
      JacPDupwindNthSymm3dgt122 = 
        kmadd(J13L,PDupwindNthSymm1dgt122,kmadd(J23L,PDupwindNthSymm2dgt122,kmul(J33L,PDupwindNthSymm3dgt122)));
      
      JacPDupwindNthSymm3dgt123 = 
        kmadd(J13L,PDupwindNthSymm1dgt123,kmadd(J23L,PDupwindNthSymm2dgt123,kmul(J33L,PDupwindNthSymm3dgt123)));
      
      JacPDupwindNthSymm3dgt133 = 
        kmadd(J13L,PDupwindNthSymm1dgt133,kmadd(J23L,PDupwindNthSymm2dgt133,kmul(J33L,PDupwindNthSymm3dgt133)));
      
      JacPDupwindNthSymm3dgt211 = 
        kmadd(J13L,PDupwindNthSymm1dgt211,kmadd(J23L,PDupwindNthSymm2dgt211,kmul(J33L,PDupwindNthSymm3dgt211)));
      
      JacPDupwindNthSymm3dgt212 = 
        kmadd(J13L,PDupwindNthSymm1dgt212,kmadd(J23L,PDupwindNthSymm2dgt212,kmul(J33L,PDupwindNthSymm3dgt212)));
      
      JacPDupwindNthSymm3dgt213 = 
        kmadd(J13L,PDupwindNthSymm1dgt213,kmadd(J23L,PDupwindNthSymm2dgt213,kmul(J33L,PDupwindNthSymm3dgt213)));
      
      JacPDupwindNthSymm3dgt222 = 
        kmadd(J13L,PDupwindNthSymm1dgt222,kmadd(J23L,PDupwindNthSymm2dgt222,kmul(J33L,PDupwindNthSymm3dgt222)));
      
      JacPDupwindNthSymm3dgt223 = 
        kmadd(J13L,PDupwindNthSymm1dgt223,kmadd(J23L,PDupwindNthSymm2dgt223,kmul(J33L,PDupwindNthSymm3dgt223)));
      
      JacPDupwindNthSymm3dgt233 = 
        kmadd(J13L,PDupwindNthSymm1dgt233,kmadd(J23L,PDupwindNthSymm2dgt233,kmul(J33L,PDupwindNthSymm3dgt233)));
      
      JacPDupwindNthSymm3dgt311 = 
        kmadd(J13L,PDupwindNthSymm1dgt311,kmadd(J23L,PDupwindNthSymm2dgt311,kmul(J33L,PDupwindNthSymm3dgt311)));
      
      JacPDupwindNthSymm3dgt312 = 
        kmadd(J13L,PDupwindNthSymm1dgt312,kmadd(J23L,PDupwindNthSymm2dgt312,kmul(J33L,PDupwindNthSymm3dgt312)));
      
      JacPDupwindNthSymm3dgt313 = 
        kmadd(J13L,PDupwindNthSymm1dgt313,kmadd(J23L,PDupwindNthSymm2dgt313,kmul(J33L,PDupwindNthSymm3dgt313)));
      
      JacPDupwindNthSymm3dgt322 = 
        kmadd(J13L,PDupwindNthSymm1dgt322,kmadd(J23L,PDupwindNthSymm2dgt322,kmul(J33L,PDupwindNthSymm3dgt322)));
      
      JacPDupwindNthSymm3dgt323 = 
        kmadd(J13L,PDupwindNthSymm1dgt323,kmadd(J23L,PDupwindNthSymm2dgt323,kmul(J33L,PDupwindNthSymm3dgt323)));
      
      JacPDupwindNthSymm3dgt333 = 
        kmadd(J13L,PDupwindNthSymm1dgt333,kmadd(J23L,PDupwindNthSymm2dgt333,kmul(J33L,PDupwindNthSymm3dgt333)));
      
      JacPDupwindNthSymm3Xt1 = 
        kmadd(J13L,PDupwindNthSymm1Xt1,kmadd(J23L,PDupwindNthSymm2Xt1,kmul(J33L,PDupwindNthSymm3Xt1)));
      
      JacPDupwindNthSymm3Xt2 = 
        kmadd(J13L,PDupwindNthSymm1Xt2,kmadd(J23L,PDupwindNthSymm2Xt2,kmul(J33L,PDupwindNthSymm3Xt2)));
      
      JacPDupwindNthSymm3Xt3 = 
        kmadd(J13L,PDupwindNthSymm1Xt3,kmadd(J23L,PDupwindNthSymm2Xt3,kmul(J33L,PDupwindNthSymm3Xt3)));
    }
    else
    {
      JacPDstandardNth1At11 = PDstandardNth1At11;
      
      JacPDstandardNth1At12 = PDstandardNth1At12;
      
      JacPDstandardNth1At13 = PDstandardNth1At13;
      
      JacPDstandardNth1At22 = PDstandardNth1At22;
      
      JacPDstandardNth1At23 = PDstandardNth1At23;
      
      JacPDstandardNth1At33 = PDstandardNth1At33;
      
      JacPDstandardNth1dbeta11 = PDstandardNth1dbeta11;
      
      JacPDstandardNth1dbeta12 = PDstandardNth1dbeta12;
      
      JacPDstandardNth1dbeta13 = PDstandardNth1dbeta13;
      
      JacPDstandardNth1dbeta21 = PDstandardNth1dbeta21;
      
      JacPDstandardNth1dbeta22 = PDstandardNth1dbeta22;
      
      JacPDstandardNth1dbeta23 = PDstandardNth1dbeta23;
      
      JacPDstandardNth1dbeta31 = PDstandardNth1dbeta31;
      
      JacPDstandardNth1dbeta32 = PDstandardNth1dbeta32;
      
      JacPDstandardNth1dbeta33 = PDstandardNth1dbeta33;
      
      JacPDstandardNth1gt11 = PDstandardNth1gt11;
      
      JacPDstandardNth1gt12 = PDstandardNth1gt12;
      
      JacPDstandardNth1gt13 = PDstandardNth1gt13;
      
      JacPDstandardNth1gt22 = PDstandardNth1gt22;
      
      JacPDstandardNth1gt23 = PDstandardNth1gt23;
      
      JacPDstandardNth1gt33 = PDstandardNth1gt33;
      
      JacPDstandardNth1trK = PDstandardNth1trK;
      
      JacPDstandardNth2At11 = PDstandardNth2At11;
      
      JacPDstandardNth2At12 = PDstandardNth2At12;
      
      JacPDstandardNth2At13 = PDstandardNth2At13;
      
      JacPDstandardNth2At22 = PDstandardNth2At22;
      
      JacPDstandardNth2At23 = PDstandardNth2At23;
      
      JacPDstandardNth2At33 = PDstandardNth2At33;
      
      JacPDstandardNth2dbeta11 = PDstandardNth2dbeta11;
      
      JacPDstandardNth2dbeta12 = PDstandardNth2dbeta12;
      
      JacPDstandardNth2dbeta13 = PDstandardNth2dbeta13;
      
      JacPDstandardNth2dbeta21 = PDstandardNth2dbeta21;
      
      JacPDstandardNth2dbeta22 = PDstandardNth2dbeta22;
      
      JacPDstandardNth2dbeta23 = PDstandardNth2dbeta23;
      
      JacPDstandardNth2dbeta31 = PDstandardNth2dbeta31;
      
      JacPDstandardNth2dbeta32 = PDstandardNth2dbeta32;
      
      JacPDstandardNth2dbeta33 = PDstandardNth2dbeta33;
      
      JacPDstandardNth2gt11 = PDstandardNth2gt11;
      
      JacPDstandardNth2gt12 = PDstandardNth2gt12;
      
      JacPDstandardNth2gt13 = PDstandardNth2gt13;
      
      JacPDstandardNth2gt22 = PDstandardNth2gt22;
      
      JacPDstandardNth2gt23 = PDstandardNth2gt23;
      
      JacPDstandardNth2gt33 = PDstandardNth2gt33;
      
      JacPDstandardNth2trK = PDstandardNth2trK;
      
      JacPDstandardNth3At11 = PDstandardNth3At11;
      
      JacPDstandardNth3At12 = PDstandardNth3At12;
      
      JacPDstandardNth3At13 = PDstandardNth3At13;
      
      JacPDstandardNth3At22 = PDstandardNth3At22;
      
      JacPDstandardNth3At23 = PDstandardNth3At23;
      
      JacPDstandardNth3At33 = PDstandardNth3At33;
      
      JacPDstandardNth3dbeta11 = PDstandardNth3dbeta11;
      
      JacPDstandardNth3dbeta12 = PDstandardNth3dbeta12;
      
      JacPDstandardNth3dbeta13 = PDstandardNth3dbeta13;
      
      JacPDstandardNth3dbeta21 = PDstandardNth3dbeta21;
      
      JacPDstandardNth3dbeta22 = PDstandardNth3dbeta22;
      
      JacPDstandardNth3dbeta23 = PDstandardNth3dbeta23;
      
      JacPDstandardNth3dbeta31 = PDstandardNth3dbeta31;
      
      JacPDstandardNth3dbeta32 = PDstandardNth3dbeta32;
      
      JacPDstandardNth3dbeta33 = PDstandardNth3dbeta33;
      
      JacPDstandardNth3gt11 = PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = PDstandardNth3gt33;
      
      JacPDstandardNth3trK = PDstandardNth3trK;
      
      JacPDupwindNthAnti1dgt111 = PDupwindNthAnti1dgt111;
      
      JacPDupwindNthAnti1dgt112 = PDupwindNthAnti1dgt112;
      
      JacPDupwindNthAnti1dgt113 = PDupwindNthAnti1dgt113;
      
      JacPDupwindNthAnti1dgt122 = PDupwindNthAnti1dgt122;
      
      JacPDupwindNthAnti1dgt123 = PDupwindNthAnti1dgt123;
      
      JacPDupwindNthAnti1dgt133 = PDupwindNthAnti1dgt133;
      
      JacPDupwindNthAnti1dgt211 = PDupwindNthAnti1dgt211;
      
      JacPDupwindNthAnti1dgt212 = PDupwindNthAnti1dgt212;
      
      JacPDupwindNthAnti1dgt213 = PDupwindNthAnti1dgt213;
      
      JacPDupwindNthAnti1dgt222 = PDupwindNthAnti1dgt222;
      
      JacPDupwindNthAnti1dgt223 = PDupwindNthAnti1dgt223;
      
      JacPDupwindNthAnti1dgt233 = PDupwindNthAnti1dgt233;
      
      JacPDupwindNthAnti1dgt311 = PDupwindNthAnti1dgt311;
      
      JacPDupwindNthAnti1dgt312 = PDupwindNthAnti1dgt312;
      
      JacPDupwindNthAnti1dgt313 = PDupwindNthAnti1dgt313;
      
      JacPDupwindNthAnti1dgt322 = PDupwindNthAnti1dgt322;
      
      JacPDupwindNthAnti1dgt323 = PDupwindNthAnti1dgt323;
      
      JacPDupwindNthAnti1dgt333 = PDupwindNthAnti1dgt333;
      
      JacPDupwindNthAnti1Xt1 = PDupwindNthAnti1Xt1;
      
      JacPDupwindNthAnti1Xt2 = PDupwindNthAnti1Xt2;
      
      JacPDupwindNthAnti1Xt3 = PDupwindNthAnti1Xt3;
      
      JacPDupwindNthSymm1dgt111 = PDupwindNthSymm1dgt111;
      
      JacPDupwindNthSymm1dgt112 = PDupwindNthSymm1dgt112;
      
      JacPDupwindNthSymm1dgt113 = PDupwindNthSymm1dgt113;
      
      JacPDupwindNthSymm1dgt122 = PDupwindNthSymm1dgt122;
      
      JacPDupwindNthSymm1dgt123 = PDupwindNthSymm1dgt123;
      
      JacPDupwindNthSymm1dgt133 = PDupwindNthSymm1dgt133;
      
      JacPDupwindNthSymm1dgt211 = PDupwindNthSymm1dgt211;
      
      JacPDupwindNthSymm1dgt212 = PDupwindNthSymm1dgt212;
      
      JacPDupwindNthSymm1dgt213 = PDupwindNthSymm1dgt213;
      
      JacPDupwindNthSymm1dgt222 = PDupwindNthSymm1dgt222;
      
      JacPDupwindNthSymm1dgt223 = PDupwindNthSymm1dgt223;
      
      JacPDupwindNthSymm1dgt233 = PDupwindNthSymm1dgt233;
      
      JacPDupwindNthSymm1dgt311 = PDupwindNthSymm1dgt311;
      
      JacPDupwindNthSymm1dgt312 = PDupwindNthSymm1dgt312;
      
      JacPDupwindNthSymm1dgt313 = PDupwindNthSymm1dgt313;
      
      JacPDupwindNthSymm1dgt322 = PDupwindNthSymm1dgt322;
      
      JacPDupwindNthSymm1dgt323 = PDupwindNthSymm1dgt323;
      
      JacPDupwindNthSymm1dgt333 = PDupwindNthSymm1dgt333;
      
      JacPDupwindNthSymm1Xt1 = PDupwindNthSymm1Xt1;
      
      JacPDupwindNthSymm1Xt2 = PDupwindNthSymm1Xt2;
      
      JacPDupwindNthSymm1Xt3 = PDupwindNthSymm1Xt3;
      
      JacPDupwindNthAnti2dgt111 = PDupwindNthAnti2dgt111;
      
      JacPDupwindNthAnti2dgt112 = PDupwindNthAnti2dgt112;
      
      JacPDupwindNthAnti2dgt113 = PDupwindNthAnti2dgt113;
      
      JacPDupwindNthAnti2dgt122 = PDupwindNthAnti2dgt122;
      
      JacPDupwindNthAnti2dgt123 = PDupwindNthAnti2dgt123;
      
      JacPDupwindNthAnti2dgt133 = PDupwindNthAnti2dgt133;
      
      JacPDupwindNthAnti2dgt211 = PDupwindNthAnti2dgt211;
      
      JacPDupwindNthAnti2dgt212 = PDupwindNthAnti2dgt212;
      
      JacPDupwindNthAnti2dgt213 = PDupwindNthAnti2dgt213;
      
      JacPDupwindNthAnti2dgt222 = PDupwindNthAnti2dgt222;
      
      JacPDupwindNthAnti2dgt223 = PDupwindNthAnti2dgt223;
      
      JacPDupwindNthAnti2dgt233 = PDupwindNthAnti2dgt233;
      
      JacPDupwindNthAnti2dgt311 = PDupwindNthAnti2dgt311;
      
      JacPDupwindNthAnti2dgt312 = PDupwindNthAnti2dgt312;
      
      JacPDupwindNthAnti2dgt313 = PDupwindNthAnti2dgt313;
      
      JacPDupwindNthAnti2dgt322 = PDupwindNthAnti2dgt322;
      
      JacPDupwindNthAnti2dgt323 = PDupwindNthAnti2dgt323;
      
      JacPDupwindNthAnti2dgt333 = PDupwindNthAnti2dgt333;
      
      JacPDupwindNthAnti2Xt1 = PDupwindNthAnti2Xt1;
      
      JacPDupwindNthAnti2Xt2 = PDupwindNthAnti2Xt2;
      
      JacPDupwindNthAnti2Xt3 = PDupwindNthAnti2Xt3;
      
      JacPDupwindNthSymm2dgt111 = PDupwindNthSymm2dgt111;
      
      JacPDupwindNthSymm2dgt112 = PDupwindNthSymm2dgt112;
      
      JacPDupwindNthSymm2dgt113 = PDupwindNthSymm2dgt113;
      
      JacPDupwindNthSymm2dgt122 = PDupwindNthSymm2dgt122;
      
      JacPDupwindNthSymm2dgt123 = PDupwindNthSymm2dgt123;
      
      JacPDupwindNthSymm2dgt133 = PDupwindNthSymm2dgt133;
      
      JacPDupwindNthSymm2dgt211 = PDupwindNthSymm2dgt211;
      
      JacPDupwindNthSymm2dgt212 = PDupwindNthSymm2dgt212;
      
      JacPDupwindNthSymm2dgt213 = PDupwindNthSymm2dgt213;
      
      JacPDupwindNthSymm2dgt222 = PDupwindNthSymm2dgt222;
      
      JacPDupwindNthSymm2dgt223 = PDupwindNthSymm2dgt223;
      
      JacPDupwindNthSymm2dgt233 = PDupwindNthSymm2dgt233;
      
      JacPDupwindNthSymm2dgt311 = PDupwindNthSymm2dgt311;
      
      JacPDupwindNthSymm2dgt312 = PDupwindNthSymm2dgt312;
      
      JacPDupwindNthSymm2dgt313 = PDupwindNthSymm2dgt313;
      
      JacPDupwindNthSymm2dgt322 = PDupwindNthSymm2dgt322;
      
      JacPDupwindNthSymm2dgt323 = PDupwindNthSymm2dgt323;
      
      JacPDupwindNthSymm2dgt333 = PDupwindNthSymm2dgt333;
      
      JacPDupwindNthSymm2Xt1 = PDupwindNthSymm2Xt1;
      
      JacPDupwindNthSymm2Xt2 = PDupwindNthSymm2Xt2;
      
      JacPDupwindNthSymm2Xt3 = PDupwindNthSymm2Xt3;
      
      JacPDupwindNthAnti3dgt111 = PDupwindNthAnti3dgt111;
      
      JacPDupwindNthAnti3dgt112 = PDupwindNthAnti3dgt112;
      
      JacPDupwindNthAnti3dgt113 = PDupwindNthAnti3dgt113;
      
      JacPDupwindNthAnti3dgt122 = PDupwindNthAnti3dgt122;
      
      JacPDupwindNthAnti3dgt123 = PDupwindNthAnti3dgt123;
      
      JacPDupwindNthAnti3dgt133 = PDupwindNthAnti3dgt133;
      
      JacPDupwindNthAnti3dgt211 = PDupwindNthAnti3dgt211;
      
      JacPDupwindNthAnti3dgt212 = PDupwindNthAnti3dgt212;
      
      JacPDupwindNthAnti3dgt213 = PDupwindNthAnti3dgt213;
      
      JacPDupwindNthAnti3dgt222 = PDupwindNthAnti3dgt222;
      
      JacPDupwindNthAnti3dgt223 = PDupwindNthAnti3dgt223;
      
      JacPDupwindNthAnti3dgt233 = PDupwindNthAnti3dgt233;
      
      JacPDupwindNthAnti3dgt311 = PDupwindNthAnti3dgt311;
      
      JacPDupwindNthAnti3dgt312 = PDupwindNthAnti3dgt312;
      
      JacPDupwindNthAnti3dgt313 = PDupwindNthAnti3dgt313;
      
      JacPDupwindNthAnti3dgt322 = PDupwindNthAnti3dgt322;
      
      JacPDupwindNthAnti3dgt323 = PDupwindNthAnti3dgt323;
      
      JacPDupwindNthAnti3dgt333 = PDupwindNthAnti3dgt333;
      
      JacPDupwindNthAnti3Xt1 = PDupwindNthAnti3Xt1;
      
      JacPDupwindNthAnti3Xt2 = PDupwindNthAnti3Xt2;
      
      JacPDupwindNthAnti3Xt3 = PDupwindNthAnti3Xt3;
      
      JacPDupwindNthSymm3dgt111 = PDupwindNthSymm3dgt111;
      
      JacPDupwindNthSymm3dgt112 = PDupwindNthSymm3dgt112;
      
      JacPDupwindNthSymm3dgt113 = PDupwindNthSymm3dgt113;
      
      JacPDupwindNthSymm3dgt122 = PDupwindNthSymm3dgt122;
      
      JacPDupwindNthSymm3dgt123 = PDupwindNthSymm3dgt123;
      
      JacPDupwindNthSymm3dgt133 = PDupwindNthSymm3dgt133;
      
      JacPDupwindNthSymm3dgt211 = PDupwindNthSymm3dgt211;
      
      JacPDupwindNthSymm3dgt212 = PDupwindNthSymm3dgt212;
      
      JacPDupwindNthSymm3dgt213 = PDupwindNthSymm3dgt213;
      
      JacPDupwindNthSymm3dgt222 = PDupwindNthSymm3dgt222;
      
      JacPDupwindNthSymm3dgt223 = PDupwindNthSymm3dgt223;
      
      JacPDupwindNthSymm3dgt233 = PDupwindNthSymm3dgt233;
      
      JacPDupwindNthSymm3dgt311 = PDupwindNthSymm3dgt311;
      
      JacPDupwindNthSymm3dgt312 = PDupwindNthSymm3dgt312;
      
      JacPDupwindNthSymm3dgt313 = PDupwindNthSymm3dgt313;
      
      JacPDupwindNthSymm3dgt322 = PDupwindNthSymm3dgt322;
      
      JacPDupwindNthSymm3dgt323 = PDupwindNthSymm3dgt323;
      
      JacPDupwindNthSymm3dgt333 = PDupwindNthSymm3dgt333;
      
      JacPDupwindNthSymm3Xt1 = PDupwindNthSymm3Xt1;
      
      JacPDupwindNthSymm3Xt2 = PDupwindNthSymm3Xt2;
      
      JacPDupwindNthSymm3Xt3 = PDupwindNthSymm3Xt3;
    }
    
    CCTK_REAL_VEC detgt CCTK_ATTRIBUTE_UNUSED = ToReal(1);
    
    CCTK_REAL_VEC gtu11 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt22L,gt33L,kmul(gt23L,gt23L)),detgt);
    
    CCTK_REAL_VEC gtu21 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt13L,gt23L,kmul(gt12L,gt33L)),detgt);
    
    CCTK_REAL_VEC gtu31 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt12L,gt23L,kmul(gt13L,gt22L)),detgt);
    
    CCTK_REAL_VEC gtu22 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt11L,gt33L,kmul(gt13L,gt13L)),detgt);
    
    CCTK_REAL_VEC gtu32 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt12L,gt13L,kmul(gt11L,gt23L)),detgt);
    
    CCTK_REAL_VEC gtu33 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt11L,gt22L,kmul(gt12L,gt12L)),detgt);
    
    CCTK_REAL_VEC Gt111 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(dgt111L,gtu11,kmadd(kmsub(ToReal(2),dgt112L,dgt211L),gtu21,kmul(kmsub(ToReal(2),dgt113L,dgt311L),gtu31))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt211 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(dgt111L,gtu21,kmadd(kmsub(ToReal(2),dgt112L,dgt211L),gtu22,kmul(kmsub(ToReal(2),dgt113L,dgt311L),gtu32))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt311 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(dgt111L,gtu31,kmadd(kmsub(ToReal(2),dgt112L,dgt211L),gtu32,kmul(kmsub(ToReal(2),dgt113L,dgt311L),gtu33))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt112 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(dgt211L,gtu11,kmadd(dgt122L,gtu21,kmul(kadd(dgt123L,ksub(dgt213L,dgt312L)),gtu31))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt212 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(dgt211L,gtu21,kmadd(dgt122L,gtu22,kmul(kadd(dgt123L,ksub(dgt213L,dgt312L)),gtu32))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt312 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(dgt211L,gtu31,kmadd(dgt122L,gtu32,kmul(kadd(dgt123L,ksub(dgt213L,dgt312L)),gtu33))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt113 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(dgt311L,gtu11,kmadd(kadd(dgt123L,ksub(dgt312L,dgt213L)),gtu21,kmul(dgt133L,gtu31))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt213 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(dgt311L,gtu21,kmadd(kadd(dgt123L,ksub(dgt312L,dgt213L)),gtu22,kmul(dgt133L,gtu32))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt313 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(dgt311L,gtu31,kmadd(kadd(dgt123L,ksub(dgt312L,dgt213L)),gtu32,kmul(dgt133L,gtu33))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt122 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(kmsub(ToReal(2),dgt212L,dgt122L),gtu11,kmadd(dgt222L,gtu21,kmul(kmsub(ToReal(2),dgt223L,dgt322L),gtu31))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt222 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(kmsub(ToReal(2),dgt212L,dgt122L),gtu21,kmadd(dgt222L,gtu22,kmul(kmsub(ToReal(2),dgt223L,dgt322L),gtu32))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt322 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(kmsub(ToReal(2),dgt212L,dgt122L),gtu31,kmadd(dgt222L,gtu32,kmul(kmsub(ToReal(2),dgt223L,dgt322L),gtu33))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt123 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ksub(kadd(dgt213L,dgt312L),dgt123L),gtu11,kmadd(dgt322L,gtu21,kmul(dgt233L,gtu31))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt223 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ksub(kadd(dgt213L,dgt312L),dgt123L),gtu21,kmadd(dgt322L,gtu22,kmul(dgt233L,gtu32))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt323 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ksub(kadd(dgt213L,dgt312L),dgt123L),gtu31,kmadd(dgt322L,gtu32,kmul(dgt233L,gtu33))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt133 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(kmsub(ToReal(2),dgt313L,dgt133L),gtu11,kmadd(kmsub(ToReal(2),dgt323L,dgt233L),gtu21,kmul(dgt333L,gtu31))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt233 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(kmsub(ToReal(2),dgt313L,dgt133L),gtu21,kmadd(kmsub(ToReal(2),dgt323L,dgt233L),gtu22,kmul(dgt333L,gtu32))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt333 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(kmsub(ToReal(2),dgt313L,dgt133L),gtu31,kmadd(kmsub(ToReal(2),dgt323L,dgt233L),gtu32,kmul(dgt333L,gtu33))),ToReal(0.5));
    
    CCTK_REAL_VEC Xtn1 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt111,gtu11,kmadd(Gt122,gtu22,kmadd(ToReal(2),kmadd(Gt112,gtu21,kmadd(Gt113,gtu31,kmul(Gt123,gtu32))),kmul(Gt133,gtu33))));
    
    CCTK_REAL_VEC Xtn2 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt211,gtu11,kmadd(Gt222,gtu22,kmadd(ToReal(2),kmadd(Gt212,gtu21,kmadd(Gt213,gtu31,kmul(Gt223,gtu32))),kmul(Gt233,gtu33))));
    
    CCTK_REAL_VEC Xtn3 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmadd(ToReal(2),kmadd(Gt312,gtu21,kmadd(Gt313,gtu31,kmul(Gt323,gtu32))),kmul(Gt333,gtu33))));
    
    CCTK_REAL_VEC Atm11 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At11L,gtu11,kmadd(At12L,gtu21,kmul(At13L,gtu31)));
    
    CCTK_REAL_VEC Atm21 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At11L,gtu21,kmadd(At12L,gtu22,kmul(At13L,gtu32)));
    
    CCTK_REAL_VEC Atm31 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At11L,gtu31,kmadd(At12L,gtu32,kmul(At13L,gtu33)));
    
    CCTK_REAL_VEC Atm12 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At12L,gtu11,kmadd(At22L,gtu21,kmul(At23L,gtu31)));
    
    CCTK_REAL_VEC Atm22 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At12L,gtu21,kmadd(At22L,gtu22,kmul(At23L,gtu32)));
    
    CCTK_REAL_VEC Atm32 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At12L,gtu31,kmadd(At22L,gtu32,kmul(At23L,gtu33)));
    
    CCTK_REAL_VEC Atm13 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At13L,gtu11,kmadd(At23L,gtu21,kmul(At33L,gtu31)));
    
    CCTK_REAL_VEC Atm23 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At13L,gtu21,kmadd(At23L,gtu22,kmul(At33L,gtu32)));
    
    CCTK_REAL_VEC Atm33 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At13L,gtu31,kmadd(At23L,gtu32,kmul(At33L,gtu33)));
    
    CCTK_REAL_VEC Atu11 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm11,gtu11,kmadd(Atm12,gtu21,kmul(Atm13,gtu31)));
    
    CCTK_REAL_VEC Atu21 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm11,gtu21,kmadd(Atm12,gtu22,kmul(Atm13,gtu32)));
    
    CCTK_REAL_VEC Atu31 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm11,gtu31,kmadd(Atm12,gtu32,kmul(Atm13,gtu33)));
    
    CCTK_REAL_VEC Atu22 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm21,gtu21,kmadd(Atm22,gtu22,kmul(Atm23,gtu32)));
    
    CCTK_REAL_VEC Atu32 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm21,gtu31,kmadd(Atm22,gtu32,kmul(Atm23,gtu33)));
    
    CCTK_REAL_VEC Atu33 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm31,gtu31,kmadd(Atm32,gtu32,kmul(Atm33,gtu33)));
    
    CCTK_REAL_VEC S1 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(ksub(kmadd(beta1L,eTxxL,kmadd(beta2L,eTxyL,kmul(beta3L,eTxzL))),eTtxL),alphaL);
    
    CCTK_REAL_VEC S2 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(ksub(kmadd(beta1L,eTxyL,kmadd(beta2L,eTyyL,kmul(beta3L,eTyzL))),eTtyL),alphaL);
    
    CCTK_REAL_VEC S3 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(ksub(kmadd(beta1L,eTxzL,kmadd(beta2L,eTyzL,kmul(beta3L,eTzzL))),eTtzL),alphaL);
    
    cdgt111L = ksub(dgt111L,JacPDstandardNth1gt11);
    
    cdgt211L = ksub(dgt211L,JacPDstandardNth2gt11);
    
    cdgt311L = ksub(dgt311L,JacPDstandardNth3gt11);
    
    cdgt112L = ksub(dgt112L,JacPDstandardNth1gt12);
    
    cdgt212L = ksub(dgt212L,JacPDstandardNth2gt12);
    
    cdgt312L = ksub(dgt312L,JacPDstandardNth3gt12);
    
    cdgt113L = ksub(dgt113L,JacPDstandardNth1gt13);
    
    cdgt213L = ksub(dgt213L,JacPDstandardNth2gt13);
    
    cdgt313L = ksub(dgt313L,JacPDstandardNth3gt13);
    
    cdgt122L = ksub(dgt122L,JacPDstandardNth1gt22);
    
    cdgt222L = ksub(dgt222L,JacPDstandardNth2gt22);
    
    cdgt322L = ksub(dgt322L,JacPDstandardNth3gt22);
    
    cdgt123L = ksub(dgt123L,JacPDstandardNth1gt23);
    
    cdgt223L = ksub(dgt223L,JacPDstandardNth2gt23);
    
    cdgt323L = ksub(dgt323L,JacPDstandardNth3gt23);
    
    cdgt133L = ksub(dgt133L,JacPDstandardNth1gt33);
    
    cdgt233L = ksub(dgt233L,JacPDstandardNth2gt33);
    
    cdgt333L = ksub(dgt333L,JacPDstandardNth3gt33);
    
    CCTK_REAL_VEC dgts111 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gt11L,JacPDstandardNth1dbeta11,kmadd(gt12L,JacPDstandardNth1dbeta12,kmul(gt13L,JacPDstandardNth1dbeta13))),ToReal(2));
    
    CCTK_REAL_VEC dgts211 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gt11L,JacPDstandardNth2dbeta11,kmadd(gt12L,JacPDstandardNth2dbeta12,kmul(gt13L,JacPDstandardNth2dbeta13))),ToReal(2));
    
    CCTK_REAL_VEC dgts311 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gt11L,JacPDstandardNth3dbeta11,kmadd(gt12L,JacPDstandardNth3dbeta12,kmul(gt13L,JacPDstandardNth3dbeta13))),ToReal(2));
    
    CCTK_REAL_VEC dgts112 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gt22L,JacPDstandardNth1dbeta12,kmadd(gt23L,JacPDstandardNth1dbeta13,kmadd(gt11L,JacPDstandardNth1dbeta21,kmadd(gt12L,kadd(JacPDstandardNth1dbeta11,JacPDstandardNth1dbeta22),kmul(gt13L,JacPDstandardNth1dbeta23)))));
    
    CCTK_REAL_VEC dgts212 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gt22L,JacPDstandardNth2dbeta12,kmadd(gt23L,JacPDstandardNth2dbeta13,kmadd(gt11L,JacPDstandardNth2dbeta21,kmadd(gt12L,kadd(JacPDstandardNth2dbeta11,JacPDstandardNth2dbeta22),kmul(gt13L,JacPDstandardNth2dbeta23)))));
    
    CCTK_REAL_VEC dgts312 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gt22L,JacPDstandardNth3dbeta12,kmadd(gt23L,JacPDstandardNth3dbeta13,kmadd(gt11L,JacPDstandardNth3dbeta21,kmadd(gt12L,kadd(JacPDstandardNth3dbeta11,JacPDstandardNth3dbeta22),kmul(gt13L,JacPDstandardNth3dbeta23)))));
    
    CCTK_REAL_VEC dgts113 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gt23L,JacPDstandardNth1dbeta12,kmadd(gt33L,JacPDstandardNth1dbeta13,kmadd(gt11L,JacPDstandardNth1dbeta31,kmadd(gt12L,JacPDstandardNth1dbeta32,kmul(gt13L,kadd(JacPDstandardNth1dbeta11,JacPDstandardNth1dbeta33))))));
    
    CCTK_REAL_VEC dgts213 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gt23L,JacPDstandardNth2dbeta12,kmadd(gt33L,JacPDstandardNth2dbeta13,kmadd(gt11L,JacPDstandardNth2dbeta31,kmadd(gt12L,JacPDstandardNth2dbeta32,kmul(gt13L,kadd(JacPDstandardNth2dbeta11,JacPDstandardNth2dbeta33))))));
    
    CCTK_REAL_VEC dgts313 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gt23L,JacPDstandardNth3dbeta12,kmadd(gt33L,JacPDstandardNth3dbeta13,kmadd(gt11L,JacPDstandardNth3dbeta31,kmadd(gt12L,JacPDstandardNth3dbeta32,kmul(gt13L,kadd(JacPDstandardNth3dbeta11,JacPDstandardNth3dbeta33))))));
    
    CCTK_REAL_VEC dgts122 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gt12L,JacPDstandardNth1dbeta21,kmadd(gt22L,JacPDstandardNth1dbeta22,kmul(gt23L,JacPDstandardNth1dbeta23))),ToReal(2));
    
    CCTK_REAL_VEC dgts222 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gt12L,JacPDstandardNth2dbeta21,kmadd(gt22L,JacPDstandardNth2dbeta22,kmul(gt23L,JacPDstandardNth2dbeta23))),ToReal(2));
    
    CCTK_REAL_VEC dgts322 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gt12L,JacPDstandardNth3dbeta21,kmadd(gt22L,JacPDstandardNth3dbeta22,kmul(gt23L,JacPDstandardNth3dbeta23))),ToReal(2));
    
    CCTK_REAL_VEC dgts123 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gt13L,JacPDstandardNth1dbeta21,kmadd(gt33L,JacPDstandardNth1dbeta23,kmadd(gt12L,JacPDstandardNth1dbeta31,kmadd(gt22L,JacPDstandardNth1dbeta32,kmul(gt23L,kadd(JacPDstandardNth1dbeta22,JacPDstandardNth1dbeta33))))));
    
    CCTK_REAL_VEC dgts223 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gt13L,JacPDstandardNth2dbeta21,kmadd(gt33L,JacPDstandardNth2dbeta23,kmadd(gt12L,JacPDstandardNth2dbeta31,kmadd(gt22L,JacPDstandardNth2dbeta32,kmul(gt23L,kadd(JacPDstandardNth2dbeta22,JacPDstandardNth2dbeta33))))));
    
    CCTK_REAL_VEC dgts323 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gt13L,JacPDstandardNth3dbeta21,kmadd(gt33L,JacPDstandardNth3dbeta23,kmadd(gt12L,JacPDstandardNth3dbeta31,kmadd(gt22L,JacPDstandardNth3dbeta32,kmul(gt23L,kadd(JacPDstandardNth3dbeta22,JacPDstandardNth3dbeta33))))));
    
    CCTK_REAL_VEC dgts133 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gt13L,JacPDstandardNth1dbeta31,kmadd(gt23L,JacPDstandardNth1dbeta32,kmul(gt33L,JacPDstandardNth1dbeta33))),ToReal(2));
    
    CCTK_REAL_VEC dgts233 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gt13L,JacPDstandardNth2dbeta31,kmadd(gt23L,JacPDstandardNth2dbeta32,kmul(gt33L,JacPDstandardNth2dbeta33))),ToReal(2));
    
    CCTK_REAL_VEC dgts333 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gt13L,JacPDstandardNth3dbeta31,kmadd(gt23L,JacPDstandardNth3dbeta32,kmul(gt33L,JacPDstandardNth3dbeta33))),ToReal(2));
    
    CCTK_REAL_VEC trdgts1 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dgts111,gtu11,kmadd(dgts122,gtu22,kmadd(ToReal(2),kmadd(dgts112,gtu21,kmadd(dgts113,gtu31,kmul(dgts123,gtu32))),kmul(dgts133,gtu33))));
    
    CCTK_REAL_VEC trdgts2 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dgts211,gtu11,kmadd(dgts222,gtu22,kmadd(ToReal(2),kmadd(dgts212,gtu21,kmadd(dgts213,gtu31,kmul(dgts223,gtu32))),kmul(dgts233,gtu33))));
    
    CCTK_REAL_VEC trdgts3 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dgts311,gtu11,kmadd(dgts322,gtu22,kmadd(ToReal(2),kmadd(dgts312,gtu21,kmadd(dgts313,gtu31,kmul(dgts323,gtu32))),kmul(dgts333,gtu33))));
    
    CCTK_REAL_VEC trcdgt1 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(cdgt111L,gtu11,kmadd(cdgt122L,gtu22,kmadd(ToReal(2),kmadd(cdgt112L,gtu21,kmadd(cdgt113L,gtu31,kmul(cdgt123L,gtu32))),kmul(cdgt133L,gtu33))));
    
    CCTK_REAL_VEC trcdgt2 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(cdgt211L,gtu11,kmadd(cdgt222L,gtu22,kmadd(ToReal(2),kmadd(cdgt212L,gtu21,kmadd(cdgt213L,gtu31,kmul(cdgt223L,gtu32))),kmul(cdgt233L,gtu33))));
    
    CCTK_REAL_VEC trcdgt3 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(cdgt311L,gtu11,kmadd(cdgt322L,gtu22,kmadd(ToReal(2),kmadd(cdgt312L,gtu21,kmadd(cdgt313L,gtu31,kmul(cdgt323L,gtu32))),kmul(cdgt333L,gtu33))));
    
    CCTK_REAL_VEC dgt111rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kmadd(ToReal(3),dbeta11L,kmul(kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),ToReal(-0.666666666666666666666666666667))),dgt111L,kmadd(ToReal(2),kmadd(dbeta12L,dgt112L,kmul(dbeta13L,dgt113L)),kmadd(dbeta12L,dgt211L,kmadd(dbeta13L,dgt311L,kadd(dgts111,kmadd(ToReal(-2),kmadd(At11L,dalpha1L,kmul(alphaL,JacPDstandardNth1At11)),kmadd(beta1L,JacPDupwindNthAnti1dgt111,kmadd(beta2L,JacPDupwindNthAnti2dgt111,kmadd(beta3L,JacPDupwindNthAnti3dgt111,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt11L,trdgts1),kmadd(JacPDupwindNthSymm1dgt111,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt111,kfabs(beta2L),kmsub(JacPDupwindNthSymm3dgt111,kfabs(beta3L),kmul(kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt11L,trcdgt1),cdgt111L),ToReal(DgtDriver)))))))))))))));
    
    CCTK_REAL_VEC dgt211rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dbeta21L,dgt111L,kmadd(kmadd(ToReal(2),dbeta11L,kmadd(ToReal(-0.666666666666666666666666666667),kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),dbeta22L)),dgt211L,kmadd(ToReal(2),kmadd(dbeta12L,dgt212L,kmul(dbeta13L,dgt213L)),kmadd(dbeta23L,dgt311L,kadd(dgts211,kmadd(ToReal(-2),kmadd(At11L,dalpha2L,kmul(alphaL,JacPDstandardNth2At11)),kmadd(beta1L,JacPDupwindNthAnti1dgt211,kmadd(beta2L,JacPDupwindNthAnti2dgt211,kmadd(beta3L,JacPDupwindNthAnti3dgt211,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt11L,trdgts2),kmadd(JacPDupwindNthSymm1dgt211,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt211,kfabs(beta2L),kmsub(JacPDupwindNthSymm3dgt211,kfabs(beta3L),kmul(kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt11L,trcdgt2),cdgt211L),ToReal(DgtDriver)))))))))))))));
    
    CCTK_REAL_VEC dgt311rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dbeta31L,dgt111L,kmadd(dbeta32L,dgt211L,kmadd(kmadd(ToReal(2),dbeta11L,kmadd(ToReal(-0.666666666666666666666666666667),kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),dbeta33L)),dgt311L,kmadd(ToReal(2),kmadd(dbeta12L,dgt312L,kmul(dbeta13L,dgt313L)),kadd(dgts311,kmadd(ToReal(-2),kmadd(At11L,dalpha3L,kmul(alphaL,JacPDstandardNth3At11)),kmadd(beta1L,JacPDupwindNthAnti1dgt311,kmadd(beta2L,JacPDupwindNthAnti2dgt311,kmadd(beta3L,JacPDupwindNthAnti3dgt311,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt11L,trdgts3),kmadd(JacPDupwindNthSymm1dgt311,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt311,kfabs(beta2L),kmsub(JacPDupwindNthSymm3dgt311,kfabs(beta3L),kmul(kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt11L,trcdgt3),cdgt311L),ToReal(DgtDriver)))))))))))))));
    
    CCTK_REAL_VEC dgt112rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dbeta21L,dgt111L,kmadd(kmadd(ToReal(2),dbeta11L,kmadd(ToReal(-0.666666666666666666666666666667),kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),dbeta22L)),dgt112L,kmadd(dbeta23L,dgt113L,kmadd(dbeta12L,kadd(dgt122L,dgt212L),kmadd(dbeta13L,kadd(dgt123L,dgt312L),kadd(dgts112,kmadd(ToReal(-2),kmadd(At12L,dalpha1L,kmul(alphaL,JacPDstandardNth1At12)),kmadd(beta1L,JacPDupwindNthAnti1dgt112,kmadd(beta2L,JacPDupwindNthAnti2dgt112,kmadd(beta3L,JacPDupwindNthAnti3dgt112,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt12L,trdgts1),kmadd(JacPDupwindNthSymm1dgt112,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt112,kfabs(beta2L),kmsub(JacPDupwindNthSymm3dgt112,kfabs(beta3L),kmul(kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt12L,trcdgt1),cdgt112L),ToReal(DgtDriver))))))))))))))));
    
    CCTK_REAL_VEC dgt212rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dbeta21L,kadd(dgt112L,dgt211L),kmadd(kadd(dbeta11L,kmadd(ToReal(2),dbeta22L,kmul(kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),ToReal(-0.666666666666666666666666666667)))),dgt212L,kmadd(dbeta12L,dgt222L,kmadd(dbeta13L,dgt223L,kmadd(dbeta23L,kadd(dgt213L,dgt312L),kadd(dgts212,kmadd(ToReal(-2),kmadd(At12L,dalpha2L,kmul(alphaL,JacPDstandardNth2At12)),kmadd(beta1L,JacPDupwindNthAnti1dgt212,kmadd(beta2L,JacPDupwindNthAnti2dgt212,kmadd(beta3L,JacPDupwindNthAnti3dgt212,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt12L,trdgts2),kmadd(JacPDupwindNthSymm1dgt212,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt212,kfabs(beta2L),kmsub(JacPDupwindNthSymm3dgt212,kfabs(beta3L),kmul(kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt12L,trcdgt2),cdgt212L),ToReal(DgtDriver))))))))))))))));
    
    CCTK_REAL_VEC dgt312rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dbeta31L,dgt112L,kmadd(dbeta32L,dgt212L,kmadd(dbeta21L,dgt311L,kmadd(dbeta23L,dgt313L,kmadd(dbeta12L,dgt322L,kmadd(dbeta13L,dgt323L,kadd(dgts312,kmadd(ToReal(-2),kmadd(At12L,dalpha3L,kmul(alphaL,JacPDstandardNth3At12)),kmadd(beta1L,JacPDupwindNthAnti1dgt312,kmadd(beta2L,JacPDupwindNthAnti2dgt312,kmadd(beta3L,JacPDupwindNthAnti3dgt312,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt12L,trdgts3),kmadd(JacPDupwindNthSymm1dgt312,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt312,kfabs(beta2L),kmadd(JacPDupwindNthSymm3dgt312,kfabs(beta3L),kmsub(ToReal(0.333333333333333333333333333333),kmadd(kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),dgt312L,kmul(kmul(gt12L,trcdgt3),ToReal(DgtDriver))),kmul(cdgt312L,ToReal(DgtDriver))))))))))))))))));
    
    CCTK_REAL_VEC dgt113rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dbeta31L,dgt111L,kmadd(dbeta32L,dgt112L,kmadd(kmadd(ToReal(2),dbeta11L,kmadd(ToReal(-0.666666666666666666666666666667),kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),dbeta33L)),dgt113L,kmadd(dbeta12L,kadd(dgt123L,dgt213L),kmadd(dbeta13L,kadd(dgt133L,dgt313L),kadd(dgts113,kmadd(ToReal(-2),kmadd(At13L,dalpha1L,kmul(alphaL,JacPDstandardNth1At13)),kmadd(beta1L,JacPDupwindNthAnti1dgt113,kmadd(beta2L,JacPDupwindNthAnti2dgt113,kmadd(beta3L,JacPDupwindNthAnti3dgt113,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt13L,trdgts1),kmadd(JacPDupwindNthSymm1dgt113,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt113,kfabs(beta2L),kmsub(JacPDupwindNthSymm3dgt113,kfabs(beta3L),kmul(kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt13L,trcdgt1),cdgt113L),ToReal(DgtDriver))))))))))))))));
    
    CCTK_REAL_VEC dgt213rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dbeta21L,dgt113L,kmadd(dbeta31L,dgt211L,kmadd(dbeta32L,dgt212L,kmadd(dbeta12L,dgt223L,kmadd(dbeta13L,dgt233L,kmadd(dbeta23L,dgt313L,kadd(dgts213,kmadd(ToReal(-2),kmadd(At13L,dalpha2L,kmul(alphaL,JacPDstandardNth2At13)),kmadd(beta1L,JacPDupwindNthAnti1dgt213,kmadd(beta2L,JacPDupwindNthAnti2dgt213,kmadd(beta3L,JacPDupwindNthAnti3dgt213,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt13L,trdgts2),kmadd(JacPDupwindNthSymm1dgt213,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt213,kfabs(beta2L),kmadd(JacPDupwindNthSymm3dgt213,kfabs(beta3L),kmsub(ToReal(0.333333333333333333333333333333),kmadd(kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),dgt213L,kmul(kmul(gt13L,trcdgt2),ToReal(DgtDriver))),kmul(cdgt213L,ToReal(DgtDriver))))))))))))))))));
    
    CCTK_REAL_VEC dgt313rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dbeta31L,kadd(dgt113L,dgt311L),kmadd(dbeta32L,kadd(dgt213L,dgt312L),kmadd(kadd(dbeta11L,kmadd(ToReal(2),dbeta33L,kmul(kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),ToReal(-0.666666666666666666666666666667)))),dgt313L,kmadd(dbeta12L,dgt323L,kmadd(dbeta13L,dgt333L,kadd(dgts313,kmadd(ToReal(-2),kmadd(At13L,dalpha3L,kmul(alphaL,JacPDstandardNth3At13)),kmadd(beta1L,JacPDupwindNthAnti1dgt313,kmadd(beta2L,JacPDupwindNthAnti2dgt313,kmadd(beta3L,JacPDupwindNthAnti3dgt313,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt13L,trdgts3),kmadd(JacPDupwindNthSymm1dgt313,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt313,kfabs(beta2L),kmsub(JacPDupwindNthSymm3dgt313,kfabs(beta3L),kmul(kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt13L,trcdgt3),cdgt313L),ToReal(DgtDriver))))))))))))))));
    
    CCTK_REAL_VEC dgt122rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kmadd(ToReal(-0.666666666666666666666666666667),kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),dbeta11L),dgt122L,kmadd(ToReal(2),kmadd(dbeta21L,dgt112L,kmadd(dbeta22L,dgt122L,kmul(dbeta23L,dgt123L))),kmadd(dbeta12L,dgt222L,kmadd(dbeta13L,dgt322L,kadd(dgts122,kmadd(ToReal(-2),kmadd(At22L,dalpha1L,kmul(alphaL,JacPDstandardNth1At22)),kmadd(beta1L,JacPDupwindNthAnti1dgt122,kmadd(beta2L,JacPDupwindNthAnti2dgt122,kmadd(beta3L,JacPDupwindNthAnti3dgt122,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt22L,trdgts1),kmadd(JacPDupwindNthSymm1dgt122,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt122,kfabs(beta2L),kmsub(JacPDupwindNthSymm3dgt122,kfabs(beta3L),kmul(kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt22L,trcdgt1),cdgt122L),ToReal(DgtDriver)))))))))))))));
    
    CCTK_REAL_VEC dgt222rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dbeta21L,kmadd(ToReal(2),dgt212L,dgt122L),kmadd(kmadd(ToReal(3),dbeta22L,kmul(kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),ToReal(-0.666666666666666666666666666667))),dgt222L,kmadd(dbeta23L,kmadd(ToReal(2),dgt223L,dgt322L),kadd(dgts222,kmadd(ToReal(-2),kmadd(At22L,dalpha2L,kmul(alphaL,JacPDstandardNth2At22)),kmadd(beta1L,JacPDupwindNthAnti1dgt222,kmadd(beta2L,JacPDupwindNthAnti2dgt222,kmadd(beta3L,JacPDupwindNthAnti3dgt222,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt22L,trdgts2),kmadd(JacPDupwindNthSymm1dgt222,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt222,kfabs(beta2L),kmsub(JacPDupwindNthSymm3dgt222,kfabs(beta3L),kmul(kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt22L,trcdgt2),cdgt222L),ToReal(DgtDriver))))))))))))));
    
    CCTK_REAL_VEC dgt322rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dbeta31L,dgt122L,kmadd(dbeta32L,dgt222L,kmadd(kmadd(ToReal(-0.666666666666666666666666666667),kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),dbeta33L),dgt322L,kmadd(ToReal(2),kmadd(dbeta21L,dgt312L,kmadd(dbeta22L,dgt322L,kmul(dbeta23L,dgt323L))),kadd(dgts322,kmadd(ToReal(-2),kmadd(At22L,dalpha3L,kmul(alphaL,JacPDstandardNth3At22)),kmadd(beta1L,JacPDupwindNthAnti1dgt322,kmadd(beta2L,JacPDupwindNthAnti2dgt322,kmadd(beta3L,JacPDupwindNthAnti3dgt322,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt22L,trdgts3),kmadd(JacPDupwindNthSymm1dgt322,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt322,kfabs(beta2L),kmsub(JacPDupwindNthSymm3dgt322,kfabs(beta3L),kmul(kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt22L,trcdgt3),cdgt322L),ToReal(DgtDriver)))))))))))))));
    
    CCTK_REAL_VEC dgt123rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dbeta31L,dgt112L,kmadd(dbeta21L,dgt113L,kmadd(dbeta32L,dgt122L,kmadd(dbeta23L,dgt133L,kmadd(dbeta12L,dgt223L,kmadd(dbeta13L,dgt323L,kadd(dgts123,kmadd(ToReal(-2),kmadd(At23L,dalpha1L,kmul(alphaL,JacPDstandardNth1At23)),kmadd(beta1L,JacPDupwindNthAnti1dgt123,kmadd(beta2L,JacPDupwindNthAnti2dgt123,kmadd(beta3L,JacPDupwindNthAnti3dgt123,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt23L,trdgts1),kmadd(JacPDupwindNthSymm1dgt123,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt123,kfabs(beta2L),kmadd(JacPDupwindNthSymm3dgt123,kfabs(beta3L),kmsub(ToReal(0.333333333333333333333333333333),kmadd(kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),dgt123L,kmul(kmul(gt23L,trcdgt1),ToReal(DgtDriver))),kmul(cdgt123L,ToReal(DgtDriver))))))))))))))))));
    
    CCTK_REAL_VEC dgt223rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dbeta31L,dgt212L,kmadd(dbeta21L,kadd(dgt123L,dgt213L),kmadd(dbeta32L,dgt222L,kmadd(kmadd(ToReal(2),dbeta22L,kmadd(ToReal(-0.666666666666666666666666666667),kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),dbeta33L)),dgt223L,kmadd(dbeta23L,kadd(dgt233L,dgt323L),kadd(dgts223,kmadd(ToReal(-2),kmadd(At23L,dalpha2L,kmul(alphaL,JacPDstandardNth2At23)),kmadd(beta1L,JacPDupwindNthAnti1dgt223,kmadd(beta2L,JacPDupwindNthAnti2dgt223,kmadd(beta3L,JacPDupwindNthAnti3dgt223,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt23L,trdgts2),kmadd(JacPDupwindNthSymm1dgt223,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt223,kfabs(beta2L),kmsub(JacPDupwindNthSymm3dgt223,kfabs(beta3L),kmul(kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt23L,trcdgt2),cdgt223L),ToReal(DgtDriver))))))))))))))));
    
    CCTK_REAL_VEC dgt323rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dbeta31L,kadd(dgt123L,dgt312L),kmadd(dbeta21L,dgt313L,kmadd(dbeta32L,kadd(dgt223L,dgt322L),kmadd(kadd(dbeta22L,kmadd(ToReal(2),dbeta33L,kmul(kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),ToReal(-0.666666666666666666666666666667)))),dgt323L,kmadd(dbeta23L,dgt333L,kadd(dgts323,kmadd(ToReal(-2),kmadd(At23L,dalpha3L,kmul(alphaL,JacPDstandardNth3At23)),kmadd(beta1L,JacPDupwindNthAnti1dgt323,kmadd(beta2L,JacPDupwindNthAnti2dgt323,kmadd(beta3L,JacPDupwindNthAnti3dgt323,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt23L,trdgts3),kmadd(JacPDupwindNthSymm1dgt323,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt323,kfabs(beta2L),kmsub(JacPDupwindNthSymm3dgt323,kfabs(beta3L),kmul(kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt23L,trcdgt3),cdgt323L),ToReal(DgtDriver))))))))))))))));
    
    CCTK_REAL_VEC dgt133rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kmadd(ToReal(-0.666666666666666666666666666667),kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),dbeta11L),dgt133L,kmadd(ToReal(2),kmadd(dbeta31L,dgt113L,kmadd(dbeta32L,dgt123L,kmul(dbeta33L,dgt133L))),kmadd(dbeta12L,dgt233L,kmadd(dbeta13L,dgt333L,kadd(dgts133,kmadd(ToReal(-2),kmadd(At33L,dalpha1L,kmul(alphaL,JacPDstandardNth1At33)),kmadd(beta1L,JacPDupwindNthAnti1dgt133,kmadd(beta2L,JacPDupwindNthAnti2dgt133,kmadd(beta3L,JacPDupwindNthAnti3dgt133,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt33L,trdgts1),kmadd(JacPDupwindNthSymm1dgt133,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt133,kfabs(beta2L),kmsub(JacPDupwindNthSymm3dgt133,kfabs(beta3L),kmul(kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt33L,trcdgt1),cdgt133L),ToReal(DgtDriver)))))))))))))));
    
    CCTK_REAL_VEC dgt233rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dbeta21L,dgt133L,kmadd(kmadd(ToReal(-0.666666666666666666666666666667),kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),dbeta22L),dgt233L,kmadd(ToReal(2),kmadd(dbeta31L,dgt213L,kmadd(dbeta32L,dgt223L,kmul(dbeta33L,dgt233L))),kmadd(dbeta23L,dgt333L,kadd(dgts233,kmadd(ToReal(-2),kmadd(At33L,dalpha2L,kmul(alphaL,JacPDstandardNth2At33)),kmadd(beta1L,JacPDupwindNthAnti1dgt233,kmadd(beta2L,JacPDupwindNthAnti2dgt233,kmadd(beta3L,JacPDupwindNthAnti3dgt233,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt33L,trdgts2),kmadd(JacPDupwindNthSymm1dgt233,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt233,kfabs(beta2L),kmsub(JacPDupwindNthSymm3dgt233,kfabs(beta3L),kmul(kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt33L,trcdgt2),cdgt233L),ToReal(DgtDriver)))))))))))))));
    
    CCTK_REAL_VEC dgt333rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(dbeta31L,kmadd(ToReal(2),dgt313L,dgt133L),kmadd(dbeta32L,kmadd(ToReal(2),dgt323L,dgt233L),kmadd(kmadd(ToReal(3),dbeta33L,kmul(kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),ToReal(-0.666666666666666666666666666667))),dgt333L,kadd(dgts333,kmadd(ToReal(-2),kmadd(At33L,dalpha3L,kmul(alphaL,JacPDstandardNth3At33)),kmadd(beta1L,JacPDupwindNthAnti1dgt333,kmadd(beta2L,JacPDupwindNthAnti2dgt333,kmadd(beta3L,JacPDupwindNthAnti3dgt333,kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt33L,trdgts3),kmadd(JacPDupwindNthSymm1dgt333,kfabs(beta1L),kmadd(JacPDupwindNthSymm2dgt333,kfabs(beta2L),kmsub(JacPDupwindNthSymm3dgt333,kfabs(beta3L),kmul(kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt33L,trcdgt3),cdgt333L),ToReal(DgtDriver))))))))))))));
    
    CCTK_REAL_VEC Xt1rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmadd(dalpha1L,Atu11,kmadd(dalpha2L,Atu21,kmul(dalpha3L,Atu31))),kmadd(gtu11,JacPDstandardNth1dbeta11,kmadd(gtu21,kadd(JacPDstandardNth1dbeta21,JacPDstandardNth2dbeta11),kmadd(gtu22,JacPDstandardNth2dbeta21,kmadd(gtu31,kadd(JacPDstandardNth1dbeta31,JacPDstandardNth3dbeta11),kmadd(gtu32,kadd(JacPDstandardNth2dbeta31,JacPDstandardNth3dbeta21),kmadd(gtu33,JacPDstandardNth3dbeta31,kmadd(ToReal(0.333333333333333333333333333333),kmadd(gtu11,kadd(JacPDstandardNth1dbeta11,kadd(JacPDstandardNth1dbeta22,JacPDstandardNth1dbeta33)),kmadd(gtu21,kadd(JacPDstandardNth2dbeta11,kadd(JacPDstandardNth2dbeta22,JacPDstandardNth2dbeta33)),kmul(gtu31,kadd(JacPDstandardNth3dbeta11,kadd(JacPDstandardNth3dbeta22,JacPDstandardNth3dbeta33))))),kmadd(beta1L,JacPDupwindNthAnti1Xt1,kmadd(beta2L,JacPDupwindNthAnti2Xt1,kmadd(beta3L,JacPDupwindNthAnti3Xt1,kmadd(alphaL,kmadd(dphi3L,Atu31,kmadd(Atu11,kmadd(ToReal(2),Gt111,dphi1L),kmadd(Atu21,kmadd(ToReal(4),Gt112,dphi2L),kmadd(ToReal(4),kmul(Atu31,Gt113),kmadd(ToReal(2),kmul(Atu22,Gt122),kmadd(ToReal(4),kmul(Atu32,Gt123),kmadd(ToReal(2),kmul(Atu33,Gt133),kmadd(ToReal(-1.33333333333333333333333333333),kmadd(gtu11,JacPDstandardNth1trK,kmadd(gtu21,JacPDstandardNth2trK,kmul(gtu31,JacPDstandardNth3trK))),kmul(kmadd(gtu11,S1,kmadd(gtu21,S2,kmul(gtu31,S3))),ToReal(-50.26548245743669181540229413247204614715)))))))))),kmadd(kmsub(ToReal(0.666666666666666666666666666667),kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),dbeta11L),Xtn1,knmsub(dbeta21L,Xtn2,knmsub(dbeta31L,Xtn3,kmadd(JacPDupwindNthSymm1Xt1,kfabs(beta1L),kmadd(JacPDupwindNthSymm2Xt1,kfabs(beta2L),kmadd(JacPDupwindNthSymm3Xt1,kfabs(beta3L),kmul(kmadd(gtu11,kadd(JacPDstandardNth1dbeta22,ksub(JacPDstandardNth1dbeta33,kadd(JacPDstandardNth3dbeta13,JacPDstandardNth2dbeta12))),kmadd(gtu31,ksub(ksub(kadd(JacPDstandardNth3dbeta11,JacPDstandardNth3dbeta22),JacPDstandardNth2dbeta32),JacPDstandardNth1dbeta31),kmul(gtu21,ksub(kadd(JacPDstandardNth2dbeta11,ksub(JacPDstandardNth2dbeta33,JacPDstandardNth3dbeta23)),JacPDstandardNth1dbeta21)))),ToReal(sigma))))))))))))))))))));
    
    CCTK_REAL_VEC Xt2rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmadd(dalpha1L,Atu21,kmadd(dalpha2L,Atu22,kmul(dalpha3L,Atu32))),kmadd(gtu11,JacPDstandardNth1dbeta12,kmadd(gtu21,kadd(JacPDstandardNth1dbeta22,JacPDstandardNth2dbeta12),kmadd(gtu22,JacPDstandardNth2dbeta22,kmadd(gtu31,kadd(JacPDstandardNth1dbeta32,JacPDstandardNth3dbeta12),kmadd(gtu32,kadd(JacPDstandardNth2dbeta32,JacPDstandardNth3dbeta22),kmadd(gtu33,JacPDstandardNth3dbeta32,kmadd(ToReal(0.333333333333333333333333333333),kmadd(gtu21,kadd(JacPDstandardNth1dbeta11,kadd(JacPDstandardNth1dbeta22,JacPDstandardNth1dbeta33)),kmadd(gtu22,kadd(JacPDstandardNth2dbeta11,kadd(JacPDstandardNth2dbeta22,JacPDstandardNth2dbeta33)),kmul(gtu32,kadd(JacPDstandardNth3dbeta11,kadd(JacPDstandardNth3dbeta22,JacPDstandardNth3dbeta33))))),kmadd(beta1L,JacPDupwindNthAnti1Xt2,kmadd(beta2L,JacPDupwindNthAnti2Xt2,kmadd(beta3L,JacPDupwindNthAnti3Xt2,kmadd(alphaL,kmadd(dphi3L,Atu32,kmadd(ToReal(2),kmul(Atu11,Gt211),kmadd(Atu21,kmadd(ToReal(4),Gt212,dphi1L),kmadd(ToReal(4),kmul(Atu31,Gt213),kmadd(Atu22,kmadd(ToReal(2),Gt222,dphi2L),kmadd(ToReal(4),kmul(Atu32,Gt223),kmadd(ToReal(2),kmul(Atu33,Gt233),kmadd(ToReal(-1.33333333333333333333333333333),kmadd(gtu21,JacPDstandardNth1trK,kmadd(gtu22,JacPDstandardNth2trK,kmul(gtu32,JacPDstandardNth3trK))),kmul(kmadd(gtu21,S1,kmadd(gtu22,S2,kmul(gtu32,S3))),ToReal(-50.26548245743669181540229413247204614715)))))))))),knmsub(dbeta12L,Xtn1,kmadd(kmsub(ToReal(0.666666666666666666666666666667),kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),dbeta22L),Xtn2,knmsub(dbeta32L,Xtn3,kmadd(JacPDupwindNthSymm1Xt2,kfabs(beta1L),kmadd(JacPDupwindNthSymm2Xt2,kfabs(beta2L),kmadd(JacPDupwindNthSymm3Xt2,kfabs(beta3L),kmul(kmadd(gtu21,kadd(JacPDstandardNth1dbeta22,ksub(JacPDstandardNth1dbeta33,kadd(JacPDstandardNth3dbeta13,JacPDstandardNth2dbeta12))),kmadd(gtu32,ksub(ksub(kadd(JacPDstandardNth3dbeta11,JacPDstandardNth3dbeta22),JacPDstandardNth2dbeta32),JacPDstandardNth1dbeta31),kmul(gtu22,ksub(kadd(JacPDstandardNth2dbeta11,ksub(JacPDstandardNth2dbeta33,JacPDstandardNth3dbeta23)),JacPDstandardNth1dbeta21)))),ToReal(sigma))))))))))))))))))));
    
    CCTK_REAL_VEC Xt3rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmadd(dalpha1L,Atu31,kmadd(dalpha2L,Atu32,kmul(dalpha3L,Atu33))),kmadd(gtu11,JacPDstandardNth1dbeta13,kmadd(gtu21,kadd(JacPDstandardNth1dbeta23,JacPDstandardNth2dbeta13),kmadd(gtu22,JacPDstandardNth2dbeta23,kmadd(gtu31,kadd(JacPDstandardNth1dbeta33,JacPDstandardNth3dbeta13),kmadd(gtu32,kadd(JacPDstandardNth2dbeta33,JacPDstandardNth3dbeta23),kmadd(gtu33,JacPDstandardNth3dbeta33,kmadd(ToReal(0.333333333333333333333333333333),kmadd(gtu31,kadd(JacPDstandardNth1dbeta11,kadd(JacPDstandardNth1dbeta22,JacPDstandardNth1dbeta33)),kmadd(gtu32,kadd(JacPDstandardNth2dbeta11,kadd(JacPDstandardNth2dbeta22,JacPDstandardNth2dbeta33)),kmul(gtu33,kadd(JacPDstandardNth3dbeta11,kadd(JacPDstandardNth3dbeta22,JacPDstandardNth3dbeta33))))),kmadd(beta1L,JacPDupwindNthAnti1Xt3,kmadd(beta2L,JacPDupwindNthAnti2Xt3,kmadd(beta3L,JacPDupwindNthAnti3Xt3,kmadd(alphaL,kmadd(dphi3L,Atu33,kmadd(ToReal(2),kmul(Atu11,Gt311),kmadd(ToReal(4),kmul(Atu21,Gt312),kmadd(Atu31,kmadd(ToReal(4),Gt313,dphi1L),kmadd(ToReal(2),kmul(Atu22,Gt322),kmadd(Atu32,kmadd(ToReal(4),Gt323,dphi2L),kmadd(ToReal(2),kmul(Atu33,Gt333),kmadd(ToReal(-1.33333333333333333333333333333),kmadd(gtu31,JacPDstandardNth1trK,kmadd(gtu32,JacPDstandardNth2trK,kmul(gtu33,JacPDstandardNth3trK))),kmul(kmadd(gtu31,S1,kmadd(gtu32,S2,kmul(gtu33,S3))),ToReal(-50.26548245743669181540229413247204614715)))))))))),knmsub(dbeta13L,Xtn1,knmsub(dbeta23L,Xtn2,kmadd(kmsub(ToReal(0.666666666666666666666666666667),kadd(dbeta11L,kadd(dbeta22L,dbeta33L)),dbeta33L),Xtn3,kmadd(JacPDupwindNthSymm1Xt3,kfabs(beta1L),kmadd(JacPDupwindNthSymm2Xt3,kfabs(beta2L),kmadd(JacPDupwindNthSymm3Xt3,kfabs(beta3L),kmul(kmadd(gtu31,kadd(JacPDstandardNth1dbeta22,ksub(JacPDstandardNth1dbeta33,kadd(JacPDstandardNth3dbeta13,JacPDstandardNth2dbeta12))),kmadd(gtu33,ksub(ksub(kadd(JacPDstandardNth3dbeta11,JacPDstandardNth3dbeta22),JacPDstandardNth2dbeta32),JacPDstandardNth1dbeta31),kmul(gtu32,ksub(kadd(JacPDstandardNth2dbeta11,ksub(JacPDstandardNth2dbeta33,JacPDstandardNth3dbeta23)),JacPDstandardNth1dbeta21)))),ToReal(sigma))))))))))))))))))));
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(cdgt111[index],cdgt111L);
    vec_store_nta_partial(cdgt112[index],cdgt112L);
    vec_store_nta_partial(cdgt113[index],cdgt113L);
    vec_store_nta_partial(cdgt122[index],cdgt122L);
    vec_store_nta_partial(cdgt123[index],cdgt123L);
    vec_store_nta_partial(cdgt133[index],cdgt133L);
    vec_store_nta_partial(cdgt211[index],cdgt211L);
    vec_store_nta_partial(cdgt212[index],cdgt212L);
    vec_store_nta_partial(cdgt213[index],cdgt213L);
    vec_store_nta_partial(cdgt222[index],cdgt222L);
    vec_store_nta_partial(cdgt223[index],cdgt223L);
    vec_store_nta_partial(cdgt233[index],cdgt233L);
    vec_store_nta_partial(cdgt311[index],cdgt311L);
    vec_store_nta_partial(cdgt312[index],cdgt312L);
    vec_store_nta_partial(cdgt313[index],cdgt313L);
    vec_store_nta_partial(cdgt322[index],cdgt322L);
    vec_store_nta_partial(cdgt323[index],cdgt323L);
    vec_store_nta_partial(cdgt333[index],cdgt333L);
    vec_store_nta_partial(dgt111rhs[index],dgt111rhsL);
    vec_store_nta_partial(dgt112rhs[index],dgt112rhsL);
    vec_store_nta_partial(dgt113rhs[index],dgt113rhsL);
    vec_store_nta_partial(dgt122rhs[index],dgt122rhsL);
    vec_store_nta_partial(dgt123rhs[index],dgt123rhsL);
    vec_store_nta_partial(dgt133rhs[index],dgt133rhsL);
    vec_store_nta_partial(dgt211rhs[index],dgt211rhsL);
    vec_store_nta_partial(dgt212rhs[index],dgt212rhsL);
    vec_store_nta_partial(dgt213rhs[index],dgt213rhsL);
    vec_store_nta_partial(dgt222rhs[index],dgt222rhsL);
    vec_store_nta_partial(dgt223rhs[index],dgt223rhsL);
    vec_store_nta_partial(dgt233rhs[index],dgt233rhsL);
    vec_store_nta_partial(dgt311rhs[index],dgt311rhsL);
    vec_store_nta_partial(dgt312rhs[index],dgt312rhsL);
    vec_store_nta_partial(dgt313rhs[index],dgt313rhsL);
    vec_store_nta_partial(dgt322rhs[index],dgt322rhsL);
    vec_store_nta_partial(dgt323rhs[index],dgt323rhsL);
    vec_store_nta_partial(dgt333rhs[index],dgt333rhsL);
    vec_store_nta_partial(Xt1rhs[index],Xt1rhsL);
    vec_store_nta_partial(Xt2rhs[index],Xt2rhsL);
    vec_store_nta_partial(Xt3rhs[index],Xt3rhsL);
  }
  CCTK_ENDLOOP3STR(CL_BSSN_RHS1);
}
extern "C" void CL_BSSN_RHS1(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CL_BSSN_RHS1
  DECLARE_CCTK_ARGUMENTS_CHECKED(CL_BSSN_RHS1);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering CL_BSSN_RHS1_Body");
  }
  if (cctk_iteration % CL_BSSN_RHS1_calc_every != CL_BSSN_RHS1_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "CL_BSSN::CL_cons_dmetric",
    "CL_BSSN::CL_curv",
    "CL_BSSN::CL_dlapse",
    "CL_BSSN::CL_dlog_confac",
    "CL_BSSN::CL_dmetric",
    "CL_BSSN::CL_dmetricrhs",
    "CL_BSSN::CL_dshift",
    "CL_BSSN::CL_Gamma",
    "CL_BSSN::CL_Gammarhs",
    "CL_BSSN::CL_lapse",
    "CL_BSSN::CL_metric",
    "CL_BSSN::CL_shift",
    "CL_BSSN::CL_trace_curv"};
  AssertGroupStorage(cctkGH, "CL_BSSN_RHS1", 13, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "CL_BSSN_RHS1", 2, 2, 2);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "CL_BSSN_RHS1", 3, 3, 3);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "CL_BSSN_RHS1", 4, 4, 4);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "CL_BSSN_RHS1", 5, 5, 5);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, CL_BSSN_RHS1_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving CL_BSSN_RHS1_Body");
  }
}

} // namespace CL_BSSN
