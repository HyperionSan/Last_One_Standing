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

extern "C" void CL_BSSN_constraints_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CL_BSSN_constraints_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(CL_BSSN_constraints_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % CL_BSSN_constraints_calc_every != CL_BSSN_constraints_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CL_BSSN::CL_cons_detg","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CL_BSSN::CL_cons_detg.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CL_BSSN::CL_cons_Gamma","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CL_BSSN::CL_cons_Gamma.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CL_BSSN::CL_cons_traceA","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CL_BSSN::CL_cons_traceA.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CL_BSSN::CL_Ham","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CL_BSSN::CL_Ham.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CL_BSSN::CL_mom","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CL_BSSN::CL_mom.");
  return;
}

static void CL_BSSN_constraints_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3STR(CL_BSSN_constraints,
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
    CCTK_REAL_VEC phiL CCTK_ATTRIBUTE_UNUSED = vec_load(phi[index]);
    CCTK_REAL_VEC trKL CCTK_ATTRIBUTE_UNUSED = vec_load(trK[index]);
    CCTK_REAL_VEC Xt1L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt1[index]);
    CCTK_REAL_VEC Xt2L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt2[index]);
    CCTK_REAL_VEC Xt3L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt3[index]);
    
    CCTK_REAL_VEC eTttL, eTtxL, eTtyL, eTtzL, eTxxL, eTxyL, eTxzL, eTyyL, eTyzL, eTzzL CCTK_ATTRIBUTE_UNUSED ;
    
    if (assume_stress_energy_state>=0 ? assume_stress_energy_state : *stress_energy_state)
    {
      eTttL = vec_load(eTtt[index]);
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
      eTttL = ToReal(0.);
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
    CCTK_REAL_VEC PDstandardNth1dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dphi1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dphi1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dphi1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dphi2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dphi2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dphi2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1dphi3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2dphi3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3dphi3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3Xt3 CCTK_ATTRIBUTE_UNUSED;
    
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
        PDstandardNth1dgt111 = PDstandardNthfdOrder21(&dgt111[index]);
        PDstandardNth2dgt111 = PDstandardNthfdOrder22(&dgt111[index]);
        PDstandardNth3dgt111 = PDstandardNthfdOrder23(&dgt111[index]);
        PDstandardNth1dgt112 = PDstandardNthfdOrder21(&dgt112[index]);
        PDstandardNth2dgt112 = PDstandardNthfdOrder22(&dgt112[index]);
        PDstandardNth3dgt112 = PDstandardNthfdOrder23(&dgt112[index]);
        PDstandardNth1dgt113 = PDstandardNthfdOrder21(&dgt113[index]);
        PDstandardNth2dgt113 = PDstandardNthfdOrder22(&dgt113[index]);
        PDstandardNth3dgt113 = PDstandardNthfdOrder23(&dgt113[index]);
        PDstandardNth1dgt122 = PDstandardNthfdOrder21(&dgt122[index]);
        PDstandardNth2dgt122 = PDstandardNthfdOrder22(&dgt122[index]);
        PDstandardNth3dgt122 = PDstandardNthfdOrder23(&dgt122[index]);
        PDstandardNth1dgt123 = PDstandardNthfdOrder21(&dgt123[index]);
        PDstandardNth2dgt123 = PDstandardNthfdOrder22(&dgt123[index]);
        PDstandardNth3dgt123 = PDstandardNthfdOrder23(&dgt123[index]);
        PDstandardNth1dgt133 = PDstandardNthfdOrder21(&dgt133[index]);
        PDstandardNth2dgt133 = PDstandardNthfdOrder22(&dgt133[index]);
        PDstandardNth3dgt133 = PDstandardNthfdOrder23(&dgt133[index]);
        PDstandardNth1dgt211 = PDstandardNthfdOrder21(&dgt211[index]);
        PDstandardNth2dgt211 = PDstandardNthfdOrder22(&dgt211[index]);
        PDstandardNth3dgt211 = PDstandardNthfdOrder23(&dgt211[index]);
        PDstandardNth1dgt212 = PDstandardNthfdOrder21(&dgt212[index]);
        PDstandardNth2dgt212 = PDstandardNthfdOrder22(&dgt212[index]);
        PDstandardNth3dgt212 = PDstandardNthfdOrder23(&dgt212[index]);
        PDstandardNth1dgt213 = PDstandardNthfdOrder21(&dgt213[index]);
        PDstandardNth2dgt213 = PDstandardNthfdOrder22(&dgt213[index]);
        PDstandardNth3dgt213 = PDstandardNthfdOrder23(&dgt213[index]);
        PDstandardNth1dgt222 = PDstandardNthfdOrder21(&dgt222[index]);
        PDstandardNth2dgt222 = PDstandardNthfdOrder22(&dgt222[index]);
        PDstandardNth3dgt222 = PDstandardNthfdOrder23(&dgt222[index]);
        PDstandardNth1dgt223 = PDstandardNthfdOrder21(&dgt223[index]);
        PDstandardNth2dgt223 = PDstandardNthfdOrder22(&dgt223[index]);
        PDstandardNth3dgt223 = PDstandardNthfdOrder23(&dgt223[index]);
        PDstandardNth1dgt233 = PDstandardNthfdOrder21(&dgt233[index]);
        PDstandardNth2dgt233 = PDstandardNthfdOrder22(&dgt233[index]);
        PDstandardNth3dgt233 = PDstandardNthfdOrder23(&dgt233[index]);
        PDstandardNth1dgt311 = PDstandardNthfdOrder21(&dgt311[index]);
        PDstandardNth2dgt311 = PDstandardNthfdOrder22(&dgt311[index]);
        PDstandardNth3dgt311 = PDstandardNthfdOrder23(&dgt311[index]);
        PDstandardNth1dgt312 = PDstandardNthfdOrder21(&dgt312[index]);
        PDstandardNth2dgt312 = PDstandardNthfdOrder22(&dgt312[index]);
        PDstandardNth3dgt312 = PDstandardNthfdOrder23(&dgt312[index]);
        PDstandardNth1dgt313 = PDstandardNthfdOrder21(&dgt313[index]);
        PDstandardNth2dgt313 = PDstandardNthfdOrder22(&dgt313[index]);
        PDstandardNth3dgt313 = PDstandardNthfdOrder23(&dgt313[index]);
        PDstandardNth1dgt322 = PDstandardNthfdOrder21(&dgt322[index]);
        PDstandardNth2dgt322 = PDstandardNthfdOrder22(&dgt322[index]);
        PDstandardNth3dgt322 = PDstandardNthfdOrder23(&dgt322[index]);
        PDstandardNth1dgt323 = PDstandardNthfdOrder21(&dgt323[index]);
        PDstandardNth2dgt323 = PDstandardNthfdOrder22(&dgt323[index]);
        PDstandardNth3dgt323 = PDstandardNthfdOrder23(&dgt323[index]);
        PDstandardNth1dgt333 = PDstandardNthfdOrder21(&dgt333[index]);
        PDstandardNth2dgt333 = PDstandardNthfdOrder22(&dgt333[index]);
        PDstandardNth3dgt333 = PDstandardNthfdOrder23(&dgt333[index]);
        PDstandardNth1dphi1 = PDstandardNthfdOrder21(&dphi1[index]);
        PDstandardNth2dphi1 = PDstandardNthfdOrder22(&dphi1[index]);
        PDstandardNth3dphi1 = PDstandardNthfdOrder23(&dphi1[index]);
        PDstandardNth1dphi2 = PDstandardNthfdOrder21(&dphi2[index]);
        PDstandardNth2dphi2 = PDstandardNthfdOrder22(&dphi2[index]);
        PDstandardNth3dphi2 = PDstandardNthfdOrder23(&dphi2[index]);
        PDstandardNth1dphi3 = PDstandardNthfdOrder21(&dphi3[index]);
        PDstandardNth2dphi3 = PDstandardNthfdOrder22(&dphi3[index]);
        PDstandardNth3dphi3 = PDstandardNthfdOrder23(&dphi3[index]);
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
        PDstandardNth1dgt111 = PDstandardNthfdOrder41(&dgt111[index]);
        PDstandardNth2dgt111 = PDstandardNthfdOrder42(&dgt111[index]);
        PDstandardNth3dgt111 = PDstandardNthfdOrder43(&dgt111[index]);
        PDstandardNth1dgt112 = PDstandardNthfdOrder41(&dgt112[index]);
        PDstandardNth2dgt112 = PDstandardNthfdOrder42(&dgt112[index]);
        PDstandardNth3dgt112 = PDstandardNthfdOrder43(&dgt112[index]);
        PDstandardNth1dgt113 = PDstandardNthfdOrder41(&dgt113[index]);
        PDstandardNth2dgt113 = PDstandardNthfdOrder42(&dgt113[index]);
        PDstandardNth3dgt113 = PDstandardNthfdOrder43(&dgt113[index]);
        PDstandardNth1dgt122 = PDstandardNthfdOrder41(&dgt122[index]);
        PDstandardNth2dgt122 = PDstandardNthfdOrder42(&dgt122[index]);
        PDstandardNth3dgt122 = PDstandardNthfdOrder43(&dgt122[index]);
        PDstandardNth1dgt123 = PDstandardNthfdOrder41(&dgt123[index]);
        PDstandardNth2dgt123 = PDstandardNthfdOrder42(&dgt123[index]);
        PDstandardNth3dgt123 = PDstandardNthfdOrder43(&dgt123[index]);
        PDstandardNth1dgt133 = PDstandardNthfdOrder41(&dgt133[index]);
        PDstandardNth2dgt133 = PDstandardNthfdOrder42(&dgt133[index]);
        PDstandardNth3dgt133 = PDstandardNthfdOrder43(&dgt133[index]);
        PDstandardNth1dgt211 = PDstandardNthfdOrder41(&dgt211[index]);
        PDstandardNth2dgt211 = PDstandardNthfdOrder42(&dgt211[index]);
        PDstandardNth3dgt211 = PDstandardNthfdOrder43(&dgt211[index]);
        PDstandardNth1dgt212 = PDstandardNthfdOrder41(&dgt212[index]);
        PDstandardNth2dgt212 = PDstandardNthfdOrder42(&dgt212[index]);
        PDstandardNth3dgt212 = PDstandardNthfdOrder43(&dgt212[index]);
        PDstandardNth1dgt213 = PDstandardNthfdOrder41(&dgt213[index]);
        PDstandardNth2dgt213 = PDstandardNthfdOrder42(&dgt213[index]);
        PDstandardNth3dgt213 = PDstandardNthfdOrder43(&dgt213[index]);
        PDstandardNth1dgt222 = PDstandardNthfdOrder41(&dgt222[index]);
        PDstandardNth2dgt222 = PDstandardNthfdOrder42(&dgt222[index]);
        PDstandardNth3dgt222 = PDstandardNthfdOrder43(&dgt222[index]);
        PDstandardNth1dgt223 = PDstandardNthfdOrder41(&dgt223[index]);
        PDstandardNth2dgt223 = PDstandardNthfdOrder42(&dgt223[index]);
        PDstandardNth3dgt223 = PDstandardNthfdOrder43(&dgt223[index]);
        PDstandardNth1dgt233 = PDstandardNthfdOrder41(&dgt233[index]);
        PDstandardNth2dgt233 = PDstandardNthfdOrder42(&dgt233[index]);
        PDstandardNth3dgt233 = PDstandardNthfdOrder43(&dgt233[index]);
        PDstandardNth1dgt311 = PDstandardNthfdOrder41(&dgt311[index]);
        PDstandardNth2dgt311 = PDstandardNthfdOrder42(&dgt311[index]);
        PDstandardNth3dgt311 = PDstandardNthfdOrder43(&dgt311[index]);
        PDstandardNth1dgt312 = PDstandardNthfdOrder41(&dgt312[index]);
        PDstandardNth2dgt312 = PDstandardNthfdOrder42(&dgt312[index]);
        PDstandardNth3dgt312 = PDstandardNthfdOrder43(&dgt312[index]);
        PDstandardNth1dgt313 = PDstandardNthfdOrder41(&dgt313[index]);
        PDstandardNth2dgt313 = PDstandardNthfdOrder42(&dgt313[index]);
        PDstandardNth3dgt313 = PDstandardNthfdOrder43(&dgt313[index]);
        PDstandardNth1dgt322 = PDstandardNthfdOrder41(&dgt322[index]);
        PDstandardNth2dgt322 = PDstandardNthfdOrder42(&dgt322[index]);
        PDstandardNth3dgt322 = PDstandardNthfdOrder43(&dgt322[index]);
        PDstandardNth1dgt323 = PDstandardNthfdOrder41(&dgt323[index]);
        PDstandardNth2dgt323 = PDstandardNthfdOrder42(&dgt323[index]);
        PDstandardNth3dgt323 = PDstandardNthfdOrder43(&dgt323[index]);
        PDstandardNth1dgt333 = PDstandardNthfdOrder41(&dgt333[index]);
        PDstandardNth2dgt333 = PDstandardNthfdOrder42(&dgt333[index]);
        PDstandardNth3dgt333 = PDstandardNthfdOrder43(&dgt333[index]);
        PDstandardNth1dphi1 = PDstandardNthfdOrder41(&dphi1[index]);
        PDstandardNth2dphi1 = PDstandardNthfdOrder42(&dphi1[index]);
        PDstandardNth3dphi1 = PDstandardNthfdOrder43(&dphi1[index]);
        PDstandardNth1dphi2 = PDstandardNthfdOrder41(&dphi2[index]);
        PDstandardNth2dphi2 = PDstandardNthfdOrder42(&dphi2[index]);
        PDstandardNth3dphi2 = PDstandardNthfdOrder43(&dphi2[index]);
        PDstandardNth1dphi3 = PDstandardNthfdOrder41(&dphi3[index]);
        PDstandardNth2dphi3 = PDstandardNthfdOrder42(&dphi3[index]);
        PDstandardNth3dphi3 = PDstandardNthfdOrder43(&dphi3[index]);
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
        PDstandardNth1dgt111 = PDstandardNthfdOrder61(&dgt111[index]);
        PDstandardNth2dgt111 = PDstandardNthfdOrder62(&dgt111[index]);
        PDstandardNth3dgt111 = PDstandardNthfdOrder63(&dgt111[index]);
        PDstandardNth1dgt112 = PDstandardNthfdOrder61(&dgt112[index]);
        PDstandardNth2dgt112 = PDstandardNthfdOrder62(&dgt112[index]);
        PDstandardNth3dgt112 = PDstandardNthfdOrder63(&dgt112[index]);
        PDstandardNth1dgt113 = PDstandardNthfdOrder61(&dgt113[index]);
        PDstandardNth2dgt113 = PDstandardNthfdOrder62(&dgt113[index]);
        PDstandardNth3dgt113 = PDstandardNthfdOrder63(&dgt113[index]);
        PDstandardNth1dgt122 = PDstandardNthfdOrder61(&dgt122[index]);
        PDstandardNth2dgt122 = PDstandardNthfdOrder62(&dgt122[index]);
        PDstandardNth3dgt122 = PDstandardNthfdOrder63(&dgt122[index]);
        PDstandardNth1dgt123 = PDstandardNthfdOrder61(&dgt123[index]);
        PDstandardNth2dgt123 = PDstandardNthfdOrder62(&dgt123[index]);
        PDstandardNth3dgt123 = PDstandardNthfdOrder63(&dgt123[index]);
        PDstandardNth1dgt133 = PDstandardNthfdOrder61(&dgt133[index]);
        PDstandardNth2dgt133 = PDstandardNthfdOrder62(&dgt133[index]);
        PDstandardNth3dgt133 = PDstandardNthfdOrder63(&dgt133[index]);
        PDstandardNth1dgt211 = PDstandardNthfdOrder61(&dgt211[index]);
        PDstandardNth2dgt211 = PDstandardNthfdOrder62(&dgt211[index]);
        PDstandardNth3dgt211 = PDstandardNthfdOrder63(&dgt211[index]);
        PDstandardNth1dgt212 = PDstandardNthfdOrder61(&dgt212[index]);
        PDstandardNth2dgt212 = PDstandardNthfdOrder62(&dgt212[index]);
        PDstandardNth3dgt212 = PDstandardNthfdOrder63(&dgt212[index]);
        PDstandardNth1dgt213 = PDstandardNthfdOrder61(&dgt213[index]);
        PDstandardNth2dgt213 = PDstandardNthfdOrder62(&dgt213[index]);
        PDstandardNth3dgt213 = PDstandardNthfdOrder63(&dgt213[index]);
        PDstandardNth1dgt222 = PDstandardNthfdOrder61(&dgt222[index]);
        PDstandardNth2dgt222 = PDstandardNthfdOrder62(&dgt222[index]);
        PDstandardNth3dgt222 = PDstandardNthfdOrder63(&dgt222[index]);
        PDstandardNth1dgt223 = PDstandardNthfdOrder61(&dgt223[index]);
        PDstandardNth2dgt223 = PDstandardNthfdOrder62(&dgt223[index]);
        PDstandardNth3dgt223 = PDstandardNthfdOrder63(&dgt223[index]);
        PDstandardNth1dgt233 = PDstandardNthfdOrder61(&dgt233[index]);
        PDstandardNth2dgt233 = PDstandardNthfdOrder62(&dgt233[index]);
        PDstandardNth3dgt233 = PDstandardNthfdOrder63(&dgt233[index]);
        PDstandardNth1dgt311 = PDstandardNthfdOrder61(&dgt311[index]);
        PDstandardNth2dgt311 = PDstandardNthfdOrder62(&dgt311[index]);
        PDstandardNth3dgt311 = PDstandardNthfdOrder63(&dgt311[index]);
        PDstandardNth1dgt312 = PDstandardNthfdOrder61(&dgt312[index]);
        PDstandardNth2dgt312 = PDstandardNthfdOrder62(&dgt312[index]);
        PDstandardNth3dgt312 = PDstandardNthfdOrder63(&dgt312[index]);
        PDstandardNth1dgt313 = PDstandardNthfdOrder61(&dgt313[index]);
        PDstandardNth2dgt313 = PDstandardNthfdOrder62(&dgt313[index]);
        PDstandardNth3dgt313 = PDstandardNthfdOrder63(&dgt313[index]);
        PDstandardNth1dgt322 = PDstandardNthfdOrder61(&dgt322[index]);
        PDstandardNth2dgt322 = PDstandardNthfdOrder62(&dgt322[index]);
        PDstandardNth3dgt322 = PDstandardNthfdOrder63(&dgt322[index]);
        PDstandardNth1dgt323 = PDstandardNthfdOrder61(&dgt323[index]);
        PDstandardNth2dgt323 = PDstandardNthfdOrder62(&dgt323[index]);
        PDstandardNth3dgt323 = PDstandardNthfdOrder63(&dgt323[index]);
        PDstandardNth1dgt333 = PDstandardNthfdOrder61(&dgt333[index]);
        PDstandardNth2dgt333 = PDstandardNthfdOrder62(&dgt333[index]);
        PDstandardNth3dgt333 = PDstandardNthfdOrder63(&dgt333[index]);
        PDstandardNth1dphi1 = PDstandardNthfdOrder61(&dphi1[index]);
        PDstandardNth2dphi1 = PDstandardNthfdOrder62(&dphi1[index]);
        PDstandardNth3dphi1 = PDstandardNthfdOrder63(&dphi1[index]);
        PDstandardNth1dphi2 = PDstandardNthfdOrder61(&dphi2[index]);
        PDstandardNth2dphi2 = PDstandardNthfdOrder62(&dphi2[index]);
        PDstandardNth3dphi2 = PDstandardNthfdOrder63(&dphi2[index]);
        PDstandardNth1dphi3 = PDstandardNthfdOrder61(&dphi3[index]);
        PDstandardNth2dphi3 = PDstandardNthfdOrder62(&dphi3[index]);
        PDstandardNth3dphi3 = PDstandardNthfdOrder63(&dphi3[index]);
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
        PDstandardNth1dgt111 = PDstandardNthfdOrder81(&dgt111[index]);
        PDstandardNth2dgt111 = PDstandardNthfdOrder82(&dgt111[index]);
        PDstandardNth3dgt111 = PDstandardNthfdOrder83(&dgt111[index]);
        PDstandardNth1dgt112 = PDstandardNthfdOrder81(&dgt112[index]);
        PDstandardNth2dgt112 = PDstandardNthfdOrder82(&dgt112[index]);
        PDstandardNth3dgt112 = PDstandardNthfdOrder83(&dgt112[index]);
        PDstandardNth1dgt113 = PDstandardNthfdOrder81(&dgt113[index]);
        PDstandardNth2dgt113 = PDstandardNthfdOrder82(&dgt113[index]);
        PDstandardNth3dgt113 = PDstandardNthfdOrder83(&dgt113[index]);
        PDstandardNth1dgt122 = PDstandardNthfdOrder81(&dgt122[index]);
        PDstandardNth2dgt122 = PDstandardNthfdOrder82(&dgt122[index]);
        PDstandardNth3dgt122 = PDstandardNthfdOrder83(&dgt122[index]);
        PDstandardNth1dgt123 = PDstandardNthfdOrder81(&dgt123[index]);
        PDstandardNth2dgt123 = PDstandardNthfdOrder82(&dgt123[index]);
        PDstandardNth3dgt123 = PDstandardNthfdOrder83(&dgt123[index]);
        PDstandardNth1dgt133 = PDstandardNthfdOrder81(&dgt133[index]);
        PDstandardNth2dgt133 = PDstandardNthfdOrder82(&dgt133[index]);
        PDstandardNth3dgt133 = PDstandardNthfdOrder83(&dgt133[index]);
        PDstandardNth1dgt211 = PDstandardNthfdOrder81(&dgt211[index]);
        PDstandardNth2dgt211 = PDstandardNthfdOrder82(&dgt211[index]);
        PDstandardNth3dgt211 = PDstandardNthfdOrder83(&dgt211[index]);
        PDstandardNth1dgt212 = PDstandardNthfdOrder81(&dgt212[index]);
        PDstandardNth2dgt212 = PDstandardNthfdOrder82(&dgt212[index]);
        PDstandardNth3dgt212 = PDstandardNthfdOrder83(&dgt212[index]);
        PDstandardNth1dgt213 = PDstandardNthfdOrder81(&dgt213[index]);
        PDstandardNth2dgt213 = PDstandardNthfdOrder82(&dgt213[index]);
        PDstandardNth3dgt213 = PDstandardNthfdOrder83(&dgt213[index]);
        PDstandardNth1dgt222 = PDstandardNthfdOrder81(&dgt222[index]);
        PDstandardNth2dgt222 = PDstandardNthfdOrder82(&dgt222[index]);
        PDstandardNth3dgt222 = PDstandardNthfdOrder83(&dgt222[index]);
        PDstandardNth1dgt223 = PDstandardNthfdOrder81(&dgt223[index]);
        PDstandardNth2dgt223 = PDstandardNthfdOrder82(&dgt223[index]);
        PDstandardNth3dgt223 = PDstandardNthfdOrder83(&dgt223[index]);
        PDstandardNth1dgt233 = PDstandardNthfdOrder81(&dgt233[index]);
        PDstandardNth2dgt233 = PDstandardNthfdOrder82(&dgt233[index]);
        PDstandardNth3dgt233 = PDstandardNthfdOrder83(&dgt233[index]);
        PDstandardNth1dgt311 = PDstandardNthfdOrder81(&dgt311[index]);
        PDstandardNth2dgt311 = PDstandardNthfdOrder82(&dgt311[index]);
        PDstandardNth3dgt311 = PDstandardNthfdOrder83(&dgt311[index]);
        PDstandardNth1dgt312 = PDstandardNthfdOrder81(&dgt312[index]);
        PDstandardNth2dgt312 = PDstandardNthfdOrder82(&dgt312[index]);
        PDstandardNth3dgt312 = PDstandardNthfdOrder83(&dgt312[index]);
        PDstandardNth1dgt313 = PDstandardNthfdOrder81(&dgt313[index]);
        PDstandardNth2dgt313 = PDstandardNthfdOrder82(&dgt313[index]);
        PDstandardNth3dgt313 = PDstandardNthfdOrder83(&dgt313[index]);
        PDstandardNth1dgt322 = PDstandardNthfdOrder81(&dgt322[index]);
        PDstandardNth2dgt322 = PDstandardNthfdOrder82(&dgt322[index]);
        PDstandardNth3dgt322 = PDstandardNthfdOrder83(&dgt322[index]);
        PDstandardNth1dgt323 = PDstandardNthfdOrder81(&dgt323[index]);
        PDstandardNth2dgt323 = PDstandardNthfdOrder82(&dgt323[index]);
        PDstandardNth3dgt323 = PDstandardNthfdOrder83(&dgt323[index]);
        PDstandardNth1dgt333 = PDstandardNthfdOrder81(&dgt333[index]);
        PDstandardNth2dgt333 = PDstandardNthfdOrder82(&dgt333[index]);
        PDstandardNth3dgt333 = PDstandardNthfdOrder83(&dgt333[index]);
        PDstandardNth1dphi1 = PDstandardNthfdOrder81(&dphi1[index]);
        PDstandardNth2dphi1 = PDstandardNthfdOrder82(&dphi1[index]);
        PDstandardNth3dphi1 = PDstandardNthfdOrder83(&dphi1[index]);
        PDstandardNth1dphi2 = PDstandardNthfdOrder81(&dphi2[index]);
        PDstandardNth2dphi2 = PDstandardNthfdOrder82(&dphi2[index]);
        PDstandardNth3dphi2 = PDstandardNthfdOrder83(&dphi2[index]);
        PDstandardNth1dphi3 = PDstandardNthfdOrder81(&dphi3[index]);
        PDstandardNth2dphi3 = PDstandardNthfdOrder82(&dphi3[index]);
        PDstandardNth3dphi3 = PDstandardNthfdOrder83(&dphi3[index]);
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
    CCTK_REAL_VEC JacPDstandardNth1At11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1At12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1At13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1At22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1At23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1At33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dphi1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dphi2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1dphi3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2At11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2At12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2At13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2At22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2At23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2At33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dphi1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dphi2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2dphi3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3At11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3At12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3At13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3At22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3At23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3At33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt111 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt112 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt113 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt122 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt123 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt133 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt211 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt212 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt213 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt222 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt223 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt233 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt311 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt312 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt313 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt322 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt323 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dgt333 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dphi1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dphi2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3dphi3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3Xt3 CCTK_ATTRIBUTE_UNUSED;
    
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
      
      JacPDstandardNth1dgt111 = 
        kmadd(J11L,PDstandardNth1dgt111,kmadd(J21L,PDstandardNth2dgt111,kmul(J31L,PDstandardNth3dgt111)));
      
      JacPDstandardNth1dgt112 = 
        kmadd(J11L,PDstandardNth1dgt112,kmadd(J21L,PDstandardNth2dgt112,kmul(J31L,PDstandardNth3dgt112)));
      
      JacPDstandardNth1dgt113 = 
        kmadd(J11L,PDstandardNth1dgt113,kmadd(J21L,PDstandardNth2dgt113,kmul(J31L,PDstandardNth3dgt113)));
      
      JacPDstandardNth1dgt122 = 
        kmadd(J11L,PDstandardNth1dgt122,kmadd(J21L,PDstandardNth2dgt122,kmul(J31L,PDstandardNth3dgt122)));
      
      JacPDstandardNth1dgt123 = 
        kmadd(J11L,PDstandardNth1dgt123,kmadd(J21L,PDstandardNth2dgt123,kmul(J31L,PDstandardNth3dgt123)));
      
      JacPDstandardNth1dgt133 = 
        kmadd(J11L,PDstandardNth1dgt133,kmadd(J21L,PDstandardNth2dgt133,kmul(J31L,PDstandardNth3dgt133)));
      
      JacPDstandardNth1dgt211 = 
        kmadd(J11L,PDstandardNth1dgt211,kmadd(J21L,PDstandardNth2dgt211,kmul(J31L,PDstandardNth3dgt211)));
      
      JacPDstandardNth1dgt212 = 
        kmadd(J11L,PDstandardNth1dgt212,kmadd(J21L,PDstandardNth2dgt212,kmul(J31L,PDstandardNth3dgt212)));
      
      JacPDstandardNth1dgt213 = 
        kmadd(J11L,PDstandardNth1dgt213,kmadd(J21L,PDstandardNth2dgt213,kmul(J31L,PDstandardNth3dgt213)));
      
      JacPDstandardNth1dgt222 = 
        kmadd(J11L,PDstandardNth1dgt222,kmadd(J21L,PDstandardNth2dgt222,kmul(J31L,PDstandardNth3dgt222)));
      
      JacPDstandardNth1dgt223 = 
        kmadd(J11L,PDstandardNth1dgt223,kmadd(J21L,PDstandardNth2dgt223,kmul(J31L,PDstandardNth3dgt223)));
      
      JacPDstandardNth1dgt233 = 
        kmadd(J11L,PDstandardNth1dgt233,kmadd(J21L,PDstandardNth2dgt233,kmul(J31L,PDstandardNth3dgt233)));
      
      JacPDstandardNth1dgt311 = 
        kmadd(J11L,PDstandardNth1dgt311,kmadd(J21L,PDstandardNth2dgt311,kmul(J31L,PDstandardNth3dgt311)));
      
      JacPDstandardNth1dgt312 = 
        kmadd(J11L,PDstandardNth1dgt312,kmadd(J21L,PDstandardNth2dgt312,kmul(J31L,PDstandardNth3dgt312)));
      
      JacPDstandardNth1dgt313 = 
        kmadd(J11L,PDstandardNth1dgt313,kmadd(J21L,PDstandardNth2dgt313,kmul(J31L,PDstandardNth3dgt313)));
      
      JacPDstandardNth1dgt322 = 
        kmadd(J11L,PDstandardNth1dgt322,kmadd(J21L,PDstandardNth2dgt322,kmul(J31L,PDstandardNth3dgt322)));
      
      JacPDstandardNth1dgt323 = 
        kmadd(J11L,PDstandardNth1dgt323,kmadd(J21L,PDstandardNth2dgt323,kmul(J31L,PDstandardNth3dgt323)));
      
      JacPDstandardNth1dgt333 = 
        kmadd(J11L,PDstandardNth1dgt333,kmadd(J21L,PDstandardNth2dgt333,kmul(J31L,PDstandardNth3dgt333)));
      
      JacPDstandardNth1dphi1 = 
        kmadd(J11L,PDstandardNth1dphi1,kmadd(J21L,PDstandardNth2dphi1,kmul(J31L,PDstandardNth3dphi1)));
      
      JacPDstandardNth1dphi2 = 
        kmadd(J11L,PDstandardNth1dphi2,kmadd(J21L,PDstandardNth2dphi2,kmul(J31L,PDstandardNth3dphi2)));
      
      JacPDstandardNth1dphi3 = 
        kmadd(J11L,PDstandardNth1dphi3,kmadd(J21L,PDstandardNth2dphi3,kmul(J31L,PDstandardNth3dphi3)));
      
      JacPDstandardNth1trK = 
        kmadd(J11L,PDstandardNth1trK,kmadd(J21L,PDstandardNth2trK,kmul(J31L,PDstandardNth3trK)));
      
      JacPDstandardNth1Xt1 = 
        kmadd(J11L,PDstandardNth1Xt1,kmadd(J21L,PDstandardNth2Xt1,kmul(J31L,PDstandardNth3Xt1)));
      
      JacPDstandardNth1Xt2 = 
        kmadd(J11L,PDstandardNth1Xt2,kmadd(J21L,PDstandardNth2Xt2,kmul(J31L,PDstandardNth3Xt2)));
      
      JacPDstandardNth1Xt3 = 
        kmadd(J11L,PDstandardNth1Xt3,kmadd(J21L,PDstandardNth2Xt3,kmul(J31L,PDstandardNth3Xt3)));
      
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
      
      JacPDstandardNth2dgt111 = 
        kmadd(J12L,PDstandardNth1dgt111,kmadd(J22L,PDstandardNth2dgt111,kmul(J32L,PDstandardNth3dgt111)));
      
      JacPDstandardNth2dgt112 = 
        kmadd(J12L,PDstandardNth1dgt112,kmadd(J22L,PDstandardNth2dgt112,kmul(J32L,PDstandardNth3dgt112)));
      
      JacPDstandardNth2dgt113 = 
        kmadd(J12L,PDstandardNth1dgt113,kmadd(J22L,PDstandardNth2dgt113,kmul(J32L,PDstandardNth3dgt113)));
      
      JacPDstandardNth2dgt122 = 
        kmadd(J12L,PDstandardNth1dgt122,kmadd(J22L,PDstandardNth2dgt122,kmul(J32L,PDstandardNth3dgt122)));
      
      JacPDstandardNth2dgt123 = 
        kmadd(J12L,PDstandardNth1dgt123,kmadd(J22L,PDstandardNth2dgt123,kmul(J32L,PDstandardNth3dgt123)));
      
      JacPDstandardNth2dgt133 = 
        kmadd(J12L,PDstandardNth1dgt133,kmadd(J22L,PDstandardNth2dgt133,kmul(J32L,PDstandardNth3dgt133)));
      
      JacPDstandardNth2dgt211 = 
        kmadd(J12L,PDstandardNth1dgt211,kmadd(J22L,PDstandardNth2dgt211,kmul(J32L,PDstandardNth3dgt211)));
      
      JacPDstandardNth2dgt212 = 
        kmadd(J12L,PDstandardNth1dgt212,kmadd(J22L,PDstandardNth2dgt212,kmul(J32L,PDstandardNth3dgt212)));
      
      JacPDstandardNth2dgt213 = 
        kmadd(J12L,PDstandardNth1dgt213,kmadd(J22L,PDstandardNth2dgt213,kmul(J32L,PDstandardNth3dgt213)));
      
      JacPDstandardNth2dgt222 = 
        kmadd(J12L,PDstandardNth1dgt222,kmadd(J22L,PDstandardNth2dgt222,kmul(J32L,PDstandardNth3dgt222)));
      
      JacPDstandardNth2dgt223 = 
        kmadd(J12L,PDstandardNth1dgt223,kmadd(J22L,PDstandardNth2dgt223,kmul(J32L,PDstandardNth3dgt223)));
      
      JacPDstandardNth2dgt233 = 
        kmadd(J12L,PDstandardNth1dgt233,kmadd(J22L,PDstandardNth2dgt233,kmul(J32L,PDstandardNth3dgt233)));
      
      JacPDstandardNth2dgt311 = 
        kmadd(J12L,PDstandardNth1dgt311,kmadd(J22L,PDstandardNth2dgt311,kmul(J32L,PDstandardNth3dgt311)));
      
      JacPDstandardNth2dgt312 = 
        kmadd(J12L,PDstandardNth1dgt312,kmadd(J22L,PDstandardNth2dgt312,kmul(J32L,PDstandardNth3dgt312)));
      
      JacPDstandardNth2dgt313 = 
        kmadd(J12L,PDstandardNth1dgt313,kmadd(J22L,PDstandardNth2dgt313,kmul(J32L,PDstandardNth3dgt313)));
      
      JacPDstandardNth2dgt322 = 
        kmadd(J12L,PDstandardNth1dgt322,kmadd(J22L,PDstandardNth2dgt322,kmul(J32L,PDstandardNth3dgt322)));
      
      JacPDstandardNth2dgt323 = 
        kmadd(J12L,PDstandardNth1dgt323,kmadd(J22L,PDstandardNth2dgt323,kmul(J32L,PDstandardNth3dgt323)));
      
      JacPDstandardNth2dgt333 = 
        kmadd(J12L,PDstandardNth1dgt333,kmadd(J22L,PDstandardNth2dgt333,kmul(J32L,PDstandardNth3dgt333)));
      
      JacPDstandardNth2dphi1 = 
        kmadd(J12L,PDstandardNth1dphi1,kmadd(J22L,PDstandardNth2dphi1,kmul(J32L,PDstandardNth3dphi1)));
      
      JacPDstandardNth2dphi2 = 
        kmadd(J12L,PDstandardNth1dphi2,kmadd(J22L,PDstandardNth2dphi2,kmul(J32L,PDstandardNth3dphi2)));
      
      JacPDstandardNth2dphi3 = 
        kmadd(J12L,PDstandardNth1dphi3,kmadd(J22L,PDstandardNth2dphi3,kmul(J32L,PDstandardNth3dphi3)));
      
      JacPDstandardNth2trK = 
        kmadd(J12L,PDstandardNth1trK,kmadd(J22L,PDstandardNth2trK,kmul(J32L,PDstandardNth3trK)));
      
      JacPDstandardNth2Xt1 = 
        kmadd(J12L,PDstandardNth1Xt1,kmadd(J22L,PDstandardNth2Xt1,kmul(J32L,PDstandardNth3Xt1)));
      
      JacPDstandardNth2Xt2 = 
        kmadd(J12L,PDstandardNth1Xt2,kmadd(J22L,PDstandardNth2Xt2,kmul(J32L,PDstandardNth3Xt2)));
      
      JacPDstandardNth2Xt3 = 
        kmadd(J12L,PDstandardNth1Xt3,kmadd(J22L,PDstandardNth2Xt3,kmul(J32L,PDstandardNth3Xt3)));
      
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
      
      JacPDstandardNth3dgt111 = 
        kmadd(J13L,PDstandardNth1dgt111,kmadd(J23L,PDstandardNth2dgt111,kmul(J33L,PDstandardNth3dgt111)));
      
      JacPDstandardNth3dgt112 = 
        kmadd(J13L,PDstandardNth1dgt112,kmadd(J23L,PDstandardNth2dgt112,kmul(J33L,PDstandardNth3dgt112)));
      
      JacPDstandardNth3dgt113 = 
        kmadd(J13L,PDstandardNth1dgt113,kmadd(J23L,PDstandardNth2dgt113,kmul(J33L,PDstandardNth3dgt113)));
      
      JacPDstandardNth3dgt122 = 
        kmadd(J13L,PDstandardNth1dgt122,kmadd(J23L,PDstandardNth2dgt122,kmul(J33L,PDstandardNth3dgt122)));
      
      JacPDstandardNth3dgt123 = 
        kmadd(J13L,PDstandardNth1dgt123,kmadd(J23L,PDstandardNth2dgt123,kmul(J33L,PDstandardNth3dgt123)));
      
      JacPDstandardNth3dgt133 = 
        kmadd(J13L,PDstandardNth1dgt133,kmadd(J23L,PDstandardNth2dgt133,kmul(J33L,PDstandardNth3dgt133)));
      
      JacPDstandardNth3dgt211 = 
        kmadd(J13L,PDstandardNth1dgt211,kmadd(J23L,PDstandardNth2dgt211,kmul(J33L,PDstandardNth3dgt211)));
      
      JacPDstandardNth3dgt212 = 
        kmadd(J13L,PDstandardNth1dgt212,kmadd(J23L,PDstandardNth2dgt212,kmul(J33L,PDstandardNth3dgt212)));
      
      JacPDstandardNth3dgt213 = 
        kmadd(J13L,PDstandardNth1dgt213,kmadd(J23L,PDstandardNth2dgt213,kmul(J33L,PDstandardNth3dgt213)));
      
      JacPDstandardNth3dgt222 = 
        kmadd(J13L,PDstandardNth1dgt222,kmadd(J23L,PDstandardNth2dgt222,kmul(J33L,PDstandardNth3dgt222)));
      
      JacPDstandardNth3dgt223 = 
        kmadd(J13L,PDstandardNth1dgt223,kmadd(J23L,PDstandardNth2dgt223,kmul(J33L,PDstandardNth3dgt223)));
      
      JacPDstandardNth3dgt233 = 
        kmadd(J13L,PDstandardNth1dgt233,kmadd(J23L,PDstandardNth2dgt233,kmul(J33L,PDstandardNth3dgt233)));
      
      JacPDstandardNth3dgt311 = 
        kmadd(J13L,PDstandardNth1dgt311,kmadd(J23L,PDstandardNth2dgt311,kmul(J33L,PDstandardNth3dgt311)));
      
      JacPDstandardNth3dgt312 = 
        kmadd(J13L,PDstandardNth1dgt312,kmadd(J23L,PDstandardNth2dgt312,kmul(J33L,PDstandardNth3dgt312)));
      
      JacPDstandardNth3dgt313 = 
        kmadd(J13L,PDstandardNth1dgt313,kmadd(J23L,PDstandardNth2dgt313,kmul(J33L,PDstandardNth3dgt313)));
      
      JacPDstandardNth3dgt322 = 
        kmadd(J13L,PDstandardNth1dgt322,kmadd(J23L,PDstandardNth2dgt322,kmul(J33L,PDstandardNth3dgt322)));
      
      JacPDstandardNth3dgt323 = 
        kmadd(J13L,PDstandardNth1dgt323,kmadd(J23L,PDstandardNth2dgt323,kmul(J33L,PDstandardNth3dgt323)));
      
      JacPDstandardNth3dgt333 = 
        kmadd(J13L,PDstandardNth1dgt333,kmadd(J23L,PDstandardNth2dgt333,kmul(J33L,PDstandardNth3dgt333)));
      
      JacPDstandardNth3dphi1 = 
        kmadd(J13L,PDstandardNth1dphi1,kmadd(J23L,PDstandardNth2dphi1,kmul(J33L,PDstandardNth3dphi1)));
      
      JacPDstandardNth3dphi2 = 
        kmadd(J13L,PDstandardNth1dphi2,kmadd(J23L,PDstandardNth2dphi2,kmul(J33L,PDstandardNth3dphi2)));
      
      JacPDstandardNth3dphi3 = 
        kmadd(J13L,PDstandardNth1dphi3,kmadd(J23L,PDstandardNth2dphi3,kmul(J33L,PDstandardNth3dphi3)));
      
      JacPDstandardNth3trK = 
        kmadd(J13L,PDstandardNth1trK,kmadd(J23L,PDstandardNth2trK,kmul(J33L,PDstandardNth3trK)));
      
      JacPDstandardNth3Xt1 = 
        kmadd(J13L,PDstandardNth1Xt1,kmadd(J23L,PDstandardNth2Xt1,kmul(J33L,PDstandardNth3Xt1)));
      
      JacPDstandardNth3Xt2 = 
        kmadd(J13L,PDstandardNth1Xt2,kmadd(J23L,PDstandardNth2Xt2,kmul(J33L,PDstandardNth3Xt2)));
      
      JacPDstandardNth3Xt3 = 
        kmadd(J13L,PDstandardNth1Xt3,kmadd(J23L,PDstandardNth2Xt3,kmul(J33L,PDstandardNth3Xt3)));
    }
    else
    {
      JacPDstandardNth1At11 = PDstandardNth1At11;
      
      JacPDstandardNth1At12 = PDstandardNth1At12;
      
      JacPDstandardNth1At13 = PDstandardNth1At13;
      
      JacPDstandardNth1At22 = PDstandardNth1At22;
      
      JacPDstandardNth1At23 = PDstandardNth1At23;
      
      JacPDstandardNth1At33 = PDstandardNth1At33;
      
      JacPDstandardNth1dgt111 = PDstandardNth1dgt111;
      
      JacPDstandardNth1dgt112 = PDstandardNth1dgt112;
      
      JacPDstandardNth1dgt113 = PDstandardNth1dgt113;
      
      JacPDstandardNth1dgt122 = PDstandardNth1dgt122;
      
      JacPDstandardNth1dgt123 = PDstandardNth1dgt123;
      
      JacPDstandardNth1dgt133 = PDstandardNth1dgt133;
      
      JacPDstandardNth1dgt211 = PDstandardNth1dgt211;
      
      JacPDstandardNth1dgt212 = PDstandardNth1dgt212;
      
      JacPDstandardNth1dgt213 = PDstandardNth1dgt213;
      
      JacPDstandardNth1dgt222 = PDstandardNth1dgt222;
      
      JacPDstandardNth1dgt223 = PDstandardNth1dgt223;
      
      JacPDstandardNth1dgt233 = PDstandardNth1dgt233;
      
      JacPDstandardNth1dgt311 = PDstandardNth1dgt311;
      
      JacPDstandardNth1dgt312 = PDstandardNth1dgt312;
      
      JacPDstandardNth1dgt313 = PDstandardNth1dgt313;
      
      JacPDstandardNth1dgt322 = PDstandardNth1dgt322;
      
      JacPDstandardNth1dgt323 = PDstandardNth1dgt323;
      
      JacPDstandardNth1dgt333 = PDstandardNth1dgt333;
      
      JacPDstandardNth1dphi1 = PDstandardNth1dphi1;
      
      JacPDstandardNth1dphi2 = PDstandardNth1dphi2;
      
      JacPDstandardNth1dphi3 = PDstandardNth1dphi3;
      
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
      
      JacPDstandardNth2dgt111 = PDstandardNth2dgt111;
      
      JacPDstandardNth2dgt112 = PDstandardNth2dgt112;
      
      JacPDstandardNth2dgt113 = PDstandardNth2dgt113;
      
      JacPDstandardNth2dgt122 = PDstandardNth2dgt122;
      
      JacPDstandardNth2dgt123 = PDstandardNth2dgt123;
      
      JacPDstandardNth2dgt133 = PDstandardNth2dgt133;
      
      JacPDstandardNth2dgt211 = PDstandardNth2dgt211;
      
      JacPDstandardNth2dgt212 = PDstandardNth2dgt212;
      
      JacPDstandardNth2dgt213 = PDstandardNth2dgt213;
      
      JacPDstandardNth2dgt222 = PDstandardNth2dgt222;
      
      JacPDstandardNth2dgt223 = PDstandardNth2dgt223;
      
      JacPDstandardNth2dgt233 = PDstandardNth2dgt233;
      
      JacPDstandardNth2dgt311 = PDstandardNth2dgt311;
      
      JacPDstandardNth2dgt312 = PDstandardNth2dgt312;
      
      JacPDstandardNth2dgt313 = PDstandardNth2dgt313;
      
      JacPDstandardNth2dgt322 = PDstandardNth2dgt322;
      
      JacPDstandardNth2dgt323 = PDstandardNth2dgt323;
      
      JacPDstandardNth2dgt333 = PDstandardNth2dgt333;
      
      JacPDstandardNth2dphi1 = PDstandardNth2dphi1;
      
      JacPDstandardNth2dphi2 = PDstandardNth2dphi2;
      
      JacPDstandardNth2dphi3 = PDstandardNth2dphi3;
      
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
      
      JacPDstandardNth3dgt111 = PDstandardNth3dgt111;
      
      JacPDstandardNth3dgt112 = PDstandardNth3dgt112;
      
      JacPDstandardNth3dgt113 = PDstandardNth3dgt113;
      
      JacPDstandardNth3dgt122 = PDstandardNth3dgt122;
      
      JacPDstandardNth3dgt123 = PDstandardNth3dgt123;
      
      JacPDstandardNth3dgt133 = PDstandardNth3dgt133;
      
      JacPDstandardNth3dgt211 = PDstandardNth3dgt211;
      
      JacPDstandardNth3dgt212 = PDstandardNth3dgt212;
      
      JacPDstandardNth3dgt213 = PDstandardNth3dgt213;
      
      JacPDstandardNth3dgt222 = PDstandardNth3dgt222;
      
      JacPDstandardNth3dgt223 = PDstandardNth3dgt223;
      
      JacPDstandardNth3dgt233 = PDstandardNth3dgt233;
      
      JacPDstandardNth3dgt311 = PDstandardNth3dgt311;
      
      JacPDstandardNth3dgt312 = PDstandardNth3dgt312;
      
      JacPDstandardNth3dgt313 = PDstandardNth3dgt313;
      
      JacPDstandardNth3dgt322 = PDstandardNth3dgt322;
      
      JacPDstandardNth3dgt323 = PDstandardNth3dgt323;
      
      JacPDstandardNth3dgt333 = PDstandardNth3dgt333;
      
      JacPDstandardNth3dphi1 = PDstandardNth3dphi1;
      
      JacPDstandardNth3dphi2 = PDstandardNth3dphi2;
      
      JacPDstandardNth3dphi3 = PDstandardNth3dphi3;
      
      JacPDstandardNth3trK = PDstandardNth3trK;
      
      JacPDstandardNth3Xt1 = PDstandardNth3Xt1;
      
      JacPDstandardNth3Xt2 = PDstandardNth3Xt2;
      
      JacPDstandardNth3Xt3 = PDstandardNth3Xt3;
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
    
    CCTK_REAL_VEC Gtl111 CCTK_ATTRIBUTE_UNUSED = 
      kmul(dgt111L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl112 CCTK_ATTRIBUTE_UNUSED = 
      kmul(dgt211L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl113 CCTK_ATTRIBUTE_UNUSED = 
      kmul(dgt311L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl122 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-0.5),dgt122L,dgt212L);
    
    CCTK_REAL_VEC Gtl123 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ksub(kadd(dgt213L,dgt312L),dgt123L),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl133 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-0.5),dgt133L,dgt313L);
    
    CCTK_REAL_VEC Gtl211 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-0.5),dgt211L,dgt112L);
    
    CCTK_REAL_VEC Gtl212 CCTK_ATTRIBUTE_UNUSED = 
      kmul(dgt122L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl213 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kadd(dgt123L,ksub(dgt312L,dgt213L)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl222 CCTK_ATTRIBUTE_UNUSED = 
      kmul(dgt222L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl223 CCTK_ATTRIBUTE_UNUSED = 
      kmul(dgt322L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl233 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-0.5),dgt233L,dgt323L);
    
    CCTK_REAL_VEC Gtl311 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-0.5),dgt311L,dgt113L);
    
    CCTK_REAL_VEC Gtl312 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kadd(dgt123L,ksub(dgt213L,dgt312L)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl313 CCTK_ATTRIBUTE_UNUSED = 
      kmul(dgt133L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl322 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-0.5),dgt322L,dgt223L);
    
    CCTK_REAL_VEC Gtl323 CCTK_ATTRIBUTE_UNUSED = 
      kmul(dgt233L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl333 CCTK_ATTRIBUTE_UNUSED = 
      kmul(dgt333L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtlu111 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu11,kmadd(Gtl112,gtu21,kmul(Gtl113,gtu31)));
    
    CCTK_REAL_VEC Gtlu112 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu21,kmadd(Gtl112,gtu22,kmul(Gtl113,gtu32)));
    
    CCTK_REAL_VEC Gtlu113 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu31,kmadd(Gtl112,gtu32,kmul(Gtl113,gtu33)));
    
    CCTK_REAL_VEC Gtlu121 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu11,kmadd(Gtl122,gtu21,kmul(Gtl123,gtu31)));
    
    CCTK_REAL_VEC Gtlu122 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu21,kmadd(Gtl122,gtu22,kmul(Gtl123,gtu32)));
    
    CCTK_REAL_VEC Gtlu123 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu31,kmadd(Gtl122,gtu32,kmul(Gtl123,gtu33)));
    
    CCTK_REAL_VEC Gtlu131 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu11,kmadd(Gtl123,gtu21,kmul(Gtl133,gtu31)));
    
    CCTK_REAL_VEC Gtlu132 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu21,kmadd(Gtl123,gtu22,kmul(Gtl133,gtu32)));
    
    CCTK_REAL_VEC Gtlu133 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu31,kmadd(Gtl123,gtu32,kmul(Gtl133,gtu33)));
    
    CCTK_REAL_VEC Gtlu211 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl211,gtu11,kmadd(Gtl212,gtu21,kmul(Gtl213,gtu31)));
    
    CCTK_REAL_VEC Gtlu212 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl211,gtu21,kmadd(Gtl212,gtu22,kmul(Gtl213,gtu32)));
    
    CCTK_REAL_VEC Gtlu213 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl211,gtu31,kmadd(Gtl212,gtu32,kmul(Gtl213,gtu33)));
    
    CCTK_REAL_VEC Gtlu221 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl212,gtu11,kmadd(Gtl222,gtu21,kmul(Gtl223,gtu31)));
    
    CCTK_REAL_VEC Gtlu222 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl212,gtu21,kmadd(Gtl222,gtu22,kmul(Gtl223,gtu32)));
    
    CCTK_REAL_VEC Gtlu223 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl212,gtu31,kmadd(Gtl222,gtu32,kmul(Gtl223,gtu33)));
    
    CCTK_REAL_VEC Gtlu231 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl213,gtu11,kmadd(Gtl223,gtu21,kmul(Gtl233,gtu31)));
    
    CCTK_REAL_VEC Gtlu232 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl213,gtu21,kmadd(Gtl223,gtu22,kmul(Gtl233,gtu32)));
    
    CCTK_REAL_VEC Gtlu233 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl213,gtu31,kmadd(Gtl223,gtu32,kmul(Gtl233,gtu33)));
    
    CCTK_REAL_VEC Gtlu311 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl311,gtu11,kmadd(Gtl312,gtu21,kmul(Gtl313,gtu31)));
    
    CCTK_REAL_VEC Gtlu312 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl311,gtu21,kmadd(Gtl312,gtu22,kmul(Gtl313,gtu32)));
    
    CCTK_REAL_VEC Gtlu313 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl311,gtu31,kmadd(Gtl312,gtu32,kmul(Gtl313,gtu33)));
    
    CCTK_REAL_VEC Gtlu321 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl312,gtu11,kmadd(Gtl322,gtu21,kmul(Gtl323,gtu31)));
    
    CCTK_REAL_VEC Gtlu322 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl312,gtu21,kmadd(Gtl322,gtu22,kmul(Gtl323,gtu32)));
    
    CCTK_REAL_VEC Gtlu323 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl312,gtu31,kmadd(Gtl322,gtu32,kmul(Gtl323,gtu33)));
    
    CCTK_REAL_VEC Gtlu331 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl313,gtu11,kmadd(Gtl323,gtu21,kmul(Gtl333,gtu31)));
    
    CCTK_REAL_VEC Gtlu332 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl313,gtu21,kmadd(Gtl323,gtu22,kmul(Gtl333,gtu32)));
    
    CCTK_REAL_VEC Gtlu333 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl313,gtu31,kmadd(Gtl323,gtu32,kmul(Gtl333,gtu33)));
    
    CCTK_REAL_VEC Gt111 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu11,kmadd(Gtl211,gtu21,kmul(Gtl311,gtu31)));
    
    CCTK_REAL_VEC Gt211 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu21,kmadd(Gtl211,gtu22,kmul(Gtl311,gtu32)));
    
    CCTK_REAL_VEC Gt311 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu31,kmadd(Gtl211,gtu32,kmul(Gtl311,gtu33)));
    
    CCTK_REAL_VEC Gt112 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu11,kmadd(Gtl212,gtu21,kmul(Gtl312,gtu31)));
    
    CCTK_REAL_VEC Gt212 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu21,kmadd(Gtl212,gtu22,kmul(Gtl312,gtu32)));
    
    CCTK_REAL_VEC Gt312 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu31,kmadd(Gtl212,gtu32,kmul(Gtl312,gtu33)));
    
    CCTK_REAL_VEC Gt113 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu11,kmadd(Gtl213,gtu21,kmul(Gtl313,gtu31)));
    
    CCTK_REAL_VEC Gt213 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu21,kmadd(Gtl213,gtu22,kmul(Gtl313,gtu32)));
    
    CCTK_REAL_VEC Gt313 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu31,kmadd(Gtl213,gtu32,kmul(Gtl313,gtu33)));
    
    CCTK_REAL_VEC Gt122 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl122,gtu11,kmadd(Gtl222,gtu21,kmul(Gtl322,gtu31)));
    
    CCTK_REAL_VEC Gt222 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl122,gtu21,kmadd(Gtl222,gtu22,kmul(Gtl322,gtu32)));
    
    CCTK_REAL_VEC Gt322 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl122,gtu31,kmadd(Gtl222,gtu32,kmul(Gtl322,gtu33)));
    
    CCTK_REAL_VEC Gt123 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl123,gtu11,kmadd(Gtl223,gtu21,kmul(Gtl323,gtu31)));
    
    CCTK_REAL_VEC Gt223 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl123,gtu21,kmadd(Gtl223,gtu22,kmul(Gtl323,gtu32)));
    
    CCTK_REAL_VEC Gt323 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl123,gtu31,kmadd(Gtl223,gtu32,kmul(Gtl323,gtu33)));
    
    CCTK_REAL_VEC Gt133 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl133,gtu11,kmadd(Gtl233,gtu21,kmul(Gtl333,gtu31)));
    
    CCTK_REAL_VEC Gt233 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl133,gtu21,kmadd(Gtl233,gtu22,kmul(Gtl333,gtu32)));
    
    CCTK_REAL_VEC Gt333 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl133,gtu31,kmadd(Gtl233,gtu32,kmul(Gtl333,gtu33)));
    
    CCTK_REAL_VEC Xtn1 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt111,gtu11,kmadd(Gt122,gtu22,kmadd(ToReal(2),kmadd(Gt112,gtu21,kmadd(Gt113,gtu31,kmul(Gt123,gtu32))),kmul(Gt133,gtu33))));
    
    CCTK_REAL_VEC Xtn2 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt211,gtu11,kmadd(Gt222,gtu22,kmadd(ToReal(2),kmadd(Gt212,gtu21,kmadd(Gt213,gtu31,kmul(Gt223,gtu32))),kmul(Gt233,gtu33))));
    
    CCTK_REAL_VEC Xtn3 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmadd(ToReal(2),kmadd(Gt312,gtu21,kmadd(Gt313,gtu31,kmul(Gt323,gtu32))),kmul(Gt333,gtu33))));
    
    CCTK_REAL_VEC Rt11 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(3),kmadd(Gt111,Gtlu111,kmadd(Gt112,Gtlu112,kmul(Gt113,Gtlu113))),kmadd(ToReal(2),kmadd(Gt211,Gtlu121,kmadd(Gt212,Gtlu122,kmadd(Gt213,Gtlu123,kmadd(Gt311,Gtlu131,kmadd(Gt312,Gtlu132,kmul(Gt313,Gtlu133)))))),kmadd(Gt211,Gtlu211,kmadd(Gt212,Gtlu212,kmadd(Gt213,Gtlu213,kmadd(Gt311,Gtlu311,kmadd(Gt312,Gtlu312,kmadd(Gt313,Gtlu313,kmadd(gt11L,JacPDstandardNth1Xt1,kmadd(gt12L,JacPDstandardNth1Xt2,kmadd(gt13L,JacPDstandardNth1Xt3,knmsub(ToReal(0.5),kmadd(gtu11,JacPDstandardNth1dgt111,kmadd(gtu21,kadd(JacPDstandardNth1dgt211,JacPDstandardNth2dgt111),kmadd(gtu22,JacPDstandardNth2dgt211,kmadd(gtu31,kadd(JacPDstandardNth1dgt311,JacPDstandardNth3dgt111),kmadd(gtu32,kadd(JacPDstandardNth2dgt311,JacPDstandardNth3dgt211),kmul(gtu33,JacPDstandardNth3dgt311)))))),kmadd(Gtl111,Xtn1,kmadd(Gtl112,Xtn2,kmul(Gtl113,Xtn3)))))))))))))));
    
    CCTK_REAL_VEC Rt12 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(4),kmadd(Gt211,Gtlu221,kmadd(Gt212,Gtlu222,kmul(Gt213,Gtlu223))),kmadd(ToReal(2),kmadd(Gt122,Gtlu112,kmadd(Gt123,Gtlu113,kmadd(Gt111,Gtlu121,kmadd(Gt212,Gtlu121,kmadd(Gt222,Gtlu122,kmadd(Gt113,Gtlu123,kmadd(Gt223,Gtlu123,kmadd(Gt312,Gtlu131,kmadd(Gt322,Gtlu132,kmadd(Gt323,Gtlu133,kmadd(Gt111,Gtlu211,kmadd(Gt112,kadd(Gtlu111,kadd(Gtlu122,Gtlu212)),kmadd(Gt113,Gtlu213,kmadd(Gt311,Gtlu231,kmadd(Gt312,Gtlu232,kmadd(Gt313,Gtlu233,kmadd(Gt311,Gtlu321,kmadd(Gt312,Gtlu322,kmul(Gt313,Gtlu323))))))))))))))))))),knmsub(gtu11,JacPDstandardNth1dgt112,kmadd(gt12L,JacPDstandardNth1Xt1,kmadd(gt22L,JacPDstandardNth1Xt2,kmadd(gt23L,JacPDstandardNth1Xt3,knmsub(gtu21,kadd(JacPDstandardNth2dgt112,JacPDstandardNth1dgt212),knmsub(gtu22,JacPDstandardNth2dgt212,kmadd(gt11L,JacPDstandardNth2Xt1,kmadd(gt12L,JacPDstandardNth2Xt2,kmadd(gt13L,JacPDstandardNth2Xt3,knmsub(gtu31,kadd(JacPDstandardNth3dgt112,JacPDstandardNth1dgt312),knmsub(gtu32,kadd(JacPDstandardNth3dgt212,JacPDstandardNth2dgt312),knmsub(gtu33,JacPDstandardNth3dgt312,kmadd(Gtl112,Xtn1,kmadd(Gtl211,Xtn1,kmadd(Gtl122,Xtn2,kmadd(Gtl212,Xtn2,kmadd(Gtl123,Xtn3,kmul(Gtl213,Xtn3)))))))))))))))))))),ToReal(0.5));
    
    CCTK_REAL_VEC Rt13 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(2),kmadd(Gt123,Gtlu112,kmadd(Gt133,Gtlu113,kmadd(Gt213,Gtlu121,kmadd(Gt223,Gtlu122,kmadd(Gt233,Gtlu123,kmadd(Gt111,Gtlu131,kmadd(Gt313,Gtlu131,kmadd(Gt112,Gtlu132,kmadd(Gt323,Gtlu132,kmadd(Gt333,Gtlu133,kmadd(Gt211,Gtlu231,kmadd(Gt212,Gtlu232,kmadd(Gt213,Gtlu233,kmadd(Gt111,Gtlu311,kmadd(Gt112,Gtlu312,kmadd(Gt113,kadd(Gtlu111,kadd(Gtlu133,Gtlu313)),kmadd(Gt211,Gtlu321,kmadd(Gt212,Gtlu322,kmul(Gt213,Gtlu323))))))))))))))))))),kmadd(ToReal(4),kmadd(Gt311,Gtlu331,kmadd(Gt312,Gtlu332,kmul(Gt313,Gtlu333))),knmsub(gtu11,JacPDstandardNth1dgt113,kmadd(gt13L,JacPDstandardNth1Xt1,kmadd(gt23L,JacPDstandardNth1Xt2,kmadd(gt33L,JacPDstandardNth1Xt3,knmsub(gtu21,kadd(JacPDstandardNth2dgt113,JacPDstandardNth1dgt213),knmsub(gtu22,JacPDstandardNth2dgt213,knmsub(gtu31,kadd(JacPDstandardNth3dgt113,JacPDstandardNth1dgt313),knmsub(gtu32,kadd(JacPDstandardNth3dgt213,JacPDstandardNth2dgt313),knmsub(gtu33,JacPDstandardNth3dgt313,kmadd(gt11L,JacPDstandardNth3Xt1,kmadd(gt12L,JacPDstandardNth3Xt2,kmadd(gt13L,JacPDstandardNth3Xt3,kmadd(Gtl113,Xtn1,kmadd(Gtl311,Xtn1,kmadd(Gtl123,Xtn2,kmadd(Gtl312,Xtn2,kmadd(Gtl133,Xtn3,kmul(Gtl313,Xtn3)))))))))))))))))))),ToReal(0.5));
    
    CCTK_REAL_VEC Rt22 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt112,kmadd(ToReal(2),Gtlu211,Gtlu121),kmadd(Gt122,kmadd(ToReal(2),Gtlu212,Gtlu122),kmadd(Gt123,kmadd(ToReal(2),Gtlu213,Gtlu123),kmadd(ToReal(3),kmadd(Gt212,Gtlu221,kmadd(Gt222,Gtlu222,kmul(Gt223,Gtlu223))),kmadd(ToReal(2),kmadd(Gt312,Gtlu231,kmadd(Gt322,Gtlu232,kmul(Gt323,Gtlu233))),kmadd(Gt312,Gtlu321,kmadd(Gt322,Gtlu322,kmadd(Gt323,Gtlu323,kmadd(gt12L,JacPDstandardNth2Xt1,kmadd(gt22L,JacPDstandardNth2Xt2,kmadd(gt23L,JacPDstandardNth2Xt3,knmsub(ToReal(0.5),kmadd(gtu11,JacPDstandardNth1dgt122,kmadd(gtu21,kadd(JacPDstandardNth1dgt222,JacPDstandardNth2dgt122),kmadd(gtu22,JacPDstandardNth2dgt222,kmadd(gtu31,kadd(JacPDstandardNth1dgt322,JacPDstandardNth3dgt122),kmadd(gtu32,kadd(JacPDstandardNth2dgt322,JacPDstandardNth3dgt222),kmul(gtu33,JacPDstandardNth3dgt322)))))),kmadd(Gtl212,Xtn1,kmadd(Gtl222,Xtn2,kmul(Gtl223,Xtn3)))))))))))))));
    
    CCTK_REAL_VEC Rt23 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(2),kmadd(Gt123,Gtlu133,kmadd(Gt113,Gtlu211,kmadd(Gt123,Gtlu212,kmadd(Gt133,Gtlu213,kmadd(Gt213,Gtlu221,kmadd(Gt223,Gtlu222,kmadd(Gt233,Gtlu223,kmadd(Gt212,Gtlu231,kmadd(Gt313,Gtlu231,kmadd(Gt222,Gtlu232,kmadd(Gt323,Gtlu232,kmadd(Gt223,Gtlu233,kmadd(Gt333,Gtlu233,kmadd(Gt112,kadd(Gtlu131,Gtlu311),kmadd(Gt122,kadd(Gtlu132,Gtlu312),kmadd(Gt123,Gtlu313,kmadd(Gt212,Gtlu321,kmadd(Gt222,Gtlu322,kmul(Gt223,Gtlu323))))))))))))))))))),kmadd(ToReal(4),kmadd(Gt312,Gtlu331,kmadd(Gt322,Gtlu332,kmul(Gt323,Gtlu333))),knmsub(gtu11,JacPDstandardNth1dgt123,knmsub(gtu21,kadd(JacPDstandardNth2dgt123,JacPDstandardNth1dgt223),knmsub(gtu22,JacPDstandardNth2dgt223,kmadd(gt13L,JacPDstandardNth2Xt1,kmadd(gt23L,JacPDstandardNth2Xt2,kmadd(gt33L,JacPDstandardNth2Xt3,knmsub(gtu31,kadd(JacPDstandardNth3dgt123,JacPDstandardNth1dgt323),knmsub(gtu32,kadd(JacPDstandardNth3dgt223,JacPDstandardNth2dgt323),knmsub(gtu33,JacPDstandardNth3dgt323,kmadd(gt12L,JacPDstandardNth3Xt1,kmadd(gt22L,JacPDstandardNth3Xt2,kmadd(gt23L,JacPDstandardNth3Xt3,kmadd(Gtl213,Xtn1,kmadd(Gtl312,Xtn1,kmadd(Gtl223,Xtn2,kmadd(Gtl322,Xtn2,kmadd(Gtl233,Xtn3,kmul(Gtl323,Xtn3)))))))))))))))))))),ToReal(0.5));
    
    CCTK_REAL_VEC Rt33 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt113,kmadd(ToReal(2),Gtlu311,Gtlu131),kmadd(Gt123,kmadd(ToReal(2),Gtlu312,Gtlu132),kmadd(Gt133,kmadd(ToReal(2),Gtlu313,Gtlu133),kmadd(Gt213,kmadd(ToReal(2),Gtlu321,Gtlu231),kmadd(Gt223,kmadd(ToReal(2),Gtlu322,Gtlu232),kmadd(Gt233,kmadd(ToReal(2),Gtlu323,Gtlu233),kmadd(ToReal(3),kmadd(Gt313,Gtlu331,kmadd(Gt323,Gtlu332,kmul(Gt333,Gtlu333))),knmsub(ToReal(0.5),kmadd(gtu11,JacPDstandardNth1dgt133,kmadd(gtu21,kadd(JacPDstandardNth1dgt233,JacPDstandardNth2dgt133),kmadd(gtu22,JacPDstandardNth2dgt233,kmadd(gtu31,kadd(JacPDstandardNth1dgt333,JacPDstandardNth3dgt133),kmadd(gtu32,kadd(JacPDstandardNth2dgt333,JacPDstandardNth3dgt233),kmul(gtu33,JacPDstandardNth3dgt333)))))),kmadd(gt13L,JacPDstandardNth3Xt1,kmadd(gt23L,JacPDstandardNth3Xt2,kmadd(gt33L,JacPDstandardNth3Xt3,kmadd(Gtl313,Xtn1,kmadd(Gtl323,Xtn2,kmul(Gtl333,Xtn3))))))))))))));
    
    CCTK_REAL_VEC Rphi11 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(6),kmadd(kmadd(dphi1L,Gt111,kmul(dphi2L,Gt211)),kmadd(gt11L,gtu11,ToReal(1)),kmul(dphi3L,kmadd(gt11L,kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmul(Gt333,gtu33))),Gt311))),kmadd(ToReal(-6),kmadd(gt11L,kmadd(gtu21,kadd(JacPDstandardNth1dphi2,JacPDstandardNth2dphi1),kmul(gtu32,kadd(JacPDstandardNth2dphi3,JacPDstandardNth3dphi2))),JacPDstandardNth1dphi1),kmadd(knmsub(gt11L,gtu11,ToReal(1)),kmul(dphi1L,dphi1L),kmul(gt11L,kmadd(ToReal(12),kmadd(dphi3L,kmul(Gt312,gtu21),kmul(dphi3L,kmul(Gt313,gtu31))),kmadd(ToReal(12),kmul(dphi3L,kmul(Gt323,gtu32)),kmadd(ToReal(2),kmul(dphi1L,kmadd(kmsub(ToReal(6),Gt112,dphi2L),gtu21,kmadd(ToReal(3),kmul(Gt122,gtu22),kmadd(kmsub(ToReal(6),Gt113,dphi3L),gtu31,kmadd(ToReal(6),kmul(Gt123,gtu32),kmul(ToReal(3),kmul(Gt133,gtu33))))))),kmadd(ToReal(2),kmul(dphi2L,kmadd(ToReal(6),kmul(Gt212,gtu21),kmadd(ToReal(3),kmul(Gt222,gtu22),kmadd(ToReal(6),kmul(Gt213,gtu31),kmadd(kmsub(ToReal(6),Gt223,dphi3L),gtu32,kmul(ToReal(3),kmul(Gt233,gtu33))))))),kmadd(ToReal(-6),kmul(gtu11,JacPDstandardNth1dphi1),kmadd(ToReal(-6),kmul(gtu31,JacPDstandardNth1dphi3),kmadd(ToReal(-6),kmul(gtu31,JacPDstandardNth3dphi1),kmadd(gtu22,kmsub(ToReal(-6),JacPDstandardNth2dphi2,kmul(dphi2L,dphi2L)),kmul(gtu33,kmsub(ToReal(-6),JacPDstandardNth3dphi3,kmul(dphi3L,dphi3L))))))))))))))),ToReal(0.0277777777777777777777777777778));
    
    CCTK_REAL_VEC Rphi12 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(6),kmul(dphi3L,kmadd(gt12L,kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmul(Gt333,gtu33))),Gt312)),kmadd(ToReal(2),kmul(dphi2L,kmadd(Gt212,kmadd(ToReal(6),kmul(gt12L,gtu21),ToReal(3)),kmul(gt12L,knmsub(dphi3L,gtu32,kmadd(ToReal(6),kmadd(Gt213,gtu31,kmul(Gt223,gtu32)),kmul(ToReal(3),kmadd(Gt211,gtu11,kmadd(Gt222,gtu22,kmul(Gt233,gtu33))))))))),kmadd(dphi1L,kadd(dphi2L,kmadd(ToReal(6),kmul(Gt112,kmadd(ToReal(2),kmul(gt12L,gtu21),ToReal(1))),kmul(gt12L,kmadd(ToReal(-2),kmul(dphi2L,gtu21),kmul(ToReal(2),kmadd(kmsub(ToReal(6),Gt113,dphi3L),gtu31,kmadd(ToReal(6),kmul(Gt123,gtu32),kmul(ToReal(3),kmadd(Gt111,gtu11,kmadd(Gt122,gtu22,kmul(Gt133,gtu33))))))))))),kmadd(ToReal(-6),kmadd(gt12L,kmadd(gtu22,JacPDstandardNth2dphi2,kmadd(gtu31,JacPDstandardNth3dphi1,kmul(gtu33,JacPDstandardNth3dphi3))),JacPDstandardNth2dphi1),kmul(gt12L,kmadd(ToReal(12),kmul(dphi3L,kmul(Gt312,gtu21)),kmadd(ToReal(12),kmul(dphi3L,kmul(Gt313,gtu31)),kmadd(ToReal(12),kmul(dphi3L,kmul(Gt323,gtu32)),kmadd(ToReal(-6),kmul(gtu21,JacPDstandardNth1dphi2),kmadd(ToReal(-6),kmul(gtu31,JacPDstandardNth1dphi3),kmadd(ToReal(-6),kmul(gtu21,JacPDstandardNth2dphi1),kmadd(ToReal(-6),kmul(gtu32,JacPDstandardNth2dphi3),kmadd(ToReal(-6),kmul(gtu32,JacPDstandardNth3dphi2),kmsub(gtu11,kmsub(ToReal(-6),JacPDstandardNth1dphi1,kmul(dphi1L,dphi1L)),kmadd(gtu33,kmul(dphi3L,dphi3L),kmul(gtu22,kmul(dphi2L,dphi2L))))))))))))))))),ToReal(0.0277777777777777777777777777778));
    
    CCTK_REAL_VEC Rphi13 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(6),kmul(dphi3L,kmadd(gt13L,kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmul(Gt333,gtu33))),Gt313)),kmadd(ToReal(2),kmul(dphi2L,kmadd(Gt213,kmadd(ToReal(6),kmul(gt13L,gtu31),ToReal(3)),kmul(gt13L,knmsub(dphi3L,gtu32,kmadd(ToReal(6),kmadd(Gt212,gtu21,kmul(Gt223,gtu32)),kmul(ToReal(3),kmadd(Gt211,gtu11,kmadd(Gt222,gtu22,kmul(Gt233,gtu33))))))))),kmadd(dphi1L,kadd(dphi3L,kmadd(ToReal(6),kmul(Gt113,kmadd(ToReal(2),kmul(gt13L,gtu31),ToReal(1))),kmul(gt13L,kmadd(ToReal(-2),kmul(dphi3L,gtu31),kmul(ToReal(2),kmadd(kmsub(ToReal(6),Gt112,dphi2L),gtu21,kmadd(ToReal(6),kmul(Gt123,gtu32),kmul(ToReal(3),kmadd(Gt111,gtu11,kmadd(Gt122,gtu22,kmul(Gt133,gtu33))))))))))),kmadd(ToReal(-6),kmadd(gt13L,kmul(gtu32,JacPDstandardNth3dphi2),JacPDstandardNth3dphi1),kmul(gt13L,kmadd(ToReal(12),kmul(dphi3L,kmul(Gt312,gtu21)),kmadd(ToReal(12),kmul(dphi3L,kmul(Gt313,gtu31)),kmadd(ToReal(12),kmul(dphi3L,kmul(Gt323,gtu32)),kmadd(ToReal(-6),kmul(gtu21,JacPDstandardNth1dphi2),kmadd(ToReal(-6),kmul(gtu31,JacPDstandardNth1dphi3),kmadd(ToReal(-6),kmul(gtu21,JacPDstandardNth2dphi1),kmadd(ToReal(-6),kmul(gtu32,JacPDstandardNth2dphi3),kmadd(ToReal(-6),kmul(gtu31,JacPDstandardNth3dphi1),kmadd(gtu11,kmsub(ToReal(-6),JacPDstandardNth1dphi1,kmul(dphi1L,dphi1L)),kmadd(gtu22,kmsub(ToReal(-6),JacPDstandardNth2dphi2,kmul(dphi2L,dphi2L)),kmul(gtu33,kmsub(ToReal(-6),JacPDstandardNth3dphi3,kmul(dphi3L,dphi3L)))))))))))))))))),ToReal(0.0277777777777777777777777777778));
    
    CCTK_REAL_VEC Rphi22 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(6),kmadd(kmadd(dphi1L,Gt122,kmul(dphi2L,Gt222)),kmadd(gt22L,gtu22,ToReal(1)),kmul(dphi3L,kmadd(gt22L,kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmul(Gt333,gtu33))),Gt322))),kmadd(ToReal(-6),kmadd(gt22L,kmul(gtu32,kadd(JacPDstandardNth2dphi3,JacPDstandardNth3dphi2)),JacPDstandardNth2dphi2),kmadd(knmsub(gt22L,gtu22,ToReal(1)),kmul(dphi2L,dphi2L),kmul(gt22L,kmadd(ToReal(12),kmul(dphi3L,kmul(Gt312,gtu21)),kmadd(ToReal(12),kmul(dphi3L,kmul(Gt313,gtu31)),kmadd(ToReal(12),kmul(dphi3L,kmul(Gt323,gtu32)),kmadd(ToReal(2),kmul(dphi1L,kmadd(ToReal(3),kmul(Gt111,gtu11),kmadd(ToReal(6),kmul(Gt112,gtu21),kmadd(kmsub(ToReal(6),Gt113,dphi3L),gtu31,kmadd(ToReal(6),kmul(Gt123,gtu32),kmul(ToReal(3),kmul(Gt133,gtu33))))))),kmadd(ToReal(2),kmul(dphi2L,kmadd(ToReal(3),kmul(Gt211,gtu11),kmadd(kmsub(ToReal(6),Gt212,dphi1L),gtu21,kmadd(ToReal(6),kmul(Gt213,gtu31),kmadd(kmsub(ToReal(6),Gt223,dphi3L),gtu32,kmul(ToReal(3),kmul(Gt233,gtu33))))))),kmadd(ToReal(-6),kmul(gtu21,JacPDstandardNth1dphi2),kmadd(ToReal(-6),kmul(gtu31,JacPDstandardNth1dphi3),kmadd(ToReal(-6),kmul(gtu21,JacPDstandardNth2dphi1),kmadd(ToReal(-6),kmul(gtu22,JacPDstandardNth2dphi2),kmadd(ToReal(-6),kmul(gtu31,JacPDstandardNth3dphi1),kmadd(gtu11,kmsub(ToReal(-6),JacPDstandardNth1dphi1,kmul(dphi1L,dphi1L)),kmul(gtu33,kmsub(ToReal(-6),JacPDstandardNth3dphi3,kmul(dphi3L,dphi3L)))))))))))))))))),ToReal(0.0277777777777777777777777777778));
    
    CCTK_REAL_VEC Rphi23 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(6),kmul(dphi3L,kmadd(gt23L,kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmul(Gt333,gtu33))),Gt323)),kmadd(ToReal(2),kmul(dphi1L,kmadd(Gt123,kmadd(ToReal(6),kmul(gt23L,gtu32),ToReal(3)),kmul(gt23L,knmsub(dphi3L,gtu31,kmadd(ToReal(6),kmadd(Gt112,gtu21,kmul(Gt113,gtu31)),kmul(ToReal(3),kmadd(Gt111,gtu11,kmadd(Gt122,gtu22,kmul(Gt133,gtu33))))))))),kmadd(dphi2L,kadd(dphi3L,kmadd(ToReal(6),kmul(Gt223,kmadd(ToReal(2),kmul(gt23L,gtu32),ToReal(1))),kmul(gt23L,kmadd(ToReal(-2),kmul(dphi3L,gtu32),kmul(ToReal(2),kmadd(kmsub(ToReal(6),Gt212,dphi1L),gtu21,kmadd(ToReal(6),kmul(Gt213,gtu31),kmul(ToReal(3),kmadd(Gt211,gtu11,kmadd(Gt222,gtu22,kmul(Gt233,gtu33))))))))))),kmadd(ToReal(-6),kmadd(gt23L,kmul(gtu33,JacPDstandardNth3dphi3),JacPDstandardNth3dphi2),kmul(gt23L,kmadd(ToReal(12),kmul(dphi3L,kmul(Gt312,gtu21)),kmadd(ToReal(12),kmul(dphi3L,kmul(Gt313,gtu31)),kmadd(ToReal(12),kmul(dphi3L,kmul(Gt323,gtu32)),kmadd(ToReal(-6),kmul(gtu21,JacPDstandardNth1dphi2),kmadd(ToReal(-6),kmul(gtu31,JacPDstandardNth1dphi3),kmadd(ToReal(-6),kmul(gtu21,JacPDstandardNth2dphi1),kmadd(ToReal(-6),kmul(gtu32,JacPDstandardNth2dphi3),kmadd(ToReal(-6),kmul(gtu31,JacPDstandardNth3dphi1),kmadd(ToReal(-6),kmul(gtu32,JacPDstandardNth3dphi2),kmadd(gtu11,kmsub(ToReal(-6),JacPDstandardNth1dphi1,kmul(dphi1L,dphi1L)),kmsub(gtu22,kmsub(ToReal(-6),JacPDstandardNth2dphi2,kmul(dphi2L,dphi2L)),kmul(gtu33,kmul(dphi3L,dphi3L)))))))))))))))))),ToReal(0.0277777777777777777777777777778));
    
    CCTK_REAL_VEC Rphi33 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(6),kmadd(kmadd(dphi1L,Gt133,kmul(dphi3L,Gt333)),kmadd(gt33L,gtu33,ToReal(1)),kmul(dphi2L,kmadd(gt33L,kmadd(Gt211,gtu11,kmadd(Gt222,gtu22,kmul(Gt233,gtu33))),Gt233))),kmadd(ToReal(-6),JacPDstandardNth3dphi3,kmadd(gt33L,kmadd(ToReal(12),kmul(dphi2L,kmul(Gt212,gtu21)),kmadd(ToReal(12),kmul(dphi2L,kmul(Gt213,gtu31)),kmadd(ToReal(12),kmul(dphi2L,kmul(Gt223,gtu32)),kmadd(ToReal(2),kmul(dphi1L,kmadd(ToReal(3),kmul(Gt111,gtu11),kmadd(kmsub(ToReal(6),Gt112,dphi2L),gtu21,kmadd(ToReal(3),kmul(Gt122,gtu22),kmadd(ToReal(6),kmul(Gt113,gtu31),kmul(ToReal(6),kmul(Gt123,gtu32))))))),kmadd(ToReal(2),kmul(dphi3L,kmadd(ToReal(3),kmul(Gt311,gtu11),kmadd(ToReal(6),kmul(Gt312,gtu21),kmadd(ToReal(3),kmul(Gt322,gtu22),kmadd(kmsub(ToReal(6),Gt313,dphi1L),gtu31,kmul(kmsub(ToReal(6),Gt323,dphi2L),gtu32)))))),kmadd(ToReal(-6),kmul(gtu21,JacPDstandardNth1dphi2),kmadd(ToReal(-6),kmul(gtu31,JacPDstandardNth1dphi3),kmadd(ToReal(-6),kmul(gtu21,JacPDstandardNth2dphi1),kmadd(ToReal(-6),kmul(gtu32,JacPDstandardNth2dphi3),kmadd(ToReal(-6),kmul(gtu31,JacPDstandardNth3dphi1),kmadd(ToReal(-6),kmul(gtu32,JacPDstandardNth3dphi2),kmadd(ToReal(-6),kmul(gtu33,JacPDstandardNth3dphi3),kmadd(gtu11,kmsub(ToReal(-6),JacPDstandardNth1dphi1,kmul(dphi1L,dphi1L)),kmul(gtu22,kmsub(ToReal(-6),JacPDstandardNth2dphi2,kmul(dphi2L,dphi2L)))))))))))))))),kmul(knmsub(gt33L,gtu33,ToReal(1)),kmul(dphi3L,dphi3L))))),ToReal(0.0277777777777777777777777777778));
    
    CCTK_REAL_VEC e4phi CCTK_ATTRIBUTE_UNUSED = 
      kexp(kmul(phiL,ToReal(4)));
    
    CCTK_REAL_VEC em4phi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),e4phi);
    
    CCTK_REAL_VEC R11 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi11,Rt11);
    
    CCTK_REAL_VEC R12 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi12,Rt12);
    
    CCTK_REAL_VEC R13 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi13,Rt13);
    
    CCTK_REAL_VEC R22 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi22,Rt22);
    
    CCTK_REAL_VEC R23 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi23,Rt23);
    
    CCTK_REAL_VEC R33 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi33,Rt33);
    
    CCTK_REAL_VEC trR CCTK_ATTRIBUTE_UNUSED = 
      kmul(em4phi,kmadd(gtu11,R11,kmadd(gtu22,R22,kmadd(ToReal(2),kmadd(gtu21,R12,kmadd(gtu31,R13,kmul(gtu32,R23))),kmul(gtu33,R33)))));
    
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
    
    CCTK_REAL_VEC rho CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kadd(eTttL,kmadd(ToReal(-2),kmadd(beta2L,eTtyL,kmul(beta3L,eTtzL)),kmadd(ToReal(2),kmadd(beta1L,ksub(kmadd(beta2L,eTxyL,kmul(beta3L,eTxzL)),eTtxL),kmul(beta2L,kmul(beta3L,eTyzL))),kmadd(eTxxL,kmul(beta1L,beta1L),kmadd(eTyyL,kmul(beta2L,beta2L),kmul(eTzzL,kmul(beta3L,beta3L))))))),kmul(alphaL,alphaL));
    
    CCTK_REAL_VEC S1 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(ksub(kmadd(beta1L,eTxxL,kmadd(beta2L,eTxyL,kmul(beta3L,eTxzL))),eTtxL),alphaL);
    
    CCTK_REAL_VEC S2 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(ksub(kmadd(beta1L,eTxyL,kmadd(beta2L,eTyyL,kmul(beta3L,eTyzL))),eTtyL),alphaL);
    
    CCTK_REAL_VEC S3 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(ksub(kmadd(beta1L,eTxzL,kmadd(beta2L,eTyzL,kmul(beta3L,eTzzL))),eTtzL),alphaL);
    
    CCTK_REAL_VEC HL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2.),kmadd(Atm12,Atm21,kmadd(Atm13,Atm31,kmul(Atm23,Atm32))),kmadd(ToReal(-50.26548245743669181540229413247204614715),rho,kadd(trR,kmadd(ToReal(0.6666666666666666666666666666666666666667),kmul(trKL,trKL),kmul(kmadd(Atm11,Atm11,kmadd(Atm22,Atm22,kmul(Atm33,Atm33))),ToReal(-1.))))));
    
    CCTK_REAL_VEC M1L CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(Gt311,kmadd(ToReal(-4.),kmul(At13L,gtu11),kmul(ToReal(-2.),kmadd(At23L,gtu21,kmul(At33L,gtu31)))),kmadd(At12L,kmadd(kadd(dphi1L,kmadd(ToReal(-2.),Gt111,kmul(ToReal(-6.),Gt212))),gtu21,kmadd(kmadd(ToReal(-2.),kadd(Gt112,Gt222),dphi2L),gtu22,kmadd(ToReal(-6.),kmul(Gt213,gtu31),kmadd(kmadd(ToReal(-2.),Gt113,dphi3L),gtu32,kmadd(ToReal(-4.),kmadd(Gt211,gtu11,kmul(Gt223,gtu32)),kmul(ToReal(-2.),kmul(Gt233,gtu33))))))),kmadd(ToReal(-2.),kmadd(kmadd(At22L,Gt212,kmadd(At23L,Gt312,kmul(At13L,Gt322))),gtu22,kmul(At13L,kmadd(Gt111,gtu31,kmadd(Gt112,gtu32,kmul(Gt113,gtu33))))),kmadd(At11L,kmadd(kmadd(ToReal(-4.),Gt111,dphi1L),gtu11,kmadd(kmadd(ToReal(-6.),Gt112,dphi2L),gtu21,kmadd(kmadd(ToReal(-6.),Gt113,dphi3L),gtu31,kmadd(ToReal(-4.),kmul(Gt123,gtu32),kmul(ToReal(-2.),kmadd(Gt122,gtu22,kmul(Gt133,gtu33))))))),kmadd(gtu21,kmadd(ToReal(-2.),kmul(At22L,Gt211),kmadd(ToReal(-6.),kmul(At13L,Gt312),kmul(ToReal(2.),JacPDstandardNth1At12))),kmadd(ToReal(-1.333333333333333333333333333333333333333),JacPDstandardNth1trK,kmadd(gtu32,kmadd(ToReal(-2.),kmadd(At22L,Gt213,kmadd(At33L,Gt312,kmul(At23L,kadd(Gt212,Gt313)))),kmadd(At13L,kmadd(ToReal(-4.),Gt323,dphi2L),kmul(ToReal(2.),JacPDstandardNth2At13))),kmadd(gtu31,kmadd(ToReal(-2.),kmul(At23L,Gt211),kmadd(At13L,kmadd(ToReal(-6.),Gt313,dphi1L),kmul(ToReal(2.),kadd(JacPDstandardNth1At13,JacPDstandardNth3At11)))),kmadd(ToReal(2.),kmadd(gtu11,JacPDstandardNth1At11,kmadd(gtu21,JacPDstandardNth2At11,kmadd(gtu22,JacPDstandardNth2At12,kmul(gtu32,JacPDstandardNth3At12)))),kmadd(gtu33,kmadd(ToReal(-2.),kmadd(At23L,Gt213,kmul(At33L,Gt313)),kmadd(At13L,kmadd(ToReal(-2.),Gt333,dphi3L),kmul(ToReal(2.),JacPDstandardNth3At13))),kmul(ToReal(-50.26548245743669181540229413247204614715),S1))))))))))),ToReal(0.5));
    
    CCTK_REAL_VEC M2L CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(At12L,kmadd(kmadd(ToReal(-2.),kadd(Gt111,Gt212),dphi1L),gtu11,kmadd(kadd(dphi2L,kmadd(ToReal(-6.),Gt112,kmul(ToReal(-2.),Gt222))),gtu21,kmadd(kmadd(ToReal(-2.),Gt223,dphi3L),gtu31,kmadd(ToReal(-4.),kmadd(Gt122,gtu22,kmul(Gt113,gtu31)),kmadd(ToReal(-6.),kmul(Gt123,gtu32),kmul(ToReal(-2.),kmul(Gt133,gtu33))))))),kmadd(At22L,kmadd(kmadd(ToReal(-6.),Gt212,dphi1L),gtu21,kmadd(kmadd(ToReal(-4.),Gt222,dphi2L),gtu22,kmadd(ToReal(-4.),kmul(Gt213,gtu31),kmadd(dphi3L,gtu32,kmul(ToReal(-2.),kmul(Gt233,gtu33)))))),kmadd(At23L,kmadd(ToReal(-6.),kmul(Gt312,gtu21),kmadd(ToReal(-4.),kmul(Gt322,gtu22),kmadd(dphi1L,gtu31,kmadd(dphi2L,gtu32,kmul(kmadd(ToReal(-2.),Gt333,dphi3L),gtu33))))),kmadd(ToReal(-2.),kmadd(kmadd(At22L,Gt211,kmadd(At23L,Gt311,kmul(At13L,Gt312))),gtu11,kmadd(kmadd(At23L,Gt212,kmul(At33L,Gt312)),gtu31,kmadd(At11L,kmadd(Gt112,gtu11,kmadd(Gt122,gtu21,kmul(Gt123,gtu31))),kmadd(kmadd(At23L,Gt223,kmul(At33L,Gt323)),gtu33,kmul(At13L,kmadd(Gt322,gtu21,kmadd(Gt112,gtu31,kmadd(Gt122,gtu32,kmul(Gt123,gtu33))))))))),kmadd(gtu32,kmadd(ToReal(-2.),kmadd(At23L,Gt222,kmul(At33L,Gt322)),kmadd(ToReal(-6.),kmadd(At22L,Gt223,kmul(At23L,Gt323)),kmul(ToReal(2.),JacPDstandardNth2At23))),kmadd(ToReal(-1.333333333333333333333333333333333333333),JacPDstandardNth2trK,kmadd(gtu31,kmadd(ToReal(-4.),kmul(At23L,Gt313),kmadd(ToReal(-2.),kmul(At13L,Gt323),kmul(ToReal(2.),kadd(JacPDstandardNth1At23,JacPDstandardNth3At12)))),kmadd(ToReal(2.),kmadd(gtu11,JacPDstandardNth1At12,kmadd(gtu21,kadd(JacPDstandardNth1At22,JacPDstandardNth2At12),kmadd(gtu22,JacPDstandardNth2At22,kmadd(gtu32,JacPDstandardNth3At22,kmul(gtu33,JacPDstandardNth3At23))))),kmul(ToReal(-50.26548245743669181540229413247204614715),S2))))))))),ToReal(0.5));
    
    CCTK_REAL_VEC M3L CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(-2.),kmadd(kmadd(At23L,Gt211,kmadd(At12L,Gt213,kmul(At33L,Gt311))),gtu11,kmadd(kmadd(At22L,Gt213,kmul(At12L,kadd(Gt113,Gt223))),gtu21,kmadd(kmadd(At12L,Gt123,kmadd(At22L,Gt223,kmadd(At33L,Gt322,kmul(At23L,kadd(Gt222,Gt323))))),gtu22,kmadd(At11L,kmadd(Gt113,gtu11,kmadd(Gt123,gtu21,kmul(Gt133,gtu31))),kmul(At12L,kmadd(Gt233,gtu31,kmul(Gt133,gtu32))))))),kmadd(At13L,kmadd(kmadd(ToReal(-2.),kadd(Gt111,Gt313),dphi1L),gtu11,kmadd(kadd(dphi2L,kmadd(ToReal(-4.),Gt112,kmul(ToReal(-2.),Gt323))),gtu21,kmadd(kmadd(ToReal(-6.),Gt113,dphi3L),gtu31,kmadd(ToReal(-2.),kmadd(Gt122,gtu22,kmul(Gt333,gtu31)),kmadd(ToReal(-6.),kmul(Gt123,gtu32),kmul(ToReal(-4.),kmul(Gt133,gtu33))))))),kmadd(At23L,kmadd(kmadd(ToReal(-4.),Gt212,dphi1L),gtu21,kmadd(dphi2L,gtu22,kmadd(ToReal(-6.),kmul(Gt213,gtu31),kmadd(dphi3L,gtu32,kmul(ToReal(-4.),kmul(Gt233,gtu33)))))),kmadd(At33L,kmadd(kmadd(ToReal(-6.),Gt313,dphi1L),gtu31,kmadd(dphi2L,gtu32,kmul(kmadd(ToReal(-4.),Gt333,dphi3L),gtu33))),kmadd(gtu21,kmadd(ToReal(-4.),kmul(At33L,Gt312),kmadd(ToReal(-2.),kmul(At23L,Gt313),kmul(ToReal(2.),kadd(JacPDstandardNth1At23,JacPDstandardNth2At13)))),kmadd(gtu32,kmadd(ToReal(-6.),kmadd(At23L,Gt223,kmul(At33L,Gt323)),kmadd(ToReal(-2.),kmadd(At22L,Gt233,kmul(At23L,Gt333)),kmul(ToReal(2.),kadd(JacPDstandardNth2At33,JacPDstandardNth3At23)))),kmadd(ToReal(2.),kmadd(gtu11,JacPDstandardNth1At13,kmadd(gtu22,JacPDstandardNth2At23,kmadd(gtu31,kadd(JacPDstandardNth1At33,JacPDstandardNth3At13),kmul(gtu33,JacPDstandardNth3At33)))),kmadd(ToReal(-1.333333333333333333333333333333333333333),JacPDstandardNth3trK,kmul(ToReal(-50.26548245743669181540229413247204614715),S3))))))))),ToReal(0.5));
    
    CCTK_REAL_VEC cSL CCTK_ATTRIBUTE_UNUSED = klog(detgt);
    
    CCTK_REAL_VEC cXt1L CCTK_ATTRIBUTE_UNUSED = 
      ksub(kmadd(Gt111,gtu11,kmadd(Gt122,gtu22,kmadd(ToReal(2),kmadd(Gt112,gtu21,kmadd(Gt113,gtu31,kmul(Gt123,gtu32))),kmul(Gt133,gtu33)))),Xt1L);
    
    CCTK_REAL_VEC cXt2L CCTK_ATTRIBUTE_UNUSED = 
      ksub(kmadd(Gt211,gtu11,kmadd(Gt222,gtu22,kmadd(ToReal(2),kmadd(Gt212,gtu21,kmadd(Gt213,gtu31,kmul(Gt223,gtu32))),kmul(Gt233,gtu33)))),Xt2L);
    
    CCTK_REAL_VEC cXt3L CCTK_ATTRIBUTE_UNUSED = 
      ksub(kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmadd(ToReal(2),kmadd(Gt312,gtu21,kmadd(Gt313,gtu31,kmul(Gt323,gtu32))),kmul(Gt333,gtu33)))),Xt3L);
    
    CCTK_REAL_VEC cAL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At11L,gtu11,kmadd(At22L,gtu22,kmadd(ToReal(2),kmadd(At12L,gtu21,kmadd(At13L,gtu31,kmul(At23L,gtu32))),kmul(At33L,gtu33))));
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(cA[index],cAL);
    vec_store_nta_partial(cS[index],cSL);
    vec_store_nta_partial(cXt1[index],cXt1L);
    vec_store_nta_partial(cXt2[index],cXt2L);
    vec_store_nta_partial(cXt3[index],cXt3L);
    vec_store_nta_partial(H[index],HL);
    vec_store_nta_partial(M1[index],M1L);
    vec_store_nta_partial(M2[index],M2L);
    vec_store_nta_partial(M3[index],M3L);
  }
  CCTK_ENDLOOP3STR(CL_BSSN_constraints);
}
extern "C" void CL_BSSN_constraints(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CL_BSSN_constraints
  DECLARE_CCTK_ARGUMENTS_CHECKED(CL_BSSN_constraints);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering CL_BSSN_constraints_Body");
  }
  if (cctk_iteration % CL_BSSN_constraints_calc_every != CL_BSSN_constraints_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "CL_BSSN::CL_cons_detg",
    "CL_BSSN::CL_cons_Gamma",
    "CL_BSSN::CL_cons_traceA",
    "CL_BSSN::CL_curv",
    "CL_BSSN::CL_dlog_confac",
    "CL_BSSN::CL_dmetric",
    "CL_BSSN::CL_Gamma",
    "CL_BSSN::CL_Ham",
    "CL_BSSN::CL_lapse",
    "CL_BSSN::CL_log_confac",
    "CL_BSSN::CL_metric",
    "CL_BSSN::CL_mom",
    "CL_BSSN::CL_shift",
    "CL_BSSN::CL_trace_curv"};
  AssertGroupStorage(cctkGH, "CL_BSSN_constraints", 14, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "CL_BSSN_constraints", 1, 1, 1);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "CL_BSSN_constraints", 2, 2, 2);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "CL_BSSN_constraints", 3, 3, 3);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "CL_BSSN_constraints", 4, 4, 4);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, CL_BSSN_constraints_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving CL_BSSN_constraints_Body");
  }
}

} // namespace CL_BSSN
