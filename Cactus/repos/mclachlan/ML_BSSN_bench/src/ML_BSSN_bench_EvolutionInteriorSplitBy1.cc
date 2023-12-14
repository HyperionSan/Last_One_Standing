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

namespace ML_BSSN_bench {

extern "C" void ML_BSSN_bench_EvolutionInteriorSplitBy1_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_BSSN_bench_EvolutionInteriorSplitBy1_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_bench_EvolutionInteriorSplitBy1_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_bench_EvolutionInteriorSplitBy1_calc_every != ML_BSSN_bench_EvolutionInteriorSplitBy1_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_bench::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN_bench::ML_lapserhs.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_bench::ML_log_confacrhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN_bench::ML_log_confacrhs.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_bench::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN_bench::ML_metricrhs.");
  return;
}

static void ML_BSSN_bench_EvolutionInteriorSplitBy1_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  const CCTK_REAL_VEC p1o180dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dx,dx));
  const CCTK_REAL_VEC p1o180dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dy,dy));
  const CCTK_REAL_VEC p1o180dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dz,dz));
  const CCTK_REAL_VEC p1o24dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dx);
  const CCTK_REAL_VEC p1o24dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dy);
  const CCTK_REAL_VEC p1o24dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dz);
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
  const CCTK_REAL_VEC p1o6dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.166666666666666666666666666667),dx);
  const CCTK_REAL_VEC p1o6dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.166666666666666666666666666667),dy);
  const CCTK_REAL_VEC p1o6dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.166666666666666666666666666667),dz);
  const CCTK_REAL_VEC p1o705600dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dx,dy));
  const CCTK_REAL_VEC p1o705600dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dx,dz));
  const CCTK_REAL_VEC p1o705600dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dy,dz));
  const CCTK_REAL_VEC p1o840dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dx);
  const CCTK_REAL_VEC p1o840dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dy);
  const CCTK_REAL_VEC p1o840dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dz);
  const CCTK_REAL_VEC p1odx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dx,dx));
  const CCTK_REAL_VEC p1ody2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dy,dy));
  const CCTK_REAL_VEC p1odz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dz,dz));
  const CCTK_REAL_VEC pm1o120dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00833333333333333333333333333333),dx);
  const CCTK_REAL_VEC pm1o120dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00833333333333333333333333333333),dy);
  const CCTK_REAL_VEC pm1o120dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00833333333333333333333333333333),dz);
  const CCTK_REAL_VEC pm1o12dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),dx);
  const CCTK_REAL_VEC pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dx,dx));
  const CCTK_REAL_VEC pm1o12dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),dy);
  const CCTK_REAL_VEC pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dy,dy));
  const CCTK_REAL_VEC pm1o12dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),dz);
  const CCTK_REAL_VEC pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dz,dz));
  const CCTK_REAL_VEC pm1o16dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0625),dx);
  const CCTK_REAL_VEC pm1o16dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0625),dy);
  const CCTK_REAL_VEC pm1o16dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0625),dz);
  const CCTK_REAL_VEC pm1o256dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00390625),dx);
  const CCTK_REAL_VEC pm1o256dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00390625),dy);
  const CCTK_REAL_VEC pm1o256dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00390625),dz);
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
  CCTK_LOOP3STR(ML_BSSN_bench_EvolutionInteriorSplitBy1,
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
    CCTK_REAL_VEC gt11L CCTK_ATTRIBUTE_UNUSED = vec_load(gt11[index]);
    CCTK_REAL_VEC gt12L CCTK_ATTRIBUTE_UNUSED = vec_load(gt12[index]);
    CCTK_REAL_VEC gt13L CCTK_ATTRIBUTE_UNUSED = vec_load(gt13[index]);
    CCTK_REAL_VEC gt22L CCTK_ATTRIBUTE_UNUSED = vec_load(gt22[index]);
    CCTK_REAL_VEC gt23L CCTK_ATTRIBUTE_UNUSED = vec_load(gt23[index]);
    CCTK_REAL_VEC gt33L CCTK_ATTRIBUTE_UNUSED = vec_load(gt33[index]);
    CCTK_REAL_VEC phiL CCTK_ATTRIBUTE_UNUSED = vec_load(phi[index]);
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
    CCTK_REAL_VEC PDupwindNthSymm1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDdissipationNth3phi CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder21(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder22(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder23(&alpha[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder21(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder22(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder23(&alpha[index]);
        PDdissipationNth1alpha = PDdissipationNthfdOrder21(&alpha[index]);
        PDdissipationNth2alpha = PDdissipationNthfdOrder22(&alpha[index]);
        PDdissipationNth3alpha = PDdissipationNthfdOrder23(&alpha[index]);
        PDstandardNth1beta1 = PDstandardNthfdOrder21(&beta1[index]);
        PDstandardNth2beta1 = PDstandardNthfdOrder22(&beta1[index]);
        PDstandardNth3beta1 = PDstandardNthfdOrder23(&beta1[index]);
        PDstandardNth1beta2 = PDstandardNthfdOrder21(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder22(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder23(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder21(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder22(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder23(&beta3[index]);
        PDupwindNthSymm1gt11 = PDupwindNthSymmfdOrder21(&gt11[index]);
        PDupwindNthSymm2gt11 = PDupwindNthSymmfdOrder22(&gt11[index]);
        PDupwindNthSymm3gt11 = PDupwindNthSymmfdOrder23(&gt11[index]);
        PDupwindNthAnti1gt11 = PDupwindNthAntifdOrder21(&gt11[index]);
        PDupwindNthAnti2gt11 = PDupwindNthAntifdOrder22(&gt11[index]);
        PDupwindNthAnti3gt11 = PDupwindNthAntifdOrder23(&gt11[index]);
        PDdissipationNth1gt11 = PDdissipationNthfdOrder21(&gt11[index]);
        PDdissipationNth2gt11 = PDdissipationNthfdOrder22(&gt11[index]);
        PDdissipationNth3gt11 = PDdissipationNthfdOrder23(&gt11[index]);
        PDupwindNthSymm1gt12 = PDupwindNthSymmfdOrder21(&gt12[index]);
        PDupwindNthSymm2gt12 = PDupwindNthSymmfdOrder22(&gt12[index]);
        PDupwindNthSymm3gt12 = PDupwindNthSymmfdOrder23(&gt12[index]);
        PDupwindNthAnti1gt12 = PDupwindNthAntifdOrder21(&gt12[index]);
        PDupwindNthAnti2gt12 = PDupwindNthAntifdOrder22(&gt12[index]);
        PDupwindNthAnti3gt12 = PDupwindNthAntifdOrder23(&gt12[index]);
        PDdissipationNth1gt12 = PDdissipationNthfdOrder21(&gt12[index]);
        PDdissipationNth2gt12 = PDdissipationNthfdOrder22(&gt12[index]);
        PDdissipationNth3gt12 = PDdissipationNthfdOrder23(&gt12[index]);
        PDupwindNthSymm1gt13 = PDupwindNthSymmfdOrder21(&gt13[index]);
        PDupwindNthSymm2gt13 = PDupwindNthSymmfdOrder22(&gt13[index]);
        PDupwindNthSymm3gt13 = PDupwindNthSymmfdOrder23(&gt13[index]);
        PDupwindNthAnti1gt13 = PDupwindNthAntifdOrder21(&gt13[index]);
        PDupwindNthAnti2gt13 = PDupwindNthAntifdOrder22(&gt13[index]);
        PDupwindNthAnti3gt13 = PDupwindNthAntifdOrder23(&gt13[index]);
        PDdissipationNth1gt13 = PDdissipationNthfdOrder21(&gt13[index]);
        PDdissipationNth2gt13 = PDdissipationNthfdOrder22(&gt13[index]);
        PDdissipationNth3gt13 = PDdissipationNthfdOrder23(&gt13[index]);
        PDupwindNthSymm1gt22 = PDupwindNthSymmfdOrder21(&gt22[index]);
        PDupwindNthSymm2gt22 = PDupwindNthSymmfdOrder22(&gt22[index]);
        PDupwindNthSymm3gt22 = PDupwindNthSymmfdOrder23(&gt22[index]);
        PDupwindNthAnti1gt22 = PDupwindNthAntifdOrder21(&gt22[index]);
        PDupwindNthAnti2gt22 = PDupwindNthAntifdOrder22(&gt22[index]);
        PDupwindNthAnti3gt22 = PDupwindNthAntifdOrder23(&gt22[index]);
        PDdissipationNth1gt22 = PDdissipationNthfdOrder21(&gt22[index]);
        PDdissipationNth2gt22 = PDdissipationNthfdOrder22(&gt22[index]);
        PDdissipationNth3gt22 = PDdissipationNthfdOrder23(&gt22[index]);
        PDupwindNthSymm1gt23 = PDupwindNthSymmfdOrder21(&gt23[index]);
        PDupwindNthSymm2gt23 = PDupwindNthSymmfdOrder22(&gt23[index]);
        PDupwindNthSymm3gt23 = PDupwindNthSymmfdOrder23(&gt23[index]);
        PDupwindNthAnti1gt23 = PDupwindNthAntifdOrder21(&gt23[index]);
        PDupwindNthAnti2gt23 = PDupwindNthAntifdOrder22(&gt23[index]);
        PDupwindNthAnti3gt23 = PDupwindNthAntifdOrder23(&gt23[index]);
        PDdissipationNth1gt23 = PDdissipationNthfdOrder21(&gt23[index]);
        PDdissipationNth2gt23 = PDdissipationNthfdOrder22(&gt23[index]);
        PDdissipationNth3gt23 = PDdissipationNthfdOrder23(&gt23[index]);
        PDupwindNthSymm1gt33 = PDupwindNthSymmfdOrder21(&gt33[index]);
        PDupwindNthSymm2gt33 = PDupwindNthSymmfdOrder22(&gt33[index]);
        PDupwindNthSymm3gt33 = PDupwindNthSymmfdOrder23(&gt33[index]);
        PDupwindNthAnti1gt33 = PDupwindNthAntifdOrder21(&gt33[index]);
        PDupwindNthAnti2gt33 = PDupwindNthAntifdOrder22(&gt33[index]);
        PDupwindNthAnti3gt33 = PDupwindNthAntifdOrder23(&gt33[index]);
        PDdissipationNth1gt33 = PDdissipationNthfdOrder21(&gt33[index]);
        PDdissipationNth2gt33 = PDdissipationNthfdOrder22(&gt33[index]);
        PDdissipationNth3gt33 = PDdissipationNthfdOrder23(&gt33[index]);
        PDupwindNthSymm1phi = PDupwindNthSymmfdOrder21(&phi[index]);
        PDupwindNthSymm2phi = PDupwindNthSymmfdOrder22(&phi[index]);
        PDupwindNthSymm3phi = PDupwindNthSymmfdOrder23(&phi[index]);
        PDupwindNthAnti1phi = PDupwindNthAntifdOrder21(&phi[index]);
        PDupwindNthAnti2phi = PDupwindNthAntifdOrder22(&phi[index]);
        PDupwindNthAnti3phi = PDupwindNthAntifdOrder23(&phi[index]);
        PDdissipationNth1phi = PDdissipationNthfdOrder21(&phi[index]);
        PDdissipationNth2phi = PDdissipationNthfdOrder22(&phi[index]);
        PDdissipationNth3phi = PDdissipationNthfdOrder23(&phi[index]);
        break;
      }
      
      case 4:
      {
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder41(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder42(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder43(&alpha[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder41(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder42(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder43(&alpha[index]);
        PDdissipationNth1alpha = PDdissipationNthfdOrder41(&alpha[index]);
        PDdissipationNth2alpha = PDdissipationNthfdOrder42(&alpha[index]);
        PDdissipationNth3alpha = PDdissipationNthfdOrder43(&alpha[index]);
        PDstandardNth1beta1 = PDstandardNthfdOrder41(&beta1[index]);
        PDstandardNth2beta1 = PDstandardNthfdOrder42(&beta1[index]);
        PDstandardNth3beta1 = PDstandardNthfdOrder43(&beta1[index]);
        PDstandardNth1beta2 = PDstandardNthfdOrder41(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder42(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder43(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder41(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder42(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder43(&beta3[index]);
        PDupwindNthSymm1gt11 = PDupwindNthSymmfdOrder41(&gt11[index]);
        PDupwindNthSymm2gt11 = PDupwindNthSymmfdOrder42(&gt11[index]);
        PDupwindNthSymm3gt11 = PDupwindNthSymmfdOrder43(&gt11[index]);
        PDupwindNthAnti1gt11 = PDupwindNthAntifdOrder41(&gt11[index]);
        PDupwindNthAnti2gt11 = PDupwindNthAntifdOrder42(&gt11[index]);
        PDupwindNthAnti3gt11 = PDupwindNthAntifdOrder43(&gt11[index]);
        PDdissipationNth1gt11 = PDdissipationNthfdOrder41(&gt11[index]);
        PDdissipationNth2gt11 = PDdissipationNthfdOrder42(&gt11[index]);
        PDdissipationNth3gt11 = PDdissipationNthfdOrder43(&gt11[index]);
        PDupwindNthSymm1gt12 = PDupwindNthSymmfdOrder41(&gt12[index]);
        PDupwindNthSymm2gt12 = PDupwindNthSymmfdOrder42(&gt12[index]);
        PDupwindNthSymm3gt12 = PDupwindNthSymmfdOrder43(&gt12[index]);
        PDupwindNthAnti1gt12 = PDupwindNthAntifdOrder41(&gt12[index]);
        PDupwindNthAnti2gt12 = PDupwindNthAntifdOrder42(&gt12[index]);
        PDupwindNthAnti3gt12 = PDupwindNthAntifdOrder43(&gt12[index]);
        PDdissipationNth1gt12 = PDdissipationNthfdOrder41(&gt12[index]);
        PDdissipationNth2gt12 = PDdissipationNthfdOrder42(&gt12[index]);
        PDdissipationNth3gt12 = PDdissipationNthfdOrder43(&gt12[index]);
        PDupwindNthSymm1gt13 = PDupwindNthSymmfdOrder41(&gt13[index]);
        PDupwindNthSymm2gt13 = PDupwindNthSymmfdOrder42(&gt13[index]);
        PDupwindNthSymm3gt13 = PDupwindNthSymmfdOrder43(&gt13[index]);
        PDupwindNthAnti1gt13 = PDupwindNthAntifdOrder41(&gt13[index]);
        PDupwindNthAnti2gt13 = PDupwindNthAntifdOrder42(&gt13[index]);
        PDupwindNthAnti3gt13 = PDupwindNthAntifdOrder43(&gt13[index]);
        PDdissipationNth1gt13 = PDdissipationNthfdOrder41(&gt13[index]);
        PDdissipationNth2gt13 = PDdissipationNthfdOrder42(&gt13[index]);
        PDdissipationNth3gt13 = PDdissipationNthfdOrder43(&gt13[index]);
        PDupwindNthSymm1gt22 = PDupwindNthSymmfdOrder41(&gt22[index]);
        PDupwindNthSymm2gt22 = PDupwindNthSymmfdOrder42(&gt22[index]);
        PDupwindNthSymm3gt22 = PDupwindNthSymmfdOrder43(&gt22[index]);
        PDupwindNthAnti1gt22 = PDupwindNthAntifdOrder41(&gt22[index]);
        PDupwindNthAnti2gt22 = PDupwindNthAntifdOrder42(&gt22[index]);
        PDupwindNthAnti3gt22 = PDupwindNthAntifdOrder43(&gt22[index]);
        PDdissipationNth1gt22 = PDdissipationNthfdOrder41(&gt22[index]);
        PDdissipationNth2gt22 = PDdissipationNthfdOrder42(&gt22[index]);
        PDdissipationNth3gt22 = PDdissipationNthfdOrder43(&gt22[index]);
        PDupwindNthSymm1gt23 = PDupwindNthSymmfdOrder41(&gt23[index]);
        PDupwindNthSymm2gt23 = PDupwindNthSymmfdOrder42(&gt23[index]);
        PDupwindNthSymm3gt23 = PDupwindNthSymmfdOrder43(&gt23[index]);
        PDupwindNthAnti1gt23 = PDupwindNthAntifdOrder41(&gt23[index]);
        PDupwindNthAnti2gt23 = PDupwindNthAntifdOrder42(&gt23[index]);
        PDupwindNthAnti3gt23 = PDupwindNthAntifdOrder43(&gt23[index]);
        PDdissipationNth1gt23 = PDdissipationNthfdOrder41(&gt23[index]);
        PDdissipationNth2gt23 = PDdissipationNthfdOrder42(&gt23[index]);
        PDdissipationNth3gt23 = PDdissipationNthfdOrder43(&gt23[index]);
        PDupwindNthSymm1gt33 = PDupwindNthSymmfdOrder41(&gt33[index]);
        PDupwindNthSymm2gt33 = PDupwindNthSymmfdOrder42(&gt33[index]);
        PDupwindNthSymm3gt33 = PDupwindNthSymmfdOrder43(&gt33[index]);
        PDupwindNthAnti1gt33 = PDupwindNthAntifdOrder41(&gt33[index]);
        PDupwindNthAnti2gt33 = PDupwindNthAntifdOrder42(&gt33[index]);
        PDupwindNthAnti3gt33 = PDupwindNthAntifdOrder43(&gt33[index]);
        PDdissipationNth1gt33 = PDdissipationNthfdOrder41(&gt33[index]);
        PDdissipationNth2gt33 = PDdissipationNthfdOrder42(&gt33[index]);
        PDdissipationNth3gt33 = PDdissipationNthfdOrder43(&gt33[index]);
        PDupwindNthSymm1phi = PDupwindNthSymmfdOrder41(&phi[index]);
        PDupwindNthSymm2phi = PDupwindNthSymmfdOrder42(&phi[index]);
        PDupwindNthSymm3phi = PDupwindNthSymmfdOrder43(&phi[index]);
        PDupwindNthAnti1phi = PDupwindNthAntifdOrder41(&phi[index]);
        PDupwindNthAnti2phi = PDupwindNthAntifdOrder42(&phi[index]);
        PDupwindNthAnti3phi = PDupwindNthAntifdOrder43(&phi[index]);
        PDdissipationNth1phi = PDdissipationNthfdOrder41(&phi[index]);
        PDdissipationNth2phi = PDdissipationNthfdOrder42(&phi[index]);
        PDdissipationNth3phi = PDdissipationNthfdOrder43(&phi[index]);
        break;
      }
      
      case 6:
      {
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder61(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder62(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder63(&alpha[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder61(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder62(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder63(&alpha[index]);
        PDdissipationNth1alpha = PDdissipationNthfdOrder61(&alpha[index]);
        PDdissipationNth2alpha = PDdissipationNthfdOrder62(&alpha[index]);
        PDdissipationNth3alpha = PDdissipationNthfdOrder63(&alpha[index]);
        PDstandardNth1beta1 = PDstandardNthfdOrder61(&beta1[index]);
        PDstandardNth2beta1 = PDstandardNthfdOrder62(&beta1[index]);
        PDstandardNth3beta1 = PDstandardNthfdOrder63(&beta1[index]);
        PDstandardNth1beta2 = PDstandardNthfdOrder61(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder62(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder63(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder61(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder62(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder63(&beta3[index]);
        PDupwindNthSymm1gt11 = PDupwindNthSymmfdOrder61(&gt11[index]);
        PDupwindNthSymm2gt11 = PDupwindNthSymmfdOrder62(&gt11[index]);
        PDupwindNthSymm3gt11 = PDupwindNthSymmfdOrder63(&gt11[index]);
        PDupwindNthAnti1gt11 = PDupwindNthAntifdOrder61(&gt11[index]);
        PDupwindNthAnti2gt11 = PDupwindNthAntifdOrder62(&gt11[index]);
        PDupwindNthAnti3gt11 = PDupwindNthAntifdOrder63(&gt11[index]);
        PDdissipationNth1gt11 = PDdissipationNthfdOrder61(&gt11[index]);
        PDdissipationNth2gt11 = PDdissipationNthfdOrder62(&gt11[index]);
        PDdissipationNth3gt11 = PDdissipationNthfdOrder63(&gt11[index]);
        PDupwindNthSymm1gt12 = PDupwindNthSymmfdOrder61(&gt12[index]);
        PDupwindNthSymm2gt12 = PDupwindNthSymmfdOrder62(&gt12[index]);
        PDupwindNthSymm3gt12 = PDupwindNthSymmfdOrder63(&gt12[index]);
        PDupwindNthAnti1gt12 = PDupwindNthAntifdOrder61(&gt12[index]);
        PDupwindNthAnti2gt12 = PDupwindNthAntifdOrder62(&gt12[index]);
        PDupwindNthAnti3gt12 = PDupwindNthAntifdOrder63(&gt12[index]);
        PDdissipationNth1gt12 = PDdissipationNthfdOrder61(&gt12[index]);
        PDdissipationNth2gt12 = PDdissipationNthfdOrder62(&gt12[index]);
        PDdissipationNth3gt12 = PDdissipationNthfdOrder63(&gt12[index]);
        PDupwindNthSymm1gt13 = PDupwindNthSymmfdOrder61(&gt13[index]);
        PDupwindNthSymm2gt13 = PDupwindNthSymmfdOrder62(&gt13[index]);
        PDupwindNthSymm3gt13 = PDupwindNthSymmfdOrder63(&gt13[index]);
        PDupwindNthAnti1gt13 = PDupwindNthAntifdOrder61(&gt13[index]);
        PDupwindNthAnti2gt13 = PDupwindNthAntifdOrder62(&gt13[index]);
        PDupwindNthAnti3gt13 = PDupwindNthAntifdOrder63(&gt13[index]);
        PDdissipationNth1gt13 = PDdissipationNthfdOrder61(&gt13[index]);
        PDdissipationNth2gt13 = PDdissipationNthfdOrder62(&gt13[index]);
        PDdissipationNth3gt13 = PDdissipationNthfdOrder63(&gt13[index]);
        PDupwindNthSymm1gt22 = PDupwindNthSymmfdOrder61(&gt22[index]);
        PDupwindNthSymm2gt22 = PDupwindNthSymmfdOrder62(&gt22[index]);
        PDupwindNthSymm3gt22 = PDupwindNthSymmfdOrder63(&gt22[index]);
        PDupwindNthAnti1gt22 = PDupwindNthAntifdOrder61(&gt22[index]);
        PDupwindNthAnti2gt22 = PDupwindNthAntifdOrder62(&gt22[index]);
        PDupwindNthAnti3gt22 = PDupwindNthAntifdOrder63(&gt22[index]);
        PDdissipationNth1gt22 = PDdissipationNthfdOrder61(&gt22[index]);
        PDdissipationNth2gt22 = PDdissipationNthfdOrder62(&gt22[index]);
        PDdissipationNth3gt22 = PDdissipationNthfdOrder63(&gt22[index]);
        PDupwindNthSymm1gt23 = PDupwindNthSymmfdOrder61(&gt23[index]);
        PDupwindNthSymm2gt23 = PDupwindNthSymmfdOrder62(&gt23[index]);
        PDupwindNthSymm3gt23 = PDupwindNthSymmfdOrder63(&gt23[index]);
        PDupwindNthAnti1gt23 = PDupwindNthAntifdOrder61(&gt23[index]);
        PDupwindNthAnti2gt23 = PDupwindNthAntifdOrder62(&gt23[index]);
        PDupwindNthAnti3gt23 = PDupwindNthAntifdOrder63(&gt23[index]);
        PDdissipationNth1gt23 = PDdissipationNthfdOrder61(&gt23[index]);
        PDdissipationNth2gt23 = PDdissipationNthfdOrder62(&gt23[index]);
        PDdissipationNth3gt23 = PDdissipationNthfdOrder63(&gt23[index]);
        PDupwindNthSymm1gt33 = PDupwindNthSymmfdOrder61(&gt33[index]);
        PDupwindNthSymm2gt33 = PDupwindNthSymmfdOrder62(&gt33[index]);
        PDupwindNthSymm3gt33 = PDupwindNthSymmfdOrder63(&gt33[index]);
        PDupwindNthAnti1gt33 = PDupwindNthAntifdOrder61(&gt33[index]);
        PDupwindNthAnti2gt33 = PDupwindNthAntifdOrder62(&gt33[index]);
        PDupwindNthAnti3gt33 = PDupwindNthAntifdOrder63(&gt33[index]);
        PDdissipationNth1gt33 = PDdissipationNthfdOrder61(&gt33[index]);
        PDdissipationNth2gt33 = PDdissipationNthfdOrder62(&gt33[index]);
        PDdissipationNth3gt33 = PDdissipationNthfdOrder63(&gt33[index]);
        PDupwindNthSymm1phi = PDupwindNthSymmfdOrder61(&phi[index]);
        PDupwindNthSymm2phi = PDupwindNthSymmfdOrder62(&phi[index]);
        PDupwindNthSymm3phi = PDupwindNthSymmfdOrder63(&phi[index]);
        PDupwindNthAnti1phi = PDupwindNthAntifdOrder61(&phi[index]);
        PDupwindNthAnti2phi = PDupwindNthAntifdOrder62(&phi[index]);
        PDupwindNthAnti3phi = PDupwindNthAntifdOrder63(&phi[index]);
        PDdissipationNth1phi = PDdissipationNthfdOrder61(&phi[index]);
        PDdissipationNth2phi = PDdissipationNthfdOrder62(&phi[index]);
        PDdissipationNth3phi = PDdissipationNthfdOrder63(&phi[index]);
        break;
      }
      
      case 8:
      {
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder81(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder82(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder83(&alpha[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder81(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder82(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder83(&alpha[index]);
        PDdissipationNth1alpha = PDdissipationNthfdOrder81(&alpha[index]);
        PDdissipationNth2alpha = PDdissipationNthfdOrder82(&alpha[index]);
        PDdissipationNth3alpha = PDdissipationNthfdOrder83(&alpha[index]);
        PDstandardNth1beta1 = PDstandardNthfdOrder81(&beta1[index]);
        PDstandardNth2beta1 = PDstandardNthfdOrder82(&beta1[index]);
        PDstandardNth3beta1 = PDstandardNthfdOrder83(&beta1[index]);
        PDstandardNth1beta2 = PDstandardNthfdOrder81(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder82(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder83(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder81(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder82(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder83(&beta3[index]);
        PDupwindNthSymm1gt11 = PDupwindNthSymmfdOrder81(&gt11[index]);
        PDupwindNthSymm2gt11 = PDupwindNthSymmfdOrder82(&gt11[index]);
        PDupwindNthSymm3gt11 = PDupwindNthSymmfdOrder83(&gt11[index]);
        PDupwindNthAnti1gt11 = PDupwindNthAntifdOrder81(&gt11[index]);
        PDupwindNthAnti2gt11 = PDupwindNthAntifdOrder82(&gt11[index]);
        PDupwindNthAnti3gt11 = PDupwindNthAntifdOrder83(&gt11[index]);
        PDdissipationNth1gt11 = PDdissipationNthfdOrder81(&gt11[index]);
        PDdissipationNth2gt11 = PDdissipationNthfdOrder82(&gt11[index]);
        PDdissipationNth3gt11 = PDdissipationNthfdOrder83(&gt11[index]);
        PDupwindNthSymm1gt12 = PDupwindNthSymmfdOrder81(&gt12[index]);
        PDupwindNthSymm2gt12 = PDupwindNthSymmfdOrder82(&gt12[index]);
        PDupwindNthSymm3gt12 = PDupwindNthSymmfdOrder83(&gt12[index]);
        PDupwindNthAnti1gt12 = PDupwindNthAntifdOrder81(&gt12[index]);
        PDupwindNthAnti2gt12 = PDupwindNthAntifdOrder82(&gt12[index]);
        PDupwindNthAnti3gt12 = PDupwindNthAntifdOrder83(&gt12[index]);
        PDdissipationNth1gt12 = PDdissipationNthfdOrder81(&gt12[index]);
        PDdissipationNth2gt12 = PDdissipationNthfdOrder82(&gt12[index]);
        PDdissipationNth3gt12 = PDdissipationNthfdOrder83(&gt12[index]);
        PDupwindNthSymm1gt13 = PDupwindNthSymmfdOrder81(&gt13[index]);
        PDupwindNthSymm2gt13 = PDupwindNthSymmfdOrder82(&gt13[index]);
        PDupwindNthSymm3gt13 = PDupwindNthSymmfdOrder83(&gt13[index]);
        PDupwindNthAnti1gt13 = PDupwindNthAntifdOrder81(&gt13[index]);
        PDupwindNthAnti2gt13 = PDupwindNthAntifdOrder82(&gt13[index]);
        PDupwindNthAnti3gt13 = PDupwindNthAntifdOrder83(&gt13[index]);
        PDdissipationNth1gt13 = PDdissipationNthfdOrder81(&gt13[index]);
        PDdissipationNth2gt13 = PDdissipationNthfdOrder82(&gt13[index]);
        PDdissipationNth3gt13 = PDdissipationNthfdOrder83(&gt13[index]);
        PDupwindNthSymm1gt22 = PDupwindNthSymmfdOrder81(&gt22[index]);
        PDupwindNthSymm2gt22 = PDupwindNthSymmfdOrder82(&gt22[index]);
        PDupwindNthSymm3gt22 = PDupwindNthSymmfdOrder83(&gt22[index]);
        PDupwindNthAnti1gt22 = PDupwindNthAntifdOrder81(&gt22[index]);
        PDupwindNthAnti2gt22 = PDupwindNthAntifdOrder82(&gt22[index]);
        PDupwindNthAnti3gt22 = PDupwindNthAntifdOrder83(&gt22[index]);
        PDdissipationNth1gt22 = PDdissipationNthfdOrder81(&gt22[index]);
        PDdissipationNth2gt22 = PDdissipationNthfdOrder82(&gt22[index]);
        PDdissipationNth3gt22 = PDdissipationNthfdOrder83(&gt22[index]);
        PDupwindNthSymm1gt23 = PDupwindNthSymmfdOrder81(&gt23[index]);
        PDupwindNthSymm2gt23 = PDupwindNthSymmfdOrder82(&gt23[index]);
        PDupwindNthSymm3gt23 = PDupwindNthSymmfdOrder83(&gt23[index]);
        PDupwindNthAnti1gt23 = PDupwindNthAntifdOrder81(&gt23[index]);
        PDupwindNthAnti2gt23 = PDupwindNthAntifdOrder82(&gt23[index]);
        PDupwindNthAnti3gt23 = PDupwindNthAntifdOrder83(&gt23[index]);
        PDdissipationNth1gt23 = PDdissipationNthfdOrder81(&gt23[index]);
        PDdissipationNth2gt23 = PDdissipationNthfdOrder82(&gt23[index]);
        PDdissipationNth3gt23 = PDdissipationNthfdOrder83(&gt23[index]);
        PDupwindNthSymm1gt33 = PDupwindNthSymmfdOrder81(&gt33[index]);
        PDupwindNthSymm2gt33 = PDupwindNthSymmfdOrder82(&gt33[index]);
        PDupwindNthSymm3gt33 = PDupwindNthSymmfdOrder83(&gt33[index]);
        PDupwindNthAnti1gt33 = PDupwindNthAntifdOrder81(&gt33[index]);
        PDupwindNthAnti2gt33 = PDupwindNthAntifdOrder82(&gt33[index]);
        PDupwindNthAnti3gt33 = PDupwindNthAntifdOrder83(&gt33[index]);
        PDdissipationNth1gt33 = PDdissipationNthfdOrder81(&gt33[index]);
        PDdissipationNth2gt33 = PDdissipationNthfdOrder82(&gt33[index]);
        PDdissipationNth3gt33 = PDdissipationNthfdOrder83(&gt33[index]);
        PDupwindNthSymm1phi = PDupwindNthSymmfdOrder81(&phi[index]);
        PDupwindNthSymm2phi = PDupwindNthSymmfdOrder82(&phi[index]);
        PDupwindNthSymm3phi = PDupwindNthSymmfdOrder83(&phi[index]);
        PDupwindNthAnti1phi = PDupwindNthAntifdOrder81(&phi[index]);
        PDupwindNthAnti2phi = PDupwindNthAntifdOrder82(&phi[index]);
        PDupwindNthAnti3phi = PDupwindNthAntifdOrder83(&phi[index]);
        PDdissipationNth1phi = PDdissipationNthfdOrder81(&phi[index]);
        PDdissipationNth2phi = PDdissipationNthfdOrder82(&phi[index]);
        PDdissipationNth3phi = PDdissipationNthfdOrder83(&phi[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC JacPDdissipationNth1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDdissipationNth3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3phi CCTK_ATTRIBUTE_UNUSED;
    
    if (usejacobian)
    {
      JacPDstandardNth1beta1 = 
        kmadd(J11L,PDstandardNth1beta1,kmadd(J21L,PDstandardNth2beta1,kmul(J31L,PDstandardNth3beta1)));
      
      JacPDstandardNth1beta2 = 
        kmadd(J11L,PDstandardNth1beta2,kmadd(J21L,PDstandardNth2beta2,kmul(J31L,PDstandardNth3beta2)));
      
      JacPDstandardNth1beta3 = 
        kmadd(J11L,PDstandardNth1beta3,kmadd(J21L,PDstandardNth2beta3,kmul(J31L,PDstandardNth3beta3)));
      
      JacPDstandardNth2beta1 = 
        kmadd(J12L,PDstandardNth1beta1,kmadd(J22L,PDstandardNth2beta1,kmul(J32L,PDstandardNth3beta1)));
      
      JacPDstandardNth2beta2 = 
        kmadd(J12L,PDstandardNth1beta2,kmadd(J22L,PDstandardNth2beta2,kmul(J32L,PDstandardNth3beta2)));
      
      JacPDstandardNth2beta3 = 
        kmadd(J12L,PDstandardNth1beta3,kmadd(J22L,PDstandardNth2beta3,kmul(J32L,PDstandardNth3beta3)));
      
      JacPDstandardNth3beta1 = 
        kmadd(J13L,PDstandardNth1beta1,kmadd(J23L,PDstandardNth2beta1,kmul(J33L,PDstandardNth3beta1)));
      
      JacPDstandardNth3beta2 = 
        kmadd(J13L,PDstandardNth1beta2,kmadd(J23L,PDstandardNth2beta2,kmul(J33L,PDstandardNth3beta2)));
      
      JacPDstandardNth3beta3 = 
        kmadd(J13L,PDstandardNth1beta3,kmadd(J23L,PDstandardNth2beta3,kmul(J33L,PDstandardNth3beta3)));
      
      JacPDupwindNthSymm1alpha = 
        kmadd(J11L,PDupwindNthSymm1alpha,kmadd(J21L,PDupwindNthSymm2alpha,kmul(J31L,PDupwindNthSymm3alpha)));
      
      JacPDupwindNthSymm1gt11 = 
        kmadd(J11L,PDupwindNthSymm1gt11,kmadd(J21L,PDupwindNthSymm2gt11,kmul(J31L,PDupwindNthSymm3gt11)));
      
      JacPDupwindNthSymm1gt12 = 
        kmadd(J11L,PDupwindNthSymm1gt12,kmadd(J21L,PDupwindNthSymm2gt12,kmul(J31L,PDupwindNthSymm3gt12)));
      
      JacPDupwindNthSymm1gt13 = 
        kmadd(J11L,PDupwindNthSymm1gt13,kmadd(J21L,PDupwindNthSymm2gt13,kmul(J31L,PDupwindNthSymm3gt13)));
      
      JacPDupwindNthSymm1gt22 = 
        kmadd(J11L,PDupwindNthSymm1gt22,kmadd(J21L,PDupwindNthSymm2gt22,kmul(J31L,PDupwindNthSymm3gt22)));
      
      JacPDupwindNthSymm1gt23 = 
        kmadd(J11L,PDupwindNthSymm1gt23,kmadd(J21L,PDupwindNthSymm2gt23,kmul(J31L,PDupwindNthSymm3gt23)));
      
      JacPDupwindNthSymm1gt33 = 
        kmadd(J11L,PDupwindNthSymm1gt33,kmadd(J21L,PDupwindNthSymm2gt33,kmul(J31L,PDupwindNthSymm3gt33)));
      
      JacPDupwindNthSymm1phi = 
        kmadd(J11L,PDupwindNthSymm1phi,kmadd(J21L,PDupwindNthSymm2phi,kmul(J31L,PDupwindNthSymm3phi)));
      
      JacPDupwindNthSymm2alpha = 
        kmadd(J12L,PDupwindNthSymm1alpha,kmadd(J22L,PDupwindNthSymm2alpha,kmul(J32L,PDupwindNthSymm3alpha)));
      
      JacPDupwindNthSymm2gt11 = 
        kmadd(J12L,PDupwindNthSymm1gt11,kmadd(J22L,PDupwindNthSymm2gt11,kmul(J32L,PDupwindNthSymm3gt11)));
      
      JacPDupwindNthSymm2gt12 = 
        kmadd(J12L,PDupwindNthSymm1gt12,kmadd(J22L,PDupwindNthSymm2gt12,kmul(J32L,PDupwindNthSymm3gt12)));
      
      JacPDupwindNthSymm2gt13 = 
        kmadd(J12L,PDupwindNthSymm1gt13,kmadd(J22L,PDupwindNthSymm2gt13,kmul(J32L,PDupwindNthSymm3gt13)));
      
      JacPDupwindNthSymm2gt22 = 
        kmadd(J12L,PDupwindNthSymm1gt22,kmadd(J22L,PDupwindNthSymm2gt22,kmul(J32L,PDupwindNthSymm3gt22)));
      
      JacPDupwindNthSymm2gt23 = 
        kmadd(J12L,PDupwindNthSymm1gt23,kmadd(J22L,PDupwindNthSymm2gt23,kmul(J32L,PDupwindNthSymm3gt23)));
      
      JacPDupwindNthSymm2gt33 = 
        kmadd(J12L,PDupwindNthSymm1gt33,kmadd(J22L,PDupwindNthSymm2gt33,kmul(J32L,PDupwindNthSymm3gt33)));
      
      JacPDupwindNthSymm2phi = 
        kmadd(J12L,PDupwindNthSymm1phi,kmadd(J22L,PDupwindNthSymm2phi,kmul(J32L,PDupwindNthSymm3phi)));
      
      JacPDupwindNthSymm3alpha = 
        kmadd(J13L,PDupwindNthSymm1alpha,kmadd(J23L,PDupwindNthSymm2alpha,kmul(J33L,PDupwindNthSymm3alpha)));
      
      JacPDupwindNthSymm3gt11 = 
        kmadd(J13L,PDupwindNthSymm1gt11,kmadd(J23L,PDupwindNthSymm2gt11,kmul(J33L,PDupwindNthSymm3gt11)));
      
      JacPDupwindNthSymm3gt12 = 
        kmadd(J13L,PDupwindNthSymm1gt12,kmadd(J23L,PDupwindNthSymm2gt12,kmul(J33L,PDupwindNthSymm3gt12)));
      
      JacPDupwindNthSymm3gt13 = 
        kmadd(J13L,PDupwindNthSymm1gt13,kmadd(J23L,PDupwindNthSymm2gt13,kmul(J33L,PDupwindNthSymm3gt13)));
      
      JacPDupwindNthSymm3gt22 = 
        kmadd(J13L,PDupwindNthSymm1gt22,kmadd(J23L,PDupwindNthSymm2gt22,kmul(J33L,PDupwindNthSymm3gt22)));
      
      JacPDupwindNthSymm3gt23 = 
        kmadd(J13L,PDupwindNthSymm1gt23,kmadd(J23L,PDupwindNthSymm2gt23,kmul(J33L,PDupwindNthSymm3gt23)));
      
      JacPDupwindNthSymm3gt33 = 
        kmadd(J13L,PDupwindNthSymm1gt33,kmadd(J23L,PDupwindNthSymm2gt33,kmul(J33L,PDupwindNthSymm3gt33)));
      
      JacPDupwindNthSymm3phi = 
        kmadd(J13L,PDupwindNthSymm1phi,kmadd(J23L,PDupwindNthSymm2phi,kmul(J33L,PDupwindNthSymm3phi)));
      
      JacPDupwindNthAnti1alpha = 
        kmadd(J11L,PDupwindNthAnti1alpha,kmadd(J21L,PDupwindNthAnti2alpha,kmul(J31L,PDupwindNthAnti3alpha)));
      
      JacPDupwindNthAnti1gt11 = 
        kmadd(J11L,PDupwindNthAnti1gt11,kmadd(J21L,PDupwindNthAnti2gt11,kmul(J31L,PDupwindNthAnti3gt11)));
      
      JacPDupwindNthAnti1gt12 = 
        kmadd(J11L,PDupwindNthAnti1gt12,kmadd(J21L,PDupwindNthAnti2gt12,kmul(J31L,PDupwindNthAnti3gt12)));
      
      JacPDupwindNthAnti1gt13 = 
        kmadd(J11L,PDupwindNthAnti1gt13,kmadd(J21L,PDupwindNthAnti2gt13,kmul(J31L,PDupwindNthAnti3gt13)));
      
      JacPDupwindNthAnti1gt22 = 
        kmadd(J11L,PDupwindNthAnti1gt22,kmadd(J21L,PDupwindNthAnti2gt22,kmul(J31L,PDupwindNthAnti3gt22)));
      
      JacPDupwindNthAnti1gt23 = 
        kmadd(J11L,PDupwindNthAnti1gt23,kmadd(J21L,PDupwindNthAnti2gt23,kmul(J31L,PDupwindNthAnti3gt23)));
      
      JacPDupwindNthAnti1gt33 = 
        kmadd(J11L,PDupwindNthAnti1gt33,kmadd(J21L,PDupwindNthAnti2gt33,kmul(J31L,PDupwindNthAnti3gt33)));
      
      JacPDupwindNthAnti1phi = 
        kmadd(J11L,PDupwindNthAnti1phi,kmadd(J21L,PDupwindNthAnti2phi,kmul(J31L,PDupwindNthAnti3phi)));
      
      JacPDupwindNthAnti2alpha = 
        kmadd(J12L,PDupwindNthAnti1alpha,kmadd(J22L,PDupwindNthAnti2alpha,kmul(J32L,PDupwindNthAnti3alpha)));
      
      JacPDupwindNthAnti2gt11 = 
        kmadd(J12L,PDupwindNthAnti1gt11,kmadd(J22L,PDupwindNthAnti2gt11,kmul(J32L,PDupwindNthAnti3gt11)));
      
      JacPDupwindNthAnti2gt12 = 
        kmadd(J12L,PDupwindNthAnti1gt12,kmadd(J22L,PDupwindNthAnti2gt12,kmul(J32L,PDupwindNthAnti3gt12)));
      
      JacPDupwindNthAnti2gt13 = 
        kmadd(J12L,PDupwindNthAnti1gt13,kmadd(J22L,PDupwindNthAnti2gt13,kmul(J32L,PDupwindNthAnti3gt13)));
      
      JacPDupwindNthAnti2gt22 = 
        kmadd(J12L,PDupwindNthAnti1gt22,kmadd(J22L,PDupwindNthAnti2gt22,kmul(J32L,PDupwindNthAnti3gt22)));
      
      JacPDupwindNthAnti2gt23 = 
        kmadd(J12L,PDupwindNthAnti1gt23,kmadd(J22L,PDupwindNthAnti2gt23,kmul(J32L,PDupwindNthAnti3gt23)));
      
      JacPDupwindNthAnti2gt33 = 
        kmadd(J12L,PDupwindNthAnti1gt33,kmadd(J22L,PDupwindNthAnti2gt33,kmul(J32L,PDupwindNthAnti3gt33)));
      
      JacPDupwindNthAnti2phi = 
        kmadd(J12L,PDupwindNthAnti1phi,kmadd(J22L,PDupwindNthAnti2phi,kmul(J32L,PDupwindNthAnti3phi)));
      
      JacPDupwindNthAnti3alpha = 
        kmadd(J13L,PDupwindNthAnti1alpha,kmadd(J23L,PDupwindNthAnti2alpha,kmul(J33L,PDupwindNthAnti3alpha)));
      
      JacPDupwindNthAnti3gt11 = 
        kmadd(J13L,PDupwindNthAnti1gt11,kmadd(J23L,PDupwindNthAnti2gt11,kmul(J33L,PDupwindNthAnti3gt11)));
      
      JacPDupwindNthAnti3gt12 = 
        kmadd(J13L,PDupwindNthAnti1gt12,kmadd(J23L,PDupwindNthAnti2gt12,kmul(J33L,PDupwindNthAnti3gt12)));
      
      JacPDupwindNthAnti3gt13 = 
        kmadd(J13L,PDupwindNthAnti1gt13,kmadd(J23L,PDupwindNthAnti2gt13,kmul(J33L,PDupwindNthAnti3gt13)));
      
      JacPDupwindNthAnti3gt22 = 
        kmadd(J13L,PDupwindNthAnti1gt22,kmadd(J23L,PDupwindNthAnti2gt22,kmul(J33L,PDupwindNthAnti3gt22)));
      
      JacPDupwindNthAnti3gt23 = 
        kmadd(J13L,PDupwindNthAnti1gt23,kmadd(J23L,PDupwindNthAnti2gt23,kmul(J33L,PDupwindNthAnti3gt23)));
      
      JacPDupwindNthAnti3gt33 = 
        kmadd(J13L,PDupwindNthAnti1gt33,kmadd(J23L,PDupwindNthAnti2gt33,kmul(J33L,PDupwindNthAnti3gt33)));
      
      JacPDupwindNthAnti3phi = 
        kmadd(J13L,PDupwindNthAnti1phi,kmadd(J23L,PDupwindNthAnti2phi,kmul(J33L,PDupwindNthAnti3phi)));
      
      JacPDdissipationNth1alpha = 
        kmadd(J11L,PDdissipationNth1alpha,kmadd(J21L,PDdissipationNth2alpha,kmul(J31L,PDdissipationNth3alpha)));
      
      JacPDdissipationNth1gt11 = 
        kmadd(J11L,PDdissipationNth1gt11,kmadd(J21L,PDdissipationNth2gt11,kmul(J31L,PDdissipationNth3gt11)));
      
      JacPDdissipationNth1gt12 = 
        kmadd(J11L,PDdissipationNth1gt12,kmadd(J21L,PDdissipationNth2gt12,kmul(J31L,PDdissipationNth3gt12)));
      
      JacPDdissipationNth1gt13 = 
        kmadd(J11L,PDdissipationNth1gt13,kmadd(J21L,PDdissipationNth2gt13,kmul(J31L,PDdissipationNth3gt13)));
      
      JacPDdissipationNth1gt22 = 
        kmadd(J11L,PDdissipationNth1gt22,kmadd(J21L,PDdissipationNth2gt22,kmul(J31L,PDdissipationNth3gt22)));
      
      JacPDdissipationNth1gt23 = 
        kmadd(J11L,PDdissipationNth1gt23,kmadd(J21L,PDdissipationNth2gt23,kmul(J31L,PDdissipationNth3gt23)));
      
      JacPDdissipationNth1gt33 = 
        kmadd(J11L,PDdissipationNth1gt33,kmadd(J21L,PDdissipationNth2gt33,kmul(J31L,PDdissipationNth3gt33)));
      
      JacPDdissipationNth1phi = 
        kmadd(J11L,PDdissipationNth1phi,kmadd(J21L,PDdissipationNth2phi,kmul(J31L,PDdissipationNth3phi)));
      
      JacPDdissipationNth2alpha = 
        kmadd(J12L,PDdissipationNth1alpha,kmadd(J22L,PDdissipationNth2alpha,kmul(J32L,PDdissipationNth3alpha)));
      
      JacPDdissipationNth2gt11 = 
        kmadd(J12L,PDdissipationNth1gt11,kmadd(J22L,PDdissipationNth2gt11,kmul(J32L,PDdissipationNth3gt11)));
      
      JacPDdissipationNth2gt12 = 
        kmadd(J12L,PDdissipationNth1gt12,kmadd(J22L,PDdissipationNth2gt12,kmul(J32L,PDdissipationNth3gt12)));
      
      JacPDdissipationNth2gt13 = 
        kmadd(J12L,PDdissipationNth1gt13,kmadd(J22L,PDdissipationNth2gt13,kmul(J32L,PDdissipationNth3gt13)));
      
      JacPDdissipationNth2gt22 = 
        kmadd(J12L,PDdissipationNth1gt22,kmadd(J22L,PDdissipationNth2gt22,kmul(J32L,PDdissipationNth3gt22)));
      
      JacPDdissipationNth2gt23 = 
        kmadd(J12L,PDdissipationNth1gt23,kmadd(J22L,PDdissipationNth2gt23,kmul(J32L,PDdissipationNth3gt23)));
      
      JacPDdissipationNth2gt33 = 
        kmadd(J12L,PDdissipationNth1gt33,kmadd(J22L,PDdissipationNth2gt33,kmul(J32L,PDdissipationNth3gt33)));
      
      JacPDdissipationNth2phi = 
        kmadd(J12L,PDdissipationNth1phi,kmadd(J22L,PDdissipationNth2phi,kmul(J32L,PDdissipationNth3phi)));
      
      JacPDdissipationNth3alpha = 
        kmadd(J13L,PDdissipationNth1alpha,kmadd(J23L,PDdissipationNth2alpha,kmul(J33L,PDdissipationNth3alpha)));
      
      JacPDdissipationNth3gt11 = 
        kmadd(J13L,PDdissipationNth1gt11,kmadd(J23L,PDdissipationNth2gt11,kmul(J33L,PDdissipationNth3gt11)));
      
      JacPDdissipationNth3gt12 = 
        kmadd(J13L,PDdissipationNth1gt12,kmadd(J23L,PDdissipationNth2gt12,kmul(J33L,PDdissipationNth3gt12)));
      
      JacPDdissipationNth3gt13 = 
        kmadd(J13L,PDdissipationNth1gt13,kmadd(J23L,PDdissipationNth2gt13,kmul(J33L,PDdissipationNth3gt13)));
      
      JacPDdissipationNth3gt22 = 
        kmadd(J13L,PDdissipationNth1gt22,kmadd(J23L,PDdissipationNth2gt22,kmul(J33L,PDdissipationNth3gt22)));
      
      JacPDdissipationNth3gt23 = 
        kmadd(J13L,PDdissipationNth1gt23,kmadd(J23L,PDdissipationNth2gt23,kmul(J33L,PDdissipationNth3gt23)));
      
      JacPDdissipationNth3gt33 = 
        kmadd(J13L,PDdissipationNth1gt33,kmadd(J23L,PDdissipationNth2gt33,kmul(J33L,PDdissipationNth3gt33)));
      
      JacPDdissipationNth3phi = 
        kmadd(J13L,PDdissipationNth1phi,kmadd(J23L,PDdissipationNth2phi,kmul(J33L,PDdissipationNth3phi)));
    }
    else
    {
      JacPDstandardNth1beta1 = PDstandardNth1beta1;
      
      JacPDstandardNth1beta2 = PDstandardNth1beta2;
      
      JacPDstandardNth1beta3 = PDstandardNth1beta3;
      
      JacPDstandardNth2beta1 = PDstandardNth2beta1;
      
      JacPDstandardNth2beta2 = PDstandardNth2beta2;
      
      JacPDstandardNth2beta3 = PDstandardNth2beta3;
      
      JacPDstandardNth3beta1 = PDstandardNth3beta1;
      
      JacPDstandardNth3beta2 = PDstandardNth3beta2;
      
      JacPDstandardNth3beta3 = PDstandardNth3beta3;
      
      JacPDupwindNthSymm1alpha = PDupwindNthSymm1alpha;
      
      JacPDupwindNthSymm1gt11 = PDupwindNthSymm1gt11;
      
      JacPDupwindNthSymm1gt12 = PDupwindNthSymm1gt12;
      
      JacPDupwindNthSymm1gt13 = PDupwindNthSymm1gt13;
      
      JacPDupwindNthSymm1gt22 = PDupwindNthSymm1gt22;
      
      JacPDupwindNthSymm1gt23 = PDupwindNthSymm1gt23;
      
      JacPDupwindNthSymm1gt33 = PDupwindNthSymm1gt33;
      
      JacPDupwindNthSymm1phi = PDupwindNthSymm1phi;
      
      JacPDupwindNthSymm2alpha = PDupwindNthSymm2alpha;
      
      JacPDupwindNthSymm2gt11 = PDupwindNthSymm2gt11;
      
      JacPDupwindNthSymm2gt12 = PDupwindNthSymm2gt12;
      
      JacPDupwindNthSymm2gt13 = PDupwindNthSymm2gt13;
      
      JacPDupwindNthSymm2gt22 = PDupwindNthSymm2gt22;
      
      JacPDupwindNthSymm2gt23 = PDupwindNthSymm2gt23;
      
      JacPDupwindNthSymm2gt33 = PDupwindNthSymm2gt33;
      
      JacPDupwindNthSymm2phi = PDupwindNthSymm2phi;
      
      JacPDupwindNthSymm3alpha = PDupwindNthSymm3alpha;
      
      JacPDupwindNthSymm3gt11 = PDupwindNthSymm3gt11;
      
      JacPDupwindNthSymm3gt12 = PDupwindNthSymm3gt12;
      
      JacPDupwindNthSymm3gt13 = PDupwindNthSymm3gt13;
      
      JacPDupwindNthSymm3gt22 = PDupwindNthSymm3gt22;
      
      JacPDupwindNthSymm3gt23 = PDupwindNthSymm3gt23;
      
      JacPDupwindNthSymm3gt33 = PDupwindNthSymm3gt33;
      
      JacPDupwindNthSymm3phi = PDupwindNthSymm3phi;
      
      JacPDupwindNthAnti1alpha = PDupwindNthAnti1alpha;
      
      JacPDupwindNthAnti1gt11 = PDupwindNthAnti1gt11;
      
      JacPDupwindNthAnti1gt12 = PDupwindNthAnti1gt12;
      
      JacPDupwindNthAnti1gt13 = PDupwindNthAnti1gt13;
      
      JacPDupwindNthAnti1gt22 = PDupwindNthAnti1gt22;
      
      JacPDupwindNthAnti1gt23 = PDupwindNthAnti1gt23;
      
      JacPDupwindNthAnti1gt33 = PDupwindNthAnti1gt33;
      
      JacPDupwindNthAnti1phi = PDupwindNthAnti1phi;
      
      JacPDupwindNthAnti2alpha = PDupwindNthAnti2alpha;
      
      JacPDupwindNthAnti2gt11 = PDupwindNthAnti2gt11;
      
      JacPDupwindNthAnti2gt12 = PDupwindNthAnti2gt12;
      
      JacPDupwindNthAnti2gt13 = PDupwindNthAnti2gt13;
      
      JacPDupwindNthAnti2gt22 = PDupwindNthAnti2gt22;
      
      JacPDupwindNthAnti2gt23 = PDupwindNthAnti2gt23;
      
      JacPDupwindNthAnti2gt33 = PDupwindNthAnti2gt33;
      
      JacPDupwindNthAnti2phi = PDupwindNthAnti2phi;
      
      JacPDupwindNthAnti3alpha = PDupwindNthAnti3alpha;
      
      JacPDupwindNthAnti3gt11 = PDupwindNthAnti3gt11;
      
      JacPDupwindNthAnti3gt12 = PDupwindNthAnti3gt12;
      
      JacPDupwindNthAnti3gt13 = PDupwindNthAnti3gt13;
      
      JacPDupwindNthAnti3gt22 = PDupwindNthAnti3gt22;
      
      JacPDupwindNthAnti3gt23 = PDupwindNthAnti3gt23;
      
      JacPDupwindNthAnti3gt33 = PDupwindNthAnti3gt33;
      
      JacPDupwindNthAnti3phi = PDupwindNthAnti3phi;
      
      JacPDdissipationNth1alpha = PDdissipationNth1alpha;
      
      JacPDdissipationNth1gt11 = PDdissipationNth1gt11;
      
      JacPDdissipationNth1gt12 = PDdissipationNth1gt12;
      
      JacPDdissipationNth1gt13 = PDdissipationNth1gt13;
      
      JacPDdissipationNth1gt22 = PDdissipationNth1gt22;
      
      JacPDdissipationNth1gt23 = PDdissipationNth1gt23;
      
      JacPDdissipationNth1gt33 = PDdissipationNth1gt33;
      
      JacPDdissipationNth1phi = PDdissipationNth1phi;
      
      JacPDdissipationNth2alpha = PDdissipationNth2alpha;
      
      JacPDdissipationNth2gt11 = PDdissipationNth2gt11;
      
      JacPDdissipationNth2gt12 = PDdissipationNth2gt12;
      
      JacPDdissipationNth2gt13 = PDdissipationNth2gt13;
      
      JacPDdissipationNth2gt22 = PDdissipationNth2gt22;
      
      JacPDdissipationNth2gt23 = PDdissipationNth2gt23;
      
      JacPDdissipationNth2gt33 = PDdissipationNth2gt33;
      
      JacPDdissipationNth2phi = PDdissipationNth2phi;
      
      JacPDdissipationNth3alpha = PDdissipationNth3alpha;
      
      JacPDdissipationNth3gt11 = PDdissipationNth3gt11;
      
      JacPDdissipationNth3gt12 = PDdissipationNth3gt12;
      
      JacPDdissipationNth3gt13 = PDdissipationNth3gt13;
      
      JacPDdissipationNth3gt22 = PDdissipationNth3gt22;
      
      JacPDdissipationNth3gt23 = PDdissipationNth3gt23;
      
      JacPDdissipationNth3gt33 = PDdissipationNth3gt33;
      
      JacPDdissipationNth3phi = PDdissipationNth3phi;
    }
    
    CCTK_REAL_VEC epsdiss1 CCTK_ATTRIBUTE_UNUSED = ToReal(epsDiss);
    
    CCTK_REAL_VEC epsdiss2 CCTK_ATTRIBUTE_UNUSED = ToReal(epsDiss);
    
    CCTK_REAL_VEC epsdiss3 CCTK_ATTRIBUTE_UNUSED = ToReal(epsDiss);
    
    CCTK_REAL_VEC phirhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(epsdiss1,JacPDdissipationNth1phi,kmadd(epsdiss2,JacPDdissipationNth2phi,kmadd(epsdiss3,JacPDdissipationNth3phi,kmadd(beta1L,JacPDupwindNthAnti1phi,kmadd(beta2L,JacPDupwindNthAnti2phi,kmadd(beta3L,JacPDupwindNthAnti3phi,kmadd(JacPDupwindNthSymm1phi,kfabs(beta1L),kmadd(JacPDupwindNthSymm2phi,kfabs(beta2L),kmsub(JacPDupwindNthSymm3phi,kfabs(beta3L),kmul(IfThen(conformalMethod 
      != 
      0,kmul(ToReal(0.333333333333333333333333333333),phiL),ToReal(-0.166666666666666666666666666667)),knmsub(alphaL,trKL,kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3)))))))))))));
    
    CCTK_REAL_VEC gt11rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmul(alphaL,At11L),kmadd(epsdiss1,JacPDdissipationNth1gt11,kmadd(epsdiss2,JacPDdissipationNth2gt11,kmadd(epsdiss3,JacPDdissipationNth3gt11,kmadd(ToReal(2),kmadd(gt11L,JacPDstandardNth1beta1,kmadd(gt12L,JacPDstandardNth1beta2,kmul(gt13L,JacPDstandardNth1beta3))),kmadd(ToReal(-0.666666666666666666666666666667),kmul(gt11L,kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3))),kmadd(beta1L,JacPDupwindNthAnti1gt11,kmadd(beta2L,JacPDupwindNthAnti2gt11,kmadd(beta3L,JacPDupwindNthAnti3gt11,kmadd(JacPDupwindNthSymm1gt11,kfabs(beta1L),kmadd(JacPDupwindNthSymm2gt11,kfabs(beta2L),kmul(JacPDupwindNthSymm3gt11,kfabs(beta3L)))))))))))));
    
    CCTK_REAL_VEC gt12rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmul(alphaL,At12L),kmadd(epsdiss1,JacPDdissipationNth1gt12,kmadd(epsdiss2,JacPDdissipationNth2gt12,kmadd(epsdiss3,JacPDdissipationNth3gt12,kmadd(gt22L,JacPDstandardNth1beta2,kmadd(gt23L,JacPDstandardNth1beta3,kmadd(gt11L,JacPDstandardNth2beta1,kmadd(gt13L,JacPDstandardNth2beta3,kmadd(gt12L,kadd(JacPDstandardNth1beta1,kmadd(ToReal(-0.666666666666666666666666666667),kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3)),JacPDstandardNth2beta2)),kmadd(beta1L,JacPDupwindNthAnti1gt12,kmadd(beta2L,JacPDupwindNthAnti2gt12,kmadd(beta3L,JacPDupwindNthAnti3gt12,kmadd(JacPDupwindNthSymm1gt12,kfabs(beta1L),kmadd(JacPDupwindNthSymm2gt12,kfabs(beta2L),kmul(JacPDupwindNthSymm3gt12,kfabs(beta3L))))))))))))))));
    
    CCTK_REAL_VEC gt13rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmul(alphaL,At13L),kmadd(epsdiss1,JacPDdissipationNth1gt13,kmadd(epsdiss2,JacPDdissipationNth2gt13,kmadd(epsdiss3,JacPDdissipationNth3gt13,kmadd(gt23L,JacPDstandardNth1beta2,kmadd(gt33L,JacPDstandardNth1beta3,kmadd(gt11L,JacPDstandardNth3beta1,kmadd(gt12L,JacPDstandardNth3beta2,kmadd(gt13L,kadd(JacPDstandardNth1beta1,kmadd(ToReal(-0.666666666666666666666666666667),kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3)),JacPDstandardNth3beta3)),kmadd(beta1L,JacPDupwindNthAnti1gt13,kmadd(beta2L,JacPDupwindNthAnti2gt13,kmadd(beta3L,JacPDupwindNthAnti3gt13,kmadd(JacPDupwindNthSymm1gt13,kfabs(beta1L),kmadd(JacPDupwindNthSymm2gt13,kfabs(beta2L),kmul(JacPDupwindNthSymm3gt13,kfabs(beta3L))))))))))))))));
    
    CCTK_REAL_VEC gt22rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmul(alphaL,At22L),kmadd(epsdiss1,JacPDdissipationNth1gt22,kmadd(epsdiss2,JacPDdissipationNth2gt22,kmadd(epsdiss3,JacPDdissipationNth3gt22,kmadd(ToReal(2),kmadd(gt12L,JacPDstandardNth2beta1,kmadd(gt22L,JacPDstandardNth2beta2,kmul(gt23L,JacPDstandardNth2beta3))),kmadd(ToReal(-0.666666666666666666666666666667),kmul(gt22L,kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3))),kmadd(beta1L,JacPDupwindNthAnti1gt22,kmadd(beta2L,JacPDupwindNthAnti2gt22,kmadd(beta3L,JacPDupwindNthAnti3gt22,kmadd(JacPDupwindNthSymm1gt22,kfabs(beta1L),kmadd(JacPDupwindNthSymm2gt22,kfabs(beta2L),kmul(JacPDupwindNthSymm3gt22,kfabs(beta3L)))))))))))));
    
    CCTK_REAL_VEC gt23rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmul(alphaL,At23L),kmadd(epsdiss1,JacPDdissipationNth1gt23,kmadd(epsdiss2,JacPDdissipationNth2gt23,kmadd(epsdiss3,JacPDdissipationNth3gt23,kmadd(gt13L,JacPDstandardNth2beta1,kmadd(gt33L,JacPDstandardNth2beta3,kmadd(gt12L,JacPDstandardNth3beta1,kmadd(gt22L,JacPDstandardNth3beta2,kmadd(gt23L,kadd(JacPDstandardNth2beta2,kmadd(ToReal(-0.666666666666666666666666666667),kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3)),JacPDstandardNth3beta3)),kmadd(beta1L,JacPDupwindNthAnti1gt23,kmadd(beta2L,JacPDupwindNthAnti2gt23,kmadd(beta3L,JacPDupwindNthAnti3gt23,kmadd(JacPDupwindNthSymm1gt23,kfabs(beta1L),kmadd(JacPDupwindNthSymm2gt23,kfabs(beta2L),kmul(JacPDupwindNthSymm3gt23,kfabs(beta3L))))))))))))))));
    
    CCTK_REAL_VEC gt33rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmul(alphaL,At33L),kmadd(epsdiss1,JacPDdissipationNth1gt33,kmadd(epsdiss2,JacPDdissipationNth2gt33,kmadd(epsdiss3,JacPDdissipationNth3gt33,kmadd(ToReal(-0.666666666666666666666666666667),kmul(gt33L,kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3))),kmadd(ToReal(2),kmadd(gt13L,JacPDstandardNth3beta1,kmadd(gt23L,JacPDstandardNth3beta2,kmul(gt33L,JacPDstandardNth3beta3))),kmadd(beta1L,JacPDupwindNthAnti1gt33,kmadd(beta2L,JacPDupwindNthAnti2gt33,kmadd(beta3L,JacPDupwindNthAnti3gt33,kmadd(JacPDupwindNthSymm1gt33,kfabs(beta1L),kmadd(JacPDupwindNthSymm2gt33,kfabs(beta2L),kmul(JacPDupwindNthSymm3gt33,kfabs(beta3L)))))))))))));
    
    CCTK_REAL_VEC dotalpha CCTK_ATTRIBUTE_UNUSED = 
      kneg(kmul(kmul(kmadd(ToReal(alphaDriver),kadd(ToReal(-1),alphaL),trKL),kpow(alphaL,harmonicN)),ToReal(harmonicF)));
    
    CCTK_REAL_VEC alpharhsL CCTK_ATTRIBUTE_UNUSED = 
      kadd(dotalpha,kmadd(epsdiss1,JacPDdissipationNth1alpha,kmadd(epsdiss2,JacPDdissipationNth2alpha,kmadd(epsdiss3,JacPDdissipationNth3alpha,IfThen(advectLapse 
      != 
      0,kmadd(beta1L,JacPDupwindNthAnti1alpha,kmadd(beta2L,JacPDupwindNthAnti2alpha,kmadd(beta3L,JacPDupwindNthAnti3alpha,kmadd(JacPDupwindNthSymm1alpha,kfabs(beta1L),kmadd(JacPDupwindNthSymm2alpha,kfabs(beta2L),kmul(JacPDupwindNthSymm3alpha,kfabs(beta3L))))))),ToReal(0))))));
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(alpharhs[index],alpharhsL);
    vec_store_nta_partial(gt11rhs[index],gt11rhsL);
    vec_store_nta_partial(gt12rhs[index],gt12rhsL);
    vec_store_nta_partial(gt13rhs[index],gt13rhsL);
    vec_store_nta_partial(gt22rhs[index],gt22rhsL);
    vec_store_nta_partial(gt23rhs[index],gt23rhsL);
    vec_store_nta_partial(gt33rhs[index],gt33rhsL);
    vec_store_nta_partial(phirhs[index],phirhsL);
  }
  CCTK_ENDLOOP3STR(ML_BSSN_bench_EvolutionInteriorSplitBy1);
}
extern "C" void ML_BSSN_bench_EvolutionInteriorSplitBy1(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_BSSN_bench_EvolutionInteriorSplitBy1
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_bench_EvolutionInteriorSplitBy1);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_bench_EvolutionInteriorSplitBy1_Body");
  }
  if (cctk_iteration % ML_BSSN_bench_EvolutionInteriorSplitBy1_calc_every != ML_BSSN_bench_EvolutionInteriorSplitBy1_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ML_BSSN_bench::ML_curv",
    "ML_BSSN_bench::ML_lapse",
    "ML_BSSN_bench::ML_lapserhs",
    "ML_BSSN_bench::ML_log_confac",
    "ML_BSSN_bench::ML_log_confacrhs",
    "ML_BSSN_bench::ML_metric",
    "ML_BSSN_bench::ML_metricrhs",
    "ML_BSSN_bench::ML_shift",
    "ML_BSSN_bench::ML_trace_curv"};
  AssertGroupStorage(cctkGH, "ML_BSSN_bench_EvolutionInteriorSplitBy1", 9, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_bench_EvolutionInteriorSplitBy1", 2, 2, 2);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_bench_EvolutionInteriorSplitBy1", 3, 3, 3);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_bench_EvolutionInteriorSplitBy1", 4, 4, 4);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_bench_EvolutionInteriorSplitBy1", 5, 5, 5);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, ML_BSSN_bench_EvolutionInteriorSplitBy1_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_bench_EvolutionInteriorSplitBy1_Body");
  }
}

} // namespace ML_BSSN_bench
