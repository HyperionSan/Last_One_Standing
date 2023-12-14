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

namespace ML_ADM {

extern "C" void ML_ADM_RHS_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_ADM_RHS_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_ADM_RHS_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_ADM_RHS_calc_every != ML_ADM_RHS_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_ADM::ML_curvrhs.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_ADM::ML_lapserhs.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_ADM::ML_metricrhs.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_ADM::ML_shiftrhs.");
  return;
}

static void ML_ADM_RHS_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  const CCTK_REAL_VEC p1o12dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dx);
  const CCTK_REAL_VEC p1o12dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dy);
  const CCTK_REAL_VEC p1o12dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dz);
  const CCTK_REAL_VEC p1o144dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dx,dy));
  const CCTK_REAL_VEC p1o144dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dx,dz));
  const CCTK_REAL_VEC p1o144dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dy,dz));
  const CCTK_REAL_VEC p1o180dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dx,dx));
  const CCTK_REAL_VEC p1o180dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dy,dy));
  const CCTK_REAL_VEC p1o180dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dz,dz));
  const CCTK_REAL_VEC p1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dx);
  const CCTK_REAL_VEC p1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dy);
  const CCTK_REAL_VEC p1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dz);
  const CCTK_REAL_VEC p1o3600dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dx,dy));
  const CCTK_REAL_VEC p1o3600dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dx,dz));
  const CCTK_REAL_VEC p1o3600dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dy,dz));
  const CCTK_REAL_VEC p1o4dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dx,dy));
  const CCTK_REAL_VEC p1o4dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dx,dz));
  const CCTK_REAL_VEC p1o4dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dy,dz));
  const CCTK_REAL_VEC p1o5040dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dx,dx));
  const CCTK_REAL_VEC p1o5040dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dy,dy));
  const CCTK_REAL_VEC p1o5040dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dz,dz));
  const CCTK_REAL_VEC p1o60dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0166666666666666666666666666667),dx);
  const CCTK_REAL_VEC p1o60dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0166666666666666666666666666667),dy);
  const CCTK_REAL_VEC p1o60dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0166666666666666666666666666667),dz);
  const CCTK_REAL_VEC p1o705600dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dx,dy));
  const CCTK_REAL_VEC p1o705600dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dx,dz));
  const CCTK_REAL_VEC p1o705600dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dy,dz));
  const CCTK_REAL_VEC p1o840dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dx);
  const CCTK_REAL_VEC p1o840dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dy);
  const CCTK_REAL_VEC p1o840dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dz);
  const CCTK_REAL_VEC p1odx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dx,dx));
  const CCTK_REAL_VEC p1ody2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dy,dy));
  const CCTK_REAL_VEC p1odz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dz,dz));
  const CCTK_REAL_VEC pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dx,dx));
  const CCTK_REAL_VEC pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dy,dy));
  const CCTK_REAL_VEC pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dz,dz));
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
  CCTK_LOOP3STR(ML_ADM_RHS,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alphaL CCTK_ATTRIBUTE_UNUSED = vec_load(alpha[index]);
    CCTK_REAL_VEC beta1L CCTK_ATTRIBUTE_UNUSED = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L CCTK_ATTRIBUTE_UNUSED = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L CCTK_ATTRIBUTE_UNUSED = vec_load(beta3[index]);
    CCTK_REAL_VEC g11L CCTK_ATTRIBUTE_UNUSED = vec_load(g11[index]);
    CCTK_REAL_VEC g12L CCTK_ATTRIBUTE_UNUSED = vec_load(g12[index]);
    CCTK_REAL_VEC g13L CCTK_ATTRIBUTE_UNUSED = vec_load(g13[index]);
    CCTK_REAL_VEC g22L CCTK_ATTRIBUTE_UNUSED = vec_load(g22[index]);
    CCTK_REAL_VEC g23L CCTK_ATTRIBUTE_UNUSED = vec_load(g23[index]);
    CCTK_REAL_VEC g33L CCTK_ATTRIBUTE_UNUSED = vec_load(g33[index]);
    CCTK_REAL_VEC K11L CCTK_ATTRIBUTE_UNUSED = vec_load(K11[index]);
    CCTK_REAL_VEC K12L CCTK_ATTRIBUTE_UNUSED = vec_load(K12[index]);
    CCTK_REAL_VEC K13L CCTK_ATTRIBUTE_UNUSED = vec_load(K13[index]);
    CCTK_REAL_VEC K22L CCTK_ATTRIBUTE_UNUSED = vec_load(K22[index]);
    CCTK_REAL_VEC K23L CCTK_ATTRIBUTE_UNUSED = vec_load(K23[index]);
    CCTK_REAL_VEC K33L CCTK_ATTRIBUTE_UNUSED = vec_load(K33[index]);
    
    
    CCTK_REAL_VEC dJ111L, dJ112L, dJ113L, dJ122L, dJ123L, dJ133L, dJ211L, dJ212L, dJ213L, dJ222L, dJ223L, dJ233L, dJ311L, dJ312L, dJ313L, dJ322L, dJ323L, dJ333L, J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L CCTK_ATTRIBUTE_UNUSED ;
    
    if (usejacobian)
    {
      dJ111L = vec_load(dJ111[index]);
      dJ112L = vec_load(dJ112[index]);
      dJ113L = vec_load(dJ113[index]);
      dJ122L = vec_load(dJ122[index]);
      dJ123L = vec_load(dJ123[index]);
      dJ133L = vec_load(dJ133[index]);
      dJ211L = vec_load(dJ211[index]);
      dJ212L = vec_load(dJ212[index]);
      dJ213L = vec_load(dJ213[index]);
      dJ222L = vec_load(dJ222[index]);
      dJ223L = vec_load(dJ223[index]);
      dJ233L = vec_load(dJ233[index]);
      dJ311L = vec_load(dJ311[index]);
      dJ312L = vec_load(dJ312[index]);
      dJ313L = vec_load(dJ313[index]);
      dJ322L = vec_load(dJ322[index]);
      dJ323L = vec_load(dJ323[index]);
      dJ333L = vec_load(dJ333[index]);
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
    CCTK_REAL_VEC PDstandardNth1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1K11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2K11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3K11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1K12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2K12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3K12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1K13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2K13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3K13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1K22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2K22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3K22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1K23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2K23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3K23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1K33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2K33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3K33 CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDstandardNth1alpha = PDstandardNthfdOrder21(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder22(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder23(&alpha[index]);
        PDstandardNth11alpha = PDstandardNthfdOrder211(&alpha[index]);
        PDstandardNth22alpha = PDstandardNthfdOrder222(&alpha[index]);
        PDstandardNth33alpha = PDstandardNthfdOrder233(&alpha[index]);
        PDstandardNth12alpha = PDstandardNthfdOrder212(&alpha[index]);
        PDstandardNth13alpha = PDstandardNthfdOrder213(&alpha[index]);
        PDstandardNth23alpha = PDstandardNthfdOrder223(&alpha[index]);
        PDstandardNth1beta1 = PDstandardNthfdOrder21(&beta1[index]);
        PDstandardNth2beta1 = PDstandardNthfdOrder22(&beta1[index]);
        PDstandardNth3beta1 = PDstandardNthfdOrder23(&beta1[index]);
        PDstandardNth1beta2 = PDstandardNthfdOrder21(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder22(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder23(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder21(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder22(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder23(&beta3[index]);
        PDstandardNth1g11 = PDstandardNthfdOrder21(&g11[index]);
        PDstandardNth2g11 = PDstandardNthfdOrder22(&g11[index]);
        PDstandardNth3g11 = PDstandardNthfdOrder23(&g11[index]);
        PDstandardNth11g11 = PDstandardNthfdOrder211(&g11[index]);
        PDstandardNth22g11 = PDstandardNthfdOrder222(&g11[index]);
        PDstandardNth33g11 = PDstandardNthfdOrder233(&g11[index]);
        PDstandardNth12g11 = PDstandardNthfdOrder212(&g11[index]);
        PDstandardNth13g11 = PDstandardNthfdOrder213(&g11[index]);
        PDstandardNth23g11 = PDstandardNthfdOrder223(&g11[index]);
        PDstandardNth1g12 = PDstandardNthfdOrder21(&g12[index]);
        PDstandardNth2g12 = PDstandardNthfdOrder22(&g12[index]);
        PDstandardNth3g12 = PDstandardNthfdOrder23(&g12[index]);
        PDstandardNth11g12 = PDstandardNthfdOrder211(&g12[index]);
        PDstandardNth22g12 = PDstandardNthfdOrder222(&g12[index]);
        PDstandardNth33g12 = PDstandardNthfdOrder233(&g12[index]);
        PDstandardNth12g12 = PDstandardNthfdOrder212(&g12[index]);
        PDstandardNth13g12 = PDstandardNthfdOrder213(&g12[index]);
        PDstandardNth23g12 = PDstandardNthfdOrder223(&g12[index]);
        PDstandardNth1g13 = PDstandardNthfdOrder21(&g13[index]);
        PDstandardNth2g13 = PDstandardNthfdOrder22(&g13[index]);
        PDstandardNth3g13 = PDstandardNthfdOrder23(&g13[index]);
        PDstandardNth11g13 = PDstandardNthfdOrder211(&g13[index]);
        PDstandardNth22g13 = PDstandardNthfdOrder222(&g13[index]);
        PDstandardNth33g13 = PDstandardNthfdOrder233(&g13[index]);
        PDstandardNth12g13 = PDstandardNthfdOrder212(&g13[index]);
        PDstandardNth13g13 = PDstandardNthfdOrder213(&g13[index]);
        PDstandardNth23g13 = PDstandardNthfdOrder223(&g13[index]);
        PDstandardNth1g22 = PDstandardNthfdOrder21(&g22[index]);
        PDstandardNth2g22 = PDstandardNthfdOrder22(&g22[index]);
        PDstandardNth3g22 = PDstandardNthfdOrder23(&g22[index]);
        PDstandardNth11g22 = PDstandardNthfdOrder211(&g22[index]);
        PDstandardNth22g22 = PDstandardNthfdOrder222(&g22[index]);
        PDstandardNth33g22 = PDstandardNthfdOrder233(&g22[index]);
        PDstandardNth12g22 = PDstandardNthfdOrder212(&g22[index]);
        PDstandardNth13g22 = PDstandardNthfdOrder213(&g22[index]);
        PDstandardNth23g22 = PDstandardNthfdOrder223(&g22[index]);
        PDstandardNth1g23 = PDstandardNthfdOrder21(&g23[index]);
        PDstandardNth2g23 = PDstandardNthfdOrder22(&g23[index]);
        PDstandardNth3g23 = PDstandardNthfdOrder23(&g23[index]);
        PDstandardNth11g23 = PDstandardNthfdOrder211(&g23[index]);
        PDstandardNth22g23 = PDstandardNthfdOrder222(&g23[index]);
        PDstandardNth33g23 = PDstandardNthfdOrder233(&g23[index]);
        PDstandardNth12g23 = PDstandardNthfdOrder212(&g23[index]);
        PDstandardNth13g23 = PDstandardNthfdOrder213(&g23[index]);
        PDstandardNth23g23 = PDstandardNthfdOrder223(&g23[index]);
        PDstandardNth1g33 = PDstandardNthfdOrder21(&g33[index]);
        PDstandardNth2g33 = PDstandardNthfdOrder22(&g33[index]);
        PDstandardNth3g33 = PDstandardNthfdOrder23(&g33[index]);
        PDstandardNth11g33 = PDstandardNthfdOrder211(&g33[index]);
        PDstandardNth22g33 = PDstandardNthfdOrder222(&g33[index]);
        PDstandardNth33g33 = PDstandardNthfdOrder233(&g33[index]);
        PDstandardNth12g33 = PDstandardNthfdOrder212(&g33[index]);
        PDstandardNth13g33 = PDstandardNthfdOrder213(&g33[index]);
        PDstandardNth23g33 = PDstandardNthfdOrder223(&g33[index]);
        PDstandardNth1K11 = PDstandardNthfdOrder21(&K11[index]);
        PDstandardNth2K11 = PDstandardNthfdOrder22(&K11[index]);
        PDstandardNth3K11 = PDstandardNthfdOrder23(&K11[index]);
        PDstandardNth1K12 = PDstandardNthfdOrder21(&K12[index]);
        PDstandardNth2K12 = PDstandardNthfdOrder22(&K12[index]);
        PDstandardNth3K12 = PDstandardNthfdOrder23(&K12[index]);
        PDstandardNth1K13 = PDstandardNthfdOrder21(&K13[index]);
        PDstandardNth2K13 = PDstandardNthfdOrder22(&K13[index]);
        PDstandardNth3K13 = PDstandardNthfdOrder23(&K13[index]);
        PDstandardNth1K22 = PDstandardNthfdOrder21(&K22[index]);
        PDstandardNth2K22 = PDstandardNthfdOrder22(&K22[index]);
        PDstandardNth3K22 = PDstandardNthfdOrder23(&K22[index]);
        PDstandardNth1K23 = PDstandardNthfdOrder21(&K23[index]);
        PDstandardNth2K23 = PDstandardNthfdOrder22(&K23[index]);
        PDstandardNth3K23 = PDstandardNthfdOrder23(&K23[index]);
        PDstandardNth1K33 = PDstandardNthfdOrder21(&K33[index]);
        PDstandardNth2K33 = PDstandardNthfdOrder22(&K33[index]);
        PDstandardNth3K33 = PDstandardNthfdOrder23(&K33[index]);
        break;
      }
      
      case 4:
      {
        PDstandardNth1alpha = PDstandardNthfdOrder41(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder42(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder43(&alpha[index]);
        PDstandardNth11alpha = PDstandardNthfdOrder411(&alpha[index]);
        PDstandardNth22alpha = PDstandardNthfdOrder422(&alpha[index]);
        PDstandardNth33alpha = PDstandardNthfdOrder433(&alpha[index]);
        PDstandardNth12alpha = PDstandardNthfdOrder412(&alpha[index]);
        PDstandardNth13alpha = PDstandardNthfdOrder413(&alpha[index]);
        PDstandardNth23alpha = PDstandardNthfdOrder423(&alpha[index]);
        PDstandardNth1beta1 = PDstandardNthfdOrder41(&beta1[index]);
        PDstandardNth2beta1 = PDstandardNthfdOrder42(&beta1[index]);
        PDstandardNth3beta1 = PDstandardNthfdOrder43(&beta1[index]);
        PDstandardNth1beta2 = PDstandardNthfdOrder41(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder42(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder43(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder41(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder42(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder43(&beta3[index]);
        PDstandardNth1g11 = PDstandardNthfdOrder41(&g11[index]);
        PDstandardNth2g11 = PDstandardNthfdOrder42(&g11[index]);
        PDstandardNth3g11 = PDstandardNthfdOrder43(&g11[index]);
        PDstandardNth11g11 = PDstandardNthfdOrder411(&g11[index]);
        PDstandardNth22g11 = PDstandardNthfdOrder422(&g11[index]);
        PDstandardNth33g11 = PDstandardNthfdOrder433(&g11[index]);
        PDstandardNth12g11 = PDstandardNthfdOrder412(&g11[index]);
        PDstandardNth13g11 = PDstandardNthfdOrder413(&g11[index]);
        PDstandardNth23g11 = PDstandardNthfdOrder423(&g11[index]);
        PDstandardNth1g12 = PDstandardNthfdOrder41(&g12[index]);
        PDstandardNth2g12 = PDstandardNthfdOrder42(&g12[index]);
        PDstandardNth3g12 = PDstandardNthfdOrder43(&g12[index]);
        PDstandardNth11g12 = PDstandardNthfdOrder411(&g12[index]);
        PDstandardNth22g12 = PDstandardNthfdOrder422(&g12[index]);
        PDstandardNth33g12 = PDstandardNthfdOrder433(&g12[index]);
        PDstandardNth12g12 = PDstandardNthfdOrder412(&g12[index]);
        PDstandardNth13g12 = PDstandardNthfdOrder413(&g12[index]);
        PDstandardNth23g12 = PDstandardNthfdOrder423(&g12[index]);
        PDstandardNth1g13 = PDstandardNthfdOrder41(&g13[index]);
        PDstandardNth2g13 = PDstandardNthfdOrder42(&g13[index]);
        PDstandardNth3g13 = PDstandardNthfdOrder43(&g13[index]);
        PDstandardNth11g13 = PDstandardNthfdOrder411(&g13[index]);
        PDstandardNth22g13 = PDstandardNthfdOrder422(&g13[index]);
        PDstandardNth33g13 = PDstandardNthfdOrder433(&g13[index]);
        PDstandardNth12g13 = PDstandardNthfdOrder412(&g13[index]);
        PDstandardNth13g13 = PDstandardNthfdOrder413(&g13[index]);
        PDstandardNth23g13 = PDstandardNthfdOrder423(&g13[index]);
        PDstandardNth1g22 = PDstandardNthfdOrder41(&g22[index]);
        PDstandardNth2g22 = PDstandardNthfdOrder42(&g22[index]);
        PDstandardNth3g22 = PDstandardNthfdOrder43(&g22[index]);
        PDstandardNth11g22 = PDstandardNthfdOrder411(&g22[index]);
        PDstandardNth22g22 = PDstandardNthfdOrder422(&g22[index]);
        PDstandardNth33g22 = PDstandardNthfdOrder433(&g22[index]);
        PDstandardNth12g22 = PDstandardNthfdOrder412(&g22[index]);
        PDstandardNth13g22 = PDstandardNthfdOrder413(&g22[index]);
        PDstandardNth23g22 = PDstandardNthfdOrder423(&g22[index]);
        PDstandardNth1g23 = PDstandardNthfdOrder41(&g23[index]);
        PDstandardNth2g23 = PDstandardNthfdOrder42(&g23[index]);
        PDstandardNth3g23 = PDstandardNthfdOrder43(&g23[index]);
        PDstandardNth11g23 = PDstandardNthfdOrder411(&g23[index]);
        PDstandardNth22g23 = PDstandardNthfdOrder422(&g23[index]);
        PDstandardNth33g23 = PDstandardNthfdOrder433(&g23[index]);
        PDstandardNth12g23 = PDstandardNthfdOrder412(&g23[index]);
        PDstandardNth13g23 = PDstandardNthfdOrder413(&g23[index]);
        PDstandardNth23g23 = PDstandardNthfdOrder423(&g23[index]);
        PDstandardNth1g33 = PDstandardNthfdOrder41(&g33[index]);
        PDstandardNth2g33 = PDstandardNthfdOrder42(&g33[index]);
        PDstandardNth3g33 = PDstandardNthfdOrder43(&g33[index]);
        PDstandardNth11g33 = PDstandardNthfdOrder411(&g33[index]);
        PDstandardNth22g33 = PDstandardNthfdOrder422(&g33[index]);
        PDstandardNth33g33 = PDstandardNthfdOrder433(&g33[index]);
        PDstandardNth12g33 = PDstandardNthfdOrder412(&g33[index]);
        PDstandardNth13g33 = PDstandardNthfdOrder413(&g33[index]);
        PDstandardNth23g33 = PDstandardNthfdOrder423(&g33[index]);
        PDstandardNth1K11 = PDstandardNthfdOrder41(&K11[index]);
        PDstandardNth2K11 = PDstandardNthfdOrder42(&K11[index]);
        PDstandardNth3K11 = PDstandardNthfdOrder43(&K11[index]);
        PDstandardNth1K12 = PDstandardNthfdOrder41(&K12[index]);
        PDstandardNth2K12 = PDstandardNthfdOrder42(&K12[index]);
        PDstandardNth3K12 = PDstandardNthfdOrder43(&K12[index]);
        PDstandardNth1K13 = PDstandardNthfdOrder41(&K13[index]);
        PDstandardNth2K13 = PDstandardNthfdOrder42(&K13[index]);
        PDstandardNth3K13 = PDstandardNthfdOrder43(&K13[index]);
        PDstandardNth1K22 = PDstandardNthfdOrder41(&K22[index]);
        PDstandardNth2K22 = PDstandardNthfdOrder42(&K22[index]);
        PDstandardNth3K22 = PDstandardNthfdOrder43(&K22[index]);
        PDstandardNth1K23 = PDstandardNthfdOrder41(&K23[index]);
        PDstandardNth2K23 = PDstandardNthfdOrder42(&K23[index]);
        PDstandardNth3K23 = PDstandardNthfdOrder43(&K23[index]);
        PDstandardNth1K33 = PDstandardNthfdOrder41(&K33[index]);
        PDstandardNth2K33 = PDstandardNthfdOrder42(&K33[index]);
        PDstandardNth3K33 = PDstandardNthfdOrder43(&K33[index]);
        break;
      }
      
      case 6:
      {
        PDstandardNth1alpha = PDstandardNthfdOrder61(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder62(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder63(&alpha[index]);
        PDstandardNth11alpha = PDstandardNthfdOrder611(&alpha[index]);
        PDstandardNth22alpha = PDstandardNthfdOrder622(&alpha[index]);
        PDstandardNth33alpha = PDstandardNthfdOrder633(&alpha[index]);
        PDstandardNth12alpha = PDstandardNthfdOrder612(&alpha[index]);
        PDstandardNth13alpha = PDstandardNthfdOrder613(&alpha[index]);
        PDstandardNth23alpha = PDstandardNthfdOrder623(&alpha[index]);
        PDstandardNth1beta1 = PDstandardNthfdOrder61(&beta1[index]);
        PDstandardNth2beta1 = PDstandardNthfdOrder62(&beta1[index]);
        PDstandardNth3beta1 = PDstandardNthfdOrder63(&beta1[index]);
        PDstandardNth1beta2 = PDstandardNthfdOrder61(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder62(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder63(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder61(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder62(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder63(&beta3[index]);
        PDstandardNth1g11 = PDstandardNthfdOrder61(&g11[index]);
        PDstandardNth2g11 = PDstandardNthfdOrder62(&g11[index]);
        PDstandardNth3g11 = PDstandardNthfdOrder63(&g11[index]);
        PDstandardNth11g11 = PDstandardNthfdOrder611(&g11[index]);
        PDstandardNth22g11 = PDstandardNthfdOrder622(&g11[index]);
        PDstandardNth33g11 = PDstandardNthfdOrder633(&g11[index]);
        PDstandardNth12g11 = PDstandardNthfdOrder612(&g11[index]);
        PDstandardNth13g11 = PDstandardNthfdOrder613(&g11[index]);
        PDstandardNth23g11 = PDstandardNthfdOrder623(&g11[index]);
        PDstandardNth1g12 = PDstandardNthfdOrder61(&g12[index]);
        PDstandardNth2g12 = PDstandardNthfdOrder62(&g12[index]);
        PDstandardNth3g12 = PDstandardNthfdOrder63(&g12[index]);
        PDstandardNth11g12 = PDstandardNthfdOrder611(&g12[index]);
        PDstandardNth22g12 = PDstandardNthfdOrder622(&g12[index]);
        PDstandardNth33g12 = PDstandardNthfdOrder633(&g12[index]);
        PDstandardNth12g12 = PDstandardNthfdOrder612(&g12[index]);
        PDstandardNth13g12 = PDstandardNthfdOrder613(&g12[index]);
        PDstandardNth23g12 = PDstandardNthfdOrder623(&g12[index]);
        PDstandardNth1g13 = PDstandardNthfdOrder61(&g13[index]);
        PDstandardNth2g13 = PDstandardNthfdOrder62(&g13[index]);
        PDstandardNth3g13 = PDstandardNthfdOrder63(&g13[index]);
        PDstandardNth11g13 = PDstandardNthfdOrder611(&g13[index]);
        PDstandardNth22g13 = PDstandardNthfdOrder622(&g13[index]);
        PDstandardNth33g13 = PDstandardNthfdOrder633(&g13[index]);
        PDstandardNth12g13 = PDstandardNthfdOrder612(&g13[index]);
        PDstandardNth13g13 = PDstandardNthfdOrder613(&g13[index]);
        PDstandardNth23g13 = PDstandardNthfdOrder623(&g13[index]);
        PDstandardNth1g22 = PDstandardNthfdOrder61(&g22[index]);
        PDstandardNth2g22 = PDstandardNthfdOrder62(&g22[index]);
        PDstandardNth3g22 = PDstandardNthfdOrder63(&g22[index]);
        PDstandardNth11g22 = PDstandardNthfdOrder611(&g22[index]);
        PDstandardNth22g22 = PDstandardNthfdOrder622(&g22[index]);
        PDstandardNth33g22 = PDstandardNthfdOrder633(&g22[index]);
        PDstandardNth12g22 = PDstandardNthfdOrder612(&g22[index]);
        PDstandardNth13g22 = PDstandardNthfdOrder613(&g22[index]);
        PDstandardNth23g22 = PDstandardNthfdOrder623(&g22[index]);
        PDstandardNth1g23 = PDstandardNthfdOrder61(&g23[index]);
        PDstandardNth2g23 = PDstandardNthfdOrder62(&g23[index]);
        PDstandardNth3g23 = PDstandardNthfdOrder63(&g23[index]);
        PDstandardNth11g23 = PDstandardNthfdOrder611(&g23[index]);
        PDstandardNth22g23 = PDstandardNthfdOrder622(&g23[index]);
        PDstandardNth33g23 = PDstandardNthfdOrder633(&g23[index]);
        PDstandardNth12g23 = PDstandardNthfdOrder612(&g23[index]);
        PDstandardNth13g23 = PDstandardNthfdOrder613(&g23[index]);
        PDstandardNth23g23 = PDstandardNthfdOrder623(&g23[index]);
        PDstandardNth1g33 = PDstandardNthfdOrder61(&g33[index]);
        PDstandardNth2g33 = PDstandardNthfdOrder62(&g33[index]);
        PDstandardNth3g33 = PDstandardNthfdOrder63(&g33[index]);
        PDstandardNth11g33 = PDstandardNthfdOrder611(&g33[index]);
        PDstandardNth22g33 = PDstandardNthfdOrder622(&g33[index]);
        PDstandardNth33g33 = PDstandardNthfdOrder633(&g33[index]);
        PDstandardNth12g33 = PDstandardNthfdOrder612(&g33[index]);
        PDstandardNth13g33 = PDstandardNthfdOrder613(&g33[index]);
        PDstandardNth23g33 = PDstandardNthfdOrder623(&g33[index]);
        PDstandardNth1K11 = PDstandardNthfdOrder61(&K11[index]);
        PDstandardNth2K11 = PDstandardNthfdOrder62(&K11[index]);
        PDstandardNth3K11 = PDstandardNthfdOrder63(&K11[index]);
        PDstandardNth1K12 = PDstandardNthfdOrder61(&K12[index]);
        PDstandardNth2K12 = PDstandardNthfdOrder62(&K12[index]);
        PDstandardNth3K12 = PDstandardNthfdOrder63(&K12[index]);
        PDstandardNth1K13 = PDstandardNthfdOrder61(&K13[index]);
        PDstandardNth2K13 = PDstandardNthfdOrder62(&K13[index]);
        PDstandardNth3K13 = PDstandardNthfdOrder63(&K13[index]);
        PDstandardNth1K22 = PDstandardNthfdOrder61(&K22[index]);
        PDstandardNth2K22 = PDstandardNthfdOrder62(&K22[index]);
        PDstandardNth3K22 = PDstandardNthfdOrder63(&K22[index]);
        PDstandardNth1K23 = PDstandardNthfdOrder61(&K23[index]);
        PDstandardNth2K23 = PDstandardNthfdOrder62(&K23[index]);
        PDstandardNth3K23 = PDstandardNthfdOrder63(&K23[index]);
        PDstandardNth1K33 = PDstandardNthfdOrder61(&K33[index]);
        PDstandardNth2K33 = PDstandardNthfdOrder62(&K33[index]);
        PDstandardNth3K33 = PDstandardNthfdOrder63(&K33[index]);
        break;
      }
      
      case 8:
      {
        PDstandardNth1alpha = PDstandardNthfdOrder81(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder82(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder83(&alpha[index]);
        PDstandardNth11alpha = PDstandardNthfdOrder811(&alpha[index]);
        PDstandardNth22alpha = PDstandardNthfdOrder822(&alpha[index]);
        PDstandardNth33alpha = PDstandardNthfdOrder833(&alpha[index]);
        PDstandardNth12alpha = PDstandardNthfdOrder812(&alpha[index]);
        PDstandardNth13alpha = PDstandardNthfdOrder813(&alpha[index]);
        PDstandardNth23alpha = PDstandardNthfdOrder823(&alpha[index]);
        PDstandardNth1beta1 = PDstandardNthfdOrder81(&beta1[index]);
        PDstandardNth2beta1 = PDstandardNthfdOrder82(&beta1[index]);
        PDstandardNth3beta1 = PDstandardNthfdOrder83(&beta1[index]);
        PDstandardNth1beta2 = PDstandardNthfdOrder81(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder82(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder83(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder81(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder82(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder83(&beta3[index]);
        PDstandardNth1g11 = PDstandardNthfdOrder81(&g11[index]);
        PDstandardNth2g11 = PDstandardNthfdOrder82(&g11[index]);
        PDstandardNth3g11 = PDstandardNthfdOrder83(&g11[index]);
        PDstandardNth11g11 = PDstandardNthfdOrder811(&g11[index]);
        PDstandardNth22g11 = PDstandardNthfdOrder822(&g11[index]);
        PDstandardNth33g11 = PDstandardNthfdOrder833(&g11[index]);
        PDstandardNth12g11 = PDstandardNthfdOrder812(&g11[index]);
        PDstandardNth13g11 = PDstandardNthfdOrder813(&g11[index]);
        PDstandardNth23g11 = PDstandardNthfdOrder823(&g11[index]);
        PDstandardNth1g12 = PDstandardNthfdOrder81(&g12[index]);
        PDstandardNth2g12 = PDstandardNthfdOrder82(&g12[index]);
        PDstandardNth3g12 = PDstandardNthfdOrder83(&g12[index]);
        PDstandardNth11g12 = PDstandardNthfdOrder811(&g12[index]);
        PDstandardNth22g12 = PDstandardNthfdOrder822(&g12[index]);
        PDstandardNth33g12 = PDstandardNthfdOrder833(&g12[index]);
        PDstandardNth12g12 = PDstandardNthfdOrder812(&g12[index]);
        PDstandardNth13g12 = PDstandardNthfdOrder813(&g12[index]);
        PDstandardNth23g12 = PDstandardNthfdOrder823(&g12[index]);
        PDstandardNth1g13 = PDstandardNthfdOrder81(&g13[index]);
        PDstandardNth2g13 = PDstandardNthfdOrder82(&g13[index]);
        PDstandardNth3g13 = PDstandardNthfdOrder83(&g13[index]);
        PDstandardNth11g13 = PDstandardNthfdOrder811(&g13[index]);
        PDstandardNth22g13 = PDstandardNthfdOrder822(&g13[index]);
        PDstandardNth33g13 = PDstandardNthfdOrder833(&g13[index]);
        PDstandardNth12g13 = PDstandardNthfdOrder812(&g13[index]);
        PDstandardNth13g13 = PDstandardNthfdOrder813(&g13[index]);
        PDstandardNth23g13 = PDstandardNthfdOrder823(&g13[index]);
        PDstandardNth1g22 = PDstandardNthfdOrder81(&g22[index]);
        PDstandardNth2g22 = PDstandardNthfdOrder82(&g22[index]);
        PDstandardNth3g22 = PDstandardNthfdOrder83(&g22[index]);
        PDstandardNth11g22 = PDstandardNthfdOrder811(&g22[index]);
        PDstandardNth22g22 = PDstandardNthfdOrder822(&g22[index]);
        PDstandardNth33g22 = PDstandardNthfdOrder833(&g22[index]);
        PDstandardNth12g22 = PDstandardNthfdOrder812(&g22[index]);
        PDstandardNth13g22 = PDstandardNthfdOrder813(&g22[index]);
        PDstandardNth23g22 = PDstandardNthfdOrder823(&g22[index]);
        PDstandardNth1g23 = PDstandardNthfdOrder81(&g23[index]);
        PDstandardNth2g23 = PDstandardNthfdOrder82(&g23[index]);
        PDstandardNth3g23 = PDstandardNthfdOrder83(&g23[index]);
        PDstandardNth11g23 = PDstandardNthfdOrder811(&g23[index]);
        PDstandardNth22g23 = PDstandardNthfdOrder822(&g23[index]);
        PDstandardNth33g23 = PDstandardNthfdOrder833(&g23[index]);
        PDstandardNth12g23 = PDstandardNthfdOrder812(&g23[index]);
        PDstandardNth13g23 = PDstandardNthfdOrder813(&g23[index]);
        PDstandardNth23g23 = PDstandardNthfdOrder823(&g23[index]);
        PDstandardNth1g33 = PDstandardNthfdOrder81(&g33[index]);
        PDstandardNth2g33 = PDstandardNthfdOrder82(&g33[index]);
        PDstandardNth3g33 = PDstandardNthfdOrder83(&g33[index]);
        PDstandardNth11g33 = PDstandardNthfdOrder811(&g33[index]);
        PDstandardNth22g33 = PDstandardNthfdOrder822(&g33[index]);
        PDstandardNth33g33 = PDstandardNthfdOrder833(&g33[index]);
        PDstandardNth12g33 = PDstandardNthfdOrder812(&g33[index]);
        PDstandardNth13g33 = PDstandardNthfdOrder813(&g33[index]);
        PDstandardNth23g33 = PDstandardNthfdOrder823(&g33[index]);
        PDstandardNth1K11 = PDstandardNthfdOrder81(&K11[index]);
        PDstandardNth2K11 = PDstandardNthfdOrder82(&K11[index]);
        PDstandardNth3K11 = PDstandardNthfdOrder83(&K11[index]);
        PDstandardNth1K12 = PDstandardNthfdOrder81(&K12[index]);
        PDstandardNth2K12 = PDstandardNthfdOrder82(&K12[index]);
        PDstandardNth3K12 = PDstandardNthfdOrder83(&K12[index]);
        PDstandardNth1K13 = PDstandardNthfdOrder81(&K13[index]);
        PDstandardNth2K13 = PDstandardNthfdOrder82(&K13[index]);
        PDstandardNth3K13 = PDstandardNthfdOrder83(&K13[index]);
        PDstandardNth1K22 = PDstandardNthfdOrder81(&K22[index]);
        PDstandardNth2K22 = PDstandardNthfdOrder82(&K22[index]);
        PDstandardNth3K22 = PDstandardNthfdOrder83(&K22[index]);
        PDstandardNth1K23 = PDstandardNthfdOrder81(&K23[index]);
        PDstandardNth2K23 = PDstandardNthfdOrder82(&K23[index]);
        PDstandardNth3K23 = PDstandardNthfdOrder83(&K23[index]);
        PDstandardNth1K33 = PDstandardNthfdOrder81(&K33[index]);
        PDstandardNth2K33 = PDstandardNthfdOrder82(&K33[index]);
        PDstandardNth3K33 = PDstandardNthfdOrder83(&K33[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC JacPDstandardNth11alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth11g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth11g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth11g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1K11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1K12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1K13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1K22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1K23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1K33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2K11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2K12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2K13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2K22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2K23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2K33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3g11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3g12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3g13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3g22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3g23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3g33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3K11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3K12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3K13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3K22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3K23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3K33 CCTK_ATTRIBUTE_UNUSED;
    
    if (usejacobian)
    {
      JacPDstandardNth1alpha = 
        kmadd(J11L,PDstandardNth1alpha,kmadd(J21L,PDstandardNth2alpha,kmul(J31L,PDstandardNth3alpha)));
      
      JacPDstandardNth1beta1 = 
        kmadd(J11L,PDstandardNth1beta1,kmadd(J21L,PDstandardNth2beta1,kmul(J31L,PDstandardNth3beta1)));
      
      JacPDstandardNth1beta2 = 
        kmadd(J11L,PDstandardNth1beta2,kmadd(J21L,PDstandardNth2beta2,kmul(J31L,PDstandardNth3beta2)));
      
      JacPDstandardNth1beta3 = 
        kmadd(J11L,PDstandardNth1beta3,kmadd(J21L,PDstandardNth2beta3,kmul(J31L,PDstandardNth3beta3)));
      
      JacPDstandardNth1g11 = 
        kmadd(J11L,PDstandardNth1g11,kmadd(J21L,PDstandardNth2g11,kmul(J31L,PDstandardNth3g11)));
      
      JacPDstandardNth1g12 = 
        kmadd(J11L,PDstandardNth1g12,kmadd(J21L,PDstandardNth2g12,kmul(J31L,PDstandardNth3g12)));
      
      JacPDstandardNth1g13 = 
        kmadd(J11L,PDstandardNth1g13,kmadd(J21L,PDstandardNth2g13,kmul(J31L,PDstandardNth3g13)));
      
      JacPDstandardNth1g22 = 
        kmadd(J11L,PDstandardNth1g22,kmadd(J21L,PDstandardNth2g22,kmul(J31L,PDstandardNth3g22)));
      
      JacPDstandardNth1g23 = 
        kmadd(J11L,PDstandardNth1g23,kmadd(J21L,PDstandardNth2g23,kmul(J31L,PDstandardNth3g23)));
      
      JacPDstandardNth1g33 = 
        kmadd(J11L,PDstandardNth1g33,kmadd(J21L,PDstandardNth2g33,kmul(J31L,PDstandardNth3g33)));
      
      JacPDstandardNth1K11 = 
        kmadd(J11L,PDstandardNth1K11,kmadd(J21L,PDstandardNth2K11,kmul(J31L,PDstandardNth3K11)));
      
      JacPDstandardNth1K12 = 
        kmadd(J11L,PDstandardNth1K12,kmadd(J21L,PDstandardNth2K12,kmul(J31L,PDstandardNth3K12)));
      
      JacPDstandardNth1K13 = 
        kmadd(J11L,PDstandardNth1K13,kmadd(J21L,PDstandardNth2K13,kmul(J31L,PDstandardNth3K13)));
      
      JacPDstandardNth1K22 = 
        kmadd(J11L,PDstandardNth1K22,kmadd(J21L,PDstandardNth2K22,kmul(J31L,PDstandardNth3K22)));
      
      JacPDstandardNth1K23 = 
        kmadd(J11L,PDstandardNth1K23,kmadd(J21L,PDstandardNth2K23,kmul(J31L,PDstandardNth3K23)));
      
      JacPDstandardNth1K33 = 
        kmadd(J11L,PDstandardNth1K33,kmadd(J21L,PDstandardNth2K33,kmul(J31L,PDstandardNth3K33)));
      
      JacPDstandardNth2alpha = 
        kmadd(J12L,PDstandardNth1alpha,kmadd(J22L,PDstandardNth2alpha,kmul(J32L,PDstandardNth3alpha)));
      
      JacPDstandardNth2beta1 = 
        kmadd(J12L,PDstandardNth1beta1,kmadd(J22L,PDstandardNth2beta1,kmul(J32L,PDstandardNth3beta1)));
      
      JacPDstandardNth2beta2 = 
        kmadd(J12L,PDstandardNth1beta2,kmadd(J22L,PDstandardNth2beta2,kmul(J32L,PDstandardNth3beta2)));
      
      JacPDstandardNth2beta3 = 
        kmadd(J12L,PDstandardNth1beta3,kmadd(J22L,PDstandardNth2beta3,kmul(J32L,PDstandardNth3beta3)));
      
      JacPDstandardNth2g11 = 
        kmadd(J12L,PDstandardNth1g11,kmadd(J22L,PDstandardNth2g11,kmul(J32L,PDstandardNth3g11)));
      
      JacPDstandardNth2g12 = 
        kmadd(J12L,PDstandardNth1g12,kmadd(J22L,PDstandardNth2g12,kmul(J32L,PDstandardNth3g12)));
      
      JacPDstandardNth2g13 = 
        kmadd(J12L,PDstandardNth1g13,kmadd(J22L,PDstandardNth2g13,kmul(J32L,PDstandardNth3g13)));
      
      JacPDstandardNth2g22 = 
        kmadd(J12L,PDstandardNth1g22,kmadd(J22L,PDstandardNth2g22,kmul(J32L,PDstandardNth3g22)));
      
      JacPDstandardNth2g23 = 
        kmadd(J12L,PDstandardNth1g23,kmadd(J22L,PDstandardNth2g23,kmul(J32L,PDstandardNth3g23)));
      
      JacPDstandardNth2g33 = 
        kmadd(J12L,PDstandardNth1g33,kmadd(J22L,PDstandardNth2g33,kmul(J32L,PDstandardNth3g33)));
      
      JacPDstandardNth2K11 = 
        kmadd(J12L,PDstandardNth1K11,kmadd(J22L,PDstandardNth2K11,kmul(J32L,PDstandardNth3K11)));
      
      JacPDstandardNth2K12 = 
        kmadd(J12L,PDstandardNth1K12,kmadd(J22L,PDstandardNth2K12,kmul(J32L,PDstandardNth3K12)));
      
      JacPDstandardNth2K13 = 
        kmadd(J12L,PDstandardNth1K13,kmadd(J22L,PDstandardNth2K13,kmul(J32L,PDstandardNth3K13)));
      
      JacPDstandardNth2K22 = 
        kmadd(J12L,PDstandardNth1K22,kmadd(J22L,PDstandardNth2K22,kmul(J32L,PDstandardNth3K22)));
      
      JacPDstandardNth2K23 = 
        kmadd(J12L,PDstandardNth1K23,kmadd(J22L,PDstandardNth2K23,kmul(J32L,PDstandardNth3K23)));
      
      JacPDstandardNth2K33 = 
        kmadd(J12L,PDstandardNth1K33,kmadd(J22L,PDstandardNth2K33,kmul(J32L,PDstandardNth3K33)));
      
      JacPDstandardNth3alpha = 
        kmadd(J13L,PDstandardNth1alpha,kmadd(J23L,PDstandardNth2alpha,kmul(J33L,PDstandardNth3alpha)));
      
      JacPDstandardNth3beta1 = 
        kmadd(J13L,PDstandardNth1beta1,kmadd(J23L,PDstandardNth2beta1,kmul(J33L,PDstandardNth3beta1)));
      
      JacPDstandardNth3beta2 = 
        kmadd(J13L,PDstandardNth1beta2,kmadd(J23L,PDstandardNth2beta2,kmul(J33L,PDstandardNth3beta2)));
      
      JacPDstandardNth3beta3 = 
        kmadd(J13L,PDstandardNth1beta3,kmadd(J23L,PDstandardNth2beta3,kmul(J33L,PDstandardNth3beta3)));
      
      JacPDstandardNth3g11 = 
        kmadd(J13L,PDstandardNth1g11,kmadd(J23L,PDstandardNth2g11,kmul(J33L,PDstandardNth3g11)));
      
      JacPDstandardNth3g12 = 
        kmadd(J13L,PDstandardNth1g12,kmadd(J23L,PDstandardNth2g12,kmul(J33L,PDstandardNth3g12)));
      
      JacPDstandardNth3g13 = 
        kmadd(J13L,PDstandardNth1g13,kmadd(J23L,PDstandardNth2g13,kmul(J33L,PDstandardNth3g13)));
      
      JacPDstandardNth3g22 = 
        kmadd(J13L,PDstandardNth1g22,kmadd(J23L,PDstandardNth2g22,kmul(J33L,PDstandardNth3g22)));
      
      JacPDstandardNth3g23 = 
        kmadd(J13L,PDstandardNth1g23,kmadd(J23L,PDstandardNth2g23,kmul(J33L,PDstandardNth3g23)));
      
      JacPDstandardNth3g33 = 
        kmadd(J13L,PDstandardNth1g33,kmadd(J23L,PDstandardNth2g33,kmul(J33L,PDstandardNth3g33)));
      
      JacPDstandardNth3K11 = 
        kmadd(J13L,PDstandardNth1K11,kmadd(J23L,PDstandardNth2K11,kmul(J33L,PDstandardNth3K11)));
      
      JacPDstandardNth3K12 = 
        kmadd(J13L,PDstandardNth1K12,kmadd(J23L,PDstandardNth2K12,kmul(J33L,PDstandardNth3K12)));
      
      JacPDstandardNth3K13 = 
        kmadd(J13L,PDstandardNth1K13,kmadd(J23L,PDstandardNth2K13,kmul(J33L,PDstandardNth3K13)));
      
      JacPDstandardNth3K22 = 
        kmadd(J13L,PDstandardNth1K22,kmadd(J23L,PDstandardNth2K22,kmul(J33L,PDstandardNth3K22)));
      
      JacPDstandardNth3K23 = 
        kmadd(J13L,PDstandardNth1K23,kmadd(J23L,PDstandardNth2K23,kmul(J33L,PDstandardNth3K23)));
      
      JacPDstandardNth3K33 = 
        kmadd(J13L,PDstandardNth1K33,kmadd(J23L,PDstandardNth2K33,kmul(J33L,PDstandardNth3K33)));
      
      JacPDstandardNth11alpha = 
        kmadd(dJ111L,PDstandardNth1alpha,kmadd(ToReal(2),kmadd(J11L,kmadd(J21L,PDstandardNth12alpha,kmul(J31L,PDstandardNth13alpha)),kmul(J21L,kmul(J31L,PDstandardNth23alpha))),kmadd(dJ211L,PDstandardNth2alpha,kmadd(dJ311L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,kmul(J11L,J11L),kmadd(PDstandardNth22alpha,kmul(J21L,J21L),kmul(PDstandardNth33alpha,kmul(J31L,J31L))))))));
      
      JacPDstandardNth11g22 = 
        kmadd(dJ111L,PDstandardNth1g22,kmadd(ToReal(2),kmadd(J11L,kmadd(J21L,PDstandardNth12g22,kmul(J31L,PDstandardNth13g22)),kmul(J21L,kmul(J31L,PDstandardNth23g22))),kmadd(dJ211L,PDstandardNth2g22,kmadd(dJ311L,PDstandardNth3g22,kmadd(PDstandardNth11g22,kmul(J11L,J11L),kmadd(PDstandardNth22g22,kmul(J21L,J21L),kmul(PDstandardNth33g22,kmul(J31L,J31L))))))));
      
      JacPDstandardNth11g23 = 
        kmadd(dJ111L,PDstandardNth1g23,kmadd(ToReal(2),kmadd(J11L,kmadd(J21L,PDstandardNth12g23,kmul(J31L,PDstandardNth13g23)),kmul(J21L,kmul(J31L,PDstandardNth23g23))),kmadd(dJ211L,PDstandardNth2g23,kmadd(dJ311L,PDstandardNth3g23,kmadd(PDstandardNth11g23,kmul(J11L,J11L),kmadd(PDstandardNth22g23,kmul(J21L,J21L),kmul(PDstandardNth33g23,kmul(J31L,J31L))))))));
      
      JacPDstandardNth11g33 = 
        kmadd(dJ111L,PDstandardNth1g33,kmadd(ToReal(2),kmadd(J11L,kmadd(J21L,PDstandardNth12g33,kmul(J31L,PDstandardNth13g33)),kmul(J21L,kmul(J31L,PDstandardNth23g33))),kmadd(dJ211L,PDstandardNth2g33,kmadd(dJ311L,PDstandardNth3g33,kmadd(PDstandardNth11g33,kmul(J11L,J11L),kmadd(PDstandardNth22g33,kmul(J21L,J21L),kmul(PDstandardNth33g33,kmul(J31L,J31L))))))));
      
      JacPDstandardNth22alpha = 
        kmadd(dJ122L,PDstandardNth1alpha,kmadd(ToReal(2),kmadd(J12L,kmadd(J22L,PDstandardNth12alpha,kmul(J32L,PDstandardNth13alpha)),kmul(J22L,kmul(J32L,PDstandardNth23alpha))),kmadd(dJ222L,PDstandardNth2alpha,kmadd(dJ322L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,kmul(J12L,J12L),kmadd(PDstandardNth22alpha,kmul(J22L,J22L),kmul(PDstandardNth33alpha,kmul(J32L,J32L))))))));
      
      JacPDstandardNth22g11 = 
        kmadd(dJ122L,PDstandardNth1g11,kmadd(ToReal(2),kmadd(J12L,kmadd(J22L,PDstandardNth12g11,kmul(J32L,PDstandardNth13g11)),kmul(J22L,kmul(J32L,PDstandardNth23g11))),kmadd(dJ222L,PDstandardNth2g11,kmadd(dJ322L,PDstandardNth3g11,kmadd(PDstandardNth11g11,kmul(J12L,J12L),kmadd(PDstandardNth22g11,kmul(J22L,J22L),kmul(PDstandardNth33g11,kmul(J32L,J32L))))))));
      
      JacPDstandardNth22g13 = 
        kmadd(dJ122L,PDstandardNth1g13,kmadd(ToReal(2),kmadd(J12L,kmadd(J22L,PDstandardNth12g13,kmul(J32L,PDstandardNth13g13)),kmul(J22L,kmul(J32L,PDstandardNth23g13))),kmadd(dJ222L,PDstandardNth2g13,kmadd(dJ322L,PDstandardNth3g13,kmadd(PDstandardNth11g13,kmul(J12L,J12L),kmadd(PDstandardNth22g13,kmul(J22L,J22L),kmul(PDstandardNth33g13,kmul(J32L,J32L))))))));
      
      JacPDstandardNth22g33 = 
        kmadd(dJ122L,PDstandardNth1g33,kmadd(ToReal(2),kmadd(J12L,kmadd(J22L,PDstandardNth12g33,kmul(J32L,PDstandardNth13g33)),kmul(J22L,kmul(J32L,PDstandardNth23g33))),kmadd(dJ222L,PDstandardNth2g33,kmadd(dJ322L,PDstandardNth3g33,kmadd(PDstandardNth11g33,kmul(J12L,J12L),kmadd(PDstandardNth22g33,kmul(J22L,J22L),kmul(PDstandardNth33g33,kmul(J32L,J32L))))))));
      
      JacPDstandardNth33alpha = 
        kmadd(dJ133L,PDstandardNth1alpha,kmadd(ToReal(2),kmadd(J13L,kmadd(J23L,PDstandardNth12alpha,kmul(J33L,PDstandardNth13alpha)),kmul(J23L,kmul(J33L,PDstandardNth23alpha))),kmadd(dJ233L,PDstandardNth2alpha,kmadd(dJ333L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,kmul(J13L,J13L),kmadd(PDstandardNth22alpha,kmul(J23L,J23L),kmul(PDstandardNth33alpha,kmul(J33L,J33L))))))));
      
      JacPDstandardNth33g11 = 
        kmadd(dJ133L,PDstandardNth1g11,kmadd(ToReal(2),kmadd(J13L,kmadd(J23L,PDstandardNth12g11,kmul(J33L,PDstandardNth13g11)),kmul(J23L,kmul(J33L,PDstandardNth23g11))),kmadd(dJ233L,PDstandardNth2g11,kmadd(dJ333L,PDstandardNth3g11,kmadd(PDstandardNth11g11,kmul(J13L,J13L),kmadd(PDstandardNth22g11,kmul(J23L,J23L),kmul(PDstandardNth33g11,kmul(J33L,J33L))))))));
      
      JacPDstandardNth33g12 = 
        kmadd(dJ133L,PDstandardNth1g12,kmadd(ToReal(2),kmadd(J13L,kmadd(J23L,PDstandardNth12g12,kmul(J33L,PDstandardNth13g12)),kmul(J23L,kmul(J33L,PDstandardNth23g12))),kmadd(dJ233L,PDstandardNth2g12,kmadd(dJ333L,PDstandardNth3g12,kmadd(PDstandardNth11g12,kmul(J13L,J13L),kmadd(PDstandardNth22g12,kmul(J23L,J23L),kmul(PDstandardNth33g12,kmul(J33L,J33L))))))));
      
      JacPDstandardNth33g22 = 
        kmadd(dJ133L,PDstandardNth1g22,kmadd(ToReal(2),kmadd(J13L,kmadd(J23L,PDstandardNth12g22,kmul(J33L,PDstandardNth13g22)),kmul(J23L,kmul(J33L,PDstandardNth23g22))),kmadd(dJ233L,PDstandardNth2g22,kmadd(dJ333L,PDstandardNth3g22,kmadd(PDstandardNth11g22,kmul(J13L,J13L),kmadd(PDstandardNth22g22,kmul(J23L,J23L),kmul(PDstandardNth33g22,kmul(J33L,J33L))))))));
      
      JacPDstandardNth12alpha = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11alpha,kmadd(J21L,PDstandardNth12alpha,kmul(J31L,PDstandardNth13alpha))),kmadd(J11L,kmadd(J22L,PDstandardNth12alpha,kmul(J32L,PDstandardNth13alpha)),kmadd(dJ112L,PDstandardNth1alpha,kmadd(J22L,kmadd(J21L,PDstandardNth22alpha,kmul(J31L,PDstandardNth23alpha)),kmadd(dJ212L,PDstandardNth2alpha,kmadd(J32L,kmadd(J21L,PDstandardNth23alpha,kmul(J31L,PDstandardNth33alpha)),kmul(dJ312L,PDstandardNth3alpha)))))));
      
      JacPDstandardNth12g11 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g11,kmadd(J21L,PDstandardNth12g11,kmul(J31L,PDstandardNth13g11))),kmadd(J11L,kmadd(J22L,PDstandardNth12g11,kmul(J32L,PDstandardNth13g11)),kmadd(dJ112L,PDstandardNth1g11,kmadd(J22L,kmadd(J21L,PDstandardNth22g11,kmul(J31L,PDstandardNth23g11)),kmadd(dJ212L,PDstandardNth2g11,kmadd(J32L,kmadd(J21L,PDstandardNth23g11,kmul(J31L,PDstandardNth33g11)),kmul(dJ312L,PDstandardNth3g11)))))));
      
      JacPDstandardNth12g12 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g12,kmadd(J21L,PDstandardNth12g12,kmul(J31L,PDstandardNth13g12))),kmadd(J11L,kmadd(J22L,PDstandardNth12g12,kmul(J32L,PDstandardNth13g12)),kmadd(dJ112L,PDstandardNth1g12,kmadd(J22L,kmadd(J21L,PDstandardNth22g12,kmul(J31L,PDstandardNth23g12)),kmadd(dJ212L,PDstandardNth2g12,kmadd(J32L,kmadd(J21L,PDstandardNth23g12,kmul(J31L,PDstandardNth33g12)),kmul(dJ312L,PDstandardNth3g12)))))));
      
      JacPDstandardNth12g13 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g13,kmadd(J21L,PDstandardNth12g13,kmul(J31L,PDstandardNth13g13))),kmadd(J11L,kmadd(J22L,PDstandardNth12g13,kmul(J32L,PDstandardNth13g13)),kmadd(dJ112L,PDstandardNth1g13,kmadd(J22L,kmadd(J21L,PDstandardNth22g13,kmul(J31L,PDstandardNth23g13)),kmadd(dJ212L,PDstandardNth2g13,kmadd(J32L,kmadd(J21L,PDstandardNth23g13,kmul(J31L,PDstandardNth33g13)),kmul(dJ312L,PDstandardNth3g13)))))));
      
      JacPDstandardNth12g22 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g22,kmadd(J21L,PDstandardNth12g22,kmul(J31L,PDstandardNth13g22))),kmadd(J11L,kmadd(J22L,PDstandardNth12g22,kmul(J32L,PDstandardNth13g22)),kmadd(dJ112L,PDstandardNth1g22,kmadd(J22L,kmadd(J21L,PDstandardNth22g22,kmul(J31L,PDstandardNth23g22)),kmadd(dJ212L,PDstandardNth2g22,kmadd(J32L,kmadd(J21L,PDstandardNth23g22,kmul(J31L,PDstandardNth33g22)),kmul(dJ312L,PDstandardNth3g22)))))));
      
      JacPDstandardNth12g23 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g23,kmadd(J21L,PDstandardNth12g23,kmul(J31L,PDstandardNth13g23))),kmadd(J11L,kmadd(J22L,PDstandardNth12g23,kmul(J32L,PDstandardNth13g23)),kmadd(dJ112L,PDstandardNth1g23,kmadd(J22L,kmadd(J21L,PDstandardNth22g23,kmul(J31L,PDstandardNth23g23)),kmadd(dJ212L,PDstandardNth2g23,kmadd(J32L,kmadd(J21L,PDstandardNth23g23,kmul(J31L,PDstandardNth33g23)),kmul(dJ312L,PDstandardNth3g23)))))));
      
      JacPDstandardNth12g33 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g33,kmadd(J21L,PDstandardNth12g33,kmul(J31L,PDstandardNth13g33))),kmadd(J11L,kmadd(J22L,PDstandardNth12g33,kmul(J32L,PDstandardNth13g33)),kmadd(dJ112L,PDstandardNth1g33,kmadd(J22L,kmadd(J21L,PDstandardNth22g33,kmul(J31L,PDstandardNth23g33)),kmadd(dJ212L,PDstandardNth2g33,kmadd(J32L,kmadd(J21L,PDstandardNth23g33,kmul(J31L,PDstandardNth33g33)),kmul(dJ312L,PDstandardNth3g33)))))));
      
      JacPDstandardNth13alpha = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11alpha,kmadd(J21L,PDstandardNth12alpha,kmul(J31L,PDstandardNth13alpha))),kmadd(J11L,kmadd(J23L,PDstandardNth12alpha,kmul(J33L,PDstandardNth13alpha)),kmadd(dJ113L,PDstandardNth1alpha,kmadd(J23L,kmadd(J21L,PDstandardNth22alpha,kmul(J31L,PDstandardNth23alpha)),kmadd(dJ213L,PDstandardNth2alpha,kmadd(J33L,kmadd(J21L,PDstandardNth23alpha,kmul(J31L,PDstandardNth33alpha)),kmul(dJ313L,PDstandardNth3alpha)))))));
      
      JacPDstandardNth13g11 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g11,kmadd(J21L,PDstandardNth12g11,kmul(J31L,PDstandardNth13g11))),kmadd(J11L,kmadd(J23L,PDstandardNth12g11,kmul(J33L,PDstandardNth13g11)),kmadd(dJ113L,PDstandardNth1g11,kmadd(J23L,kmadd(J21L,PDstandardNth22g11,kmul(J31L,PDstandardNth23g11)),kmadd(dJ213L,PDstandardNth2g11,kmadd(J33L,kmadd(J21L,PDstandardNth23g11,kmul(J31L,PDstandardNth33g11)),kmul(dJ313L,PDstandardNth3g11)))))));
      
      JacPDstandardNth13g12 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g12,kmadd(J21L,PDstandardNth12g12,kmul(J31L,PDstandardNth13g12))),kmadd(J11L,kmadd(J23L,PDstandardNth12g12,kmul(J33L,PDstandardNth13g12)),kmadd(dJ113L,PDstandardNth1g12,kmadd(J23L,kmadd(J21L,PDstandardNth22g12,kmul(J31L,PDstandardNth23g12)),kmadd(dJ213L,PDstandardNth2g12,kmadd(J33L,kmadd(J21L,PDstandardNth23g12,kmul(J31L,PDstandardNth33g12)),kmul(dJ313L,PDstandardNth3g12)))))));
      
      JacPDstandardNth13g13 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g13,kmadd(J21L,PDstandardNth12g13,kmul(J31L,PDstandardNth13g13))),kmadd(J11L,kmadd(J23L,PDstandardNth12g13,kmul(J33L,PDstandardNth13g13)),kmadd(dJ113L,PDstandardNth1g13,kmadd(J23L,kmadd(J21L,PDstandardNth22g13,kmul(J31L,PDstandardNth23g13)),kmadd(dJ213L,PDstandardNth2g13,kmadd(J33L,kmadd(J21L,PDstandardNth23g13,kmul(J31L,PDstandardNth33g13)),kmul(dJ313L,PDstandardNth3g13)))))));
      
      JacPDstandardNth13g22 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g22,kmadd(J21L,PDstandardNth12g22,kmul(J31L,PDstandardNth13g22))),kmadd(J11L,kmadd(J23L,PDstandardNth12g22,kmul(J33L,PDstandardNth13g22)),kmadd(dJ113L,PDstandardNth1g22,kmadd(J23L,kmadd(J21L,PDstandardNth22g22,kmul(J31L,PDstandardNth23g22)),kmadd(dJ213L,PDstandardNth2g22,kmadd(J33L,kmadd(J21L,PDstandardNth23g22,kmul(J31L,PDstandardNth33g22)),kmul(dJ313L,PDstandardNth3g22)))))));
      
      JacPDstandardNth13g23 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g23,kmadd(J21L,PDstandardNth12g23,kmul(J31L,PDstandardNth13g23))),kmadd(J11L,kmadd(J23L,PDstandardNth12g23,kmul(J33L,PDstandardNth13g23)),kmadd(dJ113L,PDstandardNth1g23,kmadd(J23L,kmadd(J21L,PDstandardNth22g23,kmul(J31L,PDstandardNth23g23)),kmadd(dJ213L,PDstandardNth2g23,kmadd(J33L,kmadd(J21L,PDstandardNth23g23,kmul(J31L,PDstandardNth33g23)),kmul(dJ313L,PDstandardNth3g23)))))));
      
      JacPDstandardNth13g33 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g33,kmadd(J21L,PDstandardNth12g33,kmul(J31L,PDstandardNth13g33))),kmadd(J11L,kmadd(J23L,PDstandardNth12g33,kmul(J33L,PDstandardNth13g33)),kmadd(dJ113L,PDstandardNth1g33,kmadd(J23L,kmadd(J21L,PDstandardNth22g33,kmul(J31L,PDstandardNth23g33)),kmadd(dJ213L,PDstandardNth2g33,kmadd(J33L,kmadd(J21L,PDstandardNth23g33,kmul(J31L,PDstandardNth33g33)),kmul(dJ313L,PDstandardNth3g33)))))));
      
      JacPDstandardNth21g11 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g11,kmadd(J21L,PDstandardNth12g11,kmul(J31L,PDstandardNth13g11))),kmadd(J11L,kmadd(J22L,PDstandardNth12g11,kmul(J32L,PDstandardNth13g11)),kmadd(dJ112L,PDstandardNth1g11,kmadd(J22L,kmadd(J21L,PDstandardNth22g11,kmul(J31L,PDstandardNth23g11)),kmadd(dJ212L,PDstandardNth2g11,kmadd(J32L,kmadd(J21L,PDstandardNth23g11,kmul(J31L,PDstandardNth33g11)),kmul(dJ312L,PDstandardNth3g11)))))));
      
      JacPDstandardNth21g12 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g12,kmadd(J21L,PDstandardNth12g12,kmul(J31L,PDstandardNth13g12))),kmadd(J11L,kmadd(J22L,PDstandardNth12g12,kmul(J32L,PDstandardNth13g12)),kmadd(dJ112L,PDstandardNth1g12,kmadd(J22L,kmadd(J21L,PDstandardNth22g12,kmul(J31L,PDstandardNth23g12)),kmadd(dJ212L,PDstandardNth2g12,kmadd(J32L,kmadd(J21L,PDstandardNth23g12,kmul(J31L,PDstandardNth33g12)),kmul(dJ312L,PDstandardNth3g12)))))));
      
      JacPDstandardNth21g13 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g13,kmadd(J21L,PDstandardNth12g13,kmul(J31L,PDstandardNth13g13))),kmadd(J11L,kmadd(J22L,PDstandardNth12g13,kmul(J32L,PDstandardNth13g13)),kmadd(dJ112L,PDstandardNth1g13,kmadd(J22L,kmadd(J21L,PDstandardNth22g13,kmul(J31L,PDstandardNth23g13)),kmadd(dJ212L,PDstandardNth2g13,kmadd(J32L,kmadd(J21L,PDstandardNth23g13,kmul(J31L,PDstandardNth33g13)),kmul(dJ312L,PDstandardNth3g13)))))));
      
      JacPDstandardNth21g22 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g22,kmadd(J21L,PDstandardNth12g22,kmul(J31L,PDstandardNth13g22))),kmadd(J11L,kmadd(J22L,PDstandardNth12g22,kmul(J32L,PDstandardNth13g22)),kmadd(dJ112L,PDstandardNth1g22,kmadd(J22L,kmadd(J21L,PDstandardNth22g22,kmul(J31L,PDstandardNth23g22)),kmadd(dJ212L,PDstandardNth2g22,kmadd(J32L,kmadd(J21L,PDstandardNth23g22,kmul(J31L,PDstandardNth33g22)),kmul(dJ312L,PDstandardNth3g22)))))));
      
      JacPDstandardNth21g23 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g23,kmadd(J21L,PDstandardNth12g23,kmul(J31L,PDstandardNth13g23))),kmadd(J11L,kmadd(J22L,PDstandardNth12g23,kmul(J32L,PDstandardNth13g23)),kmadd(dJ112L,PDstandardNth1g23,kmadd(J22L,kmadd(J21L,PDstandardNth22g23,kmul(J31L,PDstandardNth23g23)),kmadd(dJ212L,PDstandardNth2g23,kmadd(J32L,kmadd(J21L,PDstandardNth23g23,kmul(J31L,PDstandardNth33g23)),kmul(dJ312L,PDstandardNth3g23)))))));
      
      JacPDstandardNth21g33 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g33,kmadd(J21L,PDstandardNth12g33,kmul(J31L,PDstandardNth13g33))),kmadd(J11L,kmadd(J22L,PDstandardNth12g33,kmul(J32L,PDstandardNth13g33)),kmadd(dJ112L,PDstandardNth1g33,kmadd(J22L,kmadd(J21L,PDstandardNth22g33,kmul(J31L,PDstandardNth23g33)),kmadd(dJ212L,PDstandardNth2g33,kmadd(J32L,kmadd(J21L,PDstandardNth23g33,kmul(J31L,PDstandardNth33g33)),kmul(dJ312L,PDstandardNth3g33)))))));
      
      JacPDstandardNth23alpha = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11alpha,kmadd(J22L,PDstandardNth12alpha,kmul(J32L,PDstandardNth13alpha))),kmadd(J12L,kmadd(J23L,PDstandardNth12alpha,kmul(J33L,PDstandardNth13alpha)),kmadd(dJ123L,PDstandardNth1alpha,kmadd(J23L,kmadd(J22L,PDstandardNth22alpha,kmul(J32L,PDstandardNth23alpha)),kmadd(dJ223L,PDstandardNth2alpha,kmadd(J33L,kmadd(J22L,PDstandardNth23alpha,kmul(J32L,PDstandardNth33alpha)),kmul(dJ323L,PDstandardNth3alpha)))))));
      
      JacPDstandardNth23g11 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g11,kmadd(J22L,PDstandardNth12g11,kmul(J32L,PDstandardNth13g11))),kmadd(J12L,kmadd(J23L,PDstandardNth12g11,kmul(J33L,PDstandardNth13g11)),kmadd(dJ123L,PDstandardNth1g11,kmadd(J23L,kmadd(J22L,PDstandardNth22g11,kmul(J32L,PDstandardNth23g11)),kmadd(dJ223L,PDstandardNth2g11,kmadd(J33L,kmadd(J22L,PDstandardNth23g11,kmul(J32L,PDstandardNth33g11)),kmul(dJ323L,PDstandardNth3g11)))))));
      
      JacPDstandardNth23g12 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g12,kmadd(J22L,PDstandardNth12g12,kmul(J32L,PDstandardNth13g12))),kmadd(J12L,kmadd(J23L,PDstandardNth12g12,kmul(J33L,PDstandardNth13g12)),kmadd(dJ123L,PDstandardNth1g12,kmadd(J23L,kmadd(J22L,PDstandardNth22g12,kmul(J32L,PDstandardNth23g12)),kmadd(dJ223L,PDstandardNth2g12,kmadd(J33L,kmadd(J22L,PDstandardNth23g12,kmul(J32L,PDstandardNth33g12)),kmul(dJ323L,PDstandardNth3g12)))))));
      
      JacPDstandardNth23g13 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g13,kmadd(J22L,PDstandardNth12g13,kmul(J32L,PDstandardNth13g13))),kmadd(J12L,kmadd(J23L,PDstandardNth12g13,kmul(J33L,PDstandardNth13g13)),kmadd(dJ123L,PDstandardNth1g13,kmadd(J23L,kmadd(J22L,PDstandardNth22g13,kmul(J32L,PDstandardNth23g13)),kmadd(dJ223L,PDstandardNth2g13,kmadd(J33L,kmadd(J22L,PDstandardNth23g13,kmul(J32L,PDstandardNth33g13)),kmul(dJ323L,PDstandardNth3g13)))))));
      
      JacPDstandardNth23g22 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g22,kmadd(J22L,PDstandardNth12g22,kmul(J32L,PDstandardNth13g22))),kmadd(J12L,kmadd(J23L,PDstandardNth12g22,kmul(J33L,PDstandardNth13g22)),kmadd(dJ123L,PDstandardNth1g22,kmadd(J23L,kmadd(J22L,PDstandardNth22g22,kmul(J32L,PDstandardNth23g22)),kmadd(dJ223L,PDstandardNth2g22,kmadd(J33L,kmadd(J22L,PDstandardNth23g22,kmul(J32L,PDstandardNth33g22)),kmul(dJ323L,PDstandardNth3g22)))))));
      
      JacPDstandardNth23g23 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g23,kmadd(J22L,PDstandardNth12g23,kmul(J32L,PDstandardNth13g23))),kmadd(J12L,kmadd(J23L,PDstandardNth12g23,kmul(J33L,PDstandardNth13g23)),kmadd(dJ123L,PDstandardNth1g23,kmadd(J23L,kmadd(J22L,PDstandardNth22g23,kmul(J32L,PDstandardNth23g23)),kmadd(dJ223L,PDstandardNth2g23,kmadd(J33L,kmadd(J22L,PDstandardNth23g23,kmul(J32L,PDstandardNth33g23)),kmul(dJ323L,PDstandardNth3g23)))))));
      
      JacPDstandardNth23g33 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g33,kmadd(J22L,PDstandardNth12g33,kmul(J32L,PDstandardNth13g33))),kmadd(J12L,kmadd(J23L,PDstandardNth12g33,kmul(J33L,PDstandardNth13g33)),kmadd(dJ123L,PDstandardNth1g33,kmadd(J23L,kmadd(J22L,PDstandardNth22g33,kmul(J32L,PDstandardNth23g33)),kmadd(dJ223L,PDstandardNth2g33,kmadd(J33L,kmadd(J22L,PDstandardNth23g33,kmul(J32L,PDstandardNth33g33)),kmul(dJ323L,PDstandardNth3g33)))))));
      
      JacPDstandardNth31g11 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g11,kmadd(J21L,PDstandardNth12g11,kmul(J31L,PDstandardNth13g11))),kmadd(J11L,kmadd(J23L,PDstandardNth12g11,kmul(J33L,PDstandardNth13g11)),kmadd(dJ113L,PDstandardNth1g11,kmadd(J23L,kmadd(J21L,PDstandardNth22g11,kmul(J31L,PDstandardNth23g11)),kmadd(dJ213L,PDstandardNth2g11,kmadd(J33L,kmadd(J21L,PDstandardNth23g11,kmul(J31L,PDstandardNth33g11)),kmul(dJ313L,PDstandardNth3g11)))))));
      
      JacPDstandardNth31g12 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g12,kmadd(J21L,PDstandardNth12g12,kmul(J31L,PDstandardNth13g12))),kmadd(J11L,kmadd(J23L,PDstandardNth12g12,kmul(J33L,PDstandardNth13g12)),kmadd(dJ113L,PDstandardNth1g12,kmadd(J23L,kmadd(J21L,PDstandardNth22g12,kmul(J31L,PDstandardNth23g12)),kmadd(dJ213L,PDstandardNth2g12,kmadd(J33L,kmadd(J21L,PDstandardNth23g12,kmul(J31L,PDstandardNth33g12)),kmul(dJ313L,PDstandardNth3g12)))))));
      
      JacPDstandardNth31g13 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g13,kmadd(J21L,PDstandardNth12g13,kmul(J31L,PDstandardNth13g13))),kmadd(J11L,kmadd(J23L,PDstandardNth12g13,kmul(J33L,PDstandardNth13g13)),kmadd(dJ113L,PDstandardNth1g13,kmadd(J23L,kmadd(J21L,PDstandardNth22g13,kmul(J31L,PDstandardNth23g13)),kmadd(dJ213L,PDstandardNth2g13,kmadd(J33L,kmadd(J21L,PDstandardNth23g13,kmul(J31L,PDstandardNth33g13)),kmul(dJ313L,PDstandardNth3g13)))))));
      
      JacPDstandardNth31g22 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g22,kmadd(J21L,PDstandardNth12g22,kmul(J31L,PDstandardNth13g22))),kmadd(J11L,kmadd(J23L,PDstandardNth12g22,kmul(J33L,PDstandardNth13g22)),kmadd(dJ113L,PDstandardNth1g22,kmadd(J23L,kmadd(J21L,PDstandardNth22g22,kmul(J31L,PDstandardNth23g22)),kmadd(dJ213L,PDstandardNth2g22,kmadd(J33L,kmadd(J21L,PDstandardNth23g22,kmul(J31L,PDstandardNth33g22)),kmul(dJ313L,PDstandardNth3g22)))))));
      
      JacPDstandardNth31g23 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g23,kmadd(J21L,PDstandardNth12g23,kmul(J31L,PDstandardNth13g23))),kmadd(J11L,kmadd(J23L,PDstandardNth12g23,kmul(J33L,PDstandardNth13g23)),kmadd(dJ113L,PDstandardNth1g23,kmadd(J23L,kmadd(J21L,PDstandardNth22g23,kmul(J31L,PDstandardNth23g23)),kmadd(dJ213L,PDstandardNth2g23,kmadd(J33L,kmadd(J21L,PDstandardNth23g23,kmul(J31L,PDstandardNth33g23)),kmul(dJ313L,PDstandardNth3g23)))))));
      
      JacPDstandardNth31g33 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g33,kmadd(J21L,PDstandardNth12g33,kmul(J31L,PDstandardNth13g33))),kmadd(J11L,kmadd(J23L,PDstandardNth12g33,kmul(J33L,PDstandardNth13g33)),kmadd(dJ113L,PDstandardNth1g33,kmadd(J23L,kmadd(J21L,PDstandardNth22g33,kmul(J31L,PDstandardNth23g33)),kmadd(dJ213L,PDstandardNth2g33,kmadd(J33L,kmadd(J21L,PDstandardNth23g33,kmul(J31L,PDstandardNth33g33)),kmul(dJ313L,PDstandardNth3g33)))))));
      
      JacPDstandardNth32g11 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g11,kmadd(J22L,PDstandardNth12g11,kmul(J32L,PDstandardNth13g11))),kmadd(J12L,kmadd(J23L,PDstandardNth12g11,kmul(J33L,PDstandardNth13g11)),kmadd(dJ123L,PDstandardNth1g11,kmadd(J23L,kmadd(J22L,PDstandardNth22g11,kmul(J32L,PDstandardNth23g11)),kmadd(dJ223L,PDstandardNth2g11,kmadd(J33L,kmadd(J22L,PDstandardNth23g11,kmul(J32L,PDstandardNth33g11)),kmul(dJ323L,PDstandardNth3g11)))))));
      
      JacPDstandardNth32g12 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g12,kmadd(J22L,PDstandardNth12g12,kmul(J32L,PDstandardNth13g12))),kmadd(J12L,kmadd(J23L,PDstandardNth12g12,kmul(J33L,PDstandardNth13g12)),kmadd(dJ123L,PDstandardNth1g12,kmadd(J23L,kmadd(J22L,PDstandardNth22g12,kmul(J32L,PDstandardNth23g12)),kmadd(dJ223L,PDstandardNth2g12,kmadd(J33L,kmadd(J22L,PDstandardNth23g12,kmul(J32L,PDstandardNth33g12)),kmul(dJ323L,PDstandardNth3g12)))))));
      
      JacPDstandardNth32g13 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g13,kmadd(J22L,PDstandardNth12g13,kmul(J32L,PDstandardNth13g13))),kmadd(J12L,kmadd(J23L,PDstandardNth12g13,kmul(J33L,PDstandardNth13g13)),kmadd(dJ123L,PDstandardNth1g13,kmadd(J23L,kmadd(J22L,PDstandardNth22g13,kmul(J32L,PDstandardNth23g13)),kmadd(dJ223L,PDstandardNth2g13,kmadd(J33L,kmadd(J22L,PDstandardNth23g13,kmul(J32L,PDstandardNth33g13)),kmul(dJ323L,PDstandardNth3g13)))))));
      
      JacPDstandardNth32g22 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g22,kmadd(J22L,PDstandardNth12g22,kmul(J32L,PDstandardNth13g22))),kmadd(J12L,kmadd(J23L,PDstandardNth12g22,kmul(J33L,PDstandardNth13g22)),kmadd(dJ123L,PDstandardNth1g22,kmadd(J23L,kmadd(J22L,PDstandardNth22g22,kmul(J32L,PDstandardNth23g22)),kmadd(dJ223L,PDstandardNth2g22,kmadd(J33L,kmadd(J22L,PDstandardNth23g22,kmul(J32L,PDstandardNth33g22)),kmul(dJ323L,PDstandardNth3g22)))))));
      
      JacPDstandardNth32g23 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g23,kmadd(J22L,PDstandardNth12g23,kmul(J32L,PDstandardNth13g23))),kmadd(J12L,kmadd(J23L,PDstandardNth12g23,kmul(J33L,PDstandardNth13g23)),kmadd(dJ123L,PDstandardNth1g23,kmadd(J23L,kmadd(J22L,PDstandardNth22g23,kmul(J32L,PDstandardNth23g23)),kmadd(dJ223L,PDstandardNth2g23,kmadd(J33L,kmadd(J22L,PDstandardNth23g23,kmul(J32L,PDstandardNth33g23)),kmul(dJ323L,PDstandardNth3g23)))))));
      
      JacPDstandardNth32g33 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g33,kmadd(J22L,PDstandardNth12g33,kmul(J32L,PDstandardNth13g33))),kmadd(J12L,kmadd(J23L,PDstandardNth12g33,kmul(J33L,PDstandardNth13g33)),kmadd(dJ123L,PDstandardNth1g33,kmadd(J23L,kmadd(J22L,PDstandardNth22g33,kmul(J32L,PDstandardNth23g33)),kmadd(dJ223L,PDstandardNth2g33,kmadd(J33L,kmadd(J22L,PDstandardNth23g33,kmul(J32L,PDstandardNth33g33)),kmul(dJ323L,PDstandardNth3g33)))))));
    }
    else
    {
      JacPDstandardNth1alpha = PDstandardNth1alpha;
      
      JacPDstandardNth1beta1 = PDstandardNth1beta1;
      
      JacPDstandardNth1beta2 = PDstandardNth1beta2;
      
      JacPDstandardNth1beta3 = PDstandardNth1beta3;
      
      JacPDstandardNth1g11 = PDstandardNth1g11;
      
      JacPDstandardNth1g12 = PDstandardNth1g12;
      
      JacPDstandardNth1g13 = PDstandardNth1g13;
      
      JacPDstandardNth1g22 = PDstandardNth1g22;
      
      JacPDstandardNth1g23 = PDstandardNth1g23;
      
      JacPDstandardNth1g33 = PDstandardNth1g33;
      
      JacPDstandardNth1K11 = PDstandardNth1K11;
      
      JacPDstandardNth1K12 = PDstandardNth1K12;
      
      JacPDstandardNth1K13 = PDstandardNth1K13;
      
      JacPDstandardNth1K22 = PDstandardNth1K22;
      
      JacPDstandardNth1K23 = PDstandardNth1K23;
      
      JacPDstandardNth1K33 = PDstandardNth1K33;
      
      JacPDstandardNth2alpha = PDstandardNth2alpha;
      
      JacPDstandardNth2beta1 = PDstandardNth2beta1;
      
      JacPDstandardNth2beta2 = PDstandardNth2beta2;
      
      JacPDstandardNth2beta3 = PDstandardNth2beta3;
      
      JacPDstandardNth2g11 = PDstandardNth2g11;
      
      JacPDstandardNth2g12 = PDstandardNth2g12;
      
      JacPDstandardNth2g13 = PDstandardNth2g13;
      
      JacPDstandardNth2g22 = PDstandardNth2g22;
      
      JacPDstandardNth2g23 = PDstandardNth2g23;
      
      JacPDstandardNth2g33 = PDstandardNth2g33;
      
      JacPDstandardNth2K11 = PDstandardNth2K11;
      
      JacPDstandardNth2K12 = PDstandardNth2K12;
      
      JacPDstandardNth2K13 = PDstandardNth2K13;
      
      JacPDstandardNth2K22 = PDstandardNth2K22;
      
      JacPDstandardNth2K23 = PDstandardNth2K23;
      
      JacPDstandardNth2K33 = PDstandardNth2K33;
      
      JacPDstandardNth3alpha = PDstandardNth3alpha;
      
      JacPDstandardNth3beta1 = PDstandardNth3beta1;
      
      JacPDstandardNth3beta2 = PDstandardNth3beta2;
      
      JacPDstandardNth3beta3 = PDstandardNth3beta3;
      
      JacPDstandardNth3g11 = PDstandardNth3g11;
      
      JacPDstandardNth3g12 = PDstandardNth3g12;
      
      JacPDstandardNth3g13 = PDstandardNth3g13;
      
      JacPDstandardNth3g22 = PDstandardNth3g22;
      
      JacPDstandardNth3g23 = PDstandardNth3g23;
      
      JacPDstandardNth3g33 = PDstandardNth3g33;
      
      JacPDstandardNth3K11 = PDstandardNth3K11;
      
      JacPDstandardNth3K12 = PDstandardNth3K12;
      
      JacPDstandardNth3K13 = PDstandardNth3K13;
      
      JacPDstandardNth3K22 = PDstandardNth3K22;
      
      JacPDstandardNth3K23 = PDstandardNth3K23;
      
      JacPDstandardNth3K33 = PDstandardNth3K33;
      
      JacPDstandardNth11alpha = PDstandardNth11alpha;
      
      JacPDstandardNth11g22 = PDstandardNth11g22;
      
      JacPDstandardNth11g23 = PDstandardNth11g23;
      
      JacPDstandardNth11g33 = PDstandardNth11g33;
      
      JacPDstandardNth22alpha = PDstandardNth22alpha;
      
      JacPDstandardNth22g11 = PDstandardNth22g11;
      
      JacPDstandardNth22g13 = PDstandardNth22g13;
      
      JacPDstandardNth22g33 = PDstandardNth22g33;
      
      JacPDstandardNth33alpha = PDstandardNth33alpha;
      
      JacPDstandardNth33g11 = PDstandardNth33g11;
      
      JacPDstandardNth33g12 = PDstandardNth33g12;
      
      JacPDstandardNth33g22 = PDstandardNth33g22;
      
      JacPDstandardNth12alpha = PDstandardNth12alpha;
      
      JacPDstandardNth12g11 = PDstandardNth12g11;
      
      JacPDstandardNth12g12 = PDstandardNth12g12;
      
      JacPDstandardNth12g13 = PDstandardNth12g13;
      
      JacPDstandardNth12g22 = PDstandardNth12g22;
      
      JacPDstandardNth12g23 = PDstandardNth12g23;
      
      JacPDstandardNth12g33 = PDstandardNth12g33;
      
      JacPDstandardNth13alpha = PDstandardNth13alpha;
      
      JacPDstandardNth13g11 = PDstandardNth13g11;
      
      JacPDstandardNth13g12 = PDstandardNth13g12;
      
      JacPDstandardNth13g13 = PDstandardNth13g13;
      
      JacPDstandardNth13g22 = PDstandardNth13g22;
      
      JacPDstandardNth13g23 = PDstandardNth13g23;
      
      JacPDstandardNth13g33 = PDstandardNth13g33;
      
      JacPDstandardNth21g11 = PDstandardNth12g11;
      
      JacPDstandardNth21g12 = PDstandardNth12g12;
      
      JacPDstandardNth21g13 = PDstandardNth12g13;
      
      JacPDstandardNth21g22 = PDstandardNth12g22;
      
      JacPDstandardNth21g23 = PDstandardNth12g23;
      
      JacPDstandardNth21g33 = PDstandardNth12g33;
      
      JacPDstandardNth23alpha = PDstandardNth23alpha;
      
      JacPDstandardNth23g11 = PDstandardNth23g11;
      
      JacPDstandardNth23g12 = PDstandardNth23g12;
      
      JacPDstandardNth23g13 = PDstandardNth23g13;
      
      JacPDstandardNth23g22 = PDstandardNth23g22;
      
      JacPDstandardNth23g23 = PDstandardNth23g23;
      
      JacPDstandardNth23g33 = PDstandardNth23g33;
      
      JacPDstandardNth31g11 = PDstandardNth13g11;
      
      JacPDstandardNth31g12 = PDstandardNth13g12;
      
      JacPDstandardNth31g13 = PDstandardNth13g13;
      
      JacPDstandardNth31g22 = PDstandardNth13g22;
      
      JacPDstandardNth31g23 = PDstandardNth13g23;
      
      JacPDstandardNth31g33 = PDstandardNth13g33;
      
      JacPDstandardNth32g11 = PDstandardNth23g11;
      
      JacPDstandardNth32g12 = PDstandardNth23g12;
      
      JacPDstandardNth32g13 = PDstandardNth23g13;
      
      JacPDstandardNth32g22 = PDstandardNth23g22;
      
      JacPDstandardNth32g23 = PDstandardNth23g23;
      
      JacPDstandardNth32g33 = PDstandardNth23g33;
    }
    
    CCTK_REAL_VEC detg CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(2),kmul(g12L,kmul(g13L,g23L)),knmsub(g33L,kmul(g12L,g12L),kmsub(g22L,kmsub(g11L,g33L,kmul(g13L,g13L)),kmul(g11L,kmul(g23L,g23L)))));
    
    CCTK_REAL_VEC gu11 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g22L,g33L,kmul(g23L,g23L)),detg);
    
    CCTK_REAL_VEC gu12 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g13L,g23L,kmul(g12L,g33L)),detg);
    
    CCTK_REAL_VEC gu13 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g12L,g23L,kmul(g13L,g22L)),detg);
    
    CCTK_REAL_VEC gu21 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g13L,g23L,kmul(g12L,g33L)),detg);
    
    CCTK_REAL_VEC gu22 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g11L,g33L,kmul(g13L,g13L)),detg);
    
    CCTK_REAL_VEC gu23 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g12L,g13L,kmul(g11L,g23L)),detg);
    
    CCTK_REAL_VEC gu31 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g12L,g23L,kmul(g13L,g22L)),detg);
    
    CCTK_REAL_VEC gu32 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g12L,g13L,kmul(g11L,g23L)),detg);
    
    CCTK_REAL_VEC gu33 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g11L,g22L,kmul(g12L,g12L)),detg);
    
    CCTK_REAL_VEC G111 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu11,JacPDstandardNth1g11,kmsub(ToReal(2),kmadd(gu12,JacPDstandardNth1g12,kmul(gu13,JacPDstandardNth1g13)),kmadd(gu13,JacPDstandardNth3g11,kmul(gu12,JacPDstandardNth2g11)))),ToReal(0.5));
    
    CCTK_REAL_VEC G211 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu21,JacPDstandardNth1g11,kmsub(ToReal(2),kmadd(gu22,JacPDstandardNth1g12,kmul(gu23,JacPDstandardNth1g13)),kmadd(gu23,JacPDstandardNth3g11,kmul(gu22,JacPDstandardNth2g11)))),ToReal(0.5));
    
    CCTK_REAL_VEC G311 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu31,JacPDstandardNth1g11,kmsub(ToReal(2),kmadd(gu32,JacPDstandardNth1g12,kmul(gu33,JacPDstandardNth1g13)),kmadd(gu33,JacPDstandardNth3g11,kmul(gu32,JacPDstandardNth2g11)))),ToReal(0.5));
    
    CCTK_REAL_VEC G112 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu12,JacPDstandardNth1g22,kmadd(gu11,JacPDstandardNth2g11,kmul(gu13,kadd(JacPDstandardNth1g23,ksub(JacPDstandardNth2g13,JacPDstandardNth3g12))))),ToReal(0.5));
    
    CCTK_REAL_VEC G212 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu22,JacPDstandardNth1g22,kmadd(gu21,JacPDstandardNth2g11,kmul(gu23,kadd(JacPDstandardNth1g23,ksub(JacPDstandardNth2g13,JacPDstandardNth3g12))))),ToReal(0.5));
    
    CCTK_REAL_VEC G312 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu32,JacPDstandardNth1g22,kmadd(gu31,JacPDstandardNth2g11,kmul(gu33,kadd(JacPDstandardNth1g23,ksub(JacPDstandardNth2g13,JacPDstandardNth3g12))))),ToReal(0.5));
    
    CCTK_REAL_VEC G113 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu13,JacPDstandardNth1g33,kmadd(gu11,JacPDstandardNth3g11,kmul(gu12,kadd(JacPDstandardNth1g23,ksub(JacPDstandardNth3g12,JacPDstandardNth2g13))))),ToReal(0.5));
    
    CCTK_REAL_VEC G213 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu23,JacPDstandardNth1g33,kmadd(gu21,JacPDstandardNth3g11,kmul(gu22,kadd(JacPDstandardNth1g23,ksub(JacPDstandardNth3g12,JacPDstandardNth2g13))))),ToReal(0.5));
    
    CCTK_REAL_VEC G313 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu33,JacPDstandardNth1g33,kmadd(gu31,JacPDstandardNth3g11,kmul(gu32,kadd(JacPDstandardNth1g23,ksub(JacPDstandardNth3g12,JacPDstandardNth2g13))))),ToReal(0.5));
    
    CCTK_REAL_VEC G122 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu11,kmsub(ToReal(2),JacPDstandardNth2g12,JacPDstandardNth1g22),kmadd(gu12,JacPDstandardNth2g22,kmul(gu13,kmsub(ToReal(2),JacPDstandardNth2g23,JacPDstandardNth3g22)))),ToReal(0.5));
    
    CCTK_REAL_VEC G222 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu21,kmsub(ToReal(2),JacPDstandardNth2g12,JacPDstandardNth1g22),kmadd(gu22,JacPDstandardNth2g22,kmul(gu23,kmsub(ToReal(2),JacPDstandardNth2g23,JacPDstandardNth3g22)))),ToReal(0.5));
    
    CCTK_REAL_VEC G322 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu31,kmsub(ToReal(2),JacPDstandardNth2g12,JacPDstandardNth1g22),kmadd(gu32,JacPDstandardNth2g22,kmul(gu33,kmsub(ToReal(2),JacPDstandardNth2g23,JacPDstandardNth3g22)))),ToReal(0.5));
    
    CCTK_REAL_VEC G123 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu13,JacPDstandardNth2g33,kmadd(gu11,ksub(kadd(JacPDstandardNth2g13,JacPDstandardNth3g12),JacPDstandardNth1g23),kmul(gu12,JacPDstandardNth3g22))),ToReal(0.5));
    
    CCTK_REAL_VEC G223 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu23,JacPDstandardNth2g33,kmadd(gu21,ksub(kadd(JacPDstandardNth2g13,JacPDstandardNth3g12),JacPDstandardNth1g23),kmul(gu22,JacPDstandardNth3g22))),ToReal(0.5));
    
    CCTK_REAL_VEC G323 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu33,JacPDstandardNth2g33,kmadd(gu31,ksub(kadd(JacPDstandardNth2g13,JacPDstandardNth3g12),JacPDstandardNth1g23),kmul(gu32,JacPDstandardNth3g22))),ToReal(0.5));
    
    CCTK_REAL_VEC G133 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu11,kmsub(ToReal(2),JacPDstandardNth3g13,JacPDstandardNth1g33),kmadd(gu12,kmsub(ToReal(2),JacPDstandardNth3g23,JacPDstandardNth2g33),kmul(gu13,JacPDstandardNth3g33))),ToReal(0.5));
    
    CCTK_REAL_VEC G233 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu21,kmsub(ToReal(2),JacPDstandardNth3g13,JacPDstandardNth1g33),kmadd(gu22,kmsub(ToReal(2),JacPDstandardNth3g23,JacPDstandardNth2g33),kmul(gu23,JacPDstandardNth3g33))),ToReal(0.5));
    
    CCTK_REAL_VEC G333 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gu31,kmsub(ToReal(2),JacPDstandardNth3g13,JacPDstandardNth1g33),kmadd(gu32,kmsub(ToReal(2),JacPDstandardNth3g23,JacPDstandardNth2g33),kmul(gu33,JacPDstandardNth3g33))),ToReal(0.5));
    
    CCTK_REAL_VEC R11 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(2),kmadd(G112,G211,kmadd(G113,G311,kmul(G213,G312))),knmsub(G111,kadd(G111,kadd(G212,G313)),knmsub(G211,kadd(G112,kadd(G222,G323)),knmsub(G311,kadd(G113,kadd(G223,G333)),kmadd(ToReal(0.5),kmadd(gu12,ksub(JacPDstandardNth21g11,JacPDstandardNth12g11),knmsub(gu22,kadd(JacPDstandardNth11g22,kmadd(ToReal(-2),JacPDstandardNth21g12,JacPDstandardNth22g11)),kmadd(gu13,ksub(JacPDstandardNth31g11,JacPDstandardNth13g11),kmadd(gu23,ksub(kadd(JacPDstandardNth21g13,ksub(JacPDstandardNth31g12,JacPDstandardNth23g11)),JacPDstandardNth11g23),kmsub(gu32,ksub(kadd(JacPDstandardNth21g13,ksub(JacPDstandardNth31g12,JacPDstandardNth32g11)),JacPDstandardNth11g23),kmul(gu33,kadd(JacPDstandardNth11g33,kmadd(ToReal(-2),JacPDstandardNth31g13,JacPDstandardNth33g11)))))))),kmadd(G111,G111,kmadd(G212,G212,kmul(G313,G313))))))));
    
    CCTK_REAL_VEC R12 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(2),kmadd(G122,G211,kmadd(G123,G311,kmadd(G213,G322,kmul(G313,G323)))),kmadd(ToReal(-2),kmadd(G112,kadd(G212,G313),kmadd(G212,G323,kmadd(G312,G333,kmul(gu12,JacPDstandardNth12g12)))),knmsub(kadd(gu32,gu23),JacPDstandardNth12g23,kmadd(gu13,ksub(JacPDstandardNth11g23,kadd(JacPDstandardNth13g12,JacPDstandardNth12g13)),kmadd(gu22,ksub(JacPDstandardNth21g22,JacPDstandardNth12g22),kmadd(gu12,kadd(JacPDstandardNth11g22,JacPDstandardNth22g11),kmadd(gu23,ksub(JacPDstandardNth21g23,JacPDstandardNth23g12),kmadd(gu32,JacPDstandardNth31g22,kmadd(gu13,JacPDstandardNth32g11,kmadd(gu32,ksub(JacPDstandardNth22g13,JacPDstandardNth32g12),kmadd(gu23,JacPDstandardNth32g12,kmadd(gu33,JacPDstandardNth32g13,kmul(gu33,ksub(ksub(JacPDstandardNth31g23,JacPDstandardNth33g12),JacPDstandardNth12g33)))))))))))))),ToReal(0.5));
    
    CCTK_REAL_VEC R13 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(2),kmadd(G123,G211,kmadd(G212,G223,kmadd(G133,G311,kmul(G233,G312)))),kmadd(ToReal(-2),kmadd(G213,G222,kmadd(G223,G313,kmadd(G113,kadd(G212,G313),kmul(gu13,JacPDstandardNth13g13)))),knmsub(kadd(gu32,gu23),JacPDstandardNth13g23,kmadd(gu22,ksub(ksub(JacPDstandardNth21g23,JacPDstandardNth22g13),JacPDstandardNth13g22),kmadd(gu12,kadd(JacPDstandardNth11g23,ksub(ksub(JacPDstandardNth23g11,JacPDstandardNth13g12),JacPDstandardNth12g13)),kmadd(gu22,JacPDstandardNth23g12,kmadd(gu23,ksub(JacPDstandardNth21g33,JacPDstandardNth23g13),kmadd(gu32,JacPDstandardNth31g23,kmadd(gu33,ksub(JacPDstandardNth31g33,JacPDstandardNth13g33),kmadd(gu32,ksub(JacPDstandardNth23g13,JacPDstandardNth32g13),kmadd(gu13,kadd(JacPDstandardNth11g33,JacPDstandardNth33g11),kmul(gu23,JacPDstandardNth33g12)))))))))))),ToReal(0.5));
    
    CCTK_REAL_VEC R22 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(G122,kadd(G111,kadd(G212,G313)),kmadd(ToReal(2),kmadd(G122,G212,kmadd(G123,G312,kmul(G223,G322))),knmsub(G222,kadd(G112,kadd(G222,G323)),knmsub(G322,kadd(G113,kadd(G223,G333)),kmadd(ToReal(0.5),kmadd(gu21,ksub(JacPDstandardNth12g22,JacPDstandardNth21g22),knmsub(gu11,kadd(JacPDstandardNth11g22,kmadd(ToReal(-2),JacPDstandardNth12g12,JacPDstandardNth22g11)),kmadd(gu13,kadd(JacPDstandardNth12g23,ksub(ksub(JacPDstandardNth32g12,JacPDstandardNth22g13),JacPDstandardNth13g22)),kmadd(gu31,kadd(JacPDstandardNth12g23,ksub(ksub(JacPDstandardNth32g12,JacPDstandardNth31g22),JacPDstandardNth22g13)),kmsub(gu23,ksub(JacPDstandardNth32g22,JacPDstandardNth23g22),kmul(gu33,kadd(JacPDstandardNth22g33,kmadd(ToReal(-2),JacPDstandardNth32g23,JacPDstandardNth33g22)))))))),kmadd(G112,G112,kmadd(G222,G222,kmul(G323,G323))))))));
    
    CCTK_REAL_VEC R23 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(2),kmadd(G122,G213,kmadd(G112,ksub(G113,G223),kmadd(G133,G312,kmul(G233,G322)))),kmadd(gu11,ksub(kadd(JacPDstandardNth12g13,ksub(JacPDstandardNth13g12,JacPDstandardNth23g11)),JacPDstandardNth11g23),kmadd(gu21,kadd(JacPDstandardNth13g22,ksub(ksub(JacPDstandardNth22g13,JacPDstandardNth23g12),JacPDstandardNth21g23)),kmadd(ToReal(-2),kmadd(G111,G123,kmadd(kadd(G113,G223),G323,kmul(gu23,JacPDstandardNth23g23))),kmadd(gu31,kadd(JacPDstandardNth13g23,ksub(ksub(JacPDstandardNth32g13,JacPDstandardNth31g23),JacPDstandardNth23g13)),kmadd(gu33,ksub(JacPDstandardNth32g33,JacPDstandardNth23g33),kmadd(gu13,kadd(JacPDstandardNth12g33,ksub(ksub(JacPDstandardNth33g12,JacPDstandardNth23g13),JacPDstandardNth13g23)),kmul(gu23,kadd(JacPDstandardNth22g33,JacPDstandardNth33g22))))))))),ToReal(0.5));
    
    CCTK_REAL_VEC R33 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(G133,kadd(G111,kadd(G212,G313)),kmadd(ToReal(2),kmadd(G123,G213,kmul(G133,G313)),kmadd(G233,ksub(ksub(G323,G222),G112),knmsub(G333,kadd(G113,kadd(G223,G333)),kmadd(ToReal(0.5),kmadd(gu31,ksub(JacPDstandardNth13g33,JacPDstandardNth31g33),kmadd(gu32,ksub(JacPDstandardNth23g33,JacPDstandardNth32g33),knmsub(gu11,kadd(JacPDstandardNth11g33,kmadd(ToReal(-2),JacPDstandardNth13g13,JacPDstandardNth33g11)),kmadd(gu12,ksub(kadd(JacPDstandardNth13g23,ksub(JacPDstandardNth23g13,JacPDstandardNth33g12)),JacPDstandardNth12g33),kmsub(gu21,kadd(JacPDstandardNth13g23,ksub(ksub(JacPDstandardNth23g13,JacPDstandardNth33g12),JacPDstandardNth21g33)),kmul(gu22,kadd(JacPDstandardNth22g33,kmadd(ToReal(-2),JacPDstandardNth23g23,JacPDstandardNth33g22)))))))),kmadd(G113,G113,kmadd(G223,G223,kmul(G333,G333))))))));
    
    CCTK_REAL_VEC Km11 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(K11L,gu11,kmadd(K12L,gu12,kmul(K13L,gu13)));
    
    CCTK_REAL_VEC Km21 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(K11L,gu21,kmadd(K12L,gu22,kmul(K13L,gu23)));
    
    CCTK_REAL_VEC Km31 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(K11L,gu31,kmadd(K12L,gu32,kmul(K13L,gu33)));
    
    CCTK_REAL_VEC Km12 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(K12L,gu11,kmadd(K22L,gu12,kmul(K23L,gu13)));
    
    CCTK_REAL_VEC Km22 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(K12L,gu21,kmadd(K22L,gu22,kmul(K23L,gu23)));
    
    CCTK_REAL_VEC Km32 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(K12L,gu31,kmadd(K22L,gu32,kmul(K23L,gu33)));
    
    CCTK_REAL_VEC Km13 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(K13L,gu11,kmadd(K23L,gu12,kmul(K33L,gu13)));
    
    CCTK_REAL_VEC Km23 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(K13L,gu21,kmadd(K23L,gu22,kmul(K33L,gu23)));
    
    CCTK_REAL_VEC Km33 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(K13L,gu31,kmadd(K23L,gu32,kmul(K33L,gu33)));
    
    CCTK_REAL_VEC trK CCTK_ATTRIBUTE_UNUSED = kadd(Km11,kadd(Km22,Km33));
    
    CCTK_REAL_VEC g11rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmul(alphaL,K11L),kmadd(ToReal(2),kmadd(g11L,JacPDstandardNth1beta1,kmadd(g12L,JacPDstandardNth1beta2,kmul(g13L,JacPDstandardNth1beta3))),kmadd(beta1L,JacPDstandardNth1g11,kmadd(beta2L,JacPDstandardNth2g11,kmul(beta3L,JacPDstandardNth3g11)))));
    
    CCTK_REAL_VEC g12rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmul(alphaL,K12L),kmadd(g22L,JacPDstandardNth1beta2,kmadd(g23L,JacPDstandardNth1beta3,kmadd(beta1L,JacPDstandardNth1g12,kmadd(g11L,JacPDstandardNth2beta1,kmadd(g12L,kadd(JacPDstandardNth1beta1,JacPDstandardNth2beta2),kmadd(g13L,JacPDstandardNth2beta3,kmadd(beta2L,JacPDstandardNth2g12,kmul(beta3L,JacPDstandardNth3g12)))))))));
    
    CCTK_REAL_VEC g13rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmul(alphaL,K13L),kmadd(g23L,JacPDstandardNth1beta2,kmadd(g33L,JacPDstandardNth1beta3,kmadd(beta1L,JacPDstandardNth1g13,kmadd(beta2L,JacPDstandardNth2g13,kmadd(g11L,JacPDstandardNth3beta1,kmadd(g12L,JacPDstandardNth3beta2,kmadd(g13L,kadd(JacPDstandardNth1beta1,JacPDstandardNth3beta3),kmul(beta3L,JacPDstandardNth3g13)))))))));
    
    CCTK_REAL_VEC g22rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmul(alphaL,K22L),kmadd(beta1L,JacPDstandardNth1g22,kmadd(ToReal(2),kmadd(g12L,JacPDstandardNth2beta1,kmadd(g22L,JacPDstandardNth2beta2,kmul(g23L,JacPDstandardNth2beta3))),kmadd(beta2L,JacPDstandardNth2g22,kmul(beta3L,JacPDstandardNth3g22)))));
    
    CCTK_REAL_VEC g23rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmul(alphaL,K23L),kmadd(beta1L,JacPDstandardNth1g23,kmadd(g13L,JacPDstandardNth2beta1,kmadd(g33L,JacPDstandardNth2beta3,kmadd(beta2L,JacPDstandardNth2g23,kmadd(g12L,JacPDstandardNth3beta1,kmadd(g22L,JacPDstandardNth3beta2,kmadd(g23L,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3),kmul(beta3L,JacPDstandardNth3g23)))))))));
    
    CCTK_REAL_VEC g33rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmul(alphaL,K33L),kmadd(beta1L,JacPDstandardNth1g33,kmadd(beta2L,JacPDstandardNth2g33,kmadd(ToReal(2),kmadd(g13L,JacPDstandardNth3beta1,kmadd(g23L,JacPDstandardNth3beta2,kmul(g33L,JacPDstandardNth3beta3))),kmul(beta3L,JacPDstandardNth3g33)))));
    
    CCTK_REAL_VEC K11rhsL CCTK_ATTRIBUTE_UNUSED = 
      ksub(kmadd(G111,JacPDstandardNth1alpha,kmadd(ToReal(2),kmadd(K11L,JacPDstandardNth1beta1,kmadd(K12L,JacPDstandardNth1beta2,kmul(K13L,JacPDstandardNth1beta3))),kmadd(beta1L,JacPDstandardNth1K11,kmadd(G211,JacPDstandardNth2alpha,kmadd(beta2L,JacPDstandardNth2K11,kmadd(G311,JacPDstandardNth3alpha,kmadd(beta3L,JacPDstandardNth3K11,kmul(alphaL,kmadd(ToReal(-2),kmadd(K11L,Km11,kmadd(K12L,Km21,kmul(K13L,Km31))),kmadd(K11L,trK,R11)))))))))),JacPDstandardNth11alpha);
    
    CCTK_REAL_VEC K12rhsL CCTK_ATTRIBUTE_UNUSED = 
      ksub(kmadd(G112,JacPDstandardNth1alpha,kmadd(K22L,JacPDstandardNth1beta2,kmadd(K23L,JacPDstandardNth1beta3,kmadd(beta1L,JacPDstandardNth1K12,kmadd(G212,JacPDstandardNth2alpha,kmadd(K11L,JacPDstandardNth2beta1,kmadd(K12L,kadd(JacPDstandardNth1beta1,JacPDstandardNth2beta2),kmadd(K13L,JacPDstandardNth2beta3,kmadd(beta2L,JacPDstandardNth2K12,kmadd(G312,JacPDstandardNth3alpha,kmadd(beta3L,JacPDstandardNth3K12,kmul(alphaL,kmadd(ToReal(-2),kmadd(K11L,Km12,kmadd(K12L,Km22,kmul(K13L,Km32))),kmadd(K12L,trK,R12)))))))))))))),JacPDstandardNth12alpha);
    
    CCTK_REAL_VEC K13rhsL CCTK_ATTRIBUTE_UNUSED = 
      ksub(kmadd(G113,JacPDstandardNth1alpha,kmadd(K23L,JacPDstandardNth1beta2,kmadd(K33L,JacPDstandardNth1beta3,kmadd(beta1L,JacPDstandardNth1K13,kmadd(G213,JacPDstandardNth2alpha,kmadd(beta2L,JacPDstandardNth2K13,kmadd(G313,JacPDstandardNth3alpha,kmadd(K11L,JacPDstandardNth3beta1,kmadd(K12L,JacPDstandardNth3beta2,kmadd(K13L,kadd(JacPDstandardNth1beta1,JacPDstandardNth3beta3),kmadd(beta3L,JacPDstandardNth3K13,kmul(alphaL,kmadd(ToReal(-2),kmadd(K11L,Km13,kmadd(K12L,Km23,kmul(K13L,Km33))),kmadd(K13L,trK,R13)))))))))))))),JacPDstandardNth13alpha);
    
    CCTK_REAL_VEC K22rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(G122,JacPDstandardNth1alpha,kmadd(beta1L,JacPDstandardNth1K22,ksub(kmadd(G222,JacPDstandardNth2alpha,kmadd(ToReal(2),kmadd(K12L,JacPDstandardNth2beta1,kmadd(K22L,JacPDstandardNth2beta2,kmul(K23L,JacPDstandardNth2beta3))),kmadd(beta2L,JacPDstandardNth2K22,kmadd(G322,JacPDstandardNth3alpha,kmadd(beta3L,JacPDstandardNth3K22,kmul(alphaL,kmadd(ToReal(-2),kmadd(K12L,Km12,kmadd(K22L,Km22,kmul(K23L,Km32))),kmadd(K22L,trK,R22)))))))),JacPDstandardNth22alpha)));
    
    CCTK_REAL_VEC K23rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(G123,JacPDstandardNth1alpha,kmadd(beta1L,JacPDstandardNth1K23,ksub(kmadd(G223,JacPDstandardNth2alpha,kmadd(K13L,JacPDstandardNth2beta1,kmadd(K33L,JacPDstandardNth2beta3,kmadd(beta2L,JacPDstandardNth2K23,kmadd(G323,JacPDstandardNth3alpha,kmadd(K12L,JacPDstandardNth3beta1,kmadd(K22L,JacPDstandardNth3beta2,kmadd(K23L,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3),kmadd(beta3L,JacPDstandardNth3K23,kmul(alphaL,kmadd(ToReal(-2),kmadd(K12L,Km13,kmadd(K22L,Km23,kmul(K23L,Km33))),kmadd(K23L,trK,R23)))))))))))),JacPDstandardNth23alpha)));
    
    CCTK_REAL_VEC K33rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(G133,JacPDstandardNth1alpha,kmadd(beta1L,JacPDstandardNth1K33,kmadd(G233,JacPDstandardNth2alpha,kmadd(beta2L,JacPDstandardNth2K33,ksub(kmadd(G333,JacPDstandardNth3alpha,kmadd(ToReal(2),kmadd(K13L,JacPDstandardNth3beta1,kmadd(K23L,JacPDstandardNth3beta2,kmul(K33L,JacPDstandardNth3beta3))),kmadd(beta3L,JacPDstandardNth3K33,kmul(alphaL,kmadd(ToReal(-2),kmadd(K13L,Km13,kmadd(K23L,Km23,kmul(K33L,Km33))),kmadd(K33L,trK,R33)))))),JacPDstandardNth33alpha)))));
    
    CCTK_REAL_VEC alpharhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC beta1rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC beta2rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC beta3rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(alpharhs[index],alpharhsL);
    vec_store_nta_partial(beta1rhs[index],beta1rhsL);
    vec_store_nta_partial(beta2rhs[index],beta2rhsL);
    vec_store_nta_partial(beta3rhs[index],beta3rhsL);
    vec_store_nta_partial(g11rhs[index],g11rhsL);
    vec_store_nta_partial(g12rhs[index],g12rhsL);
    vec_store_nta_partial(g13rhs[index],g13rhsL);
    vec_store_nta_partial(g22rhs[index],g22rhsL);
    vec_store_nta_partial(g23rhs[index],g23rhsL);
    vec_store_nta_partial(g33rhs[index],g33rhsL);
    vec_store_nta_partial(K11rhs[index],K11rhsL);
    vec_store_nta_partial(K12rhs[index],K12rhsL);
    vec_store_nta_partial(K13rhs[index],K13rhsL);
    vec_store_nta_partial(K22rhs[index],K22rhsL);
    vec_store_nta_partial(K23rhs[index],K23rhsL);
    vec_store_nta_partial(K33rhs[index],K33rhsL);
  }
  CCTK_ENDLOOP3STR(ML_ADM_RHS);
}
extern "C" void ML_ADM_RHS(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_ADM_RHS
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_ADM_RHS);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADM_RHS_Body");
  }
  if (cctk_iteration % ML_ADM_RHS_calc_every != ML_ADM_RHS_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ML_ADM::ML_curv",
    "ML_ADM::ML_curvrhs",
    "ML_ADM::ML_lapse",
    "ML_ADM::ML_lapserhs",
    "ML_ADM::ML_metric",
    "ML_ADM::ML_metricrhs",
    "ML_ADM::ML_shift",
    "ML_ADM::ML_shiftrhs"};
  AssertGroupStorage(cctkGH, "ML_ADM_RHS", 8, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "ML_ADM_RHS", 1, 1, 1);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "ML_ADM_RHS", 2, 2, 2);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "ML_ADM_RHS", 3, 3, 3);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "ML_ADM_RHS", 4, 4, 4);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, ML_ADM_RHS_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_ADM_RHS_Body");
  }
}

} // namespace ML_ADM
