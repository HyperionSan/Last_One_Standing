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


static void CL_BSSN_Minkowski_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3STR(CL_BSSN_Minkowski,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    
    
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
    CCTK_REAL_VEC phiL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC gt11L CCTK_ATTRIBUTE_UNUSED = ToReal(1);
    
    CCTK_REAL_VEC gt12L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC gt13L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC gt22L CCTK_ATTRIBUTE_UNUSED = ToReal(1);
    
    CCTK_REAL_VEC gt23L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC gt33L CCTK_ATTRIBUTE_UNUSED = ToReal(1);
    
    CCTK_REAL_VEC trKL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At11L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At12L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At13L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At22L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At23L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At33L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Xt1L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Xt2L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Xt3L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dphi1L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dphi2L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dphi3L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt111L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt211L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt311L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt112L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt212L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt312L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt113L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt213L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt313L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt122L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt222L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt322L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt123L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt223L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt323L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt133L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt233L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dgt333L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC alphaL CCTK_ATTRIBUTE_UNUSED = ToReal(1);
    
    CCTK_REAL_VEC dalpha1L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dalpha2L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dalpha3L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC beta1L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC beta2L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC beta3L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dbeta11L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dbeta12L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dbeta13L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dbeta21L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dbeta22L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dbeta23L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dbeta31L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dbeta32L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC dbeta33L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC B1L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC B2L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC B3L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(alpha[index],alphaL);
    vec_store_nta_partial(At11[index],At11L);
    vec_store_nta_partial(At12[index],At12L);
    vec_store_nta_partial(At13[index],At13L);
    vec_store_nta_partial(At22[index],At22L);
    vec_store_nta_partial(At23[index],At23L);
    vec_store_nta_partial(At33[index],At33L);
    vec_store_nta_partial(B1[index],B1L);
    vec_store_nta_partial(B2[index],B2L);
    vec_store_nta_partial(B3[index],B3L);
    vec_store_nta_partial(beta1[index],beta1L);
    vec_store_nta_partial(beta2[index],beta2L);
    vec_store_nta_partial(beta3[index],beta3L);
    vec_store_nta_partial(dalpha1[index],dalpha1L);
    vec_store_nta_partial(dalpha2[index],dalpha2L);
    vec_store_nta_partial(dalpha3[index],dalpha3L);
    vec_store_nta_partial(dbeta11[index],dbeta11L);
    vec_store_nta_partial(dbeta12[index],dbeta12L);
    vec_store_nta_partial(dbeta13[index],dbeta13L);
    vec_store_nta_partial(dbeta21[index],dbeta21L);
    vec_store_nta_partial(dbeta22[index],dbeta22L);
    vec_store_nta_partial(dbeta23[index],dbeta23L);
    vec_store_nta_partial(dbeta31[index],dbeta31L);
    vec_store_nta_partial(dbeta32[index],dbeta32L);
    vec_store_nta_partial(dbeta33[index],dbeta33L);
    vec_store_nta_partial(dgt111[index],dgt111L);
    vec_store_nta_partial(dgt112[index],dgt112L);
    vec_store_nta_partial(dgt113[index],dgt113L);
    vec_store_nta_partial(dgt122[index],dgt122L);
    vec_store_nta_partial(dgt123[index],dgt123L);
    vec_store_nta_partial(dgt133[index],dgt133L);
    vec_store_nta_partial(dgt211[index],dgt211L);
    vec_store_nta_partial(dgt212[index],dgt212L);
    vec_store_nta_partial(dgt213[index],dgt213L);
    vec_store_nta_partial(dgt222[index],dgt222L);
    vec_store_nta_partial(dgt223[index],dgt223L);
    vec_store_nta_partial(dgt233[index],dgt233L);
    vec_store_nta_partial(dgt311[index],dgt311L);
    vec_store_nta_partial(dgt312[index],dgt312L);
    vec_store_nta_partial(dgt313[index],dgt313L);
    vec_store_nta_partial(dgt322[index],dgt322L);
    vec_store_nta_partial(dgt323[index],dgt323L);
    vec_store_nta_partial(dgt333[index],dgt333L);
    vec_store_nta_partial(dphi1[index],dphi1L);
    vec_store_nta_partial(dphi2[index],dphi2L);
    vec_store_nta_partial(dphi3[index],dphi3L);
    vec_store_nta_partial(gt11[index],gt11L);
    vec_store_nta_partial(gt12[index],gt12L);
    vec_store_nta_partial(gt13[index],gt13L);
    vec_store_nta_partial(gt22[index],gt22L);
    vec_store_nta_partial(gt23[index],gt23L);
    vec_store_nta_partial(gt33[index],gt33L);
    vec_store_nta_partial(phi[index],phiL);
    vec_store_nta_partial(trK[index],trKL);
    vec_store_nta_partial(Xt1[index],Xt1L);
    vec_store_nta_partial(Xt2[index],Xt2L);
    vec_store_nta_partial(Xt3[index],Xt3L);
  }
  CCTK_ENDLOOP3STR(CL_BSSN_Minkowski);
}
extern "C" void CL_BSSN_Minkowski(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CL_BSSN_Minkowski
  DECLARE_CCTK_ARGUMENTS_CHECKED(CL_BSSN_Minkowski);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering CL_BSSN_Minkowski_Body");
  }
  if (cctk_iteration % CL_BSSN_Minkowski_calc_every != CL_BSSN_Minkowski_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "CL_BSSN::CL_curv",
    "CL_BSSN::CL_dlapse",
    "CL_BSSN::CL_dlog_confac",
    "CL_BSSN::CL_dmetric",
    "CL_BSSN::CL_dshift",
    "CL_BSSN::CL_dtshift",
    "CL_BSSN::CL_Gamma",
    "CL_BSSN::CL_lapse",
    "CL_BSSN::CL_log_confac",
    "CL_BSSN::CL_metric",
    "CL_BSSN::CL_shift",
    "CL_BSSN::CL_trace_curv"};
  AssertGroupStorage(cctkGH, "CL_BSSN_Minkowski", 12, groups);
  
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
  
  LoopOverEverything(cctkGH, CL_BSSN_Minkowski_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving CL_BSSN_Minkowski_Body");
  }
}

} // namespace CL_BSSN
