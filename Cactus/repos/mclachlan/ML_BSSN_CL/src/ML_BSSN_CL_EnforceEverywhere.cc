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
#include "OpenCLRunTime.h"
#include "vectors.h"

namespace ML_BSSN_CL {


static void ML_BSSN_CL_EnforceEverywhere_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  const char* const source =
  "/* Include user-supplied include files */\n"
  "/* Initialise finite differencing variables */\n"
  "const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;\n"
  "const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = \n"
  "  CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);\n"
  "const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = \n"
  "  CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);\n"
  "const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;\n"
  "const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;\n"
  "const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;\n"
  "const ptrdiff_t cctkLbnd1 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[0];\n"
  "const ptrdiff_t cctkLbnd2 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[1];\n"
  "const ptrdiff_t cctkLbnd3 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[2];\n"
  "const CCTK_REAL_VEC t CCTK_ATTRIBUTE_UNUSED = ToReal(cctk_time);\n"
  "const CCTK_REAL_VEC cctkOriginSpace1 CCTK_ATTRIBUTE_UNUSED = \n"
  "  ToReal(CCTK_ORIGIN_SPACE(0));\n"
  "const CCTK_REAL_VEC cctkOriginSpace2 CCTK_ATTRIBUTE_UNUSED = \n"
  "  ToReal(CCTK_ORIGIN_SPACE(1));\n"
  "const CCTK_REAL_VEC cctkOriginSpace3 CCTK_ATTRIBUTE_UNUSED = \n"
  "  ToReal(CCTK_ORIGIN_SPACE(2));\n"
  "const CCTK_REAL_VEC dt CCTK_ATTRIBUTE_UNUSED = \n"
  "  ToReal(CCTK_DELTA_TIME);\n"
  "const CCTK_REAL_VEC dx CCTK_ATTRIBUTE_UNUSED = \n"
  "  ToReal(CCTK_DELTA_SPACE(0));\n"
  "const CCTK_REAL_VEC dy CCTK_ATTRIBUTE_UNUSED = \n"
  "  ToReal(CCTK_DELTA_SPACE(1));\n"
  "const CCTK_REAL_VEC dz CCTK_ATTRIBUTE_UNUSED = \n"
  "  ToReal(CCTK_DELTA_SPACE(2));\n"
  "const CCTK_REAL_VEC dxi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dx);\n"
  "const CCTK_REAL_VEC dyi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dy);\n"
  "const CCTK_REAL_VEC dzi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dz);\n"
  "const CCTK_REAL_VEC khalf CCTK_ATTRIBUTE_UNUSED = ToReal(0.5);\n"
  "const CCTK_REAL_VEC kthird CCTK_ATTRIBUTE_UNUSED = \n"
  "  ToReal(0.333333333333333333333333333333);\n"
  "const CCTK_REAL_VEC ktwothird CCTK_ATTRIBUTE_UNUSED = \n"
  "  ToReal(0.666666666666666666666666666667);\n"
  "const CCTK_REAL_VEC kfourthird CCTK_ATTRIBUTE_UNUSED = \n"
  "  ToReal(1.33333333333333333333333333333);\n"
  "const CCTK_REAL_VEC hdxi CCTK_ATTRIBUTE_UNUSED = \n"
  "  kmul(dxi,ToReal(0.5));\n"
  "const CCTK_REAL_VEC hdyi CCTK_ATTRIBUTE_UNUSED = \n"
  "  kmul(dyi,ToReal(0.5));\n"
  "const CCTK_REAL_VEC hdzi CCTK_ATTRIBUTE_UNUSED = \n"
  "  kmul(dzi,ToReal(0.5));\n"
  "/* Initialize predefined quantities */\n"
  "const CCTK_REAL_VEC p1o1024dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0009765625),dx);\n"
  "const CCTK_REAL_VEC p1o1024dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0009765625),dy);\n"
  "const CCTK_REAL_VEC p1o1024dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0009765625),dz);\n"
  "const CCTK_REAL_VEC p1o120dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00833333333333333333333333333333),dx);\n"
  "const CCTK_REAL_VEC p1o120dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00833333333333333333333333333333),dy);\n"
  "const CCTK_REAL_VEC p1o120dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00833333333333333333333333333333),dz);\n"
  "const CCTK_REAL_VEC p1o12dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dx);\n"
  "const CCTK_REAL_VEC p1o12dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dy);\n"
  "const CCTK_REAL_VEC p1o12dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dz);\n"
  "const CCTK_REAL_VEC p1o144dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dx,dy));\n"
  "const CCTK_REAL_VEC p1o144dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dx,dz));\n"
  "const CCTK_REAL_VEC p1o144dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dy,dz));\n"
  "const CCTK_REAL_VEC p1o1680dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000595238095238095238095238095238),dx);\n"
  "const CCTK_REAL_VEC p1o1680dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000595238095238095238095238095238),dy);\n"
  "const CCTK_REAL_VEC p1o1680dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000595238095238095238095238095238),dz);\n"
  "const CCTK_REAL_VEC p1o180dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dx,dx));\n"
  "const CCTK_REAL_VEC p1o180dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dy,dy));\n"
  "const CCTK_REAL_VEC p1o180dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dz,dz));\n"
  "const CCTK_REAL_VEC p1o24dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dx);\n"
  "const CCTK_REAL_VEC p1o24dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dy);\n"
  "const CCTK_REAL_VEC p1o24dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dz);\n"
  "const CCTK_REAL_VEC p1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dx);\n"
  "const CCTK_REAL_VEC p1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dy);\n"
  "const CCTK_REAL_VEC p1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dz);\n"
  "const CCTK_REAL_VEC p1o3600dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dx,dy));\n"
  "const CCTK_REAL_VEC p1o3600dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dx,dz));\n"
  "const CCTK_REAL_VEC p1o3600dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dy,dz));\n"
  "const CCTK_REAL_VEC p1o4dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dx);\n"
  "const CCTK_REAL_VEC p1o4dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dx,dy));\n"
  "const CCTK_REAL_VEC p1o4dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dx,dz));\n"
  "const CCTK_REAL_VEC p1o4dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dy);\n"
  "const CCTK_REAL_VEC p1o4dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dy,dz));\n"
  "const CCTK_REAL_VEC p1o4dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dz);\n"
  "const CCTK_REAL_VEC p1o5040dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dx,dx));\n"
  "const CCTK_REAL_VEC p1o5040dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dy,dy));\n"
  "const CCTK_REAL_VEC p1o5040dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dz,dz));\n"
  "const CCTK_REAL_VEC p1o560dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00178571428571428571428571428571),dx);\n"
  "const CCTK_REAL_VEC p1o560dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00178571428571428571428571428571),dy);\n"
  "const CCTK_REAL_VEC p1o560dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00178571428571428571428571428571),dz);\n"
  "const CCTK_REAL_VEC p1o60dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0166666666666666666666666666667),dx);\n"
  "const CCTK_REAL_VEC p1o60dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0166666666666666666666666666667),dy);\n"
  "const CCTK_REAL_VEC p1o60dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0166666666666666666666666666667),dz);\n"
  "const CCTK_REAL_VEC p1o64dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.015625),dx);\n"
  "const CCTK_REAL_VEC p1o64dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.015625),dy);\n"
  "const CCTK_REAL_VEC p1o64dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.015625),dz);\n"
  "const CCTK_REAL_VEC p1o6dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.166666666666666666666666666667),dx);\n"
  "const CCTK_REAL_VEC p1o6dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.166666666666666666666666666667),dy);\n"
  "const CCTK_REAL_VEC p1o6dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.166666666666666666666666666667),dz);\n"
  "const CCTK_REAL_VEC p1o705600dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dx,dy));\n"
  "const CCTK_REAL_VEC p1o705600dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dx,dz));\n"
  "const CCTK_REAL_VEC p1o705600dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dy,dz));\n"
  "const CCTK_REAL_VEC p1o840dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dx);\n"
  "const CCTK_REAL_VEC p1o840dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dy);\n"
  "const CCTK_REAL_VEC p1o840dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dz);\n"
  "const CCTK_REAL_VEC p1odx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dx,dx));\n"
  "const CCTK_REAL_VEC p1ody2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dy,dy));\n"
  "const CCTK_REAL_VEC p1odz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dz,dz));\n"
  "const CCTK_REAL_VEC pm1o120dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00833333333333333333333333333333),dx);\n"
  "const CCTK_REAL_VEC pm1o120dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00833333333333333333333333333333),dy);\n"
  "const CCTK_REAL_VEC pm1o120dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00833333333333333333333333333333),dz);\n"
  "const CCTK_REAL_VEC pm1o12dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),dx);\n"
  "const CCTK_REAL_VEC pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dx,dx));\n"
  "const CCTK_REAL_VEC pm1o12dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),dy);\n"
  "const CCTK_REAL_VEC pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dy,dy));\n"
  "const CCTK_REAL_VEC pm1o12dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),dz);\n"
  "const CCTK_REAL_VEC pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dz,dz));\n"
  "const CCTK_REAL_VEC pm1o16dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0625),dx);\n"
  "const CCTK_REAL_VEC pm1o16dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0625),dy);\n"
  "const CCTK_REAL_VEC pm1o16dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0625),dz);\n"
  "const CCTK_REAL_VEC pm1o256dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00390625),dx);\n"
  "const CCTK_REAL_VEC pm1o256dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00390625),dy);\n"
  "const CCTK_REAL_VEC pm1o256dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00390625),dz);\n"
  "const CCTK_REAL_VEC pm1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dx);\n"
  "const CCTK_REAL_VEC pm1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dy);\n"
  "const CCTK_REAL_VEC pm1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dz);\n"
  "const CCTK_REAL_VEC pm1o4dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.25),dx);\n"
  "const CCTK_REAL_VEC pm1o4dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.25),dy);\n"
  "const CCTK_REAL_VEC pm1o4dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.25),dz);\n"
  "const CCTK_REAL_VEC pm1o60dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0166666666666666666666666666667),dx);\n"
  "const CCTK_REAL_VEC pm1o60dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0166666666666666666666666666667),dy);\n"
  "const CCTK_REAL_VEC pm1o60dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0166666666666666666666666666667),dz);\n"
  "/* Jacobian variable pointers */\n"
  "const bool usejacobian1 = (!CCTK_IsFunctionAliased(\"MultiPatch_GetMap\") || MultiPatch_GetMap(cctkGH) != jacobian_identity_map)\n"
  "                      && strlen(jacobian_group) > 0;\n"
  "const bool usejacobian = assume_use_jacobian>=0 ? assume_use_jacobian : usejacobian1;\n"
  "if (usejacobian && (strlen(jacobian_derivative_group) == 0))\n"
  "{\n"
  "  CCTK_WARN(CCTK_WARN_ALERT, \"GenericFD::jacobian_group and GenericFD::jacobian_derivative_group must both be set to valid group names\");\n"
  "}\n"
  "\n"
  "const CCTK_REAL* restrict jacobian_ptrs[9];\n"
  "if (usejacobian) GroupDataPointers(cctkGH, jacobian_group,\n"
  "                                              9, jacobian_ptrs);\n"
  "\n"
  "const CCTK_REAL* restrict const J11 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[0] : 0;\n"
  "const CCTK_REAL* restrict const J12 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[1] : 0;\n"
  "const CCTK_REAL* restrict const J13 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[2] : 0;\n"
  "const CCTK_REAL* restrict const J21 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[3] : 0;\n"
  "const CCTK_REAL* restrict const J22 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[4] : 0;\n"
  "const CCTK_REAL* restrict const J23 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[5] : 0;\n"
  "const CCTK_REAL* restrict const J31 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[6] : 0;\n"
  "const CCTK_REAL* restrict const J32 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[7] : 0;\n"
  "const CCTK_REAL* restrict const J33 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[8] : 0;\n"
  "\n"
  "const CCTK_REAL* restrict jacobian_determinant_ptrs[1] CCTK_ATTRIBUTE_UNUSED;\n"
  "if (usejacobian && strlen(jacobian_determinant_group) > 0) GroupDataPointers(cctkGH, jacobian_determinant_group,\n"
  "                                              1, jacobian_determinant_ptrs);\n"
  "\n"
  "const CCTK_REAL* restrict const detJ CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_determinant_ptrs[0] : 0;\n"
  "\n"
  "const CCTK_REAL* restrict jacobian_inverse_ptrs[9] CCTK_ATTRIBUTE_UNUSED;\n"
  "if (usejacobian && strlen(jacobian_inverse_group) > 0) GroupDataPointers(cctkGH, jacobian_inverse_group,\n"
  "                                              9, jacobian_inverse_ptrs);\n"
  "\n"
  "const CCTK_REAL* restrict const iJ11 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[0] : 0;\n"
  "const CCTK_REAL* restrict const iJ12 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[1] : 0;\n"
  "const CCTK_REAL* restrict const iJ13 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[2] : 0;\n"
  "const CCTK_REAL* restrict const iJ21 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[3] : 0;\n"
  "const CCTK_REAL* restrict const iJ22 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[4] : 0;\n"
  "const CCTK_REAL* restrict const iJ23 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[5] : 0;\n"
  "const CCTK_REAL* restrict const iJ31 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[6] : 0;\n"
  "const CCTK_REAL* restrict const iJ32 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[7] : 0;\n"
  "const CCTK_REAL* restrict const iJ33 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[8] : 0;\n"
  "\n"
  "const CCTK_REAL* restrict jacobian_derivative_ptrs[18] CCTK_ATTRIBUTE_UNUSED;\n"
  "if (usejacobian) GroupDataPointers(cctkGH, jacobian_derivative_group,\n"
  "                                    18, jacobian_derivative_ptrs);\n"
  "\n"
  "const CCTK_REAL* restrict const dJ111 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[0] : 0;\n"
  "const CCTK_REAL* restrict const dJ112 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[1] : 0;\n"
  "const CCTK_REAL* restrict const dJ113 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[2] : 0;\n"
  "const CCTK_REAL* restrict const dJ122 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[3] : 0;\n"
  "const CCTK_REAL* restrict const dJ123 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[4] : 0;\n"
  "const CCTK_REAL* restrict const dJ133 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[5] : 0;\n"
  "const CCTK_REAL* restrict const dJ211 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[6] : 0;\n"
  "const CCTK_REAL* restrict const dJ212 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[7] : 0;\n"
  "const CCTK_REAL* restrict const dJ213 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[8] : 0;\n"
  "const CCTK_REAL* restrict const dJ222 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[9] : 0;\n"
  "const CCTK_REAL* restrict const dJ223 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[10] : 0;\n"
  "const CCTK_REAL* restrict const dJ233 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[11] : 0;\n"
  "const CCTK_REAL* restrict const dJ311 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[12] : 0;\n"
  "const CCTK_REAL* restrict const dJ312 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[13] : 0;\n"
  "const CCTK_REAL* restrict const dJ313 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[14] : 0;\n"
  "const CCTK_REAL* restrict const dJ322 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[15] : 0;\n"
  "const CCTK_REAL* restrict const dJ323 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[16] : 0;\n"
  "const CCTK_REAL* restrict const dJ333 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[17] : 0;\n"
  "/* Assign local copies of arrays functions */\n"
  "\n"
  "\n"
  "/* Calculate temporaries and arrays functions */\n"
  "/* Copy local copies back to grid functions */\n"
  "/* Loop over the grid points */\n"
  "const int imin0=imin[0];\n"
  "const int imin1=imin[1];\n"
  "const int imin2=imin[2];\n"
  "const int imax0=imax[0];\n"
  "const int imax1=imax[1];\n"
  "const int imax2=imax[2];\n"
  "#pragma omp parallel\n"
  "CCTK_LOOP3STR(ML_BSSN_CL_EnforceEverywhere,\n"
  "  i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,\n"
  "  cctk_ash[0],cctk_ash[1],cctk_ash[2],\n"
  "  vecimin,vecimax, CCTK_REAL_VEC_SIZE)\n"
  "{\n"
  "  const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;\n"
  "  /* Assign local copies of grid functions */\n"
  "  \n"
  "  CCTK_REAL_VEC alphaL CCTK_ATTRIBUTE_UNUSED = vec_load(alpha[index]);\n"
  "  CCTK_REAL_VEC At11L CCTK_ATTRIBUTE_UNUSED = vec_load(At11[index]);\n"
  "  CCTK_REAL_VEC At12L CCTK_ATTRIBUTE_UNUSED = vec_load(At12[index]);\n"
  "  CCTK_REAL_VEC At13L CCTK_ATTRIBUTE_UNUSED = vec_load(At13[index]);\n"
  "  CCTK_REAL_VEC At22L CCTK_ATTRIBUTE_UNUSED = vec_load(At22[index]);\n"
  "  CCTK_REAL_VEC At23L CCTK_ATTRIBUTE_UNUSED = vec_load(At23[index]);\n"
  "  CCTK_REAL_VEC At33L CCTK_ATTRIBUTE_UNUSED = vec_load(At33[index]);\n"
  "  CCTK_REAL_VEC gt11L CCTK_ATTRIBUTE_UNUSED = vec_load(gt11[index]);\n"
  "  CCTK_REAL_VEC gt12L CCTK_ATTRIBUTE_UNUSED = vec_load(gt12[index]);\n"
  "  CCTK_REAL_VEC gt13L CCTK_ATTRIBUTE_UNUSED = vec_load(gt13[index]);\n"
  "  CCTK_REAL_VEC gt22L CCTK_ATTRIBUTE_UNUSED = vec_load(gt22[index]);\n"
  "  CCTK_REAL_VEC gt23L CCTK_ATTRIBUTE_UNUSED = vec_load(gt23[index]);\n"
  "  CCTK_REAL_VEC gt33L CCTK_ATTRIBUTE_UNUSED = vec_load(gt33[index]);\n"
  "  \n"
  "  \n"
  "  /* Include user supplied include files */\n"
  "  /* Precompute derivatives */\n"
  "  \n"
  "  switch (fdOrder)\n"
  "  {\n"
  "    case 2:\n"
  "    {\n"
  "      break;\n"
  "    }\n"
  "    \n"
  "    case 4:\n"
  "    {\n"
  "      break;\n"
  "    }\n"
  "    \n"
  "    case 6:\n"
  "    {\n"
  "      break;\n"
  "    }\n"
  "    \n"
  "    case 8:\n"
  "    {\n"
  "      break;\n"
  "    }\n"
  "    default:\n"
  "      CCTK_BUILTIN_UNREACHABLE();\n"
  "  }\n"
  "  /* Calculate temporaries and grid functions */\n"
  "  CCTK_REAL_VEC detgt CCTK_ATTRIBUTE_UNUSED = ToReal(1);\n"
  "  \n"
  "  gt11L = kmul(gt11L,kpow(detgt,-0.333333333333333333333333333333));\n"
  "  \n"
  "  gt12L = kmul(gt12L,kpow(detgt,-0.333333333333333333333333333333));\n"
  "  \n"
  "  gt13L = kmul(gt13L,kpow(detgt,-0.333333333333333333333333333333));\n"
  "  \n"
  "  gt22L = kmul(gt22L,kpow(detgt,-0.333333333333333333333333333333));\n"
  "  \n"
  "  gt23L = kmul(gt23L,kpow(detgt,-0.333333333333333333333333333333));\n"
  "  \n"
  "  gt33L = kmul(gt33L,kpow(detgt,-0.333333333333333333333333333333));\n"
  "  \n"
  "  detgt = ToReal(1);\n"
  "  \n"
  "  CCTK_REAL_VEC gtu11 CCTK_ATTRIBUTE_UNUSED = \n"
  "    kdiv(kmsub(gt22L,gt33L,kmul(gt23L,gt23L)),detgt);\n"
  "  \n"
  "  CCTK_REAL_VEC gtu12 CCTK_ATTRIBUTE_UNUSED = \n"
  "    kdiv(kmsub(gt13L,gt23L,kmul(gt12L,gt33L)),detgt);\n"
  "  \n"
  "  CCTK_REAL_VEC gtu13 CCTK_ATTRIBUTE_UNUSED = \n"
  "    kdiv(kmsub(gt12L,gt23L,kmul(gt13L,gt22L)),detgt);\n"
  "  \n"
  "  CCTK_REAL_VEC gtu22 CCTK_ATTRIBUTE_UNUSED = \n"
  "    kdiv(kmsub(gt11L,gt33L,kmul(gt13L,gt13L)),detgt);\n"
  "  \n"
  "  CCTK_REAL_VEC gtu23 CCTK_ATTRIBUTE_UNUSED = \n"
  "    kdiv(kmsub(gt12L,gt13L,kmul(gt11L,gt23L)),detgt);\n"
  "  \n"
  "  CCTK_REAL_VEC gtu33 CCTK_ATTRIBUTE_UNUSED = \n"
  "    kdiv(kmsub(gt11L,gt22L,kmul(gt12L,gt12L)),detgt);\n"
  "  \n"
  "  CCTK_REAL_VEC Atm11 CCTK_ATTRIBUTE_UNUSED = \n"
  "    kmadd(At11L,gtu11,kmadd(At12L,gtu12,kmul(At13L,gtu13)));\n"
  "  \n"
  "  CCTK_REAL_VEC Atm22 CCTK_ATTRIBUTE_UNUSED = \n"
  "    kmadd(At12L,gtu12,kmadd(At22L,gtu22,kmul(At23L,gtu23)));\n"
  "  \n"
  "  CCTK_REAL_VEC Atm33 CCTK_ATTRIBUTE_UNUSED = \n"
  "    kmadd(At13L,gtu13,kmadd(At23L,gtu23,kmul(At33L,gtu33)));\n"
  "  \n"
  "  CCTK_REAL_VEC trAt CCTK_ATTRIBUTE_UNUSED = \n"
  "    kadd(Atm11,kadd(Atm22,Atm33));\n"
  "  \n"
  "  At11L = \n"
  "    kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt11L,trAt),At11L);\n"
  "  \n"
  "  At12L = \n"
  "    kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt12L,trAt),At12L);\n"
  "  \n"
  "  At13L = \n"
  "    kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt13L,trAt),At13L);\n"
  "  \n"
  "  At22L = \n"
  "    kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt22L,trAt),At22L);\n"
  "  \n"
  "  At23L = \n"
  "    kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt23L,trAt),At23L);\n"
  "  \n"
  "  At33L = \n"
  "    kmadd(ToReal(-0.333333333333333333333333333333),kmul(gt33L,trAt),At33L);\n"
  "  \n"
  "  alphaL = kfmax(alphaL,ToReal(minimumLapse));\n"
  "  /* Copy local copies back to grid functions */\n"
  "  vec_store_partial_prepare(i,lc_imin,lc_imax);\n"
  "  vec_store_nta_partial(alpha[index],alphaL);\n"
  "  vec_store_nta_partial(At11[index],At11L);\n"
  "  vec_store_nta_partial(At12[index],At12L);\n"
  "  vec_store_nta_partial(At13[index],At13L);\n"
  "  vec_store_nta_partial(At22[index],At22L);\n"
  "  vec_store_nta_partial(At23[index],At23L);\n"
  "  vec_store_nta_partial(At33[index],At33L);\n"
  "  vec_store_nta_partial(gt11[index],gt11L);\n"
  "  vec_store_nta_partial(gt12[index],gt12L);\n"
  "  vec_store_nta_partial(gt13[index],gt13L);\n"
  "  vec_store_nta_partial(gt22[index],gt22L);\n"
  "  vec_store_nta_partial(gt23[index],gt23L);\n"
  "  vec_store_nta_partial(gt33[index],gt33L);\n"
  "}\n"
  "CCTK_ENDLOOP3STR(ML_BSSN_CL_EnforceEverywhere);\n"
  ""
  ;
  
  const char* const groups[] = {
    "ML_BSSN_CL::ML_curv",
    "ML_BSSN_CL::ML_lapse",
    "ML_BSSN_CL::ML_metric",
    NULL};
  
  static struct OpenCLKernel *kernel = NULL;
  const char* const sources[] = {differencing, source, NULL};
  OpenCLRunTime_CallKernel(cctkGH, CCTK_THORNSTRING, "ML_BSSN_CL_EnforceEverywhere",
                           sources, groups, NULL, NULL, NULL, -1,
                           imin, imax, &kernel);
  
}
extern "C" void ML_BSSN_CL_EnforceEverywhere(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_BSSN_CL_EnforceEverywhere
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_CL_EnforceEverywhere);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_CL_EnforceEverywhere_Body");
  }
  if (cctk_iteration % ML_BSSN_CL_EnforceEverywhere_calc_every != ML_BSSN_CL_EnforceEverywhere_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ML_BSSN_CL::ML_curv",
    "ML_BSSN_CL::ML_lapse",
    "ML_BSSN_CL::ML_metric"};
  AssertGroupStorage(cctkGH, "ML_BSSN_CL_EnforceEverywhere", 3, groups);
  
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
  
  LoopOverEverything(cctkGH, ML_BSSN_CL_EnforceEverywhere_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_CL_EnforceEverywhere_Body");
  }
}

} // namespace ML_BSSN_CL
