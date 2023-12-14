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

namespace ML_WaveToy_CL {


static void WT_CL_Standing_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  "const CCTK_REAL_VEC p1o12dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dx);\n"
  "const CCTK_REAL_VEC p1o12dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dy);\n"
  "const CCTK_REAL_VEC p1o12dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dz);\n"
  "const CCTK_REAL_VEC p1o144dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dx,dy));\n"
  "const CCTK_REAL_VEC p1o144dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dx,dz));\n"
  "const CCTK_REAL_VEC p1o144dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dy,dz));\n"
  "const CCTK_REAL_VEC pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dx,dx));\n"
  "const CCTK_REAL_VEC pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dy,dy));\n"
  "const CCTK_REAL_VEC pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dz,dz));\n"
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
  "CCTK_LOOP3STR(WT_CL_Standing,\n"
  "  i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,\n"
  "  cctk_ash[0],cctk_ash[1],cctk_ash[2],\n"
  "  vecimin,vecimax, CCTK_REAL_VEC_SIZE)\n"
  "{\n"
  "  const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;\n"
  "  /* Assign local copies of grid functions */\n"
  "  \n"
  "  CCTK_REAL_VEC xL CCTK_ATTRIBUTE_UNUSED = vec_load(x[index]);\n"
  "  CCTK_REAL_VEC yL CCTK_ATTRIBUTE_UNUSED = vec_load(y[index]);\n"
  "  CCTK_REAL_VEC zL CCTK_ATTRIBUTE_UNUSED = vec_load(z[index]);\n"
  "  \n"
  "  /* Include user supplied include files */\n"
  "  /* Precompute derivatives */\n"
  "  /* Calculate temporaries and grid functions */\n"
  "  CCTK_REAL_VEC kvec CCTK_ATTRIBUTE_UNUSED = \n"
  "    ToReal(3.14159265358979323846264338328*pow(width,-1));\n"
  "  \n"
  "  CCTK_REAL_VEC omega CCTK_ATTRIBUTE_UNUSED = \n"
  "    kmul(kpow(kmul(kvec,kvec),0.5),ToReal(1.732050807568877293527446341506));\n"
  "  \n"
  "  CCTK_REAL_VEC uL CCTK_ATTRIBUTE_UNUSED = \n"
  "    kmul(kmul(kcos(kmul(xL,kvec)),kmul(kcos(kmul(yL,kvec)),kmul(kcos(kmul(zL,kvec)),kcos(kmul(omega,t))))),ToReal(amplitude));\n"
  "  \n"
  "  CCTK_REAL_VEC rhoL CCTK_ATTRIBUTE_UNUSED = \n"
  "    kneg(kmul(kmul(omega,kmul(kcos(kmul(xL,kvec)),kmul(kcos(kmul(yL,kvec)),kmul(kcos(kmul(zL,kvec)),ksin(kmul(omega,t)))))),ToReal(amplitude)));\n"
  "  /* Copy local copies back to grid functions */\n"
  "  vec_store_partial_prepare(i,lc_imin,lc_imax);\n"
  "  vec_store_nta_partial(rho[index],rhoL);\n"
  "  vec_store_nta_partial(u[index],uL);\n"
  "}\n"
  "CCTK_ENDLOOP3STR(WT_CL_Standing);\n"
  ""
  ;
  
  const char* const groups[] = {
    "grid::coordinates",
    "ML_WaveToy_CL::WT_rho",
    "ML_WaveToy_CL::WT_u",
    NULL};
  
  static struct OpenCLKernel *kernel = NULL;
  const char* const sources[] = {differencing, source, NULL};
  OpenCLRunTime_CallKernel(cctkGH, CCTK_THORNSTRING, "WT_CL_Standing",
                           sources, groups, NULL, NULL, NULL, -1,
                           imin, imax, &kernel);
  
}
extern "C" void WT_CL_Standing(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_WT_CL_Standing
  DECLARE_CCTK_ARGUMENTS_CHECKED(WT_CL_Standing);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering WT_CL_Standing_Body");
  }
  if (cctk_iteration % WT_CL_Standing_calc_every != WT_CL_Standing_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "grid::coordinates",
    "ML_WaveToy_CL::WT_rho",
    "ML_WaveToy_CL::WT_u"};
  AssertGroupStorage(cctkGH, "WT_CL_Standing", 3, groups);
  
  
  LoopOverEverything(cctkGH, WT_CL_Standing_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving WT_CL_Standing_Body");
  }
}

} // namespace ML_WaveToy_CL
