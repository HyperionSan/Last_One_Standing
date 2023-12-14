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

namespace ML_hydro {

extern "C" void hydro_RHS_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_hydro_RHS_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(hydro_RHS_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % hydro_RHS_calc_every != hydro_RHS_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_hydro::ene_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_hydro::ene_grouprhs.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_hydro::mass_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_hydro::mass_grouprhs.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_hydro::mom_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_hydro::mom_grouprhs.");
  return;
}

static void hydro_RHS_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  const CCTK_REAL p1o2dx CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dx,-1);
  const CCTK_REAL p1o2dy CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dy,-1);
  const CCTK_REAL p1o2dz CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dz,-1);
  const CCTK_REAL p1o4dxdy CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o4dxdz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dydz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1odx2 CCTK_ATTRIBUTE_UNUSED = pow(dx,-2);
  const CCTK_REAL p1ody2 CCTK_ATTRIBUTE_UNUSED = pow(dy,-2);
  const CCTK_REAL p1odz2 CCTK_ATTRIBUTE_UNUSED = pow(dz,-2);
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
  CCTK_LOOP3(hydro_RHS,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL eneflux1L CCTK_ATTRIBUTE_UNUSED = eneflux1[index];
    CCTK_REAL eneflux2L CCTK_ATTRIBUTE_UNUSED = eneflux2[index];
    CCTK_REAL eneflux3L CCTK_ATTRIBUTE_UNUSED = eneflux3[index];
    CCTK_REAL massflux1L CCTK_ATTRIBUTE_UNUSED = massflux1[index];
    CCTK_REAL massflux2L CCTK_ATTRIBUTE_UNUSED = massflux2[index];
    CCTK_REAL massflux3L CCTK_ATTRIBUTE_UNUSED = massflux3[index];
    CCTK_REAL momflux11L CCTK_ATTRIBUTE_UNUSED = momflux11[index];
    CCTK_REAL momflux12L CCTK_ATTRIBUTE_UNUSED = momflux12[index];
    CCTK_REAL momflux13L CCTK_ATTRIBUTE_UNUSED = momflux13[index];
    CCTK_REAL momflux21L CCTK_ATTRIBUTE_UNUSED = momflux21[index];
    CCTK_REAL momflux22L CCTK_ATTRIBUTE_UNUSED = momflux22[index];
    CCTK_REAL momflux23L CCTK_ATTRIBUTE_UNUSED = momflux23[index];
    CCTK_REAL momflux31L CCTK_ATTRIBUTE_UNUSED = momflux31[index];
    CCTK_REAL momflux32L CCTK_ATTRIBUTE_UNUSED = momflux32[index];
    CCTK_REAL momflux33L CCTK_ATTRIBUTE_UNUSED = momflux33[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    /* Calculate temporaries and grid functions */
    CCTK_REAL massrhsL CCTK_ATTRIBUTE_UNUSED = -PD(massflux1L,1) - 
      PD(massflux2L,2) - PD(massflux3L,3);
    
    CCTK_REAL mom1rhsL CCTK_ATTRIBUTE_UNUSED = -PD(momflux11L,1) - 
      PD(momflux12L,2) - PD(momflux13L,3);
    
    CCTK_REAL mom2rhsL CCTK_ATTRIBUTE_UNUSED = -PD(momflux21L,1) - 
      PD(momflux22L,2) - PD(momflux23L,3);
    
    CCTK_REAL mom3rhsL CCTK_ATTRIBUTE_UNUSED = -PD(momflux31L,1) - 
      PD(momflux32L,2) - PD(momflux33L,3);
    
    CCTK_REAL enerhsL CCTK_ATTRIBUTE_UNUSED = -PD(eneflux1L,1) - 
      PD(eneflux2L,2) - PD(eneflux3L,3);
    /* Copy local copies back to grid functions */
    enerhs[index] = enerhsL;
    massrhs[index] = massrhsL;
    mom1rhs[index] = mom1rhsL;
    mom2rhs[index] = mom2rhsL;
    mom3rhs[index] = mom3rhsL;
  }
  CCTK_ENDLOOP3(hydro_RHS);
}
extern "C" void hydro_RHS(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_hydro_RHS
  DECLARE_CCTK_ARGUMENTS_CHECKED(hydro_RHS);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering hydro_RHS_Body");
  }
  if (cctk_iteration % hydro_RHS_calc_every != hydro_RHS_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ML_hydro::eneflux_group",
    "ML_hydro::ene_grouprhs",
    "ML_hydro::massflux_group",
    "ML_hydro::mass_grouprhs",
    "ML_hydro::momflux_group",
    "ML_hydro::mom_grouprhs"};
  AssertGroupStorage(cctkGH, "hydro_RHS", 6, groups);
  
  
  LoopOverInterior(cctkGH, hydro_RHS_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving hydro_RHS_Body");
  }
}

} // namespace ML_hydro
