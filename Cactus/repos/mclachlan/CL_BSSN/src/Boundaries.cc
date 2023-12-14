/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Faces.h"
#include "util_Table.h"
#include "Symmetry.h"
#include "Kranc.hh"


/* the boundary treatment is split into 3 steps:    */
/* 1. excision                                      */
/* 2. symmetries                                    */
/* 3. "other" boundary conditions, e.g. radiative */

/* to simplify scheduling and testing, the 3 steps  */
/* are currently applied in separate functions      */


extern "C" void CL_BSSN_CheckBoundaries(CCTK_ARGUMENTS)
{
#ifdef DECLARE_CCTK_ARGUMENTS_CL_BSSN_CheckBoundaries
  DECLARE_CCTK_ARGUMENTS_CHECKED(CL_BSSN_CheckBoundaries);
#else
  DECLARE_CCTK_ARGUMENTS;
#endif
  DECLARE_CCTK_PARAMETERS;
  
  return;
}

extern "C" void CL_BSSN_SelectBoundConds(CCTK_ARGUMENTS)
{
#ifdef DECLARE_CCTK_ARGUMENTS_CL_BSSN_SelectBoundConds
  DECLARE_CCTK_ARGUMENTS_CHECKED(CL_BSSN_SelectBoundConds);
#else
  DECLARE_CCTK_ARGUMENTS;
#endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  
  if (CCTK_EQUALS(CL_log_confac_bound, "none"  ) ||
      CCTK_EQUALS(CL_log_confac_bound, "static") ||
      CCTK_EQUALS(CL_log_confac_bound, "flat"  ) ||
      CCTK_EQUALS(CL_log_confac_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::CL_log_confac", CL_log_confac_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register CL_log_confac_bound BC for CL_BSSN::CL_log_confac!");
  }
  
  if (CCTK_EQUALS(CL_dlog_confac_bound, "none"  ) ||
      CCTK_EQUALS(CL_dlog_confac_bound, "static") ||
      CCTK_EQUALS(CL_dlog_confac_bound, "flat"  ) ||
      CCTK_EQUALS(CL_dlog_confac_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::CL_dlog_confac", CL_dlog_confac_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register CL_dlog_confac_bound BC for CL_BSSN::CL_dlog_confac!");
  }
  
  if (CCTK_EQUALS(CL_metric_bound, "none"  ) ||
      CCTK_EQUALS(CL_metric_bound, "static") ||
      CCTK_EQUALS(CL_metric_bound, "flat"  ) ||
      CCTK_EQUALS(CL_metric_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::CL_metric", CL_metric_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register CL_metric_bound BC for CL_BSSN::CL_metric!");
  }
  
  if (CCTK_EQUALS(CL_dmetric_bound, "none"  ) ||
      CCTK_EQUALS(CL_dmetric_bound, "static") ||
      CCTK_EQUALS(CL_dmetric_bound, "flat"  ) ||
      CCTK_EQUALS(CL_dmetric_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::CL_dmetric", CL_dmetric_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register CL_dmetric_bound BC for CL_BSSN::CL_dmetric!");
  }
  
  if (CCTK_EQUALS(CL_Gamma_bound, "none"  ) ||
      CCTK_EQUALS(CL_Gamma_bound, "static") ||
      CCTK_EQUALS(CL_Gamma_bound, "flat"  ) ||
      CCTK_EQUALS(CL_Gamma_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::CL_Gamma", CL_Gamma_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register CL_Gamma_bound BC for CL_BSSN::CL_Gamma!");
  }
  
  if (CCTK_EQUALS(CL_trace_curv_bound, "none"  ) ||
      CCTK_EQUALS(CL_trace_curv_bound, "static") ||
      CCTK_EQUALS(CL_trace_curv_bound, "flat"  ) ||
      CCTK_EQUALS(CL_trace_curv_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::CL_trace_curv", CL_trace_curv_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register CL_trace_curv_bound BC for CL_BSSN::CL_trace_curv!");
  }
  
  if (CCTK_EQUALS(CL_curv_bound, "none"  ) ||
      CCTK_EQUALS(CL_curv_bound, "static") ||
      CCTK_EQUALS(CL_curv_bound, "flat"  ) ||
      CCTK_EQUALS(CL_curv_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::CL_curv", CL_curv_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register CL_curv_bound BC for CL_BSSN::CL_curv!");
  }
  
  if (CCTK_EQUALS(CL_lapse_bound, "none"  ) ||
      CCTK_EQUALS(CL_lapse_bound, "static") ||
      CCTK_EQUALS(CL_lapse_bound, "flat"  ) ||
      CCTK_EQUALS(CL_lapse_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::CL_lapse", CL_lapse_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register CL_lapse_bound BC for CL_BSSN::CL_lapse!");
  }
  
  if (CCTK_EQUALS(CL_dlapse_bound, "none"  ) ||
      CCTK_EQUALS(CL_dlapse_bound, "static") ||
      CCTK_EQUALS(CL_dlapse_bound, "flat"  ) ||
      CCTK_EQUALS(CL_dlapse_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::CL_dlapse", CL_dlapse_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register CL_dlapse_bound BC for CL_BSSN::CL_dlapse!");
  }
  
  if (CCTK_EQUALS(CL_shift_bound, "none"  ) ||
      CCTK_EQUALS(CL_shift_bound, "static") ||
      CCTK_EQUALS(CL_shift_bound, "flat"  ) ||
      CCTK_EQUALS(CL_shift_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::CL_shift", CL_shift_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register CL_shift_bound BC for CL_BSSN::CL_shift!");
  }
  
  if (CCTK_EQUALS(CL_dshift_bound, "none"  ) ||
      CCTK_EQUALS(CL_dshift_bound, "static") ||
      CCTK_EQUALS(CL_dshift_bound, "flat"  ) ||
      CCTK_EQUALS(CL_dshift_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::CL_dshift", CL_dshift_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register CL_dshift_bound BC for CL_BSSN::CL_dshift!");
  }
  
  if (CCTK_EQUALS(CL_dtshift_bound, "none"  ) ||
      CCTK_EQUALS(CL_dtshift_bound, "static") ||
      CCTK_EQUALS(CL_dtshift_bound, "flat"  ) ||
      CCTK_EQUALS(CL_dtshift_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::CL_dtshift", CL_dtshift_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register CL_dtshift_bound BC for CL_BSSN::CL_dtshift!");
  }
  
  if (CCTK_EQUALS(phi_bound, "none"  ) ||
      CCTK_EQUALS(phi_bound, "static") ||
      CCTK_EQUALS(phi_bound, "flat"  ) ||
      CCTK_EQUALS(phi_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::phi", phi_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register phi_bound BC for CL_BSSN::phi!");
  }
  
  if (CCTK_EQUALS(dphi1_bound, "none"  ) ||
      CCTK_EQUALS(dphi1_bound, "static") ||
      CCTK_EQUALS(dphi1_bound, "flat"  ) ||
      CCTK_EQUALS(dphi1_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dphi1", dphi1_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dphi1_bound BC for CL_BSSN::dphi1!");
  }
  
  if (CCTK_EQUALS(dphi2_bound, "none"  ) ||
      CCTK_EQUALS(dphi2_bound, "static") ||
      CCTK_EQUALS(dphi2_bound, "flat"  ) ||
      CCTK_EQUALS(dphi2_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dphi2", dphi2_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dphi2_bound BC for CL_BSSN::dphi2!");
  }
  
  if (CCTK_EQUALS(dphi3_bound, "none"  ) ||
      CCTK_EQUALS(dphi3_bound, "static") ||
      CCTK_EQUALS(dphi3_bound, "flat"  ) ||
      CCTK_EQUALS(dphi3_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dphi3", dphi3_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dphi3_bound BC for CL_BSSN::dphi3!");
  }
  
  if (CCTK_EQUALS(gt11_bound, "none"  ) ||
      CCTK_EQUALS(gt11_bound, "static") ||
      CCTK_EQUALS(gt11_bound, "flat"  ) ||
      CCTK_EQUALS(gt11_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::gt11", gt11_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register gt11_bound BC for CL_BSSN::gt11!");
  }
  
  if (CCTK_EQUALS(gt12_bound, "none"  ) ||
      CCTK_EQUALS(gt12_bound, "static") ||
      CCTK_EQUALS(gt12_bound, "flat"  ) ||
      CCTK_EQUALS(gt12_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::gt12", gt12_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register gt12_bound BC for CL_BSSN::gt12!");
  }
  
  if (CCTK_EQUALS(gt13_bound, "none"  ) ||
      CCTK_EQUALS(gt13_bound, "static") ||
      CCTK_EQUALS(gt13_bound, "flat"  ) ||
      CCTK_EQUALS(gt13_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::gt13", gt13_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register gt13_bound BC for CL_BSSN::gt13!");
  }
  
  if (CCTK_EQUALS(gt22_bound, "none"  ) ||
      CCTK_EQUALS(gt22_bound, "static") ||
      CCTK_EQUALS(gt22_bound, "flat"  ) ||
      CCTK_EQUALS(gt22_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::gt22", gt22_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register gt22_bound BC for CL_BSSN::gt22!");
  }
  
  if (CCTK_EQUALS(gt23_bound, "none"  ) ||
      CCTK_EQUALS(gt23_bound, "static") ||
      CCTK_EQUALS(gt23_bound, "flat"  ) ||
      CCTK_EQUALS(gt23_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::gt23", gt23_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register gt23_bound BC for CL_BSSN::gt23!");
  }
  
  if (CCTK_EQUALS(gt33_bound, "none"  ) ||
      CCTK_EQUALS(gt33_bound, "static") ||
      CCTK_EQUALS(gt33_bound, "flat"  ) ||
      CCTK_EQUALS(gt33_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::gt33", gt33_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register gt33_bound BC for CL_BSSN::gt33!");
  }
  
  if (CCTK_EQUALS(dgt111_bound, "none"  ) ||
      CCTK_EQUALS(dgt111_bound, "static") ||
      CCTK_EQUALS(dgt111_bound, "flat"  ) ||
      CCTK_EQUALS(dgt111_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt111", dgt111_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt111_bound BC for CL_BSSN::dgt111!");
  }
  
  if (CCTK_EQUALS(dgt112_bound, "none"  ) ||
      CCTK_EQUALS(dgt112_bound, "static") ||
      CCTK_EQUALS(dgt112_bound, "flat"  ) ||
      CCTK_EQUALS(dgt112_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt112", dgt112_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt112_bound BC for CL_BSSN::dgt112!");
  }
  
  if (CCTK_EQUALS(dgt113_bound, "none"  ) ||
      CCTK_EQUALS(dgt113_bound, "static") ||
      CCTK_EQUALS(dgt113_bound, "flat"  ) ||
      CCTK_EQUALS(dgt113_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt113", dgt113_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt113_bound BC for CL_BSSN::dgt113!");
  }
  
  if (CCTK_EQUALS(dgt122_bound, "none"  ) ||
      CCTK_EQUALS(dgt122_bound, "static") ||
      CCTK_EQUALS(dgt122_bound, "flat"  ) ||
      CCTK_EQUALS(dgt122_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt122", dgt122_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt122_bound BC for CL_BSSN::dgt122!");
  }
  
  if (CCTK_EQUALS(dgt123_bound, "none"  ) ||
      CCTK_EQUALS(dgt123_bound, "static") ||
      CCTK_EQUALS(dgt123_bound, "flat"  ) ||
      CCTK_EQUALS(dgt123_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt123", dgt123_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt123_bound BC for CL_BSSN::dgt123!");
  }
  
  if (CCTK_EQUALS(dgt133_bound, "none"  ) ||
      CCTK_EQUALS(dgt133_bound, "static") ||
      CCTK_EQUALS(dgt133_bound, "flat"  ) ||
      CCTK_EQUALS(dgt133_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt133", dgt133_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt133_bound BC for CL_BSSN::dgt133!");
  }
  
  if (CCTK_EQUALS(dgt211_bound, "none"  ) ||
      CCTK_EQUALS(dgt211_bound, "static") ||
      CCTK_EQUALS(dgt211_bound, "flat"  ) ||
      CCTK_EQUALS(dgt211_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt211", dgt211_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt211_bound BC for CL_BSSN::dgt211!");
  }
  
  if (CCTK_EQUALS(dgt212_bound, "none"  ) ||
      CCTK_EQUALS(dgt212_bound, "static") ||
      CCTK_EQUALS(dgt212_bound, "flat"  ) ||
      CCTK_EQUALS(dgt212_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt212", dgt212_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt212_bound BC for CL_BSSN::dgt212!");
  }
  
  if (CCTK_EQUALS(dgt213_bound, "none"  ) ||
      CCTK_EQUALS(dgt213_bound, "static") ||
      CCTK_EQUALS(dgt213_bound, "flat"  ) ||
      CCTK_EQUALS(dgt213_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt213", dgt213_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt213_bound BC for CL_BSSN::dgt213!");
  }
  
  if (CCTK_EQUALS(dgt222_bound, "none"  ) ||
      CCTK_EQUALS(dgt222_bound, "static") ||
      CCTK_EQUALS(dgt222_bound, "flat"  ) ||
      CCTK_EQUALS(dgt222_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt222", dgt222_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt222_bound BC for CL_BSSN::dgt222!");
  }
  
  if (CCTK_EQUALS(dgt223_bound, "none"  ) ||
      CCTK_EQUALS(dgt223_bound, "static") ||
      CCTK_EQUALS(dgt223_bound, "flat"  ) ||
      CCTK_EQUALS(dgt223_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt223", dgt223_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt223_bound BC for CL_BSSN::dgt223!");
  }
  
  if (CCTK_EQUALS(dgt233_bound, "none"  ) ||
      CCTK_EQUALS(dgt233_bound, "static") ||
      CCTK_EQUALS(dgt233_bound, "flat"  ) ||
      CCTK_EQUALS(dgt233_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt233", dgt233_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt233_bound BC for CL_BSSN::dgt233!");
  }
  
  if (CCTK_EQUALS(dgt311_bound, "none"  ) ||
      CCTK_EQUALS(dgt311_bound, "static") ||
      CCTK_EQUALS(dgt311_bound, "flat"  ) ||
      CCTK_EQUALS(dgt311_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt311", dgt311_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt311_bound BC for CL_BSSN::dgt311!");
  }
  
  if (CCTK_EQUALS(dgt312_bound, "none"  ) ||
      CCTK_EQUALS(dgt312_bound, "static") ||
      CCTK_EQUALS(dgt312_bound, "flat"  ) ||
      CCTK_EQUALS(dgt312_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt312", dgt312_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt312_bound BC for CL_BSSN::dgt312!");
  }
  
  if (CCTK_EQUALS(dgt313_bound, "none"  ) ||
      CCTK_EQUALS(dgt313_bound, "static") ||
      CCTK_EQUALS(dgt313_bound, "flat"  ) ||
      CCTK_EQUALS(dgt313_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt313", dgt313_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt313_bound BC for CL_BSSN::dgt313!");
  }
  
  if (CCTK_EQUALS(dgt322_bound, "none"  ) ||
      CCTK_EQUALS(dgt322_bound, "static") ||
      CCTK_EQUALS(dgt322_bound, "flat"  ) ||
      CCTK_EQUALS(dgt322_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt322", dgt322_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt322_bound BC for CL_BSSN::dgt322!");
  }
  
  if (CCTK_EQUALS(dgt323_bound, "none"  ) ||
      CCTK_EQUALS(dgt323_bound, "static") ||
      CCTK_EQUALS(dgt323_bound, "flat"  ) ||
      CCTK_EQUALS(dgt323_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt323", dgt323_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt323_bound BC for CL_BSSN::dgt323!");
  }
  
  if (CCTK_EQUALS(dgt333_bound, "none"  ) ||
      CCTK_EQUALS(dgt333_bound, "static") ||
      CCTK_EQUALS(dgt333_bound, "flat"  ) ||
      CCTK_EQUALS(dgt333_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dgt333", dgt333_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dgt333_bound BC for CL_BSSN::dgt333!");
  }
  
  if (CCTK_EQUALS(Xt1_bound, "none"  ) ||
      CCTK_EQUALS(Xt1_bound, "static") ||
      CCTK_EQUALS(Xt1_bound, "flat"  ) ||
      CCTK_EQUALS(Xt1_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::Xt1", Xt1_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register Xt1_bound BC for CL_BSSN::Xt1!");
  }
  
  if (CCTK_EQUALS(Xt2_bound, "none"  ) ||
      CCTK_EQUALS(Xt2_bound, "static") ||
      CCTK_EQUALS(Xt2_bound, "flat"  ) ||
      CCTK_EQUALS(Xt2_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::Xt2", Xt2_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register Xt2_bound BC for CL_BSSN::Xt2!");
  }
  
  if (CCTK_EQUALS(Xt3_bound, "none"  ) ||
      CCTK_EQUALS(Xt3_bound, "static") ||
      CCTK_EQUALS(Xt3_bound, "flat"  ) ||
      CCTK_EQUALS(Xt3_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::Xt3", Xt3_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register Xt3_bound BC for CL_BSSN::Xt3!");
  }
  
  if (CCTK_EQUALS(trK_bound, "none"  ) ||
      CCTK_EQUALS(trK_bound, "static") ||
      CCTK_EQUALS(trK_bound, "flat"  ) ||
      CCTK_EQUALS(trK_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::trK", trK_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register trK_bound BC for CL_BSSN::trK!");
  }
  
  if (CCTK_EQUALS(At11_bound, "none"  ) ||
      CCTK_EQUALS(At11_bound, "static") ||
      CCTK_EQUALS(At11_bound, "flat"  ) ||
      CCTK_EQUALS(At11_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::At11", At11_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register At11_bound BC for CL_BSSN::At11!");
  }
  
  if (CCTK_EQUALS(At12_bound, "none"  ) ||
      CCTK_EQUALS(At12_bound, "static") ||
      CCTK_EQUALS(At12_bound, "flat"  ) ||
      CCTK_EQUALS(At12_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::At12", At12_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register At12_bound BC for CL_BSSN::At12!");
  }
  
  if (CCTK_EQUALS(At13_bound, "none"  ) ||
      CCTK_EQUALS(At13_bound, "static") ||
      CCTK_EQUALS(At13_bound, "flat"  ) ||
      CCTK_EQUALS(At13_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::At13", At13_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register At13_bound BC for CL_BSSN::At13!");
  }
  
  if (CCTK_EQUALS(At22_bound, "none"  ) ||
      CCTK_EQUALS(At22_bound, "static") ||
      CCTK_EQUALS(At22_bound, "flat"  ) ||
      CCTK_EQUALS(At22_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::At22", At22_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register At22_bound BC for CL_BSSN::At22!");
  }
  
  if (CCTK_EQUALS(At23_bound, "none"  ) ||
      CCTK_EQUALS(At23_bound, "static") ||
      CCTK_EQUALS(At23_bound, "flat"  ) ||
      CCTK_EQUALS(At23_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::At23", At23_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register At23_bound BC for CL_BSSN::At23!");
  }
  
  if (CCTK_EQUALS(At33_bound, "none"  ) ||
      CCTK_EQUALS(At33_bound, "static") ||
      CCTK_EQUALS(At33_bound, "flat"  ) ||
      CCTK_EQUALS(At33_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::At33", At33_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register At33_bound BC for CL_BSSN::At33!");
  }
  
  if (CCTK_EQUALS(alpha_bound, "none"  ) ||
      CCTK_EQUALS(alpha_bound, "static") ||
      CCTK_EQUALS(alpha_bound, "flat"  ) ||
      CCTK_EQUALS(alpha_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::alpha", alpha_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register alpha_bound BC for CL_BSSN::alpha!");
  }
  
  if (CCTK_EQUALS(dalpha1_bound, "none"  ) ||
      CCTK_EQUALS(dalpha1_bound, "static") ||
      CCTK_EQUALS(dalpha1_bound, "flat"  ) ||
      CCTK_EQUALS(dalpha1_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dalpha1", dalpha1_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dalpha1_bound BC for CL_BSSN::dalpha1!");
  }
  
  if (CCTK_EQUALS(dalpha2_bound, "none"  ) ||
      CCTK_EQUALS(dalpha2_bound, "static") ||
      CCTK_EQUALS(dalpha2_bound, "flat"  ) ||
      CCTK_EQUALS(dalpha2_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dalpha2", dalpha2_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dalpha2_bound BC for CL_BSSN::dalpha2!");
  }
  
  if (CCTK_EQUALS(dalpha3_bound, "none"  ) ||
      CCTK_EQUALS(dalpha3_bound, "static") ||
      CCTK_EQUALS(dalpha3_bound, "flat"  ) ||
      CCTK_EQUALS(dalpha3_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dalpha3", dalpha3_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dalpha3_bound BC for CL_BSSN::dalpha3!");
  }
  
  if (CCTK_EQUALS(beta1_bound, "none"  ) ||
      CCTK_EQUALS(beta1_bound, "static") ||
      CCTK_EQUALS(beta1_bound, "flat"  ) ||
      CCTK_EQUALS(beta1_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::beta1", beta1_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register beta1_bound BC for CL_BSSN::beta1!");
  }
  
  if (CCTK_EQUALS(beta2_bound, "none"  ) ||
      CCTK_EQUALS(beta2_bound, "static") ||
      CCTK_EQUALS(beta2_bound, "flat"  ) ||
      CCTK_EQUALS(beta2_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::beta2", beta2_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register beta2_bound BC for CL_BSSN::beta2!");
  }
  
  if (CCTK_EQUALS(beta3_bound, "none"  ) ||
      CCTK_EQUALS(beta3_bound, "static") ||
      CCTK_EQUALS(beta3_bound, "flat"  ) ||
      CCTK_EQUALS(beta3_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::beta3", beta3_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register beta3_bound BC for CL_BSSN::beta3!");
  }
  
  if (CCTK_EQUALS(dbeta11_bound, "none"  ) ||
      CCTK_EQUALS(dbeta11_bound, "static") ||
      CCTK_EQUALS(dbeta11_bound, "flat"  ) ||
      CCTK_EQUALS(dbeta11_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dbeta11", dbeta11_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dbeta11_bound BC for CL_BSSN::dbeta11!");
  }
  
  if (CCTK_EQUALS(dbeta12_bound, "none"  ) ||
      CCTK_EQUALS(dbeta12_bound, "static") ||
      CCTK_EQUALS(dbeta12_bound, "flat"  ) ||
      CCTK_EQUALS(dbeta12_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dbeta12", dbeta12_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dbeta12_bound BC for CL_BSSN::dbeta12!");
  }
  
  if (CCTK_EQUALS(dbeta13_bound, "none"  ) ||
      CCTK_EQUALS(dbeta13_bound, "static") ||
      CCTK_EQUALS(dbeta13_bound, "flat"  ) ||
      CCTK_EQUALS(dbeta13_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dbeta13", dbeta13_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dbeta13_bound BC for CL_BSSN::dbeta13!");
  }
  
  if (CCTK_EQUALS(dbeta21_bound, "none"  ) ||
      CCTK_EQUALS(dbeta21_bound, "static") ||
      CCTK_EQUALS(dbeta21_bound, "flat"  ) ||
      CCTK_EQUALS(dbeta21_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dbeta21", dbeta21_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dbeta21_bound BC for CL_BSSN::dbeta21!");
  }
  
  if (CCTK_EQUALS(dbeta22_bound, "none"  ) ||
      CCTK_EQUALS(dbeta22_bound, "static") ||
      CCTK_EQUALS(dbeta22_bound, "flat"  ) ||
      CCTK_EQUALS(dbeta22_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dbeta22", dbeta22_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dbeta22_bound BC for CL_BSSN::dbeta22!");
  }
  
  if (CCTK_EQUALS(dbeta23_bound, "none"  ) ||
      CCTK_EQUALS(dbeta23_bound, "static") ||
      CCTK_EQUALS(dbeta23_bound, "flat"  ) ||
      CCTK_EQUALS(dbeta23_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dbeta23", dbeta23_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dbeta23_bound BC for CL_BSSN::dbeta23!");
  }
  
  if (CCTK_EQUALS(dbeta31_bound, "none"  ) ||
      CCTK_EQUALS(dbeta31_bound, "static") ||
      CCTK_EQUALS(dbeta31_bound, "flat"  ) ||
      CCTK_EQUALS(dbeta31_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dbeta31", dbeta31_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dbeta31_bound BC for CL_BSSN::dbeta31!");
  }
  
  if (CCTK_EQUALS(dbeta32_bound, "none"  ) ||
      CCTK_EQUALS(dbeta32_bound, "static") ||
      CCTK_EQUALS(dbeta32_bound, "flat"  ) ||
      CCTK_EQUALS(dbeta32_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dbeta32", dbeta32_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dbeta32_bound BC for CL_BSSN::dbeta32!");
  }
  
  if (CCTK_EQUALS(dbeta33_bound, "none"  ) ||
      CCTK_EQUALS(dbeta33_bound, "static") ||
      CCTK_EQUALS(dbeta33_bound, "flat"  ) ||
      CCTK_EQUALS(dbeta33_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::dbeta33", dbeta33_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register dbeta33_bound BC for CL_BSSN::dbeta33!");
  }
  
  if (CCTK_EQUALS(B1_bound, "none"  ) ||
      CCTK_EQUALS(B1_bound, "static") ||
      CCTK_EQUALS(B1_bound, "flat"  ) ||
      CCTK_EQUALS(B1_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::B1", B1_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register B1_bound BC for CL_BSSN::B1!");
  }
  
  if (CCTK_EQUALS(B2_bound, "none"  ) ||
      CCTK_EQUALS(B2_bound, "static") ||
      CCTK_EQUALS(B2_bound, "flat"  ) ||
      CCTK_EQUALS(B2_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::B2", B2_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register B2_bound BC for CL_BSSN::B2!");
  }
  
  if (CCTK_EQUALS(B3_bound, "none"  ) ||
      CCTK_EQUALS(B3_bound, "static") ||
      CCTK_EQUALS(B3_bound, "flat"  ) ||
      CCTK_EQUALS(B3_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CL_BSSN::B3", B3_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register B3_bound BC for CL_BSSN::B3!");
  }
  
  if (CCTK_EQUALS(CL_log_confac_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_CL_log_confac_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_log_confac_bound < 0) handle_CL_log_confac_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_log_confac_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_log_confac_bound , CL_log_confac_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_CL_log_confac_bound ,CL_log_confac_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_log_confac_bound, 
                      "CL_BSSN::CL_log_confac", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::CL_log_confac!");
  
  }
  
  if (CCTK_EQUALS(CL_dlog_confac_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_CL_dlog_confac_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_dlog_confac_bound < 0) handle_CL_dlog_confac_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_dlog_confac_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_dlog_confac_bound , CL_dlog_confac_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_CL_dlog_confac_bound ,CL_dlog_confac_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_dlog_confac_bound, 
                      "CL_BSSN::CL_dlog_confac", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::CL_dlog_confac!");
  
  }
  
  if (CCTK_EQUALS(CL_metric_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_CL_metric_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_metric_bound < 0) handle_CL_metric_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_metric_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_metric_bound , CL_metric_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_CL_metric_bound ,CL_metric_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_metric_bound, 
                      "CL_BSSN::CL_metric", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::CL_metric!");
  
  }
  
  if (CCTK_EQUALS(CL_dmetric_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_CL_dmetric_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_dmetric_bound < 0) handle_CL_dmetric_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_dmetric_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_dmetric_bound , CL_dmetric_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_CL_dmetric_bound ,CL_dmetric_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_dmetric_bound, 
                      "CL_BSSN::CL_dmetric", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::CL_dmetric!");
  
  }
  
  if (CCTK_EQUALS(CL_Gamma_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_CL_Gamma_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_Gamma_bound < 0) handle_CL_Gamma_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_Gamma_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_Gamma_bound , CL_Gamma_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_CL_Gamma_bound ,CL_Gamma_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_Gamma_bound, 
                      "CL_BSSN::CL_Gamma", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::CL_Gamma!");
  
  }
  
  if (CCTK_EQUALS(CL_trace_curv_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_CL_trace_curv_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_trace_curv_bound < 0) handle_CL_trace_curv_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_trace_curv_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_trace_curv_bound , CL_trace_curv_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_CL_trace_curv_bound ,CL_trace_curv_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_trace_curv_bound, 
                      "CL_BSSN::CL_trace_curv", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::CL_trace_curv!");
  
  }
  
  if (CCTK_EQUALS(CL_curv_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_CL_curv_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_curv_bound < 0) handle_CL_curv_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_curv_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_curv_bound , CL_curv_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_CL_curv_bound ,CL_curv_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_curv_bound, 
                      "CL_BSSN::CL_curv", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::CL_curv!");
  
  }
  
  if (CCTK_EQUALS(CL_lapse_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_CL_lapse_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_lapse_bound < 0) handle_CL_lapse_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_lapse_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_lapse_bound , CL_lapse_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_CL_lapse_bound ,CL_lapse_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_lapse_bound, 
                      "CL_BSSN::CL_lapse", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::CL_lapse!");
  
  }
  
  if (CCTK_EQUALS(CL_dlapse_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_CL_dlapse_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_dlapse_bound < 0) handle_CL_dlapse_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_dlapse_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_dlapse_bound , CL_dlapse_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_CL_dlapse_bound ,CL_dlapse_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_dlapse_bound, 
                      "CL_BSSN::CL_dlapse", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::CL_dlapse!");
  
  }
  
  if (CCTK_EQUALS(CL_shift_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_CL_shift_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_shift_bound < 0) handle_CL_shift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_shift_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_shift_bound , CL_shift_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_CL_shift_bound ,CL_shift_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_shift_bound, 
                      "CL_BSSN::CL_shift", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::CL_shift!");
  
  }
  
  if (CCTK_EQUALS(CL_dshift_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_CL_dshift_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_dshift_bound < 0) handle_CL_dshift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_dshift_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_dshift_bound , CL_dshift_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_CL_dshift_bound ,CL_dshift_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_dshift_bound, 
                      "CL_BSSN::CL_dshift", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::CL_dshift!");
  
  }
  
  if (CCTK_EQUALS(CL_dtshift_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_CL_dtshift_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_dtshift_bound < 0) handle_CL_dtshift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_dtshift_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_dtshift_bound , CL_dtshift_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_CL_dtshift_bound ,CL_dtshift_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_dtshift_bound, 
                      "CL_BSSN::CL_dtshift", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::CL_dtshift!");
  
  }
  
  if (CCTK_EQUALS(phi_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_phi_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_phi_bound < 0) handle_phi_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_phi_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_phi_bound , phi_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_phi_bound ,phi_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_phi_bound, 
                      "CL_BSSN::phi", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::phi!");
  
  }
  
  if (CCTK_EQUALS(dphi1_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dphi1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dphi1_bound < 0) handle_dphi1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dphi1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dphi1_bound , dphi1_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dphi1_bound ,dphi1_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dphi1_bound, 
                      "CL_BSSN::dphi1", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dphi1!");
  
  }
  
  if (CCTK_EQUALS(dphi2_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dphi2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dphi2_bound < 0) handle_dphi2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dphi2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dphi2_bound , dphi2_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dphi2_bound ,dphi2_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dphi2_bound, 
                      "CL_BSSN::dphi2", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dphi2!");
  
  }
  
  if (CCTK_EQUALS(dphi3_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dphi3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dphi3_bound < 0) handle_dphi3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dphi3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dphi3_bound , dphi3_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dphi3_bound ,dphi3_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dphi3_bound, 
                      "CL_BSSN::dphi3", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dphi3!");
  
  }
  
  if (CCTK_EQUALS(gt11_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_gt11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt11_bound < 0) handle_gt11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt11_bound , gt11_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt11_bound ,gt11_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt11_bound, 
                      "CL_BSSN::gt11", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::gt11!");
  
  }
  
  if (CCTK_EQUALS(gt12_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_gt12_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt12_bound < 0) handle_gt12_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt12_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt12_bound , gt12_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt12_bound ,gt12_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt12_bound, 
                      "CL_BSSN::gt12", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::gt12!");
  
  }
  
  if (CCTK_EQUALS(gt13_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_gt13_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt13_bound < 0) handle_gt13_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt13_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt13_bound , gt13_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt13_bound ,gt13_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt13_bound, 
                      "CL_BSSN::gt13", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::gt13!");
  
  }
  
  if (CCTK_EQUALS(gt22_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_gt22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt22_bound < 0) handle_gt22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt22_bound , gt22_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt22_bound ,gt22_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt22_bound, 
                      "CL_BSSN::gt22", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::gt22!");
  
  }
  
  if (CCTK_EQUALS(gt23_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_gt23_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt23_bound < 0) handle_gt23_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt23_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt23_bound , gt23_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt23_bound ,gt23_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt23_bound, 
                      "CL_BSSN::gt23", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::gt23!");
  
  }
  
  if (CCTK_EQUALS(gt33_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_gt33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt33_bound < 0) handle_gt33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt33_bound , gt33_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt33_bound ,gt33_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt33_bound, 
                      "CL_BSSN::gt33", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::gt33!");
  
  }
  
  if (CCTK_EQUALS(dgt111_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt111_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt111_bound < 0) handle_dgt111_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt111_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt111_bound , dgt111_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt111_bound ,dgt111_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt111_bound, 
                      "CL_BSSN::dgt111", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt111!");
  
  }
  
  if (CCTK_EQUALS(dgt112_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt112_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt112_bound < 0) handle_dgt112_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt112_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt112_bound , dgt112_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt112_bound ,dgt112_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt112_bound, 
                      "CL_BSSN::dgt112", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt112!");
  
  }
  
  if (CCTK_EQUALS(dgt113_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt113_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt113_bound < 0) handle_dgt113_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt113_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt113_bound , dgt113_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt113_bound ,dgt113_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt113_bound, 
                      "CL_BSSN::dgt113", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt113!");
  
  }
  
  if (CCTK_EQUALS(dgt122_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt122_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt122_bound < 0) handle_dgt122_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt122_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt122_bound , dgt122_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt122_bound ,dgt122_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt122_bound, 
                      "CL_BSSN::dgt122", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt122!");
  
  }
  
  if (CCTK_EQUALS(dgt123_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt123_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt123_bound < 0) handle_dgt123_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt123_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt123_bound , dgt123_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt123_bound ,dgt123_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt123_bound, 
                      "CL_BSSN::dgt123", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt123!");
  
  }
  
  if (CCTK_EQUALS(dgt133_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt133_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt133_bound < 0) handle_dgt133_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt133_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt133_bound , dgt133_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt133_bound ,dgt133_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt133_bound, 
                      "CL_BSSN::dgt133", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt133!");
  
  }
  
  if (CCTK_EQUALS(dgt211_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt211_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt211_bound < 0) handle_dgt211_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt211_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt211_bound , dgt211_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt211_bound ,dgt211_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt211_bound, 
                      "CL_BSSN::dgt211", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt211!");
  
  }
  
  if (CCTK_EQUALS(dgt212_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt212_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt212_bound < 0) handle_dgt212_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt212_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt212_bound , dgt212_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt212_bound ,dgt212_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt212_bound, 
                      "CL_BSSN::dgt212", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt212!");
  
  }
  
  if (CCTK_EQUALS(dgt213_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt213_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt213_bound < 0) handle_dgt213_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt213_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt213_bound , dgt213_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt213_bound ,dgt213_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt213_bound, 
                      "CL_BSSN::dgt213", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt213!");
  
  }
  
  if (CCTK_EQUALS(dgt222_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt222_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt222_bound < 0) handle_dgt222_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt222_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt222_bound , dgt222_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt222_bound ,dgt222_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt222_bound, 
                      "CL_BSSN::dgt222", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt222!");
  
  }
  
  if (CCTK_EQUALS(dgt223_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt223_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt223_bound < 0) handle_dgt223_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt223_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt223_bound , dgt223_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt223_bound ,dgt223_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt223_bound, 
                      "CL_BSSN::dgt223", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt223!");
  
  }
  
  if (CCTK_EQUALS(dgt233_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt233_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt233_bound < 0) handle_dgt233_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt233_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt233_bound , dgt233_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt233_bound ,dgt233_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt233_bound, 
                      "CL_BSSN::dgt233", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt233!");
  
  }
  
  if (CCTK_EQUALS(dgt311_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt311_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt311_bound < 0) handle_dgt311_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt311_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt311_bound , dgt311_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt311_bound ,dgt311_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt311_bound, 
                      "CL_BSSN::dgt311", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt311!");
  
  }
  
  if (CCTK_EQUALS(dgt312_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt312_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt312_bound < 0) handle_dgt312_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt312_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt312_bound , dgt312_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt312_bound ,dgt312_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt312_bound, 
                      "CL_BSSN::dgt312", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt312!");
  
  }
  
  if (CCTK_EQUALS(dgt313_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt313_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt313_bound < 0) handle_dgt313_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt313_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt313_bound , dgt313_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt313_bound ,dgt313_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt313_bound, 
                      "CL_BSSN::dgt313", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt313!");
  
  }
  
  if (CCTK_EQUALS(dgt322_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt322_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt322_bound < 0) handle_dgt322_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt322_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt322_bound , dgt322_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt322_bound ,dgt322_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt322_bound, 
                      "CL_BSSN::dgt322", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt322!");
  
  }
  
  if (CCTK_EQUALS(dgt323_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt323_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt323_bound < 0) handle_dgt323_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt323_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt323_bound , dgt323_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt323_bound ,dgt323_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt323_bound, 
                      "CL_BSSN::dgt323", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt323!");
  
  }
  
  if (CCTK_EQUALS(dgt333_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dgt333_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt333_bound < 0) handle_dgt333_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt333_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt333_bound , dgt333_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dgt333_bound ,dgt333_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt333_bound, 
                      "CL_BSSN::dgt333", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dgt333!");
  
  }
  
  if (CCTK_EQUALS(Xt1_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_Xt1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Xt1_bound < 0) handle_Xt1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Xt1_bound , Xt1_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_Xt1_bound ,Xt1_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt1_bound, 
                      "CL_BSSN::Xt1", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::Xt1!");
  
  }
  
  if (CCTK_EQUALS(Xt2_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_Xt2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Xt2_bound < 0) handle_Xt2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Xt2_bound , Xt2_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_Xt2_bound ,Xt2_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt2_bound, 
                      "CL_BSSN::Xt2", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::Xt2!");
  
  }
  
  if (CCTK_EQUALS(Xt3_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_Xt3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Xt3_bound < 0) handle_Xt3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Xt3_bound , Xt3_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_Xt3_bound ,Xt3_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt3_bound, 
                      "CL_BSSN::Xt3", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::Xt3!");
  
  }
  
  if (CCTK_EQUALS(trK_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_trK_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_trK_bound < 0) handle_trK_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_trK_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_trK_bound , trK_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_trK_bound ,trK_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_trK_bound, 
                      "CL_BSSN::trK", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::trK!");
  
  }
  
  if (CCTK_EQUALS(At11_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_At11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At11_bound < 0) handle_At11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At11_bound , At11_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At11_bound ,At11_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At11_bound, 
                      "CL_BSSN::At11", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::At11!");
  
  }
  
  if (CCTK_EQUALS(At12_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_At12_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At12_bound < 0) handle_At12_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At12_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At12_bound , At12_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At12_bound ,At12_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At12_bound, 
                      "CL_BSSN::At12", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::At12!");
  
  }
  
  if (CCTK_EQUALS(At13_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_At13_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At13_bound < 0) handle_At13_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At13_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At13_bound , At13_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At13_bound ,At13_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At13_bound, 
                      "CL_BSSN::At13", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::At13!");
  
  }
  
  if (CCTK_EQUALS(At22_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_At22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At22_bound < 0) handle_At22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At22_bound , At22_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At22_bound ,At22_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At22_bound, 
                      "CL_BSSN::At22", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::At22!");
  
  }
  
  if (CCTK_EQUALS(At23_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_At23_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At23_bound < 0) handle_At23_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At23_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At23_bound , At23_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At23_bound ,At23_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At23_bound, 
                      "CL_BSSN::At23", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::At23!");
  
  }
  
  if (CCTK_EQUALS(At33_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_At33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At33_bound < 0) handle_At33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At33_bound , At33_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At33_bound ,At33_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At33_bound, 
                      "CL_BSSN::At33", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::At33!");
  
  }
  
  if (CCTK_EQUALS(alpha_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_alpha_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_alpha_bound < 0) handle_alpha_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_alpha_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_alpha_bound , alpha_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_alpha_bound ,alpha_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_alpha_bound, 
                      "CL_BSSN::alpha", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::alpha!");
  
  }
  
  if (CCTK_EQUALS(dalpha1_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dalpha1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dalpha1_bound < 0) handle_dalpha1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dalpha1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dalpha1_bound , dalpha1_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dalpha1_bound ,dalpha1_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dalpha1_bound, 
                      "CL_BSSN::dalpha1", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dalpha1!");
  
  }
  
  if (CCTK_EQUALS(dalpha2_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dalpha2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dalpha2_bound < 0) handle_dalpha2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dalpha2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dalpha2_bound , dalpha2_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dalpha2_bound ,dalpha2_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dalpha2_bound, 
                      "CL_BSSN::dalpha2", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dalpha2!");
  
  }
  
  if (CCTK_EQUALS(dalpha3_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dalpha3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dalpha3_bound < 0) handle_dalpha3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dalpha3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dalpha3_bound , dalpha3_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dalpha3_bound ,dalpha3_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dalpha3_bound, 
                      "CL_BSSN::dalpha3", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dalpha3!");
  
  }
  
  if (CCTK_EQUALS(beta1_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_beta1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta1_bound < 0) handle_beta1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta1_bound , beta1_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_beta1_bound ,beta1_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta1_bound, 
                      "CL_BSSN::beta1", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::beta1!");
  
  }
  
  if (CCTK_EQUALS(beta2_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_beta2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta2_bound < 0) handle_beta2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta2_bound , beta2_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_beta2_bound ,beta2_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta2_bound, 
                      "CL_BSSN::beta2", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::beta2!");
  
  }
  
  if (CCTK_EQUALS(beta3_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_beta3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta3_bound < 0) handle_beta3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta3_bound , beta3_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_beta3_bound ,beta3_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta3_bound, 
                      "CL_BSSN::beta3", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::beta3!");
  
  }
  
  if (CCTK_EQUALS(dbeta11_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dbeta11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta11_bound < 0) handle_dbeta11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta11_bound , dbeta11_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dbeta11_bound ,dbeta11_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta11_bound, 
                      "CL_BSSN::dbeta11", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dbeta11!");
  
  }
  
  if (CCTK_EQUALS(dbeta12_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dbeta12_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta12_bound < 0) handle_dbeta12_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta12_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta12_bound , dbeta12_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dbeta12_bound ,dbeta12_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta12_bound, 
                      "CL_BSSN::dbeta12", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dbeta12!");
  
  }
  
  if (CCTK_EQUALS(dbeta13_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dbeta13_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta13_bound < 0) handle_dbeta13_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta13_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta13_bound , dbeta13_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dbeta13_bound ,dbeta13_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta13_bound, 
                      "CL_BSSN::dbeta13", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dbeta13!");
  
  }
  
  if (CCTK_EQUALS(dbeta21_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dbeta21_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta21_bound < 0) handle_dbeta21_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta21_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta21_bound , dbeta21_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dbeta21_bound ,dbeta21_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta21_bound, 
                      "CL_BSSN::dbeta21", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dbeta21!");
  
  }
  
  if (CCTK_EQUALS(dbeta22_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dbeta22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta22_bound < 0) handle_dbeta22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta22_bound , dbeta22_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dbeta22_bound ,dbeta22_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta22_bound, 
                      "CL_BSSN::dbeta22", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dbeta22!");
  
  }
  
  if (CCTK_EQUALS(dbeta23_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dbeta23_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta23_bound < 0) handle_dbeta23_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta23_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta23_bound , dbeta23_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dbeta23_bound ,dbeta23_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta23_bound, 
                      "CL_BSSN::dbeta23", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dbeta23!");
  
  }
  
  if (CCTK_EQUALS(dbeta31_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dbeta31_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta31_bound < 0) handle_dbeta31_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta31_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta31_bound , dbeta31_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dbeta31_bound ,dbeta31_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta31_bound, 
                      "CL_BSSN::dbeta31", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dbeta31!");
  
  }
  
  if (CCTK_EQUALS(dbeta32_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dbeta32_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta32_bound < 0) handle_dbeta32_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta32_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta32_bound , dbeta32_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dbeta32_bound ,dbeta32_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta32_bound, 
                      "CL_BSSN::dbeta32", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dbeta32!");
  
  }
  
  if (CCTK_EQUALS(dbeta33_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_dbeta33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta33_bound < 0) handle_dbeta33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta33_bound , dbeta33_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dbeta33_bound ,dbeta33_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta33_bound, 
                      "CL_BSSN::dbeta33", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::dbeta33!");
  
  }
  
  if (CCTK_EQUALS(B1_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_B1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B1_bound < 0) handle_B1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_B1_bound , B1_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_B1_bound ,B1_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B1_bound, 
                      "CL_BSSN::B1", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::B1!");
  
  }
  
  if (CCTK_EQUALS(B2_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_B2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B2_bound < 0) handle_B2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_B2_bound , B2_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_B2_bound ,B2_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B2_bound, 
                      "CL_BSSN::B2", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::B2!");
  
  }
  
  if (CCTK_EQUALS(B3_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_B3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B3_bound < 0) handle_B3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_B3_bound , B3_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_B3_bound ,B3_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B3_bound, 
                      "CL_BSSN::B3", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CL_BSSN::B3!");
  
  }
  
  if (CCTK_EQUALS(CL_log_confac_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_CL_log_confac_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_log_confac_bound < 0) handle_CL_log_confac_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_log_confac_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_log_confac_bound ,CL_log_confac_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_log_confac_bound, 
                      "CL_BSSN::CL_log_confac", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for CL_BSSN::CL_log_confac!");
  
  }
  
  if (CCTK_EQUALS(CL_dlog_confac_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_CL_dlog_confac_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_dlog_confac_bound < 0) handle_CL_dlog_confac_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_dlog_confac_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_dlog_confac_bound ,CL_dlog_confac_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_dlog_confac_bound, 
                      "CL_BSSN::CL_dlog_confac", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for CL_BSSN::CL_dlog_confac!");
  
  }
  
  if (CCTK_EQUALS(CL_metric_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_CL_metric_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_metric_bound < 0) handle_CL_metric_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_metric_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_metric_bound ,CL_metric_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_metric_bound, 
                      "CL_BSSN::CL_metric", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for CL_BSSN::CL_metric!");
  
  }
  
  if (CCTK_EQUALS(CL_dmetric_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_CL_dmetric_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_dmetric_bound < 0) handle_CL_dmetric_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_dmetric_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_dmetric_bound ,CL_dmetric_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_dmetric_bound, 
                      "CL_BSSN::CL_dmetric", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for CL_BSSN::CL_dmetric!");
  
  }
  
  if (CCTK_EQUALS(CL_Gamma_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_CL_Gamma_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_Gamma_bound < 0) handle_CL_Gamma_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_Gamma_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_Gamma_bound ,CL_Gamma_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_Gamma_bound, 
                      "CL_BSSN::CL_Gamma", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for CL_BSSN::CL_Gamma!");
  
  }
  
  if (CCTK_EQUALS(CL_trace_curv_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_CL_trace_curv_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_trace_curv_bound < 0) handle_CL_trace_curv_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_trace_curv_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_trace_curv_bound ,CL_trace_curv_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_trace_curv_bound, 
                      "CL_BSSN::CL_trace_curv", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for CL_BSSN::CL_trace_curv!");
  
  }
  
  if (CCTK_EQUALS(CL_curv_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_CL_curv_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_curv_bound < 0) handle_CL_curv_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_curv_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_curv_bound ,CL_curv_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_curv_bound, 
                      "CL_BSSN::CL_curv", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for CL_BSSN::CL_curv!");
  
  }
  
  if (CCTK_EQUALS(CL_lapse_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_CL_lapse_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_lapse_bound < 0) handle_CL_lapse_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_lapse_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_lapse_bound ,CL_lapse_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_lapse_bound, 
                      "CL_BSSN::CL_lapse", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for CL_BSSN::CL_lapse!");
  
  }
  
  if (CCTK_EQUALS(CL_dlapse_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_CL_dlapse_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_dlapse_bound < 0) handle_CL_dlapse_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_dlapse_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_dlapse_bound ,CL_dlapse_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_dlapse_bound, 
                      "CL_BSSN::CL_dlapse", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for CL_BSSN::CL_dlapse!");
  
  }
  
  if (CCTK_EQUALS(CL_shift_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_CL_shift_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_shift_bound < 0) handle_CL_shift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_shift_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_shift_bound ,CL_shift_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_shift_bound, 
                      "CL_BSSN::CL_shift", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for CL_BSSN::CL_shift!");
  
  }
  
  if (CCTK_EQUALS(CL_dshift_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_CL_dshift_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_dshift_bound < 0) handle_CL_dshift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_dshift_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_dshift_bound ,CL_dshift_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_dshift_bound, 
                      "CL_BSSN::CL_dshift", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for CL_BSSN::CL_dshift!");
  
  }
  
  if (CCTK_EQUALS(CL_dtshift_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_CL_dtshift_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CL_dtshift_bound < 0) handle_CL_dtshift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CL_dtshift_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CL_dtshift_bound ,CL_dtshift_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CL_dtshift_bound, 
                      "CL_BSSN::CL_dtshift", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for CL_BSSN::CL_dtshift!");
  
  }
  
  if (CCTK_EQUALS(phi_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_phi_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_phi_bound < 0) handle_phi_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_phi_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_phi_bound ,phi_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_phi_bound, 
                      "CL_BSSN::phi", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::phi!");
  
  }
  
  if (CCTK_EQUALS(dphi1_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dphi1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dphi1_bound < 0) handle_dphi1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dphi1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dphi1_bound ,dphi1_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dphi1_bound, 
                      "CL_BSSN::dphi1", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dphi1!");
  
  }
  
  if (CCTK_EQUALS(dphi2_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dphi2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dphi2_bound < 0) handle_dphi2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dphi2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dphi2_bound ,dphi2_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dphi2_bound, 
                      "CL_BSSN::dphi2", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dphi2!");
  
  }
  
  if (CCTK_EQUALS(dphi3_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dphi3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dphi3_bound < 0) handle_dphi3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dphi3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dphi3_bound ,dphi3_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dphi3_bound, 
                      "CL_BSSN::dphi3", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dphi3!");
  
  }
  
  if (CCTK_EQUALS(gt11_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_gt11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt11_bound < 0) handle_gt11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt11_bound ,gt11_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt11_bound, 
                      "CL_BSSN::gt11", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::gt11!");
  
  }
  
  if (CCTK_EQUALS(gt12_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_gt12_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt12_bound < 0) handle_gt12_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt12_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt12_bound ,gt12_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt12_bound, 
                      "CL_BSSN::gt12", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::gt12!");
  
  }
  
  if (CCTK_EQUALS(gt13_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_gt13_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt13_bound < 0) handle_gt13_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt13_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt13_bound ,gt13_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt13_bound, 
                      "CL_BSSN::gt13", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::gt13!");
  
  }
  
  if (CCTK_EQUALS(gt22_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_gt22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt22_bound < 0) handle_gt22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt22_bound ,gt22_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt22_bound, 
                      "CL_BSSN::gt22", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::gt22!");
  
  }
  
  if (CCTK_EQUALS(gt23_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_gt23_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt23_bound < 0) handle_gt23_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt23_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt23_bound ,gt23_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt23_bound, 
                      "CL_BSSN::gt23", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::gt23!");
  
  }
  
  if (CCTK_EQUALS(gt33_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_gt33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt33_bound < 0) handle_gt33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt33_bound ,gt33_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt33_bound, 
                      "CL_BSSN::gt33", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::gt33!");
  
  }
  
  if (CCTK_EQUALS(dgt111_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt111_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt111_bound < 0) handle_dgt111_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt111_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt111_bound ,dgt111_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt111_bound, 
                      "CL_BSSN::dgt111", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt111!");
  
  }
  
  if (CCTK_EQUALS(dgt112_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt112_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt112_bound < 0) handle_dgt112_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt112_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt112_bound ,dgt112_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt112_bound, 
                      "CL_BSSN::dgt112", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt112!");
  
  }
  
  if (CCTK_EQUALS(dgt113_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt113_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt113_bound < 0) handle_dgt113_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt113_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt113_bound ,dgt113_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt113_bound, 
                      "CL_BSSN::dgt113", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt113!");
  
  }
  
  if (CCTK_EQUALS(dgt122_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt122_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt122_bound < 0) handle_dgt122_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt122_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt122_bound ,dgt122_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt122_bound, 
                      "CL_BSSN::dgt122", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt122!");
  
  }
  
  if (CCTK_EQUALS(dgt123_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt123_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt123_bound < 0) handle_dgt123_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt123_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt123_bound ,dgt123_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt123_bound, 
                      "CL_BSSN::dgt123", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt123!");
  
  }
  
  if (CCTK_EQUALS(dgt133_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt133_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt133_bound < 0) handle_dgt133_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt133_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt133_bound ,dgt133_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt133_bound, 
                      "CL_BSSN::dgt133", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt133!");
  
  }
  
  if (CCTK_EQUALS(dgt211_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt211_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt211_bound < 0) handle_dgt211_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt211_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt211_bound ,dgt211_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt211_bound, 
                      "CL_BSSN::dgt211", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt211!");
  
  }
  
  if (CCTK_EQUALS(dgt212_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt212_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt212_bound < 0) handle_dgt212_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt212_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt212_bound ,dgt212_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt212_bound, 
                      "CL_BSSN::dgt212", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt212!");
  
  }
  
  if (CCTK_EQUALS(dgt213_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt213_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt213_bound < 0) handle_dgt213_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt213_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt213_bound ,dgt213_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt213_bound, 
                      "CL_BSSN::dgt213", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt213!");
  
  }
  
  if (CCTK_EQUALS(dgt222_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt222_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt222_bound < 0) handle_dgt222_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt222_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt222_bound ,dgt222_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt222_bound, 
                      "CL_BSSN::dgt222", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt222!");
  
  }
  
  if (CCTK_EQUALS(dgt223_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt223_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt223_bound < 0) handle_dgt223_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt223_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt223_bound ,dgt223_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt223_bound, 
                      "CL_BSSN::dgt223", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt223!");
  
  }
  
  if (CCTK_EQUALS(dgt233_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt233_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt233_bound < 0) handle_dgt233_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt233_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt233_bound ,dgt233_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt233_bound, 
                      "CL_BSSN::dgt233", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt233!");
  
  }
  
  if (CCTK_EQUALS(dgt311_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt311_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt311_bound < 0) handle_dgt311_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt311_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt311_bound ,dgt311_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt311_bound, 
                      "CL_BSSN::dgt311", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt311!");
  
  }
  
  if (CCTK_EQUALS(dgt312_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt312_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt312_bound < 0) handle_dgt312_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt312_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt312_bound ,dgt312_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt312_bound, 
                      "CL_BSSN::dgt312", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt312!");
  
  }
  
  if (CCTK_EQUALS(dgt313_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt313_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt313_bound < 0) handle_dgt313_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt313_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt313_bound ,dgt313_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt313_bound, 
                      "CL_BSSN::dgt313", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt313!");
  
  }
  
  if (CCTK_EQUALS(dgt322_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt322_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt322_bound < 0) handle_dgt322_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt322_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt322_bound ,dgt322_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt322_bound, 
                      "CL_BSSN::dgt322", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt322!");
  
  }
  
  if (CCTK_EQUALS(dgt323_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt323_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt323_bound < 0) handle_dgt323_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt323_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt323_bound ,dgt323_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt323_bound, 
                      "CL_BSSN::dgt323", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt323!");
  
  }
  
  if (CCTK_EQUALS(dgt333_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dgt333_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dgt333_bound < 0) handle_dgt333_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dgt333_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dgt333_bound ,dgt333_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dgt333_bound, 
                      "CL_BSSN::dgt333", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dgt333!");
  
  }
  
  if (CCTK_EQUALS(Xt1_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_Xt1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Xt1_bound < 0) handle_Xt1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Xt1_bound ,Xt1_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt1_bound, 
                      "CL_BSSN::Xt1", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::Xt1!");
  
  }
  
  if (CCTK_EQUALS(Xt2_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_Xt2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Xt2_bound < 0) handle_Xt2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Xt2_bound ,Xt2_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt2_bound, 
                      "CL_BSSN::Xt2", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::Xt2!");
  
  }
  
  if (CCTK_EQUALS(Xt3_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_Xt3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Xt3_bound < 0) handle_Xt3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Xt3_bound ,Xt3_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt3_bound, 
                      "CL_BSSN::Xt3", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::Xt3!");
  
  }
  
  if (CCTK_EQUALS(trK_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_trK_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_trK_bound < 0) handle_trK_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_trK_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_trK_bound ,trK_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_trK_bound, 
                      "CL_BSSN::trK", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::trK!");
  
  }
  
  if (CCTK_EQUALS(At11_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_At11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At11_bound < 0) handle_At11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At11_bound ,At11_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At11_bound, 
                      "CL_BSSN::At11", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::At11!");
  
  }
  
  if (CCTK_EQUALS(At12_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_At12_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At12_bound < 0) handle_At12_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At12_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At12_bound ,At12_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At12_bound, 
                      "CL_BSSN::At12", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::At12!");
  
  }
  
  if (CCTK_EQUALS(At13_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_At13_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At13_bound < 0) handle_At13_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At13_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At13_bound ,At13_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At13_bound, 
                      "CL_BSSN::At13", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::At13!");
  
  }
  
  if (CCTK_EQUALS(At22_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_At22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At22_bound < 0) handle_At22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At22_bound ,At22_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At22_bound, 
                      "CL_BSSN::At22", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::At22!");
  
  }
  
  if (CCTK_EQUALS(At23_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_At23_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At23_bound < 0) handle_At23_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At23_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At23_bound ,At23_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At23_bound, 
                      "CL_BSSN::At23", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::At23!");
  
  }
  
  if (CCTK_EQUALS(At33_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_At33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At33_bound < 0) handle_At33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At33_bound ,At33_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At33_bound, 
                      "CL_BSSN::At33", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::At33!");
  
  }
  
  if (CCTK_EQUALS(alpha_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_alpha_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_alpha_bound < 0) handle_alpha_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_alpha_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_alpha_bound ,alpha_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_alpha_bound, 
                      "CL_BSSN::alpha", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::alpha!");
  
  }
  
  if (CCTK_EQUALS(dalpha1_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dalpha1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dalpha1_bound < 0) handle_dalpha1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dalpha1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dalpha1_bound ,dalpha1_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dalpha1_bound, 
                      "CL_BSSN::dalpha1", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dalpha1!");
  
  }
  
  if (CCTK_EQUALS(dalpha2_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dalpha2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dalpha2_bound < 0) handle_dalpha2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dalpha2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dalpha2_bound ,dalpha2_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dalpha2_bound, 
                      "CL_BSSN::dalpha2", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dalpha2!");
  
  }
  
  if (CCTK_EQUALS(dalpha3_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dalpha3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dalpha3_bound < 0) handle_dalpha3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dalpha3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dalpha3_bound ,dalpha3_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dalpha3_bound, 
                      "CL_BSSN::dalpha3", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dalpha3!");
  
  }
  
  if (CCTK_EQUALS(beta1_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_beta1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta1_bound < 0) handle_beta1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta1_bound ,beta1_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta1_bound, 
                      "CL_BSSN::beta1", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::beta1!");
  
  }
  
  if (CCTK_EQUALS(beta2_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_beta2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta2_bound < 0) handle_beta2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta2_bound ,beta2_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta2_bound, 
                      "CL_BSSN::beta2", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::beta2!");
  
  }
  
  if (CCTK_EQUALS(beta3_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_beta3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta3_bound < 0) handle_beta3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta3_bound ,beta3_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta3_bound, 
                      "CL_BSSN::beta3", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::beta3!");
  
  }
  
  if (CCTK_EQUALS(dbeta11_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dbeta11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta11_bound < 0) handle_dbeta11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta11_bound ,dbeta11_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta11_bound, 
                      "CL_BSSN::dbeta11", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dbeta11!");
  
  }
  
  if (CCTK_EQUALS(dbeta12_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dbeta12_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta12_bound < 0) handle_dbeta12_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta12_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta12_bound ,dbeta12_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta12_bound, 
                      "CL_BSSN::dbeta12", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dbeta12!");
  
  }
  
  if (CCTK_EQUALS(dbeta13_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dbeta13_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta13_bound < 0) handle_dbeta13_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta13_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta13_bound ,dbeta13_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta13_bound, 
                      "CL_BSSN::dbeta13", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dbeta13!");
  
  }
  
  if (CCTK_EQUALS(dbeta21_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dbeta21_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta21_bound < 0) handle_dbeta21_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta21_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta21_bound ,dbeta21_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta21_bound, 
                      "CL_BSSN::dbeta21", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dbeta21!");
  
  }
  
  if (CCTK_EQUALS(dbeta22_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dbeta22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta22_bound < 0) handle_dbeta22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta22_bound ,dbeta22_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta22_bound, 
                      "CL_BSSN::dbeta22", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dbeta22!");
  
  }
  
  if (CCTK_EQUALS(dbeta23_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dbeta23_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta23_bound < 0) handle_dbeta23_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta23_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta23_bound ,dbeta23_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta23_bound, 
                      "CL_BSSN::dbeta23", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dbeta23!");
  
  }
  
  if (CCTK_EQUALS(dbeta31_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dbeta31_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta31_bound < 0) handle_dbeta31_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta31_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta31_bound ,dbeta31_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta31_bound, 
                      "CL_BSSN::dbeta31", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dbeta31!");
  
  }
  
  if (CCTK_EQUALS(dbeta32_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dbeta32_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta32_bound < 0) handle_dbeta32_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta32_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta32_bound ,dbeta32_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta32_bound, 
                      "CL_BSSN::dbeta32", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dbeta32!");
  
  }
  
  if (CCTK_EQUALS(dbeta33_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_dbeta33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_dbeta33_bound < 0) handle_dbeta33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dbeta33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_dbeta33_bound ,dbeta33_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dbeta33_bound, 
                      "CL_BSSN::dbeta33", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::dbeta33!");
  
  }
  
  if (CCTK_EQUALS(B1_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_B1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B1_bound < 0) handle_B1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_B1_bound ,B1_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B1_bound, 
                      "CL_BSSN::B1", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::B1!");
  
  }
  
  if (CCTK_EQUALS(B2_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_B2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B2_bound < 0) handle_B2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_B2_bound ,B2_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B2_bound, 
                      "CL_BSSN::B2", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::B2!");
  
  }
  
  if (CCTK_EQUALS(B3_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_B3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B3_bound < 0) handle_B3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_B3_bound ,B3_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B3_bound, 
                      "CL_BSSN::B3", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CL_BSSN::B3!");
  
  }
  return;
}



/* template for entries in parameter file:
#$bound$#CL_BSSN::CL_log_confac_bound       = "skip"
#$bound$#CL_BSSN::CL_log_confac_bound_speed = 1.0
#$bound$#CL_BSSN::CL_log_confac_bound_limit = 0.0
#$bound$#CL_BSSN::CL_log_confac_bound_scalar = 0.0

#$bound$#CL_BSSN::CL_dlog_confac_bound       = "skip"
#$bound$#CL_BSSN::CL_dlog_confac_bound_speed = 1.0
#$bound$#CL_BSSN::CL_dlog_confac_bound_limit = 0.0
#$bound$#CL_BSSN::CL_dlog_confac_bound_scalar = 0.0

#$bound$#CL_BSSN::CL_metric_bound       = "skip"
#$bound$#CL_BSSN::CL_metric_bound_speed = 1.0
#$bound$#CL_BSSN::CL_metric_bound_limit = 0.0
#$bound$#CL_BSSN::CL_metric_bound_scalar = 0.0

#$bound$#CL_BSSN::CL_dmetric_bound       = "skip"
#$bound$#CL_BSSN::CL_dmetric_bound_speed = 1.0
#$bound$#CL_BSSN::CL_dmetric_bound_limit = 0.0
#$bound$#CL_BSSN::CL_dmetric_bound_scalar = 0.0

#$bound$#CL_BSSN::CL_Gamma_bound       = "skip"
#$bound$#CL_BSSN::CL_Gamma_bound_speed = 1.0
#$bound$#CL_BSSN::CL_Gamma_bound_limit = 0.0
#$bound$#CL_BSSN::CL_Gamma_bound_scalar = 0.0

#$bound$#CL_BSSN::CL_trace_curv_bound       = "skip"
#$bound$#CL_BSSN::CL_trace_curv_bound_speed = 1.0
#$bound$#CL_BSSN::CL_trace_curv_bound_limit = 0.0
#$bound$#CL_BSSN::CL_trace_curv_bound_scalar = 0.0

#$bound$#CL_BSSN::CL_curv_bound       = "skip"
#$bound$#CL_BSSN::CL_curv_bound_speed = 1.0
#$bound$#CL_BSSN::CL_curv_bound_limit = 0.0
#$bound$#CL_BSSN::CL_curv_bound_scalar = 0.0

#$bound$#CL_BSSN::CL_lapse_bound       = "skip"
#$bound$#CL_BSSN::CL_lapse_bound_speed = 1.0
#$bound$#CL_BSSN::CL_lapse_bound_limit = 0.0
#$bound$#CL_BSSN::CL_lapse_bound_scalar = 0.0

#$bound$#CL_BSSN::CL_dlapse_bound       = "skip"
#$bound$#CL_BSSN::CL_dlapse_bound_speed = 1.0
#$bound$#CL_BSSN::CL_dlapse_bound_limit = 0.0
#$bound$#CL_BSSN::CL_dlapse_bound_scalar = 0.0

#$bound$#CL_BSSN::CL_shift_bound       = "skip"
#$bound$#CL_BSSN::CL_shift_bound_speed = 1.0
#$bound$#CL_BSSN::CL_shift_bound_limit = 0.0
#$bound$#CL_BSSN::CL_shift_bound_scalar = 0.0

#$bound$#CL_BSSN::CL_dshift_bound       = "skip"
#$bound$#CL_BSSN::CL_dshift_bound_speed = 1.0
#$bound$#CL_BSSN::CL_dshift_bound_limit = 0.0
#$bound$#CL_BSSN::CL_dshift_bound_scalar = 0.0

#$bound$#CL_BSSN::CL_dtshift_bound       = "skip"
#$bound$#CL_BSSN::CL_dtshift_bound_speed = 1.0
#$bound$#CL_BSSN::CL_dtshift_bound_limit = 0.0
#$bound$#CL_BSSN::CL_dtshift_bound_scalar = 0.0

#$bound$#CL_BSSN::phi_bound       = "skip"
#$bound$#CL_BSSN::phi_bound_speed = 1.0
#$bound$#CL_BSSN::phi_bound_limit = 0.0
#$bound$#CL_BSSN::phi_bound_scalar = 0.0

#$bound$#CL_BSSN::dphi1_bound       = "skip"
#$bound$#CL_BSSN::dphi1_bound_speed = 1.0
#$bound$#CL_BSSN::dphi1_bound_limit = 0.0
#$bound$#CL_BSSN::dphi1_bound_scalar = 0.0

#$bound$#CL_BSSN::dphi2_bound       = "skip"
#$bound$#CL_BSSN::dphi2_bound_speed = 1.0
#$bound$#CL_BSSN::dphi2_bound_limit = 0.0
#$bound$#CL_BSSN::dphi2_bound_scalar = 0.0

#$bound$#CL_BSSN::dphi3_bound       = "skip"
#$bound$#CL_BSSN::dphi3_bound_speed = 1.0
#$bound$#CL_BSSN::dphi3_bound_limit = 0.0
#$bound$#CL_BSSN::dphi3_bound_scalar = 0.0

#$bound$#CL_BSSN::gt11_bound       = "skip"
#$bound$#CL_BSSN::gt11_bound_speed = 1.0
#$bound$#CL_BSSN::gt11_bound_limit = 0.0
#$bound$#CL_BSSN::gt11_bound_scalar = 0.0

#$bound$#CL_BSSN::gt12_bound       = "skip"
#$bound$#CL_BSSN::gt12_bound_speed = 1.0
#$bound$#CL_BSSN::gt12_bound_limit = 0.0
#$bound$#CL_BSSN::gt12_bound_scalar = 0.0

#$bound$#CL_BSSN::gt13_bound       = "skip"
#$bound$#CL_BSSN::gt13_bound_speed = 1.0
#$bound$#CL_BSSN::gt13_bound_limit = 0.0
#$bound$#CL_BSSN::gt13_bound_scalar = 0.0

#$bound$#CL_BSSN::gt22_bound       = "skip"
#$bound$#CL_BSSN::gt22_bound_speed = 1.0
#$bound$#CL_BSSN::gt22_bound_limit = 0.0
#$bound$#CL_BSSN::gt22_bound_scalar = 0.0

#$bound$#CL_BSSN::gt23_bound       = "skip"
#$bound$#CL_BSSN::gt23_bound_speed = 1.0
#$bound$#CL_BSSN::gt23_bound_limit = 0.0
#$bound$#CL_BSSN::gt23_bound_scalar = 0.0

#$bound$#CL_BSSN::gt33_bound       = "skip"
#$bound$#CL_BSSN::gt33_bound_speed = 1.0
#$bound$#CL_BSSN::gt33_bound_limit = 0.0
#$bound$#CL_BSSN::gt33_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt111_bound       = "skip"
#$bound$#CL_BSSN::dgt111_bound_speed = 1.0
#$bound$#CL_BSSN::dgt111_bound_limit = 0.0
#$bound$#CL_BSSN::dgt111_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt112_bound       = "skip"
#$bound$#CL_BSSN::dgt112_bound_speed = 1.0
#$bound$#CL_BSSN::dgt112_bound_limit = 0.0
#$bound$#CL_BSSN::dgt112_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt113_bound       = "skip"
#$bound$#CL_BSSN::dgt113_bound_speed = 1.0
#$bound$#CL_BSSN::dgt113_bound_limit = 0.0
#$bound$#CL_BSSN::dgt113_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt122_bound       = "skip"
#$bound$#CL_BSSN::dgt122_bound_speed = 1.0
#$bound$#CL_BSSN::dgt122_bound_limit = 0.0
#$bound$#CL_BSSN::dgt122_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt123_bound       = "skip"
#$bound$#CL_BSSN::dgt123_bound_speed = 1.0
#$bound$#CL_BSSN::dgt123_bound_limit = 0.0
#$bound$#CL_BSSN::dgt123_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt133_bound       = "skip"
#$bound$#CL_BSSN::dgt133_bound_speed = 1.0
#$bound$#CL_BSSN::dgt133_bound_limit = 0.0
#$bound$#CL_BSSN::dgt133_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt211_bound       = "skip"
#$bound$#CL_BSSN::dgt211_bound_speed = 1.0
#$bound$#CL_BSSN::dgt211_bound_limit = 0.0
#$bound$#CL_BSSN::dgt211_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt212_bound       = "skip"
#$bound$#CL_BSSN::dgt212_bound_speed = 1.0
#$bound$#CL_BSSN::dgt212_bound_limit = 0.0
#$bound$#CL_BSSN::dgt212_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt213_bound       = "skip"
#$bound$#CL_BSSN::dgt213_bound_speed = 1.0
#$bound$#CL_BSSN::dgt213_bound_limit = 0.0
#$bound$#CL_BSSN::dgt213_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt222_bound       = "skip"
#$bound$#CL_BSSN::dgt222_bound_speed = 1.0
#$bound$#CL_BSSN::dgt222_bound_limit = 0.0
#$bound$#CL_BSSN::dgt222_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt223_bound       = "skip"
#$bound$#CL_BSSN::dgt223_bound_speed = 1.0
#$bound$#CL_BSSN::dgt223_bound_limit = 0.0
#$bound$#CL_BSSN::dgt223_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt233_bound       = "skip"
#$bound$#CL_BSSN::dgt233_bound_speed = 1.0
#$bound$#CL_BSSN::dgt233_bound_limit = 0.0
#$bound$#CL_BSSN::dgt233_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt311_bound       = "skip"
#$bound$#CL_BSSN::dgt311_bound_speed = 1.0
#$bound$#CL_BSSN::dgt311_bound_limit = 0.0
#$bound$#CL_BSSN::dgt311_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt312_bound       = "skip"
#$bound$#CL_BSSN::dgt312_bound_speed = 1.0
#$bound$#CL_BSSN::dgt312_bound_limit = 0.0
#$bound$#CL_BSSN::dgt312_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt313_bound       = "skip"
#$bound$#CL_BSSN::dgt313_bound_speed = 1.0
#$bound$#CL_BSSN::dgt313_bound_limit = 0.0
#$bound$#CL_BSSN::dgt313_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt322_bound       = "skip"
#$bound$#CL_BSSN::dgt322_bound_speed = 1.0
#$bound$#CL_BSSN::dgt322_bound_limit = 0.0
#$bound$#CL_BSSN::dgt322_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt323_bound       = "skip"
#$bound$#CL_BSSN::dgt323_bound_speed = 1.0
#$bound$#CL_BSSN::dgt323_bound_limit = 0.0
#$bound$#CL_BSSN::dgt323_bound_scalar = 0.0

#$bound$#CL_BSSN::dgt333_bound       = "skip"
#$bound$#CL_BSSN::dgt333_bound_speed = 1.0
#$bound$#CL_BSSN::dgt333_bound_limit = 0.0
#$bound$#CL_BSSN::dgt333_bound_scalar = 0.0

#$bound$#CL_BSSN::Xt1_bound       = "skip"
#$bound$#CL_BSSN::Xt1_bound_speed = 1.0
#$bound$#CL_BSSN::Xt1_bound_limit = 0.0
#$bound$#CL_BSSN::Xt1_bound_scalar = 0.0

#$bound$#CL_BSSN::Xt2_bound       = "skip"
#$bound$#CL_BSSN::Xt2_bound_speed = 1.0
#$bound$#CL_BSSN::Xt2_bound_limit = 0.0
#$bound$#CL_BSSN::Xt2_bound_scalar = 0.0

#$bound$#CL_BSSN::Xt3_bound       = "skip"
#$bound$#CL_BSSN::Xt3_bound_speed = 1.0
#$bound$#CL_BSSN::Xt3_bound_limit = 0.0
#$bound$#CL_BSSN::Xt3_bound_scalar = 0.0

#$bound$#CL_BSSN::trK_bound       = "skip"
#$bound$#CL_BSSN::trK_bound_speed = 1.0
#$bound$#CL_BSSN::trK_bound_limit = 0.0
#$bound$#CL_BSSN::trK_bound_scalar = 0.0

#$bound$#CL_BSSN::At11_bound       = "skip"
#$bound$#CL_BSSN::At11_bound_speed = 1.0
#$bound$#CL_BSSN::At11_bound_limit = 0.0
#$bound$#CL_BSSN::At11_bound_scalar = 0.0

#$bound$#CL_BSSN::At12_bound       = "skip"
#$bound$#CL_BSSN::At12_bound_speed = 1.0
#$bound$#CL_BSSN::At12_bound_limit = 0.0
#$bound$#CL_BSSN::At12_bound_scalar = 0.0

#$bound$#CL_BSSN::At13_bound       = "skip"
#$bound$#CL_BSSN::At13_bound_speed = 1.0
#$bound$#CL_BSSN::At13_bound_limit = 0.0
#$bound$#CL_BSSN::At13_bound_scalar = 0.0

#$bound$#CL_BSSN::At22_bound       = "skip"
#$bound$#CL_BSSN::At22_bound_speed = 1.0
#$bound$#CL_BSSN::At22_bound_limit = 0.0
#$bound$#CL_BSSN::At22_bound_scalar = 0.0

#$bound$#CL_BSSN::At23_bound       = "skip"
#$bound$#CL_BSSN::At23_bound_speed = 1.0
#$bound$#CL_BSSN::At23_bound_limit = 0.0
#$bound$#CL_BSSN::At23_bound_scalar = 0.0

#$bound$#CL_BSSN::At33_bound       = "skip"
#$bound$#CL_BSSN::At33_bound_speed = 1.0
#$bound$#CL_BSSN::At33_bound_limit = 0.0
#$bound$#CL_BSSN::At33_bound_scalar = 0.0

#$bound$#CL_BSSN::alpha_bound       = "skip"
#$bound$#CL_BSSN::alpha_bound_speed = 1.0
#$bound$#CL_BSSN::alpha_bound_limit = 0.0
#$bound$#CL_BSSN::alpha_bound_scalar = 0.0

#$bound$#CL_BSSN::dalpha1_bound       = "skip"
#$bound$#CL_BSSN::dalpha1_bound_speed = 1.0
#$bound$#CL_BSSN::dalpha1_bound_limit = 0.0
#$bound$#CL_BSSN::dalpha1_bound_scalar = 0.0

#$bound$#CL_BSSN::dalpha2_bound       = "skip"
#$bound$#CL_BSSN::dalpha2_bound_speed = 1.0
#$bound$#CL_BSSN::dalpha2_bound_limit = 0.0
#$bound$#CL_BSSN::dalpha2_bound_scalar = 0.0

#$bound$#CL_BSSN::dalpha3_bound       = "skip"
#$bound$#CL_BSSN::dalpha3_bound_speed = 1.0
#$bound$#CL_BSSN::dalpha3_bound_limit = 0.0
#$bound$#CL_BSSN::dalpha3_bound_scalar = 0.0

#$bound$#CL_BSSN::beta1_bound       = "skip"
#$bound$#CL_BSSN::beta1_bound_speed = 1.0
#$bound$#CL_BSSN::beta1_bound_limit = 0.0
#$bound$#CL_BSSN::beta1_bound_scalar = 0.0

#$bound$#CL_BSSN::beta2_bound       = "skip"
#$bound$#CL_BSSN::beta2_bound_speed = 1.0
#$bound$#CL_BSSN::beta2_bound_limit = 0.0
#$bound$#CL_BSSN::beta2_bound_scalar = 0.0

#$bound$#CL_BSSN::beta3_bound       = "skip"
#$bound$#CL_BSSN::beta3_bound_speed = 1.0
#$bound$#CL_BSSN::beta3_bound_limit = 0.0
#$bound$#CL_BSSN::beta3_bound_scalar = 0.0

#$bound$#CL_BSSN::dbeta11_bound       = "skip"
#$bound$#CL_BSSN::dbeta11_bound_speed = 1.0
#$bound$#CL_BSSN::dbeta11_bound_limit = 0.0
#$bound$#CL_BSSN::dbeta11_bound_scalar = 0.0

#$bound$#CL_BSSN::dbeta12_bound       = "skip"
#$bound$#CL_BSSN::dbeta12_bound_speed = 1.0
#$bound$#CL_BSSN::dbeta12_bound_limit = 0.0
#$bound$#CL_BSSN::dbeta12_bound_scalar = 0.0

#$bound$#CL_BSSN::dbeta13_bound       = "skip"
#$bound$#CL_BSSN::dbeta13_bound_speed = 1.0
#$bound$#CL_BSSN::dbeta13_bound_limit = 0.0
#$bound$#CL_BSSN::dbeta13_bound_scalar = 0.0

#$bound$#CL_BSSN::dbeta21_bound       = "skip"
#$bound$#CL_BSSN::dbeta21_bound_speed = 1.0
#$bound$#CL_BSSN::dbeta21_bound_limit = 0.0
#$bound$#CL_BSSN::dbeta21_bound_scalar = 0.0

#$bound$#CL_BSSN::dbeta22_bound       = "skip"
#$bound$#CL_BSSN::dbeta22_bound_speed = 1.0
#$bound$#CL_BSSN::dbeta22_bound_limit = 0.0
#$bound$#CL_BSSN::dbeta22_bound_scalar = 0.0

#$bound$#CL_BSSN::dbeta23_bound       = "skip"
#$bound$#CL_BSSN::dbeta23_bound_speed = 1.0
#$bound$#CL_BSSN::dbeta23_bound_limit = 0.0
#$bound$#CL_BSSN::dbeta23_bound_scalar = 0.0

#$bound$#CL_BSSN::dbeta31_bound       = "skip"
#$bound$#CL_BSSN::dbeta31_bound_speed = 1.0
#$bound$#CL_BSSN::dbeta31_bound_limit = 0.0
#$bound$#CL_BSSN::dbeta31_bound_scalar = 0.0

#$bound$#CL_BSSN::dbeta32_bound       = "skip"
#$bound$#CL_BSSN::dbeta32_bound_speed = 1.0
#$bound$#CL_BSSN::dbeta32_bound_limit = 0.0
#$bound$#CL_BSSN::dbeta32_bound_scalar = 0.0

#$bound$#CL_BSSN::dbeta33_bound       = "skip"
#$bound$#CL_BSSN::dbeta33_bound_speed = 1.0
#$bound$#CL_BSSN::dbeta33_bound_limit = 0.0
#$bound$#CL_BSSN::dbeta33_bound_scalar = 0.0

#$bound$#CL_BSSN::B1_bound       = "skip"
#$bound$#CL_BSSN::B1_bound_speed = 1.0
#$bound$#CL_BSSN::B1_bound_limit = 0.0
#$bound$#CL_BSSN::B1_bound_scalar = 0.0

#$bound$#CL_BSSN::B2_bound       = "skip"
#$bound$#CL_BSSN::B2_bound_speed = 1.0
#$bound$#CL_BSSN::B2_bound_limit = 0.0
#$bound$#CL_BSSN::B2_bound_scalar = 0.0

#$bound$#CL_BSSN::B3_bound       = "skip"
#$bound$#CL_BSSN::B3_bound_speed = 1.0
#$bound$#CL_BSSN::B3_bound_limit = 0.0
#$bound$#CL_BSSN::B3_bound_scalar = 0.0

*/

