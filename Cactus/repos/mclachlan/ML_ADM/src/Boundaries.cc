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


extern "C" void ML_ADM_CheckBoundaries(CCTK_ARGUMENTS)
{
#ifdef DECLARE_CCTK_ARGUMENTS_ML_ADM_CheckBoundaries
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_ADM_CheckBoundaries);
#else
  DECLARE_CCTK_ARGUMENTS;
#endif
  DECLARE_CCTK_PARAMETERS;
  
  return;
}

extern "C" void ML_ADM_SelectBoundConds(CCTK_ARGUMENTS)
{
#ifdef DECLARE_CCTK_ARGUMENTS_ML_ADM_SelectBoundConds
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_ADM_SelectBoundConds);
#else
  DECLARE_CCTK_ARGUMENTS;
#endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  
  if (CCTK_EQUALS(ML_metric_bound, "none"  ) ||
      CCTK_EQUALS(ML_metric_bound, "static") ||
      CCTK_EQUALS(ML_metric_bound, "flat"  ) ||
      CCTK_EQUALS(ML_metric_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::ML_metric", ML_metric_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register ML_metric_bound BC for ML_ADM::ML_metric!");
  }
  
  if (CCTK_EQUALS(ML_curv_bound, "none"  ) ||
      CCTK_EQUALS(ML_curv_bound, "static") ||
      CCTK_EQUALS(ML_curv_bound, "flat"  ) ||
      CCTK_EQUALS(ML_curv_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::ML_curv", ML_curv_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register ML_curv_bound BC for ML_ADM::ML_curv!");
  }
  
  if (CCTK_EQUALS(ML_lapse_bound, "none"  ) ||
      CCTK_EQUALS(ML_lapse_bound, "static") ||
      CCTK_EQUALS(ML_lapse_bound, "flat"  ) ||
      CCTK_EQUALS(ML_lapse_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::ML_lapse", ML_lapse_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register ML_lapse_bound BC for ML_ADM::ML_lapse!");
  }
  
  if (CCTK_EQUALS(ML_shift_bound, "none"  ) ||
      CCTK_EQUALS(ML_shift_bound, "static") ||
      CCTK_EQUALS(ML_shift_bound, "flat"  ) ||
      CCTK_EQUALS(ML_shift_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::ML_shift", ML_shift_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register ML_shift_bound BC for ML_ADM::ML_shift!");
  }
  
  if (CCTK_EQUALS(g11_bound, "none"  ) ||
      CCTK_EQUALS(g11_bound, "static") ||
      CCTK_EQUALS(g11_bound, "flat"  ) ||
      CCTK_EQUALS(g11_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::g11", g11_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register g11_bound BC for ML_ADM::g11!");
  }
  
  if (CCTK_EQUALS(g12_bound, "none"  ) ||
      CCTK_EQUALS(g12_bound, "static") ||
      CCTK_EQUALS(g12_bound, "flat"  ) ||
      CCTK_EQUALS(g12_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::g12", g12_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register g12_bound BC for ML_ADM::g12!");
  }
  
  if (CCTK_EQUALS(g13_bound, "none"  ) ||
      CCTK_EQUALS(g13_bound, "static") ||
      CCTK_EQUALS(g13_bound, "flat"  ) ||
      CCTK_EQUALS(g13_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::g13", g13_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register g13_bound BC for ML_ADM::g13!");
  }
  
  if (CCTK_EQUALS(g22_bound, "none"  ) ||
      CCTK_EQUALS(g22_bound, "static") ||
      CCTK_EQUALS(g22_bound, "flat"  ) ||
      CCTK_EQUALS(g22_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::g22", g22_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register g22_bound BC for ML_ADM::g22!");
  }
  
  if (CCTK_EQUALS(g23_bound, "none"  ) ||
      CCTK_EQUALS(g23_bound, "static") ||
      CCTK_EQUALS(g23_bound, "flat"  ) ||
      CCTK_EQUALS(g23_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::g23", g23_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register g23_bound BC for ML_ADM::g23!");
  }
  
  if (CCTK_EQUALS(g33_bound, "none"  ) ||
      CCTK_EQUALS(g33_bound, "static") ||
      CCTK_EQUALS(g33_bound, "flat"  ) ||
      CCTK_EQUALS(g33_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::g33", g33_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register g33_bound BC for ML_ADM::g33!");
  }
  
  if (CCTK_EQUALS(K11_bound, "none"  ) ||
      CCTK_EQUALS(K11_bound, "static") ||
      CCTK_EQUALS(K11_bound, "flat"  ) ||
      CCTK_EQUALS(K11_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::K11", K11_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register K11_bound BC for ML_ADM::K11!");
  }
  
  if (CCTK_EQUALS(K12_bound, "none"  ) ||
      CCTK_EQUALS(K12_bound, "static") ||
      CCTK_EQUALS(K12_bound, "flat"  ) ||
      CCTK_EQUALS(K12_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::K12", K12_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register K12_bound BC for ML_ADM::K12!");
  }
  
  if (CCTK_EQUALS(K13_bound, "none"  ) ||
      CCTK_EQUALS(K13_bound, "static") ||
      CCTK_EQUALS(K13_bound, "flat"  ) ||
      CCTK_EQUALS(K13_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::K13", K13_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register K13_bound BC for ML_ADM::K13!");
  }
  
  if (CCTK_EQUALS(K22_bound, "none"  ) ||
      CCTK_EQUALS(K22_bound, "static") ||
      CCTK_EQUALS(K22_bound, "flat"  ) ||
      CCTK_EQUALS(K22_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::K22", K22_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register K22_bound BC for ML_ADM::K22!");
  }
  
  if (CCTK_EQUALS(K23_bound, "none"  ) ||
      CCTK_EQUALS(K23_bound, "static") ||
      CCTK_EQUALS(K23_bound, "flat"  ) ||
      CCTK_EQUALS(K23_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::K23", K23_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register K23_bound BC for ML_ADM::K23!");
  }
  
  if (CCTK_EQUALS(K33_bound, "none"  ) ||
      CCTK_EQUALS(K33_bound, "static") ||
      CCTK_EQUALS(K33_bound, "flat"  ) ||
      CCTK_EQUALS(K33_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::K33", K33_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register K33_bound BC for ML_ADM::K33!");
  }
  
  if (CCTK_EQUALS(alpha_bound, "none"  ) ||
      CCTK_EQUALS(alpha_bound, "static") ||
      CCTK_EQUALS(alpha_bound, "flat"  ) ||
      CCTK_EQUALS(alpha_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::alpha", alpha_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register alpha_bound BC for ML_ADM::alpha!");
  }
  
  if (CCTK_EQUALS(beta1_bound, "none"  ) ||
      CCTK_EQUALS(beta1_bound, "static") ||
      CCTK_EQUALS(beta1_bound, "flat"  ) ||
      CCTK_EQUALS(beta1_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::beta1", beta1_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register beta1_bound BC for ML_ADM::beta1!");
  }
  
  if (CCTK_EQUALS(beta2_bound, "none"  ) ||
      CCTK_EQUALS(beta2_bound, "static") ||
      CCTK_EQUALS(beta2_bound, "flat"  ) ||
      CCTK_EQUALS(beta2_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::beta2", beta2_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register beta2_bound BC for ML_ADM::beta2!");
  }
  
  if (CCTK_EQUALS(beta3_bound, "none"  ) ||
      CCTK_EQUALS(beta3_bound, "static") ||
      CCTK_EQUALS(beta3_bound, "flat"  ) ||
      CCTK_EQUALS(beta3_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_ADM::beta3", beta3_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register beta3_bound BC for ML_ADM::beta3!");
  }
  
  if (CCTK_EQUALS(ML_metric_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_ML_metric_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_metric_bound < 0) handle_ML_metric_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_metric_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_metric_bound , ML_metric_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_ML_metric_bound ,ML_metric_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_metric_bound, 
                      "ML_ADM::ML_metric", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::ML_metric!");
  
  }
  
  if (CCTK_EQUALS(ML_curv_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_ML_curv_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_curv_bound < 0) handle_ML_curv_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_curv_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_curv_bound , ML_curv_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_ML_curv_bound ,ML_curv_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_curv_bound, 
                      "ML_ADM::ML_curv", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::ML_curv!");
  
  }
  
  if (CCTK_EQUALS(ML_lapse_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_ML_lapse_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_lapse_bound < 0) handle_ML_lapse_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_lapse_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_lapse_bound , ML_lapse_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_ML_lapse_bound ,ML_lapse_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_lapse_bound, 
                      "ML_ADM::ML_lapse", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::ML_lapse!");
  
  }
  
  if (CCTK_EQUALS(ML_shift_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_ML_shift_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_shift_bound < 0) handle_ML_shift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_shift_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_shift_bound , ML_shift_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_ML_shift_bound ,ML_shift_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_shift_bound, 
                      "ML_ADM::ML_shift", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::ML_shift!");
  
  }
  
  if (CCTK_EQUALS(g11_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_g11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_g11_bound < 0) handle_g11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_g11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_g11_bound , g11_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_g11_bound ,g11_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_g11_bound, 
                      "ML_ADM::g11", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::g11!");
  
  }
  
  if (CCTK_EQUALS(g12_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_g12_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_g12_bound < 0) handle_g12_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_g12_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_g12_bound , g12_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_g12_bound ,g12_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_g12_bound, 
                      "ML_ADM::g12", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::g12!");
  
  }
  
  if (CCTK_EQUALS(g13_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_g13_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_g13_bound < 0) handle_g13_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_g13_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_g13_bound , g13_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_g13_bound ,g13_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_g13_bound, 
                      "ML_ADM::g13", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::g13!");
  
  }
  
  if (CCTK_EQUALS(g22_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_g22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_g22_bound < 0) handle_g22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_g22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_g22_bound , g22_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_g22_bound ,g22_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_g22_bound, 
                      "ML_ADM::g22", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::g22!");
  
  }
  
  if (CCTK_EQUALS(g23_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_g23_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_g23_bound < 0) handle_g23_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_g23_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_g23_bound , g23_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_g23_bound ,g23_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_g23_bound, 
                      "ML_ADM::g23", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::g23!");
  
  }
  
  if (CCTK_EQUALS(g33_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_g33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_g33_bound < 0) handle_g33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_g33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_g33_bound , g33_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_g33_bound ,g33_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_g33_bound, 
                      "ML_ADM::g33", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::g33!");
  
  }
  
  if (CCTK_EQUALS(K11_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_K11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K11_bound < 0) handle_K11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K11_bound , K11_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_K11_bound ,K11_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K11_bound, 
                      "ML_ADM::K11", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::K11!");
  
  }
  
  if (CCTK_EQUALS(K12_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_K12_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K12_bound < 0) handle_K12_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K12_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K12_bound , K12_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_K12_bound ,K12_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K12_bound, 
                      "ML_ADM::K12", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::K12!");
  
  }
  
  if (CCTK_EQUALS(K13_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_K13_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K13_bound < 0) handle_K13_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K13_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K13_bound , K13_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_K13_bound ,K13_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K13_bound, 
                      "ML_ADM::K13", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::K13!");
  
  }
  
  if (CCTK_EQUALS(K22_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_K22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K22_bound < 0) handle_K22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K22_bound , K22_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_K22_bound ,K22_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K22_bound, 
                      "ML_ADM::K22", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::K22!");
  
  }
  
  if (CCTK_EQUALS(K23_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_K23_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K23_bound < 0) handle_K23_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K23_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K23_bound , K23_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_K23_bound ,K23_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K23_bound, 
                      "ML_ADM::K23", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::K23!");
  
  }
  
  if (CCTK_EQUALS(K33_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_K33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K33_bound < 0) handle_K33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K33_bound , K33_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_K33_bound ,K33_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K33_bound, 
                      "ML_ADM::K33", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::K33!");
  
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
                      "ML_ADM::alpha", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::alpha!");
  
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
                      "ML_ADM::beta1", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::beta1!");
  
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
                      "ML_ADM::beta2", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::beta2!");
  
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
                      "ML_ADM::beta3", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_ADM::beta3!");
  
  }
  
  if (CCTK_EQUALS(ML_metric_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_ML_metric_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_metric_bound < 0) handle_ML_metric_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_metric_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_metric_bound ,ML_metric_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_metric_bound, 
                      "ML_ADM::ML_metric", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_ADM::ML_metric!");
  
  }
  
  if (CCTK_EQUALS(ML_curv_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_ML_curv_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_curv_bound < 0) handle_ML_curv_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_curv_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_curv_bound ,ML_curv_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_curv_bound, 
                      "ML_ADM::ML_curv", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_ADM::ML_curv!");
  
  }
  
  if (CCTK_EQUALS(ML_lapse_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_ML_lapse_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_lapse_bound < 0) handle_ML_lapse_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_lapse_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_lapse_bound ,ML_lapse_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_lapse_bound, 
                      "ML_ADM::ML_lapse", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_ADM::ML_lapse!");
  
  }
  
  if (CCTK_EQUALS(ML_shift_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_ML_shift_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_shift_bound < 0) handle_ML_shift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_shift_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_shift_bound ,ML_shift_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_shift_bound, 
                      "ML_ADM::ML_shift", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_ADM::ML_shift!");
  
  }
  
  if (CCTK_EQUALS(g11_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_g11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_g11_bound < 0) handle_g11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_g11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_g11_bound ,g11_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_g11_bound, 
                      "ML_ADM::g11", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::g11!");
  
  }
  
  if (CCTK_EQUALS(g12_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_g12_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_g12_bound < 0) handle_g12_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_g12_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_g12_bound ,g12_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_g12_bound, 
                      "ML_ADM::g12", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::g12!");
  
  }
  
  if (CCTK_EQUALS(g13_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_g13_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_g13_bound < 0) handle_g13_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_g13_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_g13_bound ,g13_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_g13_bound, 
                      "ML_ADM::g13", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::g13!");
  
  }
  
  if (CCTK_EQUALS(g22_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_g22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_g22_bound < 0) handle_g22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_g22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_g22_bound ,g22_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_g22_bound, 
                      "ML_ADM::g22", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::g22!");
  
  }
  
  if (CCTK_EQUALS(g23_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_g23_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_g23_bound < 0) handle_g23_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_g23_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_g23_bound ,g23_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_g23_bound, 
                      "ML_ADM::g23", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::g23!");
  
  }
  
  if (CCTK_EQUALS(g33_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_g33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_g33_bound < 0) handle_g33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_g33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_g33_bound ,g33_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_g33_bound, 
                      "ML_ADM::g33", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::g33!");
  
  }
  
  if (CCTK_EQUALS(K11_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_K11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K11_bound < 0) handle_K11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K11_bound ,K11_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K11_bound, 
                      "ML_ADM::K11", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::K11!");
  
  }
  
  if (CCTK_EQUALS(K12_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_K12_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K12_bound < 0) handle_K12_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K12_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K12_bound ,K12_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K12_bound, 
                      "ML_ADM::K12", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::K12!");
  
  }
  
  if (CCTK_EQUALS(K13_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_K13_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K13_bound < 0) handle_K13_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K13_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K13_bound ,K13_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K13_bound, 
                      "ML_ADM::K13", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::K13!");
  
  }
  
  if (CCTK_EQUALS(K22_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_K22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K22_bound < 0) handle_K22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K22_bound ,K22_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K22_bound, 
                      "ML_ADM::K22", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::K22!");
  
  }
  
  if (CCTK_EQUALS(K23_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_K23_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K23_bound < 0) handle_K23_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K23_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K23_bound ,K23_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K23_bound, 
                      "ML_ADM::K23", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::K23!");
  
  }
  
  if (CCTK_EQUALS(K33_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_K33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K33_bound < 0) handle_K33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K33_bound ,K33_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K33_bound, 
                      "ML_ADM::K33", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::K33!");
  
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
                      "ML_ADM::alpha", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::alpha!");
  
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
                      "ML_ADM::beta1", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::beta1!");
  
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
                      "ML_ADM::beta2", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::beta2!");
  
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
                      "ML_ADM::beta3", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_ADM::beta3!");
  
  }
  return;
}



/* template for entries in parameter file:
#$bound$#ML_ADM::ML_metric_bound       = "skip"
#$bound$#ML_ADM::ML_metric_bound_speed = 1.0
#$bound$#ML_ADM::ML_metric_bound_limit = 0.0
#$bound$#ML_ADM::ML_metric_bound_scalar = 0.0

#$bound$#ML_ADM::ML_curv_bound       = "skip"
#$bound$#ML_ADM::ML_curv_bound_speed = 1.0
#$bound$#ML_ADM::ML_curv_bound_limit = 0.0
#$bound$#ML_ADM::ML_curv_bound_scalar = 0.0

#$bound$#ML_ADM::ML_lapse_bound       = "skip"
#$bound$#ML_ADM::ML_lapse_bound_speed = 1.0
#$bound$#ML_ADM::ML_lapse_bound_limit = 0.0
#$bound$#ML_ADM::ML_lapse_bound_scalar = 0.0

#$bound$#ML_ADM::ML_shift_bound       = "skip"
#$bound$#ML_ADM::ML_shift_bound_speed = 1.0
#$bound$#ML_ADM::ML_shift_bound_limit = 0.0
#$bound$#ML_ADM::ML_shift_bound_scalar = 0.0

#$bound$#ML_ADM::g11_bound       = "skip"
#$bound$#ML_ADM::g11_bound_speed = 1.0
#$bound$#ML_ADM::g11_bound_limit = 0.0
#$bound$#ML_ADM::g11_bound_scalar = 0.0

#$bound$#ML_ADM::g12_bound       = "skip"
#$bound$#ML_ADM::g12_bound_speed = 1.0
#$bound$#ML_ADM::g12_bound_limit = 0.0
#$bound$#ML_ADM::g12_bound_scalar = 0.0

#$bound$#ML_ADM::g13_bound       = "skip"
#$bound$#ML_ADM::g13_bound_speed = 1.0
#$bound$#ML_ADM::g13_bound_limit = 0.0
#$bound$#ML_ADM::g13_bound_scalar = 0.0

#$bound$#ML_ADM::g22_bound       = "skip"
#$bound$#ML_ADM::g22_bound_speed = 1.0
#$bound$#ML_ADM::g22_bound_limit = 0.0
#$bound$#ML_ADM::g22_bound_scalar = 0.0

#$bound$#ML_ADM::g23_bound       = "skip"
#$bound$#ML_ADM::g23_bound_speed = 1.0
#$bound$#ML_ADM::g23_bound_limit = 0.0
#$bound$#ML_ADM::g23_bound_scalar = 0.0

#$bound$#ML_ADM::g33_bound       = "skip"
#$bound$#ML_ADM::g33_bound_speed = 1.0
#$bound$#ML_ADM::g33_bound_limit = 0.0
#$bound$#ML_ADM::g33_bound_scalar = 0.0

#$bound$#ML_ADM::K11_bound       = "skip"
#$bound$#ML_ADM::K11_bound_speed = 1.0
#$bound$#ML_ADM::K11_bound_limit = 0.0
#$bound$#ML_ADM::K11_bound_scalar = 0.0

#$bound$#ML_ADM::K12_bound       = "skip"
#$bound$#ML_ADM::K12_bound_speed = 1.0
#$bound$#ML_ADM::K12_bound_limit = 0.0
#$bound$#ML_ADM::K12_bound_scalar = 0.0

#$bound$#ML_ADM::K13_bound       = "skip"
#$bound$#ML_ADM::K13_bound_speed = 1.0
#$bound$#ML_ADM::K13_bound_limit = 0.0
#$bound$#ML_ADM::K13_bound_scalar = 0.0

#$bound$#ML_ADM::K22_bound       = "skip"
#$bound$#ML_ADM::K22_bound_speed = 1.0
#$bound$#ML_ADM::K22_bound_limit = 0.0
#$bound$#ML_ADM::K22_bound_scalar = 0.0

#$bound$#ML_ADM::K23_bound       = "skip"
#$bound$#ML_ADM::K23_bound_speed = 1.0
#$bound$#ML_ADM::K23_bound_limit = 0.0
#$bound$#ML_ADM::K23_bound_scalar = 0.0

#$bound$#ML_ADM::K33_bound       = "skip"
#$bound$#ML_ADM::K33_bound_speed = 1.0
#$bound$#ML_ADM::K33_bound_limit = 0.0
#$bound$#ML_ADM::K33_bound_scalar = 0.0

#$bound$#ML_ADM::alpha_bound       = "skip"
#$bound$#ML_ADM::alpha_bound_speed = 1.0
#$bound$#ML_ADM::alpha_bound_limit = 0.0
#$bound$#ML_ADM::alpha_bound_scalar = 0.0

#$bound$#ML_ADM::beta1_bound       = "skip"
#$bound$#ML_ADM::beta1_bound_speed = 1.0
#$bound$#ML_ADM::beta1_bound_limit = 0.0
#$bound$#ML_ADM::beta1_bound_scalar = 0.0

#$bound$#ML_ADM::beta2_bound       = "skip"
#$bound$#ML_ADM::beta2_bound_speed = 1.0
#$bound$#ML_ADM::beta2_bound_limit = 0.0
#$bound$#ML_ADM::beta2_bound_scalar = 0.0

#$bound$#ML_ADM::beta3_bound       = "skip"
#$bound$#ML_ADM::beta3_bound_speed = 1.0
#$bound$#ML_ADM::beta3_bound_limit = 0.0
#$bound$#ML_ADM::beta3_bound_scalar = 0.0

*/

