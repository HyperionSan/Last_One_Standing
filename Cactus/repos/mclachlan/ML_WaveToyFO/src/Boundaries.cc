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


extern "C" void ML_WaveToyFO_CheckBoundaries(CCTK_ARGUMENTS)
{
#ifdef DECLARE_CCTK_ARGUMENTS_ML_WaveToyFO_CheckBoundaries
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_WaveToyFO_CheckBoundaries);
#else
  DECLARE_CCTK_ARGUMENTS;
#endif
  DECLARE_CCTK_PARAMETERS;
  
  return;
}

extern "C" void ML_WaveToyFO_SelectBoundConds(CCTK_ARGUMENTS)
{
#ifdef DECLARE_CCTK_ARGUMENTS_ML_WaveToyFO_SelectBoundConds
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_WaveToyFO_SelectBoundConds);
#else
  DECLARE_CCTK_ARGUMENTS;
#endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  
  if (CCTK_EQUALS(WT_u_bound, "none"  ) ||
      CCTK_EQUALS(WT_u_bound, "static") ||
      CCTK_EQUALS(WT_u_bound, "flat"  ) ||
      CCTK_EQUALS(WT_u_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_WaveToyFO::WT_u", WT_u_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register WT_u_bound BC for ML_WaveToyFO::WT_u!");
  }
  
  if (CCTK_EQUALS(WT_v_bound, "none"  ) ||
      CCTK_EQUALS(WT_v_bound, "static") ||
      CCTK_EQUALS(WT_v_bound, "flat"  ) ||
      CCTK_EQUALS(WT_v_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_WaveToyFO::WT_v", WT_v_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register WT_v_bound BC for ML_WaveToyFO::WT_v!");
  }
  
  if (CCTK_EQUALS(WT_rho_bound, "none"  ) ||
      CCTK_EQUALS(WT_rho_bound, "static") ||
      CCTK_EQUALS(WT_rho_bound, "flat"  ) ||
      CCTK_EQUALS(WT_rho_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_WaveToyFO::WT_rho", WT_rho_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register WT_rho_bound BC for ML_WaveToyFO::WT_rho!");
  }
  
  if (CCTK_EQUALS(u_bound, "none"  ) ||
      CCTK_EQUALS(u_bound, "static") ||
      CCTK_EQUALS(u_bound, "flat"  ) ||
      CCTK_EQUALS(u_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_WaveToyFO::u", u_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register u_bound BC for ML_WaveToyFO::u!");
  }
  
  if (CCTK_EQUALS(v1_bound, "none"  ) ||
      CCTK_EQUALS(v1_bound, "static") ||
      CCTK_EQUALS(v1_bound, "flat"  ) ||
      CCTK_EQUALS(v1_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_WaveToyFO::v1", v1_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register v1_bound BC for ML_WaveToyFO::v1!");
  }
  
  if (CCTK_EQUALS(v2_bound, "none"  ) ||
      CCTK_EQUALS(v2_bound, "static") ||
      CCTK_EQUALS(v2_bound, "flat"  ) ||
      CCTK_EQUALS(v2_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_WaveToyFO::v2", v2_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register v2_bound BC for ML_WaveToyFO::v2!");
  }
  
  if (CCTK_EQUALS(v3_bound, "none"  ) ||
      CCTK_EQUALS(v3_bound, "static") ||
      CCTK_EQUALS(v3_bound, "flat"  ) ||
      CCTK_EQUALS(v3_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_WaveToyFO::v3", v3_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register v3_bound BC for ML_WaveToyFO::v3!");
  }
  
  if (CCTK_EQUALS(rho_bound, "none"  ) ||
      CCTK_EQUALS(rho_bound, "static") ||
      CCTK_EQUALS(rho_bound, "flat"  ) ||
      CCTK_EQUALS(rho_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_WaveToyFO::rho", rho_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register rho_bound BC for ML_WaveToyFO::rho!");
  }
  
  if (CCTK_EQUALS(WT_u_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_WT_u_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_WT_u_bound < 0) handle_WT_u_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_WT_u_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_WT_u_bound , WT_u_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_WT_u_bound ,WT_u_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_WT_u_bound, 
                      "ML_WaveToyFO::WT_u", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_WaveToyFO::WT_u!");
  
  }
  
  if (CCTK_EQUALS(WT_v_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_WT_v_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_WT_v_bound < 0) handle_WT_v_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_WT_v_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_WT_v_bound , WT_v_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_WT_v_bound ,WT_v_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_WT_v_bound, 
                      "ML_WaveToyFO::WT_v", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_WaveToyFO::WT_v!");
  
  }
  
  if (CCTK_EQUALS(WT_rho_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_WT_rho_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_WT_rho_bound < 0) handle_WT_rho_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_WT_rho_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_WT_rho_bound , WT_rho_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_WT_rho_bound ,WT_rho_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_WT_rho_bound, 
                      "ML_WaveToyFO::WT_rho", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_WaveToyFO::WT_rho!");
  
  }
  
  if (CCTK_EQUALS(u_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_u_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_u_bound < 0) handle_u_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_u_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_u_bound , u_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_u_bound ,u_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_u_bound, 
                      "ML_WaveToyFO::u", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_WaveToyFO::u!");
  
  }
  
  if (CCTK_EQUALS(v1_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_v1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_v1_bound < 0) handle_v1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_v1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_v1_bound , v1_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_v1_bound ,v1_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_v1_bound, 
                      "ML_WaveToyFO::v1", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_WaveToyFO::v1!");
  
  }
  
  if (CCTK_EQUALS(v2_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_v2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_v2_bound < 0) handle_v2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_v2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_v2_bound , v2_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_v2_bound ,v2_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_v2_bound, 
                      "ML_WaveToyFO::v2", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_WaveToyFO::v2!");
  
  }
  
  if (CCTK_EQUALS(v3_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_v3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_v3_bound < 0) handle_v3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_v3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_v3_bound , v3_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_v3_bound ,v3_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_v3_bound, 
                      "ML_WaveToyFO::v3", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_WaveToyFO::v3!");
  
  }
  
  if (CCTK_EQUALS(rho_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_rho_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_rho_bound < 0) handle_rho_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_rho_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_rho_bound , rho_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_rho_bound ,rho_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_rho_bound, 
                      "ML_WaveToyFO::rho", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_WaveToyFO::rho!");
  
  }
  
  if (CCTK_EQUALS(WT_u_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_WT_u_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_WT_u_bound < 0) handle_WT_u_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_WT_u_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_WT_u_bound ,WT_u_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_WT_u_bound, 
                      "ML_WaveToyFO::WT_u", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_WaveToyFO::WT_u!");
  
  }
  
  if (CCTK_EQUALS(WT_v_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_WT_v_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_WT_v_bound < 0) handle_WT_v_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_WT_v_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_WT_v_bound ,WT_v_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_WT_v_bound, 
                      "ML_WaveToyFO::WT_v", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_WaveToyFO::WT_v!");
  
  }
  
  if (CCTK_EQUALS(WT_rho_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_WT_rho_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_WT_rho_bound < 0) handle_WT_rho_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_WT_rho_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_WT_rho_bound ,WT_rho_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_WT_rho_bound, 
                      "ML_WaveToyFO::WT_rho", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_WaveToyFO::WT_rho!");
  
  }
  
  if (CCTK_EQUALS(u_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_u_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_u_bound < 0) handle_u_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_u_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_u_bound ,u_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_u_bound, 
                      "ML_WaveToyFO::u", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_WaveToyFO::u!");
  
  }
  
  if (CCTK_EQUALS(v1_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_v1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_v1_bound < 0) handle_v1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_v1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_v1_bound ,v1_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_v1_bound, 
                      "ML_WaveToyFO::v1", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_WaveToyFO::v1!");
  
  }
  
  if (CCTK_EQUALS(v2_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_v2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_v2_bound < 0) handle_v2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_v2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_v2_bound ,v2_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_v2_bound, 
                      "ML_WaveToyFO::v2", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_WaveToyFO::v2!");
  
  }
  
  if (CCTK_EQUALS(v3_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_v3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_v3_bound < 0) handle_v3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_v3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_v3_bound ,v3_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_v3_bound, 
                      "ML_WaveToyFO::v3", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_WaveToyFO::v3!");
  
  }
  
  if (CCTK_EQUALS(rho_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_rho_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_rho_bound < 0) handle_rho_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_rho_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_rho_bound ,rho_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_rho_bound, 
                      "ML_WaveToyFO::rho", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_WaveToyFO::rho!");
  
  }
  return;
}



/* template for entries in parameter file:
#$bound$#ML_WaveToyFO::WT_u_bound       = "skip"
#$bound$#ML_WaveToyFO::WT_u_bound_speed = 1.0
#$bound$#ML_WaveToyFO::WT_u_bound_limit = 0.0
#$bound$#ML_WaveToyFO::WT_u_bound_scalar = 0.0

#$bound$#ML_WaveToyFO::WT_v_bound       = "skip"
#$bound$#ML_WaveToyFO::WT_v_bound_speed = 1.0
#$bound$#ML_WaveToyFO::WT_v_bound_limit = 0.0
#$bound$#ML_WaveToyFO::WT_v_bound_scalar = 0.0

#$bound$#ML_WaveToyFO::WT_rho_bound       = "skip"
#$bound$#ML_WaveToyFO::WT_rho_bound_speed = 1.0
#$bound$#ML_WaveToyFO::WT_rho_bound_limit = 0.0
#$bound$#ML_WaveToyFO::WT_rho_bound_scalar = 0.0

#$bound$#ML_WaveToyFO::u_bound       = "skip"
#$bound$#ML_WaveToyFO::u_bound_speed = 1.0
#$bound$#ML_WaveToyFO::u_bound_limit = 0.0
#$bound$#ML_WaveToyFO::u_bound_scalar = 0.0

#$bound$#ML_WaveToyFO::v1_bound       = "skip"
#$bound$#ML_WaveToyFO::v1_bound_speed = 1.0
#$bound$#ML_WaveToyFO::v1_bound_limit = 0.0
#$bound$#ML_WaveToyFO::v1_bound_scalar = 0.0

#$bound$#ML_WaveToyFO::v2_bound       = "skip"
#$bound$#ML_WaveToyFO::v2_bound_speed = 1.0
#$bound$#ML_WaveToyFO::v2_bound_limit = 0.0
#$bound$#ML_WaveToyFO::v2_bound_scalar = 0.0

#$bound$#ML_WaveToyFO::v3_bound       = "skip"
#$bound$#ML_WaveToyFO::v3_bound_speed = 1.0
#$bound$#ML_WaveToyFO::v3_bound_limit = 0.0
#$bound$#ML_WaveToyFO::v3_bound_scalar = 0.0

#$bound$#ML_WaveToyFO::rho_bound       = "skip"
#$bound$#ML_WaveToyFO::rho_bound_speed = 1.0
#$bound$#ML_WaveToyFO::rho_bound_limit = 0.0
#$bound$#ML_WaveToyFO::rho_bound_scalar = 0.0

*/

