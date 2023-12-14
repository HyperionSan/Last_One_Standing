/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_hydro_RegisterVars(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_hydro_RegisterVars
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_hydro_RegisterVars);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_hydro::ene"),  CCTK_VarIndex("ML_hydro::enerhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_hydro::mass"),  CCTK_VarIndex("ML_hydro::massrhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_hydro::mom1"),  CCTK_VarIndex("ML_hydro::mom1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_hydro::mom2"),  CCTK_VarIndex("ML_hydro::mom2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_hydro::mom3"),  CCTK_VarIndex("ML_hydro::mom3rhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
