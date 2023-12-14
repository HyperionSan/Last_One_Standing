/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_WaveToyFO_RegisterVars(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_WaveToyFO_RegisterVars
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_WaveToyFO_RegisterVars);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_WaveToyFO::rho"),  CCTK_VarIndex("ML_WaveToyFO::rhorhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_WaveToyFO::u"),  CCTK_VarIndex("ML_WaveToyFO::urhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_WaveToyFO::v1"),  CCTK_VarIndex("ML_WaveToyFO::v1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_WaveToyFO::v2"),  CCTK_VarIndex("ML_WaveToyFO::v2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_WaveToyFO::v3"),  CCTK_VarIndex("ML_WaveToyFO::v3rhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
