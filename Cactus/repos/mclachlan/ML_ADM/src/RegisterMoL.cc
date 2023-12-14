/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_ADM_RegisterVars(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_ADM_RegisterVars
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_ADM_RegisterVars);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::K11"),  CCTK_VarIndex("ML_ADM::K11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::K12"),  CCTK_VarIndex("ML_ADM::K12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::K13"),  CCTK_VarIndex("ML_ADM::K13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::K22"),  CCTK_VarIndex("ML_ADM::K22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::K23"),  CCTK_VarIndex("ML_ADM::K23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::K33"),  CCTK_VarIndex("ML_ADM::K33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::alpha"),  CCTK_VarIndex("ML_ADM::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::g11"),  CCTK_VarIndex("ML_ADM::g11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::g12"),  CCTK_VarIndex("ML_ADM::g12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::g13"),  CCTK_VarIndex("ML_ADM::g13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::g22"),  CCTK_VarIndex("ML_ADM::g22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::g23"),  CCTK_VarIndex("ML_ADM::g23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::g33"),  CCTK_VarIndex("ML_ADM::g33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::beta1"),  CCTK_VarIndex("ML_ADM::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::beta2"),  CCTK_VarIndex("ML_ADM::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::beta3"),  CCTK_VarIndex("ML_ADM::beta3rhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
