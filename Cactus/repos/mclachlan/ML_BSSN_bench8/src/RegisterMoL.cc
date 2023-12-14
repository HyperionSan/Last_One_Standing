/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_BSSN_bench8_RegisterVars(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_BSSN_bench8_RegisterVars
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_bench8_RegisterVars);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::At11"),  CCTK_VarIndex("ML_BSSN_bench8::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::At12"),  CCTK_VarIndex("ML_BSSN_bench8::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::At13"),  CCTK_VarIndex("ML_BSSN_bench8::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::At22"),  CCTK_VarIndex("ML_BSSN_bench8::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::At23"),  CCTK_VarIndex("ML_BSSN_bench8::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::At33"),  CCTK_VarIndex("ML_BSSN_bench8::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::Xt1"),  CCTK_VarIndex("ML_BSSN_bench8::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::Xt2"),  CCTK_VarIndex("ML_BSSN_bench8::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::Xt3"),  CCTK_VarIndex("ML_BSSN_bench8::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::alpha"),  CCTK_VarIndex("ML_BSSN_bench8::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::phi"),  CCTK_VarIndex("ML_BSSN_bench8::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::gt11"),  CCTK_VarIndex("ML_BSSN_bench8::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::gt12"),  CCTK_VarIndex("ML_BSSN_bench8::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::gt13"),  CCTK_VarIndex("ML_BSSN_bench8::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::gt22"),  CCTK_VarIndex("ML_BSSN_bench8::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::gt23"),  CCTK_VarIndex("ML_BSSN_bench8::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::gt33"),  CCTK_VarIndex("ML_BSSN_bench8::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::beta1"),  CCTK_VarIndex("ML_BSSN_bench8::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::beta2"),  CCTK_VarIndex("ML_BSSN_bench8::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::beta3"),  CCTK_VarIndex("ML_BSSN_bench8::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench8::trK"),  CCTK_VarIndex("ML_BSSN_bench8::trKrhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
