/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_BSSN_ND_RegisterVars(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_BSSN_ND_RegisterVars
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_ND_RegisterVars);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::At11"),  CCTK_VarIndex("ML_BSSN_ND::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::At12"),  CCTK_VarIndex("ML_BSSN_ND::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::At13"),  CCTK_VarIndex("ML_BSSN_ND::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::At22"),  CCTK_VarIndex("ML_BSSN_ND::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::At23"),  CCTK_VarIndex("ML_BSSN_ND::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::At33"),  CCTK_VarIndex("ML_BSSN_ND::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::A"),  CCTK_VarIndex("ML_BSSN_ND::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::B1"),  CCTK_VarIndex("ML_BSSN_ND::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::B2"),  CCTK_VarIndex("ML_BSSN_ND::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::B3"),  CCTK_VarIndex("ML_BSSN_ND::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::Xt1"),  CCTK_VarIndex("ML_BSSN_ND::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::Xt2"),  CCTK_VarIndex("ML_BSSN_ND::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::Xt3"),  CCTK_VarIndex("ML_BSSN_ND::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::alpha"),  CCTK_VarIndex("ML_BSSN_ND::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::phi"),  CCTK_VarIndex("ML_BSSN_ND::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::gt11"),  CCTK_VarIndex("ML_BSSN_ND::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::gt12"),  CCTK_VarIndex("ML_BSSN_ND::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::gt13"),  CCTK_VarIndex("ML_BSSN_ND::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::gt22"),  CCTK_VarIndex("ML_BSSN_ND::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::gt23"),  CCTK_VarIndex("ML_BSSN_ND::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::gt33"),  CCTK_VarIndex("ML_BSSN_ND::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::beta1"),  CCTK_VarIndex("ML_BSSN_ND::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::beta2"),  CCTK_VarIndex("ML_BSSN_ND::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::beta3"),  CCTK_VarIndex("ML_BSSN_ND::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_ND::trK"),  CCTK_VarIndex("ML_BSSN_ND::trKrhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
