/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_BSSN_NV_RegisterVars(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_BSSN_NV_RegisterVars
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_NV_RegisterVars);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::At11"),  CCTK_VarIndex("ML_BSSN_NV::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::At12"),  CCTK_VarIndex("ML_BSSN_NV::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::At13"),  CCTK_VarIndex("ML_BSSN_NV::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::At22"),  CCTK_VarIndex("ML_BSSN_NV::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::At23"),  CCTK_VarIndex("ML_BSSN_NV::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::At33"),  CCTK_VarIndex("ML_BSSN_NV::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::A"),  CCTK_VarIndex("ML_BSSN_NV::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::B1"),  CCTK_VarIndex("ML_BSSN_NV::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::B2"),  CCTK_VarIndex("ML_BSSN_NV::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::B3"),  CCTK_VarIndex("ML_BSSN_NV::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::Xt1"),  CCTK_VarIndex("ML_BSSN_NV::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::Xt2"),  CCTK_VarIndex("ML_BSSN_NV::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::Xt3"),  CCTK_VarIndex("ML_BSSN_NV::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::alpha"),  CCTK_VarIndex("ML_BSSN_NV::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::phi"),  CCTK_VarIndex("ML_BSSN_NV::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::gt11"),  CCTK_VarIndex("ML_BSSN_NV::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::gt12"),  CCTK_VarIndex("ML_BSSN_NV::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::gt13"),  CCTK_VarIndex("ML_BSSN_NV::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::gt22"),  CCTK_VarIndex("ML_BSSN_NV::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::gt23"),  CCTK_VarIndex("ML_BSSN_NV::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::gt33"),  CCTK_VarIndex("ML_BSSN_NV::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::beta1"),  CCTK_VarIndex("ML_BSSN_NV::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::beta2"),  CCTK_VarIndex("ML_BSSN_NV::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::beta3"),  CCTK_VarIndex("ML_BSSN_NV::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::Theta"),  CCTK_VarIndex("ML_BSSN_NV::Thetarhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NV::trK"),  CCTK_VarIndex("ML_BSSN_NV::trKrhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
