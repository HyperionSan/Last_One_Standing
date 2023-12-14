/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_BSSN_CL_RegisterVars(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_BSSN_CL_RegisterVars
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_CL_RegisterVars);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::At11"),  CCTK_VarIndex("ML_BSSN_CL::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::At12"),  CCTK_VarIndex("ML_BSSN_CL::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::At13"),  CCTK_VarIndex("ML_BSSN_CL::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::At22"),  CCTK_VarIndex("ML_BSSN_CL::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::At23"),  CCTK_VarIndex("ML_BSSN_CL::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::At33"),  CCTK_VarIndex("ML_BSSN_CL::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::A"),  CCTK_VarIndex("ML_BSSN_CL::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::B1"),  CCTK_VarIndex("ML_BSSN_CL::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::B2"),  CCTK_VarIndex("ML_BSSN_CL::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::B3"),  CCTK_VarIndex("ML_BSSN_CL::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::Xt1"),  CCTK_VarIndex("ML_BSSN_CL::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::Xt2"),  CCTK_VarIndex("ML_BSSN_CL::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::Xt3"),  CCTK_VarIndex("ML_BSSN_CL::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::alpha"),  CCTK_VarIndex("ML_BSSN_CL::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::phi"),  CCTK_VarIndex("ML_BSSN_CL::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::gt11"),  CCTK_VarIndex("ML_BSSN_CL::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::gt12"),  CCTK_VarIndex("ML_BSSN_CL::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::gt13"),  CCTK_VarIndex("ML_BSSN_CL::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::gt22"),  CCTK_VarIndex("ML_BSSN_CL::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::gt23"),  CCTK_VarIndex("ML_BSSN_CL::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::gt33"),  CCTK_VarIndex("ML_BSSN_CL::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::beta1"),  CCTK_VarIndex("ML_BSSN_CL::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::beta2"),  CCTK_VarIndex("ML_BSSN_CL::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::beta3"),  CCTK_VarIndex("ML_BSSN_CL::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::Theta"),  CCTK_VarIndex("ML_BSSN_CL::Thetarhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CL::trK"),  CCTK_VarIndex("ML_BSSN_CL::trKrhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
