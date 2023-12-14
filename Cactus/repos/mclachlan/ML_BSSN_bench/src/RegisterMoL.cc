/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_BSSN_bench_RegisterVars(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_BSSN_bench_RegisterVars
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_bench_RegisterVars);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::At11"),  CCTK_VarIndex("ML_BSSN_bench::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::At12"),  CCTK_VarIndex("ML_BSSN_bench::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::At13"),  CCTK_VarIndex("ML_BSSN_bench::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::At22"),  CCTK_VarIndex("ML_BSSN_bench::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::At23"),  CCTK_VarIndex("ML_BSSN_bench::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::At33"),  CCTK_VarIndex("ML_BSSN_bench::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::B1"),  CCTK_VarIndex("ML_BSSN_bench::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::B2"),  CCTK_VarIndex("ML_BSSN_bench::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::B3"),  CCTK_VarIndex("ML_BSSN_bench::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::Xt1"),  CCTK_VarIndex("ML_BSSN_bench::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::Xt2"),  CCTK_VarIndex("ML_BSSN_bench::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::Xt3"),  CCTK_VarIndex("ML_BSSN_bench::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::alpha"),  CCTK_VarIndex("ML_BSSN_bench::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::phi"),  CCTK_VarIndex("ML_BSSN_bench::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::gt11"),  CCTK_VarIndex("ML_BSSN_bench::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::gt12"),  CCTK_VarIndex("ML_BSSN_bench::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::gt13"),  CCTK_VarIndex("ML_BSSN_bench::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::gt22"),  CCTK_VarIndex("ML_BSSN_bench::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::gt23"),  CCTK_VarIndex("ML_BSSN_bench::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::gt33"),  CCTK_VarIndex("ML_BSSN_bench::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::beta1"),  CCTK_VarIndex("ML_BSSN_bench::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::beta2"),  CCTK_VarIndex("ML_BSSN_bench::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::beta3"),  CCTK_VarIndex("ML_BSSN_bench::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_bench::trK"),  CCTK_VarIndex("ML_BSSN_bench::trKrhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
