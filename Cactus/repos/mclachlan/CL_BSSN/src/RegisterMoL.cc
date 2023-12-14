/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void CL_BSSN_RegisterVars(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CL_BSSN_RegisterVars
  DECLARE_CCTK_ARGUMENTS_CHECKED(CL_BSSN_RegisterVars);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::At11"),  CCTK_VarIndex("CL_BSSN::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::At12"),  CCTK_VarIndex("CL_BSSN::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::At13"),  CCTK_VarIndex("CL_BSSN::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::At22"),  CCTK_VarIndex("CL_BSSN::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::At23"),  CCTK_VarIndex("CL_BSSN::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::At33"),  CCTK_VarIndex("CL_BSSN::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dalpha1"),  CCTK_VarIndex("CL_BSSN::dalpha1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dalpha2"),  CCTK_VarIndex("CL_BSSN::dalpha2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dalpha3"),  CCTK_VarIndex("CL_BSSN::dalpha3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dphi1"),  CCTK_VarIndex("CL_BSSN::dphi1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dphi2"),  CCTK_VarIndex("CL_BSSN::dphi2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dphi3"),  CCTK_VarIndex("CL_BSSN::dphi3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt111"),  CCTK_VarIndex("CL_BSSN::dgt111rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt112"),  CCTK_VarIndex("CL_BSSN::dgt112rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt113"),  CCTK_VarIndex("CL_BSSN::dgt113rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt122"),  CCTK_VarIndex("CL_BSSN::dgt122rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt123"),  CCTK_VarIndex("CL_BSSN::dgt123rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt133"),  CCTK_VarIndex("CL_BSSN::dgt133rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt211"),  CCTK_VarIndex("CL_BSSN::dgt211rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt212"),  CCTK_VarIndex("CL_BSSN::dgt212rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt213"),  CCTK_VarIndex("CL_BSSN::dgt213rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt222"),  CCTK_VarIndex("CL_BSSN::dgt222rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt223"),  CCTK_VarIndex("CL_BSSN::dgt223rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt233"),  CCTK_VarIndex("CL_BSSN::dgt233rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt311"),  CCTK_VarIndex("CL_BSSN::dgt311rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt312"),  CCTK_VarIndex("CL_BSSN::dgt312rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt313"),  CCTK_VarIndex("CL_BSSN::dgt313rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt322"),  CCTK_VarIndex("CL_BSSN::dgt322rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt323"),  CCTK_VarIndex("CL_BSSN::dgt323rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dgt333"),  CCTK_VarIndex("CL_BSSN::dgt333rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dbeta11"),  CCTK_VarIndex("CL_BSSN::dbeta11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dbeta12"),  CCTK_VarIndex("CL_BSSN::dbeta12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dbeta13"),  CCTK_VarIndex("CL_BSSN::dbeta13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dbeta21"),  CCTK_VarIndex("CL_BSSN::dbeta21rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dbeta22"),  CCTK_VarIndex("CL_BSSN::dbeta22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dbeta23"),  CCTK_VarIndex("CL_BSSN::dbeta23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dbeta31"),  CCTK_VarIndex("CL_BSSN::dbeta31rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dbeta32"),  CCTK_VarIndex("CL_BSSN::dbeta32rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::dbeta33"),  CCTK_VarIndex("CL_BSSN::dbeta33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::B1"),  CCTK_VarIndex("CL_BSSN::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::B2"),  CCTK_VarIndex("CL_BSSN::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::B3"),  CCTK_VarIndex("CL_BSSN::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::Xt1"),  CCTK_VarIndex("CL_BSSN::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::Xt2"),  CCTK_VarIndex("CL_BSSN::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::Xt3"),  CCTK_VarIndex("CL_BSSN::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::alpha"),  CCTK_VarIndex("CL_BSSN::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::phi"),  CCTK_VarIndex("CL_BSSN::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::gt11"),  CCTK_VarIndex("CL_BSSN::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::gt12"),  CCTK_VarIndex("CL_BSSN::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::gt13"),  CCTK_VarIndex("CL_BSSN::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::gt22"),  CCTK_VarIndex("CL_BSSN::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::gt23"),  CCTK_VarIndex("CL_BSSN::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::gt33"),  CCTK_VarIndex("CL_BSSN::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::beta1"),  CCTK_VarIndex("CL_BSSN::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::beta2"),  CCTK_VarIndex("CL_BSSN::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::beta3"),  CCTK_VarIndex("CL_BSSN::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CL_BSSN::trK"),  CCTK_VarIndex("CL_BSSN::trKrhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
