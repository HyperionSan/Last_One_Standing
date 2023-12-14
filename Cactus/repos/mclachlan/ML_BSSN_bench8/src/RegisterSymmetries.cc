/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

extern "C" void ML_BSSN_bench8_RegisterSymmetries(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_BSSN_bench8_RegisterSymmetries
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_bench8_RegisterSymmetries);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  /* array holding symmetry definitions */
  int sym[3];
  
  /* Register symmetries of grid functions */
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::phi");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::gt11");
  
  sym[0] = -1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::gt12");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::gt13");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::gt22");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::gt23");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::gt33");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::Xt1");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::Xt2");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::Xt3");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::trK");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::At11");
  
  sym[0] = -1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::At12");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::At13");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::At22");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::At23");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::At33");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::alpha");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::beta1");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::beta2");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::beta3");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::H");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::M1");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::M2");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::M3");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::cS");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::cXt1");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::cXt2");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::cXt3");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_bench8::cA");
  
}
