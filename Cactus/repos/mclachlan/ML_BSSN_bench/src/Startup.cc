/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_BSSN_bench_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ML_BSSN_bench";
  CCTK_RegisterBanner(banner);
  return 0;
}
