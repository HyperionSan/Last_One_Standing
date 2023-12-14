/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_BSSN_NV_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ML_BSSN_NV";
  CCTK_RegisterBanner(banner);
  return 0;
}
