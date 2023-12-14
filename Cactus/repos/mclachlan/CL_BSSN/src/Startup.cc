/*  File produced by Kranc */

#include "cctk.h"

extern "C" int CL_BSSN_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "CL_BSSN";
  CCTK_RegisterBanner(banner);
  return 0;
}
