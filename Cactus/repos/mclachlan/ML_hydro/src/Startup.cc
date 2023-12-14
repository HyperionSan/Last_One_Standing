/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_hydro_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ML_hydro";
  CCTK_RegisterBanner(banner);
  return 0;
}
