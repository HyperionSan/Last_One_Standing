/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_WaveToyFO_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ML_WaveToyFO";
  CCTK_RegisterBanner(banner);
  return 0;
}
