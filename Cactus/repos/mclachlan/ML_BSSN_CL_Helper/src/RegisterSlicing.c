#include <cctk.h>
#include <Slicing.h>

int ML_BSSN_CL_RegisterSlicing(void) {
  Einstein_RegisterSlicing("ML_BSSN_CL");
  return 0;
}
