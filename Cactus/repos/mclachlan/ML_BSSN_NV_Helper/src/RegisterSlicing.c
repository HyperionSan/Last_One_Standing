#include <cctk.h>
#include <Slicing.h>

int ML_BSSN_NV_RegisterSlicing(void) {
  Einstein_RegisterSlicing("ML_BSSN_NV");
  return 0;
}
