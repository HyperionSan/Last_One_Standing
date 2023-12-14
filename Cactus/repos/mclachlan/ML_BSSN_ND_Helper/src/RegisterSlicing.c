#include <cctk.h>
#include <Slicing.h>

int ML_BSSN_ND_RegisterSlicing(void) {
  Einstein_RegisterSlicing("ML_BSSN_ND");
  return 0;
}
