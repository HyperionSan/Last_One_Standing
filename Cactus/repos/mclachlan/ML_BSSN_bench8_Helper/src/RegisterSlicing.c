#include <cctk.h>
#include <Slicing.h>

int ML_BSSN_bench8_RegisterSlicing(void) {
  Einstein_RegisterSlicing("ML_BSSN_bench8");
  return 0;
}
