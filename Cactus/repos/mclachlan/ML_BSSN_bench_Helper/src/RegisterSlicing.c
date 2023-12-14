#include <cctk.h>
#include <Slicing.h>

int ML_BSSN_bench_RegisterSlicing(void) {
  Einstein_RegisterSlicing("ML_BSSN_bench");
  return 0;
}
