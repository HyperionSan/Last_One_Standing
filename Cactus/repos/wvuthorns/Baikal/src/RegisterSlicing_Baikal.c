#include "cctk.h"
#include "Slicing.h"
/*
 * Register slicing condition for NRPy+-generated thorn Baikal
 */
int RegisterSlicing_Baikal() {

  Einstein_RegisterSlicing ("Baikal");
  return 0;
}
