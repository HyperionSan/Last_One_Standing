#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Initial guess at a single point.
 */
void initial_guess_single_point(const paramstruct *restrict params,
                const REAL xx0, const REAL xx1, const REAL xx2,
                REAL *uu_exact, REAL *vv_exact) {
#include "./set_Cparameters.h"

  /*
   * NRPy+ Finite Difference Code Generation, Step 1 of 1: Evaluate SymPy expressions and write to main memory:
   */
  *uu_exact = 0;
  *vv_exact = 0;
}
