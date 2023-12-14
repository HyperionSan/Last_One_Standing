#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Functions.h>
#include <cctk_Parameters.h>

#include "util_Table.h"

#include <assert.h>
#include <stdlib.h>

#define BOUNDARY_WIDTH 1

void TestAutoSync_RegisterBCs (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  int ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                                   "TestAutoSync::Var1", "scalar");
  assert(!ierr);

  if(CCTK_EQUALS(presync_mode, "presync-only")) {
    int table_handle = Util_TableCreateFromString("scalar=1.0");
    assert(table_handle >= 0);
    ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, table_handle,
                                 "TestAutoSync::Var2", "scalar");
    assert(!ierr);
  }
}

void TestAutoSync_Init (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_TestAutoSync_Init;

  const CCTK_REAL nan = atof("nan");

  printf("ox: %g\n", cctk_origin_space[0]);
  #pragma omp parallel
  CCTK_LOOP3_ALL(TestAutoSync_Init, cctkGH, i,j,k) {
    const int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
    Var1[idx] = (i-BOUNDARY_WIDTH) * CCTK_DELTA_SPACE(0) * 100 +
                j * CCTK_DELTA_SPACE(1) * 10 * 0+
                k * CCTK_DELTA_SPACE(2) * 1 * 0;
    Var2[idx] = -10;
    TransferVar[idx] = nan;
  } CCTK_ENDLOOP3_ALL(TestAutoSync_Init);
}

// expected values are 100*i+10*j+k + cctk_iteration in the interior and 0 in
// the boundary
void TestAutoSync_Evolve1 (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_TestAutoSync_Evolve1;

  #pragma omp parallel
  CCTK_LOOP3_INT(TestAutoSync_Evolve1, cctkGH, i,j,k) {
    const int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
    Var1[idx] = Var1_p[idx] + 1.;
  } CCTK_ENDLOOP3_INT(TestAutoSync_Evolve1);
}

void TestAutoSync_Transfer (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_TestAutoSync_Transfer;

  #pragma omp parallel
  CCTK_LOOP3_ALL(TestAutoSync_Transfer, cctkGH, i,j,k) {
    const int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
    TransferVar[idx] = Var1[idx];
  } CCTK_ENDLOOP3_ALL(TestAutoSync_Transfer);
}

// expected values are 1000. + Var1 in the interior and 1 in the boundary
void TestAutoSync_Evolve2 (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  #pragma omp parallel
  CCTK_LOOP3_INT(TestAutoSync_Evolve2, cctkGH, i,j,k) {
    const int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
    Var2[idx] = TransferVar[idx] + 1000.;
  } CCTK_ENDLOOP3_INT(TestAutoSync_Evolve2);
}

void TestAutoSync_SelectBCsForVar1 (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  int ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, BOUNDARY_WIDTH, -1,
                                     "TestAutoSync::Var1", "scalar");
  assert(!ierr);
}

void TestAutoSync_SelectBCsForVar2 (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  static int table_handle = -1;
  if(table_handle == -1)
    table_handle = Util_TableCreateFromString("scalar=1.0");
  assert(table_handle >= 0);
  int ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, BOUNDARY_WIDTH, table_handle,
                                     "TestAutoSync::Var2", "scalar");
  assert(!ierr);
}

void TestAutoSync_Dummy (CCTK_ARGUMENTS)
{
  return;
}
