#include "cctk.h"
#include "cctk_Arguments.h"

#include <vector>

using namespace std;

// Multipole_Interp:
//      This function interpolates psi4 onto the sphere in cartesian 
//      coordinates as created by Multipole_CoordSetup.
void Multipole_Interp(CCTK_ARGUMENTS,
                      vector<CCTK_REAL> const &x, vector<CCTK_REAL> const &y,
                      vector<CCTK_REAL> const &z,
                      int real_idx, int imag_idx,
                      vector<CCTK_REAL> &psi4r, vector<CCTK_REAL> &psi4i);

void Multipole_InterpVar(CCTK_ARGUMENTS,
                         vector<CCTK_REAL> const &x, vector<CCTK_REAL> const &y,
                         vector<CCTK_REAL> const &z, const char *var_name,
                         vector<CCTK_REAL> &sphere_var);
