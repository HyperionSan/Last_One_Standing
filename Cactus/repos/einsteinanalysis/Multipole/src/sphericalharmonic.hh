#ifndef sphericalharmonic_h
#define sphericalharmonic_h
#include "cctk.h"

#include <vector>

using namespace std;

void Multipole_HarmonicSetup(int s, int l, int m,
                             vector<CCTK_REAL> const &th, vector<CCTK_REAL> const &ph,
                             vector<CCTK_REAL> &reY, vector<CCTK_REAL> &imY);


void Multipole_SphericalHarmonic(int s, int l, int m,
                                 CCTK_REAL th, CCTK_REAL ph,
                                 CCTK_REAL *reY, CCTK_REAL *imY);
#endif
