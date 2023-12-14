#ifndef __utils_h
#define __utils_h

#include "cctk.h"
#include <vector>
#include <string>

using namespace std;

enum mp_coord {mp_theta, mp_phi};

namespace Multipole {
  class mode_array;
  struct variable_desc;
}

void Multipole_OutputArrayToFile(CCTK_ARGUMENTS, const string &name,
                                 vector<CCTK_REAL> const &th, vector<CCTK_REAL> const &ph,
                                 vector<CCTK_REAL> const &x, vector<CCTK_REAL> const &y,
                                 vector<CCTK_REAL> const &z, vector<CCTK_REAL> const &data);

void Multipole_Output1D(CCTK_ARGUMENTS, const string &name,
                        vector<CCTK_REAL> const &th, vector<CCTK_REAL> const &ph,
                        mp_coord coord, vector<CCTK_REAL> const &data);

void Multipole_OutputComplexToFile(CCTK_ARGUMENTS, const string &name, CCTK_REAL redata, CCTK_REAL imdata);

void Multipole_OutputComplexToH5File(CCTK_ARGUMENTS, const vector<Multipole::variable_desc> &vars,
                                     const CCTK_REAL radii[],
                                     const Multipole::mode_array& modes);

void Multipole_CoordSetup(vector<CCTK_REAL> &xhat, vector<CCTK_REAL> &yhat,
                          vector<CCTK_REAL> &zhat, vector<CCTK_REAL> &th,
                          vector<CCTK_REAL> &ph);

void Multipole_ScaleCartesian(CCTK_REAL r,
                              vector<CCTK_REAL> const &xhat, vector<CCTK_REAL> const &yhat,
                              vector<CCTK_REAL> const &zhat,
                              vector<CCTK_REAL> &x, vector<CCTK_REAL> &y,
                              vector<CCTK_REAL> &z);

static inline int Multipole_Index(int it, int ip, int ntheta)
{
  return it + (ntheta+1)*ip;
}

void Multipole_Integrate(
    vector<CCTK_REAL> const &array1r, vector<CCTK_REAL> const &array1i,
    vector<CCTK_REAL> const &array2r, vector<CCTK_REAL> const &array2i,
    vector<CCTK_REAL> const &th, vector<CCTK_REAL> const &pph,
    CCTK_REAL *out_valr, CCTK_REAL *out_vali);

#endif
