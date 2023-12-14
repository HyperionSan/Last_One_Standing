#ifndef _NUSURFACER_H_
#define _NUSURFACER_H_

#include "IsoSurfacerGH.h"
void NuFindSurface(CCTK_REAL *data, 
                   int nx,int ny,int nz,
                   CCTK_REAL *xcoords, CCTK_REAL *ycoords, CCTK_REAL *zcoords,
                   CCTK_REAL isovalue, int compute_normals, polypatch *results);

#endif
