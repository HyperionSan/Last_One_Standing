
/* Copyright 2013 Peter Diener, Nils Dorband, Roland Haas, Ian Hinder,
Christian Ott, Denis Pollney, Thomas Radke, Christian Reisswig, Erik
Schnetter, Barry Wardell and Burkhard Zink

This file is part of Llama.

Llama is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 2 of the License, or (at your
option) any later version.

Llama is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with Llama.  If not, see <http://www.gnu.org/licenses/>. */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include <assert.h>
#include <iostream>

#include "register.h"
#include "spherical_slices.hh"



using namespace WorldTube;
using namespace SPS;

double norm_inf(double const * const f, int const ni, int const nj);
int get_1patch_metric_component(double g_sph[4][4], double const r, double const theta, double const phi, 
                               const double lapse, const double shiftx, const double shifty, const double shiftz,
                               const double gxx, const double gxy, const double gxz,
                               const double gyy, const double gyz,
                               const double gzz);
double calc_R(double const target_radius, 
              double const rp, double const tp, double const pp, const double g_sph[4][4], 
              double const diff_theta, double const diff_phi);

extern "C" void WorldTube_ArealSphere(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   
   int mode_dim = (lmax+1)*(lmax+1);
   
   
   for (int ii=0; ii < ntubes; ++ii)
   {
      const int si = which_slice_to_take[ii];
   
      if (seek_const_areal_radius[ii])
      {
         int it = 0;
         double res_Linf = 2.0*tolerance;
         while ((res_Linf > tolerance) && (it < max_it))
         {
            // first of all we do a collective interpolation for the current iterated radius onto the sphere
            SphericalSlice_CollectiveSyncVariables(cctkGH, &sid[si].front(), 10);
            
            
            if (CCTK_Equals(type[si], "1patch"))
            {
               vector<CCTK_REAL> R(slices_1patch(INDEX1P(ss_radius_id[si]), 0).lsh(0)[0] * slices_1patch(INDEX1P(ss_radius_id[si]), 0).lsh(0)[1], 0);
               vector<CCTK_REAL> res(slices_1patch(INDEX1P(ss_radius_id[si]), 0).lsh(0)[0] * slices_1patch(INDEX1P(ss_radius_id[si]), 0).lsh(0)[1], 0);
               
               for (spheredata_1patch<CCTK_REAL>::const_iter i=slices_1patch(INDEX1P(ss_radius_id[si]), 0).begin(); !i.done(); ++i)
               {
                  // dont't integrate in ghostzones
                  if (i.ghostzone())
                     continue;
                  
                  int const idx = i.idx().ij;
                  
                  // get metric components at current point in polar-spherical basis
                  CCTK_REAL g_sph[4][4];
                  const CCTK_REAL lapse  = slices_1patch(INDEX1P(sid[si][0]), 0)(i);
                  const CCTK_REAL shiftx = slices_1patch(INDEX1P(sid[si][1]), 0)(i);
                  const CCTK_REAL shifty = slices_1patch(INDEX1P(sid[si][2]), 0)(i);
                  const CCTK_REAL shiftz = slices_1patch(INDEX1P(sid[si][3]), 0)(i);
                  const CCTK_REAL gxx    = slices_1patch(INDEX1P(sid[si][4]), 0)(i);
                  const CCTK_REAL gxy    = slices_1patch(INDEX1P(sid[si][5]), 0)(i);
                  const CCTK_REAL gxz    = slices_1patch(INDEX1P(sid[si][6]), 0)(i);
                  const CCTK_REAL gyy    = slices_1patch(INDEX1P(sid[si][7]), 0)(i);
                  const CCTK_REAL gyz    = slices_1patch(INDEX1P(sid[si][8]), 0)(i);
                  const CCTK_REAL gzz    = slices_1patch(INDEX1P(sid[si][9]), 0)(i);
                  get_1patch_metric_component(g_sph, *i, i.idx().theta, i.idx().phi, 
                                              lapse,
                                              shiftx, shifty, shiftz,
                                              gxx, gxy, gxz,
                                              gyy, gyz,
                                              gzz);
                  
                  // get angular derivatives of radius
                  CCTK_REAL diff_theta = slices_1patch(INDEX1P(ss_radius_id[si]), 0).dx(i);
                  CCTK_REAL diff_phi   = slices_1patch(INDEX1P(ss_radius_id[si]), 0).dy(i);
                  
                  R[idx] = calc_R(radius[si], *i, i.idx().theta, i.idx().phi, g_sph, diff_theta, diff_phi);
                  res[idx] = fabs(R[idx] - r[idx]);
               }
               
               res_Linf = norm_inf(&res.front(), slices_1patch(INDEX1P(ss_radius_id[si]), 0).lsh(0)[0], 
                                                 slices_1patch(INDEX1P(ss_radius_id[si]), 0).lsh(0)[1]);
               
               for (spheredata_1patch<CCTK_REAL>::iter i=slices_1patch(INDEX1P(ss_radius_id[si]), 0).begin(); !i.done(); ++i)
                  *i = R[i.idx().ij];
               
               ++it; 
            }
            if (CCTK_Equals(type[si], "6patch"))
               CCTK_WARN(0, "Sorry: Finding surface of constant areal radius not implemented for 6-patches yet :(");
            
            
         }
      }
      
   }
}



/*
 *  Infinity norm over the slice.
 */
double
norm_inf(double const * const f, int const ni, int const nj)
{
  double norm = f[0];
  for (int idx=0; idx<ni*nj; ++idx)
    if (fabs(f[idx]) > norm)
      norm = fabs(f[idx]);
  double global_norm;
  MPI_Allreduce (&norm, &global_norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
  return global_norm;
}


/*
 *  Get the metric at a point on the sphere (1-patch coordinate basis).
 *
 *  The metric components should be given in cartesian components to
 *  fill the g_cart[] array. This is transformed to spherical polar
 *  coordinates g_sph[] using the jacobian J[].
 *
 *  In practice, the jacobian J[] will also need to recognise the possibility
 *  of alternate patch systems that are in use, eg. thornburg04 style 6-patch
 *  systems.
 *
 *  The metric components should be obtained by interpolation from the
 *  numerically generated spacetime.
 *
 */
int
get_1patch_metric_component(double g_sph[4][4], double const r, double const theta, double const phi, 
                            const double lapse, const double shiftx, const double shifty, const double shiftz,
                            const double gxx, const double gxy, const double gxz,
                            const double gyy, const double gyz,
                            const double gzz)
{
  double const J[4][4] = {
    { 1,                   0,                     0,                      0 },
    { 0, sin(theta)*cos(phi), r*cos(theta)*cos(phi), -r*sin(theta)*sin(phi) },
    { 0, sin(theta)*sin(phi), r*cos(theta)*sin(phi),  r*sin(theta)*cos(phi) },
    { 0,          cos(theta),         -r*sin(theta),                      0 }
  };

  double const dt2 =   -lapse*lapse
                   +   gxx*shiftx*shiftx
                   + 2*gxy*shiftx*shifty
                   + 2*gxz*shiftx*shiftz
                   +   gyy*shifty*shifty
                   + 2*gyz*shifty*shiftz
                   +   gzz*shiftz*shiftz;

  double const dtdx = gxx*shiftx+gxy*shifty+gxz*shiftz;
  double const dtdy = gxy*shiftx+gyy*shifty+gyz*shiftz;
  double const dtdz = gxz*shiftx+gyz*shifty+gzz*shiftz;

  double const g_cart[4][4] = {
    {  dt2, dtdx, dtdy, dtdz },
    { dtdx,  gxx,  gxy,  gxz },
    { dtdy,  gxy,  gyy,  gyz },
    { dtdz,  gxz,  gyz,  gzz } 
  };

  for (int a=0; a<4; ++a)
    for (int b=0; b<4; ++b)
      {
        g_sph[a][b] = 0.0;
        for (int ix=1; ix<4; ++ix)
          for (int jx=1; jx<4; ++jx)
            g_sph[a][b] += g_cart[ix][jx]*J[ix][a]*J[jx][b];
      }


  return 0;
}



/*
 *  Calculate the updated radius, R, at a point given an initial guess
 *  for the radius, r, and the metric on the sphere. See Nigel's notes
 *  for details on the computation.
 */
double
calc_R(double const target_radius, 
       double const rp, double const tp, double const pp, const double g_sph[4][4], 
       double const diff_theta, double const diff_phi)
{
  double const grr = g_sph[1][1];
  double const hr2 = g_sph[1][2] / rp;
  double const hr3 = g_sph[1][3] / rp;
  double const h22 = g_sph[2][2] / (rp*rp);
  double const h23 = g_sph[2][3] / (rp*rp);
  double const h33 = g_sph[3][3] / (rp*rp);
  
  double const log_r_theta = diff_theta / rp;
  double const log_r_phi = diff_phi / rp;

  double const cg = sin(tp)*sin(tp) /
    (h22*h33 - h23 * h23
     + 2.0 * log_r_theta * hr2 * log_r_phi * hr3
     + 2.0 * log_r_theta * hr2 * h33
     + h22 * log_r_phi * log_r_phi * grr
     + log_r_theta * log_r_theta * grr * h33
     - 2.0 * log_r_theta * log_r_phi * grr * h23
     - log_r_phi * log_r_phi * hr2 * hr2
     + 2.0 * h22 * log_r_phi * hr3
     - log_r_theta * log_r_theta * hr3 * hr3
     - 2.0 * log_r_theta * hr3 * h23
     - 2.0 * log_r_phi * hr2 * h23);

  return  target_radius * pow(cg, 0.25);
}





