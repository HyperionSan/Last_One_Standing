
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

namespace HDecomp {

using namespace SPS;


//void collective_sync_1patch(const cGH* const cctkGH, const int si);
//void collective_sync_6patch(const cGH* const cctkGH, const int si);


extern "C" void HarmonicDecomposition_Decompose(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   
   int mode_dim = (lmax-lmin+1)*(lmax-lmin+1);
   
   
   for (int si=0; si < nslices; ++si)
   {
      
      // first of all we do a collective interpolation onto the sphere
      SphericalSlice_CollectiveSyncVariables(cctkGH, &sid[si].front(), number_of_vars);
      /*if (is_1patch(sid[si][0]))
         collective_sync_1patch(cctkGH, si);
   
      // first of all we do a collective interpolation onto the sphere
      if (is_6patch(sid[si][0]))
         collective_sync_6patch(cctkGH, si);
      */
      for (int j=0; j < number_of_vars; ++j)
      {
         //if (!precalc_sYlms || spin_weight[j] != 0) //spin_weight[j] != spin_weight[0])
         {
            SphericalSlice_ContractVariableWithAllsYlm(sid[si][j], 0, spin_weight[j], lmin, lmax, &decomp_vars[mode_dim*(si + nslices*j)]);
         }
         /*else
         {
            int s = spin_weight[0];  // we only precalculate sYlms for variables with the same spin-weight as the first variable!
         
            if (is_1patch(sid[si][j]))
            {
               assert(INDEX1P(sid[si][j]) < slices_1patch.slice().size());
               
               if (!slices_1patch(INDEX1P(sid[si][j]), 0).has_constant_radius() && slices_1patch(INDEX1P(sid[si][j]), 0).nghosts() < 2) 
                  CCTK_WARN(0, "ghostzone-width too small!");
                           
               vector<integrator_1patch> intsYlm_re = vector<integrator_1patch>(mode_dim, slices_1patch(INDEX1P(sid[si][j]), 0));
               vector<integrator_1patch> intsYlm_im = vector<integrator_1patch>(mode_dim, slices_1patch(INDEX1P(sid[si][j]), 0));
               
               for (const_iter_1patch i=slices_1patch(INDEX1P(sid[si][j]), 0).begin(); !i.done(); ++i)
               {
                  // dont't integrate in ghostzones
                  if (i.ghostzone())
                     continue;
                  
                  // get local starting index on this processor
                  vect<int,2> lbnd = slices_1patch(INDEX1P(sid[si][j]), 0).lbnd(i.idx().p);
                  // get global number of gridpoints for this slice
                  vect<int,2> gsh = slices_1patch(INDEX1P(sid[si][j]), 0).gsh(i.idx().p);
                  
                  CCTK_REAL r = slices_1patch(INDEX1P(sid[si][j]), 0).radius(i);
                  CCTK_REAL det = slices_1patch(INDEX1P(sid[si][j]), 0).det(i) / (r*r);
                  
                  int lm = 0;
                  int lm2 = 0;
                  for (int l=lmin; l <= lmax; ++l)
                  {
                     int save_lm2 = lm2;
                     lm2 += l;
                     for (int m=-l; m < 0; ++m)
                     {
                        if (l >= abs(spin_weight[j]))
                        {
                           CCTK_REAL val_sYlm_re = slices_1patch(INDEX1P(sid_sYlm_re[si][lm2]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                           CCTK_REAL val_sYlm_im = -slices_1patch(INDEX1P(sid_sYlm_im[si][lm2]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                           intsYlm_re[lm].sum(i, det, pow(-1.0, m)*val_sYlm_re);
                           intsYlm_im[lm].sum(i, det, pow(-1.0, m)*val_sYlm_im);
                        }
                        lm++;
                        lm2--;
                     }
                     lm2 = save_lm2;
                     for (int m=0; m <= l; ++m)
                     {
                        if (l >= abs(spin_weight[j]))
                        {
                           CCTK_REAL val_sYlm_re = slices_1patch(INDEX1P(sid_sYlm_re[si][lm2]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                           CCTK_REAL val_sYlm_im = slices_1patch(INDEX1P(sid_sYlm_im[si][lm2]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                           intsYlm_re[lm].sum(i, det, val_sYlm_re);
                           intsYlm_im[lm].sum(i, det, val_sYlm_im);
                        }
                        lm++;
                        lm2++;
                     }
                  }
               }
               
               // create a collective commstack in order to do only one MPI_Allreduce call instead of 2*l^2
               
               commstack cs(MPI_SUM, MPI_COMM_WORLD);
               
               int lm = 0;
               for (int l=lmin; l <= lmax; ++l)
               {
                  for (int m=-l; m <= l; ++m)
                  {
                     if (l < abs(s))
                        decomp_vars[lm + mode_dim*(si + nslices*j)] = CCTK_Cmplx(0.0, 0.0);
                     else
                     {
                        //decomp_vars[lm + mode_dim*(si + nslices*j)] = CCTK_Cmplx(intsYlm_re[lm].finalize(), -intsYlm_im[lm].finalize());
                        intsYlm_re[lm].finalize(&cs); 
                        intsYlm_im[lm].finalize(&cs);
                     }
                     lm++;
                  }
               }
               
               // collective reduction
               cs.reduce();
               
               // write back results
               lm = 0;
               for (int l=lmin; l <= lmax; ++l)
               {
                  for (int m=-l; m <= l; ++m)
                  {
                     if (l < abs(s))
                        decomp_vars[lm + mode_dim*(si + nslices*j)] = CCTK_Cmplx(0.0, 0.0);
                     else
                     {
                        //decomp_vars[lm + mode_dim*(si + nslices*j)] = CCTK_Cmplx(intsYlm_re[lm].finalize(), -intsYlm_im[lm].finalize());
                        decomp_vars[lm + mode_dim*(si + nslices*j)] = CCTK_Cmplx(cs.reduced_val(2*(lm-s*s)), -cs.reduced_val(2*(lm-s*s)+1));
                     }
                     lm++;
                  }
               }
               
               
            }
            
            
            if (is_6patch(sid[si][j]))
            {
               assert(INDEX6P(sid[si][j]) < slices_6patch.slice().size());
               
               if (!slices_6patch(INDEX6P(sid[si][j]), 0).has_constant_radius() && slices_6patch(INDEX6P(sid[si][j]), 0).nghosts() < 2) 
                  CCTK_WARN(0, "ghostzone-width too small!");
                           
               vector<integrator_6patch> intsYlm_re = vector<integrator_6patch>(mode_dim, slices_6patch(INDEX6P(sid[si][j]), 0));
               vector<integrator_6patch> intsYlm_im = vector<integrator_6patch>(mode_dim, slices_6patch(INDEX6P(sid[si][j]), 0));
               
               for (const_iter_6patch i=slices_6patch(INDEX6P(sid[si][j]), 0).begin(); !i.done(); ++i)
               {
                  // dont't integrate in ghostzones
                  if (i.ghostzone())
                     continue;
                  
                  // get local starting index on this processor
                  vect<int,2> lbnd = slices_6patch(INDEX6P(sid[si][j]), 0).lbnd(i.idx().p);
                  
                  CCTK_REAL r = slices_6patch(INDEX6P(sid[si][j]), 0).radius(i);
                  CCTK_REAL det = slices_6patch(INDEX6P(sid[si][j]), 0).det(i) / (r*r);
                  
                  int lm = 0;
                  int lm2 = 0;
                  for (int l=lmin; l <= lmax; ++l)
                  {
                     int save_lm2 = lm2;
                     lm2 += l;
                     for (int m=-l; m < 0; ++m)
                     {
                        if (l >= abs(spin_weight[j]))
                        {
                           CCTK_REAL val_sYlm_re = slices_6patch(INDEX6P(sid_sYlm_re[si][lm2]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                           CCTK_REAL val_sYlm_im = -slices_6patch(INDEX6P(sid_sYlm_im[si][lm2]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                           intsYlm_re[lm].sum(i, det, pow(-1.0, m)*val_sYlm_re);
                           intsYlm_im[lm].sum(i, det, pow(-1.0, m)*val_sYlm_im);
                        }
                        lm++;
                        lm2--;
                     }
                     lm2 = save_lm2;
                     for (int m=0; m <= l; ++m)
                     {
                        if (l >= abs(spin_weight[j]))
                        {
                           CCTK_REAL val_sYlm_re = slices_6patch(INDEX6P(sid_sYlm_re[si][lm2]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                           CCTK_REAL val_sYlm_im = slices_6patch(INDEX6P(sid_sYlm_im[si][lm2]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                           intsYlm_re[lm].sum(i, det, val_sYlm_re);
                           intsYlm_im[lm].sum(i, det, val_sYlm_im);
                        }
                        lm++;
                        lm2++;
                     }
                  }
               }
               
               
               // create a collective commstack in order to do only one MPI_Allreduce call instead of 2*l^2
               
               commstack cs(MPI_SUM, MPI_COMM_WORLD);
               
               int lm = 0;
               for (int l=lmin; l <= lmax; ++l)
               {
                  for (int m=-l; m <= l; ++m)
                  {
                     if (l < abs(s))
                        decomp_vars[lm + mode_dim*(si + nslices*j)] = CCTK_Cmplx(0.0, 0.0);
                     else
                     {
                        //decomp_vars[lm + mode_dim*(si + nslices*j)] = CCTK_Cmplx(intsYlm_re[lm].finalize(), -intsYlm_im[lm].finalize());
                        intsYlm_re[lm].finalize(&cs); 
                        intsYlm_im[lm].finalize(&cs);
                     }
                     lm++;
                  }
               }
               
               // collective reduction
               cs.reduce();
               
               // write back results
               lm = 0;
               for (int l=lmin; l <= lmax; ++l)
               {
                  for (int m=-l; m <= l; ++m)
                  {
                     if (l < abs(s))
                        decomp_vars[lm + mode_dim*(si + nslices*j)] = CCTK_Cmplx(0.0, 0.0);
                     else
                     {
                        //decomp_vars[lm + mode_dim*(si + nslices*j)] = CCTK_Cmplx(intsYlm_re[lm].finalize(), -intsYlm_im[lm].finalize());
                        decomp_vars[lm + mode_dim*(si + nslices*j)] = CCTK_Cmplx(cs.reduced_val(2*(lm-s*s)), -cs.reduced_val(2*(lm-s*s)+1));
                     }
                     lm++;
                  }
               }
            }*/
            
         /*}*/
      }
   }
}


/*
/// interpolate from Cactus gridfunctions onto sphere
void collective_sync_1patch(const cGH* const cctkGH, const int si) 
{
   DECLARE_CCTK_PARAMETERS

   // get the input variable's index
   vector<int> varindices (number_of_vars, -1);
   
   for (int i=0; i < number_of_vars; ++i)
   {
      varindices[i] = CCTK_VarIndex(slices_1patch(INDEX1P(sid[si][i]), 0).varname().c_str());
      if (varindices[i] < 0)
         CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "couldn't get index of slice variable '%s'", slices_1patch(INDEX1P(sid[si][i]), 0).varname().c_str());
   }
   
   vector<CCTK_REAL*> values (number_of_vars, static_cast<CCTK_REAL*>(NULL));
   
   for (int i=0; i < number_of_vars; ++i)
      values[i] = (CCTK_REAL*) slices_1patch(INDEX1P(sid[si][i]), 0).data_pointer();
   
   // set up the interpolation
   assert(interp_setups.size() > slices_1patch(INDEX1P(sid[si][0]), 0).ID());
   interp_setup_t* &interp_setup = interp_setups[slices_1patch(INDEX1P(sid[si][0]), 0).ID()];
   if (not interp_setup) {
      interp_setup = new interp_setup_t(slices_1patch(INDEX1P(sid[si][0]), 0).lsh(0)[0] * slices_1patch(INDEX1P(sid[si][0]), 0).lsh(0)[1]);

      // allocate storage for coordinates
      fasterp_glocs_t locations (interp_setup->npoints);

      // get Cartesian coordinate values of gridpoints on spherical surface
      for (int j=0; j < slices_1patch(INDEX1P(sid[si][0]), 0).lsh(0)[1]; ++j) {
      for (int k=0; k < slices_1patch(INDEX1P(sid[si][0]), 0).lsh(0)[0]; ++k) {
         const int l = k + slices_1patch(INDEX1P(sid[si][0]), 0).lsh(0)[0]*j;
         locations.coords[0][l] = slices_1patch(INDEX1P(sid[si][0]), 0).cart_x(0, k, j);
         locations.coords[1][l] = slices_1patch(INDEX1P(sid[si][0]), 0).cart_y(0, k, j);
         locations.coords[2][l] = slices_1patch(INDEX1P(sid[si][0]), 0).cart_z(0, k, j);
      }
      }

      interp_setup->fasterp_setup =
      new fasterp_setup_t(cctkGH, locations, interpolator_order);
   }

   // do the interpolation
   assert(interp_setup->fasterp_setup);
   assert(interp_setup->npoints == slices_1patch(INDEX1P(sid[si][0]), 0).lsh(0)[0] * slices_1patch(INDEX1P(sid[si][0]), 0).lsh(0)[1]);
   interp_setup->fasterp_setup->interpolate (cctkGH, varindices, values);

   if (not slices_1patch(INDEX1P(sid[si][0]), 0).has_constant_radius()) {
      delete interp_setup;
      interp_setup = NULL;
   }

}










void collective_sync_6patch(const cGH* const cctkGH, const int si) 
{
   DECLARE_CCTK_PARAMETERS

   const void*   interp_coords[N_DIMS];
   vector<int>   input_array_indices(number_of_vars);
   vector<void*> output_arrays(number_of_vars);
   const vector<CCTK_INT> input_array_type_codes(number_of_vars, CCTK_VARIABLE_REAL);
            //= { CCTK_VARIABLE_REAL };
   const vector<CCTK_INT> output_array_type_codes(number_of_vars, CCTK_VARIABLE_REAL);
            //= { CCTK_VARIABLE_REAL };
   
   // set interpolation handles
   int operator_handle = CCTK_InterpHandle("Lagrange polynomial interpolation");
   if (operator_handle < 0)
      CCTK_WARN(0, "can’t get interpolation handle!");
   int param_table_handle = Util_TableCreateFromString("order=4");
   if (param_table_handle < 0)
      CCTK_WARN(0, "can’t create parameter table!");
   int coordsys_handle = CCTK_CoordSystemHandle("cart3d");
   if (coordsys_handle < 0)
      CCTK_WARN(0, "can’t create coordsys handle! Forgot to activate a Coordinate-Thorn?");
   
   for (int i=0; i < number_of_vars; ++i)
   {
      input_array_indices[i] = CCTK_VarIndex(slices_6patch(INDEX6P(sid[si][i]), 0).varname().c_str());
      if (input_array_indices[i] < 0)
         CCTK_WARN(0, "error getting VarIndex of variable that shall be sliced");
   }
   
   for (int m=0; m < 6; ++m)
   {
      // project variable onto sphere
      // this means we have to 3d-interpolate every non-coinciding gridpoint of
      // the spherical grid from neighboring gridpoints of the Cartesian grid.
      // for this, we first calculate the Cartesian coordinate of the jk-th spherical
      // surface gridpoint
      int n_interp_points = slices_6patch(INDEX6P(sid[si][0]), 0).lsh(m)[0]
                           *slices_6patch(INDEX6P(sid[si][0]), 0).lsh(m)[1];
                           //lsh(m)[0]*lsh(m)[1];
      
      // Arrays of Cartesian coordinates of the surface points onto which we want to interpolate
      vector<CCTK_REAL> interp_x(n_interp_points, POISON_VAL);
      vector<CCTK_REAL> interp_y(n_interp_points, POISON_VAL);
      vector<CCTK_REAL> interp_z(n_interp_points, POISON_VAL);
      
      // get Cartesian coordinate values of gridpoints on spherical surface
      for (int j=0; j < slices_6patch(INDEX6P(sid[si][0]), 0).lsh(0)[1]; ++j)
         for (int k=0; k < slices_6patch(INDEX6P(sid[si][0]), 0).lsh(0)[0]; ++k)
         {
            const int l = k + slices_6patch(INDEX6P(sid[si][0]), 0).lsh(m)[0]*j;
            interp_x[l] = slices_6patch(INDEX6P(sid[si][0]), 0).cart_x(m, k, j);
            interp_y[l] = slices_6patch(INDEX6P(sid[si][0]), 0).cart_y(m, k, j);
            interp_z[l] = slices_6patch(INDEX6P(sid[si][0]), 0).cart_z(m, k, j);
         }
      
      interp_coords[0] = (const void*) &interp_x.front();
      interp_coords[1] = (const void*) &interp_y.front();
      interp_coords[2] = (const void*) &interp_z.front();
      
      for (int i=0; i < number_of_vars; ++i)
         output_arrays[i] = slices_6patch(INDEX6P(sid[si][i]), 0).data_pointer(m);//(void*) &((vector<vector<CCTK_REAL> >*) slices_6patch(INDEX6P(sid[si][i]), 0).data_pointer())[m].front();
      //output_arrays[0] = (void*) &data[m].front();

      // Do the actual interpolation.
      // Only those processes interpolate that contain data
      // All other processes simply send their
      // grid-data.
      //if (data[m].size() == 0)
      //   n_interp_points = 0;      // all other processes shall not interpolate
      
      if (CCTK_InterpGridArrays(cctkGH,
                                 N_DIMS,
                                 operator_handle, param_table_handle,
                                 coordsys_handle,
                                 n_interp_points,
                                    CCTK_VARIABLE_REAL,
                                    interp_coords,
                                 number_of_vars,
                                    &input_array_indices.front(),
                                 number_of_vars,
                                    &output_array_type_codes.front(),
                                    &output_arrays.front()) < 0)
         CCTK_WARN(1, "error return from interpolator!");
         
   }
}

*/







}
