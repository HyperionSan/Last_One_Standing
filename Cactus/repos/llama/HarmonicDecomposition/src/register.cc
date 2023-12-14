
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

#include "register.h"

#include "spherical_slices.hh" 


namespace HDecomp {

using namespace SPS;

// the slice ids for each slice and variable
vector<vector<int> > sid;

// the slice ids for each slice and precalculated sYlm
//vector<vector<int> > sid_sYlm_re;
//vector<vector<int> > sid_sYlm_im;

extern "C" void HarmonicDecomposition_Register(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS

   sid = vector<vector<int> >(nslices);
   
   for (int i=0; i < nslices; ++i)
   {
      sid[i] = vector<int>(number_of_vars, 0);
      for (int j=0; j < number_of_vars; ++j)
      {
         sid[i][j] = SphericalSlice_RegisterVariable(vars[j], which_slice_to_take[i], 1, "split");
      }
   }
   
   /*
   // precalculate and store sYlm's
   if (precalc_sYlms)
   {
      int s = spin_weight[0];  // we only precalculate sYlms for variables with the same spin-weight as the first variable!
      
      int mode_dim = (lmax-lmin+1)*(lmax-lmin+1);
      
      sid_sYlm_re = vector<vector<int> >(nslices);
      sid_sYlm_im = vector<vector<int> >(nslices);
      
      for (int si=0; si < nslices; ++si)
      {
         sid_sYlm_re[si] = vector<int>(mode_dim, 0);
         sid_sYlm_im[si] = vector<int>(mode_dim, 0);
      
         int lm = 0;
         for (int l=lmin; l <= lmax; ++l)
            for (int m=0; m <= l; ++m)  // we only need to store +m modes because -m can be recovered from them via CC!
            {
               // register
               sid_sYlm_re[si][lm] = SphericalSlice_RegisterVariable("sYlm_re", which_slice_to_take[si], 1, "const");
               sid_sYlm_im[si][lm] = SphericalSlice_RegisterVariable("sYlm_im", which_slice_to_take[si], 1, "const");
               
               // precalculate
               if (is_1patch(sid_sYlm_re[si][lm]))
               {
                  
                  for (iter_1patch i=slices_1patch(INDEX1P(sid_sYlm_re[si][lm]), 0).begin(); !i.done(); ++i)
                  {
                     double sYlm_re, sYlm_im;
	             sYlm(s,l,m, i.idx().theta, i.idx().phi,
			   &sYlm_re, &sYlm_im);
                     *i = sYlm_re;
                  }
                  for (iter_1patch i=slices_1patch(INDEX1P(sid_sYlm_im[si][lm]), 0).begin(); !i.done(); ++i)
                  {
                     double sYlm_re, sYlm_im;
	             sYlm(s,l,m, i.idx().theta, i.idx().phi,
			   &sYlm_re, &sYlm_im);
                     *i = sYlm_im;
                  }
               }
               
               if (is_6patch(sid_sYlm_re[si][lm]))
               {
                  
                  for (iter_6patch i=slices_6patch(INDEX6P(sid_sYlm_re[si][lm]), 0).begin(); !i.done(); ++i)
                  {
                     double sYlm_re, sYlm_im;
	             sYlm(s,l,m, slices_6patch(INDEX6P(sid_sYlm_re[si][lm]), 0).coord_spherical(i)[0],
                                 slices_6patch(INDEX6P(sid_sYlm_re[si][lm]), 0).coord_spherical(i)[1],
			         &sYlm_re, &sYlm_im);
                     *i = sYlm_re;
                  }
                  for (iter_6patch i=slices_6patch(INDEX6P(sid_sYlm_im[si][lm]), 0).begin(); !i.done(); ++i)
                  {
                     double sYlm_re, sYlm_im;
	             sYlm(s,l,m, slices_6patch(INDEX6P(sid_sYlm_re[si][lm]), 0).coord_spherical(i)[0],
                                 slices_6patch(INDEX6P(sid_sYlm_re[si][lm]), 0).coord_spherical(i)[1],
			         &sYlm_re, &sYlm_im);
                     *i = sYlm_im;
                  }
               }
               
               lm++;
            }
         
      }
   }*/

}


}
