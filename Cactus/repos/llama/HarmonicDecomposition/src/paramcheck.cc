
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
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"



extern "C" void HarmonicDecomposition_ParamCheck(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   
   // do this only on processor 0
   if (CCTK_MyProc(cctkGH) != 0)
      return;
   
   if (lmin > lmax)
      CCTK_PARAMWARN("ERROR: lmin > lmax!");
   
   for (int i=0; i < nslices; ++i)
   {
      if (which_slice_to_take[i] >= ss_nslices)
         CCTK_PARAMWARN("which_slice_to_take[i]: Slice number is greater than available slices!");
      
      for (int j=0; j < nslices; ++j)
         if (i != j && which_slice_to_take[i] == which_slice_to_take[j])
            CCTK_PARAMWARN("which_slice_to_take[i]: You have selected the same slice twice!");
      
      if (!CCTK_Equals(grid_type, type[which_slice_to_take[i]]))
         CCTK_WARN(1, "The requested grid-type of the slice doesn't match the type as provided by SphericalSlice!\n This is not a problem.\n"\
	              "This warning is here only to notify you that I will use the type as set by SphericalSlice.\n" \
		      "To get rid of this warning you can choose a different slice or change the slice's grid-type.");

      // check, if the number of points required for the (lmax, mmax)-mode is met 
      // Eg: l=2, m=2 has 4 maxima -> we need at least twice as many points in order to resolve that
      int required_points = lmax*4;
      int recommended_points = lmax*8;

      if (CCTK_Equals(type[which_slice_to_take[i]], "1patch"))
      {
         if (nphi[which_slice_to_take[i]] < required_points || 2*ntheta[which_slice_to_take[i]] < required_points)
            CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "For mode lmax=%d you need at least nphi=2*ntheta=%d to resolve it.", int(lmax), required_points);
         if (nphi[which_slice_to_take[i]] < recommended_points || 2*ntheta[which_slice_to_take[i]] < recommended_points)
            CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "For mode lmax=%d I recommend nphi=2*ntheta=%d points.", int(lmax), recommended_points);
      }
      if (CCTK_Equals(type[which_slice_to_take[i]], "6patch"))
      {
         if (6*nphi[which_slice_to_take[i]] < required_points || 6*ntheta[which_slice_to_take[i]] < required_points)
            CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "For mode lmax=%d you need at least 6*nphi=6*ntheta=%d to resolve it.", int(lmax), required_points);
         if (6*nphi[which_slice_to_take[i]] < recommended_points || 6*ntheta[which_slice_to_take[i]] < recommended_points)
            CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "For mode lmax=%d I recommend 6*nphi=6*ntheta=%d points.", int(lmax), recommended_points);
      }
   }
   
   for (int i=0; i < number_of_vars; ++i)
   {
      if (CCTK_VarIndex(vars[i]) < 0)
         CCTK_PARAMWARN("vars[i]: grid-variable not known!");
      
      //if (precalc_sYlms && spin_weight[i] != 0) //spin_weight[i] != spin_weight[0])
      //   CCTK_WARN(1, "Current limitation: I will only use precalculated sYlm for spin-weight[i]==0!");
         //I will only use precalculated sYlm for those variables, whose spin-weight[i] is equal to the spin-weight[0] of the first variable!");
   }
}
