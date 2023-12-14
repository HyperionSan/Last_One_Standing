
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





namespace WorldTube {


extern "C" void WorldTube_Decompose(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   
   int mode_dim = (lmax+1)*(lmax+1);
   
   // set all coeffs to zero initially
   for (int i=0; i < 50*ntubes*(lmax+1)*(lmax+1); ++i)
      extraction_vars[i] = CCTK_Cmplx(0,0);
   
   for (int si=0; si < ntubes; ++si)
   {
      
      // first of all we do a collective interpolation onto the sphere
      SphericalSlice_CollectiveSyncVariables(cctkGH, &sid[si].front(), extr_vars.size());

      for (int j=0; j < extr_vars.size(); ++j)
      {
         SphericalSlice_ContractVariableWithAllsYlm(sid[si][j], 0, 0, 0, lmax, &extraction_vars[mode_dim*(si + ntubes*j)]);
      }
   }
}


}



