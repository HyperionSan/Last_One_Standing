
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

#include "register.h"

namespace WorldTube {


using namespace std;


// the slice ids for each slice and variable
vector<vector<int> > sid;
// all Cauchy extracted variables
vector<string> extr_vars;



extern "C" void WorldTube_RegisterSlices(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   
   if (CCTK_Equals(boundary_behavior, "extract metric"))
   {
      extr_vars = vector<string>(10);
      extr_vars[0] = "ADMBase::alp";
      extr_vars[1] = "ADMBase::betax";
      extr_vars[2] = "ADMBase::betay";
      extr_vars[3] = "ADMBase::betaz";
      extr_vars[4] = "ADMBase::gxx";
      extr_vars[5] = "ADMBase::gxy";
      extr_vars[6] = "ADMBase::gxz";
      extr_vars[7] = "ADMBase::gyy";
      extr_vars[8] = "ADMBase::gyz";
      extr_vars[9] = "ADMBase::gzz";
   }
   if (CCTK_Equals(boundary_behavior, "CCE"))
   {
      extr_vars = vector<string>(30);
      extr_vars[0] = "ADMBase::alp";
      extr_vars[1] = "ADMBase::betax";
      extr_vars[2] = "ADMBase::betay";
      extr_vars[3] = "ADMBase::betaz";
      extr_vars[4] = "ADMBase::gxx";
      extr_vars[5] = "ADMBase::gxy";
      extr_vars[6] = "ADMBase::gxz";
      extr_vars[7] = "ADMBase::gyy";
      extr_vars[8] = "ADMBase::gyz";
      extr_vars[9] = "ADMBase::gzz";
      
      // radial derivatives
      
      extr_vars[10] = "ADMDerivatives::alp_dr";
      extr_vars[11] = "ADMDerivatives::betax_dr";
      extr_vars[12] = "ADMDerivatives::betay_dr";
      extr_vars[13] = "ADMDerivatives::betaz_dr";
      extr_vars[14] = "ADMDerivatives::gxx_dr";
      extr_vars[15] = "ADMDerivatives::gxy_dr";
      extr_vars[16] = "ADMDerivatives::gxz_dr";
      extr_vars[17] = "ADMDerivatives::gyy_dr";
      extr_vars[18] = "ADMDerivatives::gyz_dr";
      extr_vars[19] = "ADMDerivatives::gzz_dr";
      
      // time derivatives
      
      // this can only be checked here (and not in parametercheck) because this is initialized in basegrid!
      if (*dtlapse_state == 0 || *dtshift_state == 0)
         CCTK_WARN(0, "You need to activate storage for ADMBase::dtalp and ADMBase::dtshift for CCE!");
      
      extr_vars[20] = "ADMBase::dtalp";
      extr_vars[21] = "ADMBase::dtbetax";
      extr_vars[22] = "ADMBase::dtbetay";
      extr_vars[23] = "ADMBase::dtbetaz";
      extr_vars[24] = "ADMDerivatives::gxx_dt";
      extr_vars[25] = "ADMDerivatives::gxy_dt";
      extr_vars[26] = "ADMDerivatives::gxz_dt";
      extr_vars[27] = "ADMDerivatives::gyy_dt";
      extr_vars[28] = "ADMDerivatives::gyz_dt";
      extr_vars[29] = "ADMDerivatives::gzz_dt";
   }
   
   if (CCTK_Equals(boundary_behavior, "CCE_cart"))
   {
      extr_vars = vector<string>(50);
      extr_vars[0] = "ADMBase::alp";
      extr_vars[1] = "ADMBase::betax";
      extr_vars[2] = "ADMBase::betay";
      extr_vars[3] = "ADMBase::betaz";
      extr_vars[4] = "ADMBase::gxx";
      extr_vars[5] = "ADMBase::gxy";
      extr_vars[6] = "ADMBase::gxz";
      extr_vars[7] = "ADMBase::gyy";
      extr_vars[8] = "ADMBase::gyz";
      extr_vars[9] = "ADMBase::gzz";
      
      // cartesian derivatives
      
      extr_vars[10] = "ADMDerivatives::alp_dx";
      extr_vars[11] = "ADMDerivatives::betax_dx";
      extr_vars[12] = "ADMDerivatives::betay_dx";
      extr_vars[13] = "ADMDerivatives::betaz_dx";
      extr_vars[14] = "ADMDerivatives::gxx_dx";
      extr_vars[15] = "ADMDerivatives::gxy_dx";
      extr_vars[16] = "ADMDerivatives::gxz_dx";
      extr_vars[17] = "ADMDerivatives::gyy_dx";
      extr_vars[18] = "ADMDerivatives::gyz_dx";
      extr_vars[19] = "ADMDerivatives::gzz_dx";
      
      extr_vars[20] = "ADMDerivatives::alp_dy";
      extr_vars[21] = "ADMDerivatives::betax_dy";
      extr_vars[22] = "ADMDerivatives::betay_dy";
      extr_vars[23] = "ADMDerivatives::betaz_dy";
      extr_vars[24] = "ADMDerivatives::gxx_dy";
      extr_vars[25] = "ADMDerivatives::gxy_dy";
      extr_vars[26] = "ADMDerivatives::gxz_dy";
      extr_vars[27] = "ADMDerivatives::gyy_dy";
      extr_vars[28] = "ADMDerivatives::gyz_dy";
      extr_vars[29] = "ADMDerivatives::gzz_dy";
      
      extr_vars[30] = "ADMDerivatives::alp_dz";
      extr_vars[31] = "ADMDerivatives::betax_dz";
      extr_vars[32] = "ADMDerivatives::betay_dz";
      extr_vars[33] = "ADMDerivatives::betaz_dz";
      extr_vars[34] = "ADMDerivatives::gxx_dz";
      extr_vars[35] = "ADMDerivatives::gxy_dz";
      extr_vars[36] = "ADMDerivatives::gxz_dz";
      extr_vars[37] = "ADMDerivatives::gyy_dz";
      extr_vars[38] = "ADMDerivatives::gyz_dz";
      extr_vars[39] = "ADMDerivatives::gzz_dz";
      
      // time derivatives
      
      // this can only be checked here (and not in parametercheck) because this is initialized in basegrid!
      if (*dtlapse_state == 0 || *dtshift_state == 0)
         CCTK_WARN(0, "You need to activate storage for ADMBase::dtalp and ADMBase::dtshift for CCE!");
      
      extr_vars[40] = "ADMBase::dtalp";
      extr_vars[41] = "ADMBase::dtbetax";
      extr_vars[42] = "ADMBase::dtbetay";
      extr_vars[43] = "ADMBase::dtbetaz";
      extr_vars[44] = "ADMDerivatives::gxx_dt";
      extr_vars[45] = "ADMDerivatives::gxy_dt";
      extr_vars[46] = "ADMDerivatives::gxz_dt";
      extr_vars[47] = "ADMDerivatives::gyy_dt";
      extr_vars[48] = "ADMDerivatives::gyz_dt";
      extr_vars[49] = "ADMDerivatives::gzz_dt";
   }
   
   
   sid = vector<vector<int> >(ntubes);
   
   for (int i=0; i < ntubes; ++i)
   {
      sid[i] = vector<int>(extr_vars.size(), 0);
      for (int j=0; j < extr_vars.size(); ++j)
      {
         sid[i][j] = SphericalSlice_RegisterVariable(extr_vars[j].c_str(), which_slice_to_take[i], 1, "split");
      }
   }
}



}

