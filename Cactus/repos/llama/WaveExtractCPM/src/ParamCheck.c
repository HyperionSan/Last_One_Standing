 /*@@
   @file      ParamCheck.c
   @date      April 26 2002
   @author    Gabrielle Allen
   @desc 
              Check parameters for WaveExtractCPM
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <stdio.h>

void WavExtrCPM_ParamCheck(CCTK_ARGUMENTS);

 /*@@
   @routine    WavExtrCPM_ParamCheck
   @date       April 26 2002
   @author     Gabrielle Allen
   @desc 
               Check parameters for WaveExtractCPM
   @enddesc 
   @calls     
@@*/
void WavExtrCPM_ParamCheck(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (verbose>3)
    CCTK_INFO("Checking parameters");

  if(! (CCTK_EQUALS(metric_type, "physical") || 
        CCTK_EQUALS(metric_type, "static conformal")))
  {
    CCTK_PARAMWARN("WaveExtractCPM only works currently with metric_type \"static conformal\" or \"physical\"");
  }
  
  int store_radial_derivs = *(CCTK_INT const*) CCTK_ParameterGet("store_radial_derivatives", "ADMDerivatives", NULL);
  if (!store_radial_derivs)
     CCTK_PARAMWARN("Radial derivatives need to be calculated in ADMDerivatives!");

  if (origin_x != 0.0 || origin_y != 0.0 || origin_z !=0)
     CCTK_PARAMWARN("Radial derivatives as calculated in ADMDerivatives assume origin(x,y,z) = 0.0!");
}
