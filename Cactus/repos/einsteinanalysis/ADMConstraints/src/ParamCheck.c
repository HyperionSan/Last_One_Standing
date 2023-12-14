 /*@@
   @file      ParamCheck.c
   @date      Fri Apr 26 11:40:36 2002
   @author    Tom Goodale
   @desc 
   Parameter checking stuff for ADMConstraints
   @enddesc
   @version $Header$
 @@*/

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_ADMConstraints_ParamCheck_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void ADMConstraints_ParamCheck(CCTK_ARGUMENTS);
void ADMConstraints_ConformalCheck(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    ADMConstraints_ParamCheck
   @date       Fri Apr 26 11:40:36 2002
   @author     Tom Goodale
   @desc 
   Scheduled routine to detect invalid parameter settings.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void ADMConstraints_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(! CCTK_EQUALS(metric_type, "physical"))
  {
    CCTK_PARAMWARN("Unknown ADMBase::metric_type - known type is \"physical\"");
  }

}

