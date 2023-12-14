 /*@@
   @file      Startup.c
   @date      Fri May 21 1999
   @author    Thomas Radke
   @desc 
   Startup routines for IsoSurfacer.
   @enddesc 
   @history
   @endhistory
 @@*/

#include <stdio.h>

#include "cctk_Flesh.h"
#include "cctk_GHExtensions.h"
#include "cctk_IOMethods.h"
#include "IsoSurfacerInit.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Id: Startup.c 49 2001-06-14 12:57:23Z tradke $";
CCTK_FILEVERSION(CactusPUGHIO_IsoSurfacer_Startup_c)

int IsoSurfacer_Startup(void);
void IsoSurfacer_Worker(CCTK_ARGUMENTS);
/*int Iso_SetupServer(int dataport, int clientport, int queue_size, int hunt);
Iso_SetupServer(7051, 7050, 5, 1);  needs to move into InitGH */

int IsoSurfacer_Startup(void)
{
  int IOMethod;
  int IsoSurfacer_GHExtension;

  IsoSurfacer_GHExtension = CCTK_RegisterGHExtension ("IsoSurfacer");
  CCTK_RegisterGHExtensionSetupGH (IsoSurfacer_GHExtension,IsoSurfacer_SetupGH);
  CCTK_RegisterGHExtensionInitGH (IsoSurfacer_GHExtension, IsoSurfacer_InitGH);

  /* Register the IsoSurfacer routine as output method */
  IOMethod = CCTK_RegisterIOMethod ("IsoSurfacer");
  CCTK_RegisterIOMethodOutputGH (IOMethod, IsoSurfacer);
  CCTK_RegisterIOMethodTimeToOutput (IOMethod, IsoSurfacer_TimeForOutput);
  CCTK_RegisterIOMethodTriggerOutput (IOMethod, IsoSurfacer_TriggerOutput);

  return (0);
}

int Iso_Poll(cGH *cctkGH, long sec, long usec);

void IsoSurfacer_Worker(CCTK_ARGUMENTS)
{
  Iso_Poll(cctkGH, 0, 0);
}
