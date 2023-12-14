 /*@@
   @file      Startup.c
   @date      Sat 12 June 2004
   @author    Thomas Radke
   @desc
              Startup routines for IOSDF.
   @enddesc
   @version   $Id$
 @@*/


#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_IOMethods.h"
#include "cctk_Parameters.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"
#include "ioSDFGH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Id$";
CCTK_FILEVERSION(CactusIO_IOSDF_Startup_c)


/********************************************************************
 ********************    Macro Definitions   ************************
 ********************************************************************/
#define CREATE_OUTDIR(dim, dir)                                               \
        {                                                                     \
          const char *_dir = dir;                                             \
                                                                              \
                                                                              \
          /* check whether "dir" was set; if not default to "IO::out_dir" */  \
          if (*_dir == 0)                                                     \
          {                                                                   \
            _dir = out_dir;                                                   \
          }                                                                   \
                                                                              \
          /* omit the directory name if it's the current working dir */       \
          if (strcmp (_dir, ".") == 0)                                        \
          {                                                                   \
            myGH->dir = strdup ("");                                          \
          }                                                                   \
          else                                                                \
          {                                                                   \
            myGH->dir = malloc (strlen (_dir) + 2);                           \
            sprintf (myGH->dir, "%s/", _dir);                                 \
          }                                                                   \
                                                                              \
          /* create the directory */                                          \
          i = IOUtil_CreateDirectory (GH, myGH->dir, 0, 0);                   \
          if (i < 0)                                                          \
          {                                                                   \
            CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,              \
                        "IOSDF_SetupGH: Problem creating directory '%s' "   \
                        "for %dD output", myGH->dir, dim);                    \
          }                                                                   \
          else if (i >= 0 && CCTK_Equals (verbose, "full"))                   \
          {                                                                   \
            CCTK_VInfo (CCTK_THORNSTRING, "%dD output to directory '%s'",     \
                        dim, myGH->dir);                                      \
          }                                                                   \
        }


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static void *IOSDF_SetupGH (tFleshConfig *config, int conv_level, cGH *GH);


 /*@@
   @routine   IOSDF_Startup
   @date      Sat 12 June 2004
   @author    Thomas Radke
   @desc
              The startup registration routine for IOSDF.
              Registers the GH extensions needed for IOSDF
              along with its setup routine.
   @enddesc
   @calls     CCTK_RegisterGHExtension
              CCTK_RegisterGHExtensionSetupGH
@@*/
void IOSDF_Startup (void)
{
  CCTK_RegisterGHExtensionSetupGH (CCTK_RegisterGHExtension ("IOSDF"),
                                   IOSDF_SetupGH);
}


 /*@@
   @routine   IOSDF_Terminate
   @date      Sat 12 June 2004
   @author    Thomas Radke
   @desc
              The termination routine for IOSDF.
              Closes all open SDF output files.
   @enddesc
@@*/
void IOSDF_Terminate (void)
{
  gft_close_all ();
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
 /*@@
   @routine   IOSDF_SetupGH
   @date      Sat 12 June 2004
   @author    Thomas Radke
   @desc
              Allocates and sets up IOSDF's GH extension structure
   @enddesc

   @calls     CCTK_RegisterIOMethod
              CCTK_RegisterIOMethodOutputGH
              CCTK_RegisterIOMethodOutputVarAs
              CCTK_RegisterIOMethodTimeToOutput
              CCTK_RegisterIOMethodTriggerOutput

   @var       config
   @vdesc     the CCTK configuration as provided by the flesh
   @vtype     tFleshConfig *
   @vio       unused
   @endvar
   @var       conv_level
   @vdesc     the convergence level
   @vtype     int
   @vio       unused
   @endvar
   @var       GH
   @vdesc     Pointer to CCTK grid hierarchy
   @vtype     cGH *
   @vio       in
   @endvar

   @returntype void *
   @returndesc
               pointer to the allocated GH extension structure
   @endreturndesc
@@*/
static void *IOSDF_SetupGH (tFleshConfig *config, int conv_level, cGH *GH)
{
  int i, maxdim, numvars;
  ioSDFGH *myGH;
  DECLARE_CCTK_PARAMETERS


  /* suppress compiler warnings about unused variables */
  (void) (config + 0);
  (void) (conv_level + 0);
  (void) (GH + 0);

  /* allocate the GH extension and its components */
  myGH = malloc (sizeof (ioSDFGH));
  if (! myGH)
  {
    CCTK_WARN (0, "IOSDF_SetupGH: Unable to allocate memory for GH");
  }

  /* register the IOSDF routines as output methods  */
  i = CCTK_RegisterIOMethod ("IOSDF_1D");
  CCTK_RegisterIOMethodOutputGH (i, IOSDF_Output1DGH);
  CCTK_RegisterIOMethodOutputVarAs (i, IOSDF_Output1DVarAs);
  CCTK_RegisterIOMethodTimeToOutput (i, IOSDF_TimeFor1D);
  CCTK_RegisterIOMethodTriggerOutput (i, IOSDF_TriggerOutput1D);

  /* only register N-D IOSDF I/O methods
     if at least N-dimensional grid variables are defined by thorns */
  maxdim = CCTK_MaxDim ();
  if (maxdim >= 2)
  {
    i = CCTK_RegisterIOMethod ("IOSDF_2D");
    CCTK_RegisterIOMethodOutputGH (i, IOSDF_Output2DGH);
    CCTK_RegisterIOMethodOutputVarAs (i, IOSDF_Output2DVarAs);
    CCTK_RegisterIOMethodTimeToOutput (i, IOSDF_TimeFor2D);
    CCTK_RegisterIOMethodTriggerOutput (i, IOSDF_TriggerOutput2D);
  }

  if (maxdim >= 3)
  {
    i = CCTK_RegisterIOMethod ("IOSDF_3D");
    CCTK_RegisterIOMethodOutputGH (i, IOSDF_Output3DGH);
    CCTK_RegisterIOMethodOutputVarAs (i, IOSDF_Output3DVarAs);
    CCTK_RegisterIOMethodTimeToOutput (i, IOSDF_TimeFor3D);
    CCTK_RegisterIOMethodTriggerOutput (i, IOSDF_TriggerOutput3D);
  }

  numvars = CCTK_NumVars ();
  myGH->out1D_every = malloc (numvars * sizeof (CCTK_INT));
  myGH->out2D_every = malloc (numvars * sizeof (CCTK_INT));
  myGH->out3D_every = malloc (numvars * sizeof (CCTK_INT));
  myGH->out1D_last = malloc (numvars * sizeof (int));
  myGH->out2D_last = malloc (numvars * sizeof (int));
  myGH->out3D_last = malloc (numvars * sizeof (int));

  for (i = 0; i < numvars; i++)
  {
    myGH->out1D_last[i] = -1;
    myGH->out2D_last[i] = -1;
    myGH->out3D_last[i] = -1;
  }

  myGH->out1D_vars = strdup ("");
  myGH->out2D_vars = strdup ("");
  myGH->out3D_vars = strdup ("");
  myGH->out1D_every_default = out1D_every - 1;
  myGH->out2D_every_default = out2D_every - 1;
  myGH->out3D_every_default = out3D_every - 1;

  myGH->fileList_1D = NULL;
  myGH->fileList_2D = NULL;
  myGH->fileList_3D = NULL;

  myGH->stop_on_parse_errors = strict_io_parameter_check;
  if (! CCTK_Equals (verbose, "none"))
  {
    CCTK_INFO ("I/O Method 'IOSDF_1D' registered: output of 1D lines of grid "
               "functions/arrays to SDF files");
  }
  IOSDF_CheckSteerableParameters1D (myGH);
  if (maxdim >= 2)
  {
    if (! CCTK_Equals (verbose, "none"))
    {
      CCTK_INFO ("I/O Method 'IOSDF_2D' registered: output of 2D planes of "
                 "grid functions/arrays to SDF files");
    }
    IOSDF_CheckSteerableParameters2D (myGH);
  }
  if (maxdim >= 3)
  {
    if (! CCTK_Equals (verbose, "none"))
    {
      CCTK_INFO ("I/O Method 'IOSDF_3D' registered: output of 3D grid "
                 "functions/arrays to SDF files");
    }
    IOSDF_CheckSteerableParameters3D (myGH);
  }
  myGH->stop_on_parse_errors = 0;

  /* make sure all output directories exist */
  CREATE_OUTDIR (1, out1D_dir);
  CREATE_OUTDIR (2, out2D_dir);
  CREATE_OUTDIR (3, out3D_dir);

  return (myGH);
}
