 /*@@
   @file      Output3D.c
   @date      Sat 12 June 2004
   @author    Thomas Radke
   @desc
              Functions to deal with 3D SDF output of variables
   @enddesc
   @version   $Id$
 @@*/

#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "util_String.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"
#include "ioSDFGH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Id$";
CCTK_FILEVERSION(CactusIO_IOSDF_Output3D_c)


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static int CheckOutputVar (int vindex);
static void SetOutputFlag (int vindex, const char *optstring, void *arg);


/*@@
   @routine    IOSDF_Output3DGH
   @date       Sat 12 June 2004
   @author     Thomas Radke
   @desc
               Loops over all variables and outputs them if necessary
   @enddesc
   @calls      IOSDF_TimeFor3D
               IOSDF_Write3D

   @var        GH
   @vdesc      pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               the number of variables which were output at this iteration
               (or 0 if it wasn't time to output yet)
   @endreturndesc
@@*/
int IOSDF_Output3DGH (const cGH *GH)
{
  int vindex, retval;
  const ioSDFGH *myGH;


  retval = 0;
  myGH = CCTK_GHExtension (GH, "IOSDF");

  /* loop over all variables */
  for (vindex = CCTK_NumVars () - 1; vindex >= 0; vindex--)
  {
    if (IOSDF_TimeFor3D (GH, vindex) &&
        IOSDF_Write3D (GH, vindex, CCTK_VarName (vindex)) == 0)
    {
      /* register variable as having 3D output this iteration */
      myGH->out3D_last[vindex] = GH->cctk_iteration;
      retval++;
    }
  }

  return (retval);
}


/*@@
   @routine    IOSDF_Output3DVarAs
   @date       Sat 12 June 2004
   @author     Thomas Radke
   @desc
               Unconditional output of a variable
               using the IOSDF 3D output method
   @enddesc
   @calls      IOSDF_Write3D

   @var        GH
   @vdesc      pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        fullname
   @vdesc      complete name of variable to output
   @vtype      const char *
   @vio        in
   @endvar
   @var        alias
   @vdesc      alias name of variable to output
               (used to generate output filename)
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine IOSDF_Write3D, or<BR>
               -1 if variable cannot be output by IOSDF_3D
   @endreturndesc
@@*/
int IOSDF_Output3DVarAs (const cGH *GH, const char *fullname, const char *alias)
{
  int vindex, retval;


  retval = -1;
  vindex = CCTK_VarIndex (fullname);
  if (vindex<0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "There is no such variable %s.  Ignoring IOSDF 3D output.",
                fullname);
  } else if (CheckOutputVar (vindex) == 0)
  {
    retval = IOSDF_Write3D (GH, vindex, alias);
  }

  return (retval);
}


 /*@@
   @routine    IOSDF_TimeFor3D
   @date       Sat 12 June 2004
   @author     Thomas Radke
   @desc
               Decides if it is time to output a variable
               using the IOSDF 3D output method
   @enddesc

   @var        GH
   @vdesc      pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        vindex
   @vdesc      index of variable
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               1 if output should take place at this iteration, or<BR>
               0 if not
   @endreturndesc
@@*/
int IOSDF_TimeFor3D (const cGH *GH, int vindex)
{
  int retval;
  char *fullname;
  ioSDFGH *myGH;


  myGH = CCTK_GHExtension (GH, "IOSDF");
  IOSDF_CheckSteerableParameters3D (myGH);

  /* check if this variable should be output */
  retval = myGH->out3D_every[vindex] > 0 &&
           GH->cctk_iteration % myGH->out3D_every[vindex] == 0;
  if (retval)
  {
    /* check if variable wasn't already output this iteration */
    if (myGH->out3D_last[vindex] == GH->cctk_iteration)
    {
      fullname = CCTK_FullName (vindex);
      CCTK_VWarn (5, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Already done IOSDF 3D output for '%s' in current "
                  "iteration (probably via triggers)", fullname);
      free (fullname);
      retval = 0;
    }
  }

  return (retval);
}


/*@@
   @routine    IOSDF_TriggerOutput3D
   @date       Sat 12 June 2004
   @author     Thomas Radke
   @desc
               Triggers the output of a variable
               using the IOSDF 3D output method
   @enddesc
   @calls      IOSDF_Write3D

   @var        GH
   @vdesc      pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        vindex
   @vdesc      index of variable to output
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine IOSDF_Write3D
   @endreturndesc
@@*/
int IOSDF_TriggerOutput3D (const cGH *GH, int vindex)
{
  int retval;
  ioSDFGH *myGH;


  /* do the 3D output */
  retval = IOSDF_Write3D (GH, vindex, CCTK_VarName (vindex));
  if (retval == 0)
  {
    /* register variables as having 3D output this iteration */
    myGH = CCTK_GHExtension (GH, "IOSDF");
    myGH->out3D_last[vindex] = GH->cctk_iteration;
  }

  return (retval);
}


/*@@
   @routine    IOSDF_CheckSteerableParameters3D
   @date       Sat 12 June 2004
   @author     Thomas Radke
   @desc
               Check if steerable IOSDF 3D parameters have changed
   @enddesc
   @calls      CCTK_TraverseString

   @var        myGH
   @vdesc      Pointer to IOSDF GH
   @vtype      ioSDFGH *
   @vio        inout
   @endvar
@@*/
void IOSDF_CheckSteerableParameters3D (ioSDFGH *myGH)
{
  int i, num_vars;
  char *fullname, *msg;
  DECLARE_CCTK_PARAMETERS


  /* how often to output */
  i = myGH->out3D_every_default;
  myGH->out3D_every_default = out3D_every >= 0 ? out3D_every : out_every;

  /* report if frequency changed */
  if (myGH->out3D_every_default != i && ! CCTK_Equals (verbose, "none"))
  {
    if (myGH->out3D_every_default > 0)
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Periodic 3D output every %d "
                  "iterations", myGH->out3D_every_default);
    }
    else
    {
      CCTK_INFO ("Periodic 3D output turned off");
    }
  }

  /* re-parse the 'out3D_vars' parameter if it was changed */
  if (strcmp (out3D_vars, myGH->out3D_vars) || myGH->out3D_every_default != i)
  {
    num_vars = CCTK_NumVars ();
    memset (myGH->out3D_every, 0, num_vars * sizeof (CCTK_INT));
    if (CCTK_TraverseString (out3D_vars, SetOutputFlag, myGH,
                             CCTK_GROUP_OR_VAR) < 0)
    {
      CCTK_WARN (myGH->stop_on_parse_errors ? 0 : 1,
                 "error while parsing parameter 'IOSDF::out3D_vars'");
    }

    if (myGH->out3D_every_default == i || ! CCTK_Equals (verbose, "none"))
    {
      msg = NULL;
      for (i = 0; i < num_vars; i++)
      {
        if (myGH->out3D_every[i] > 0)
        {
          fullname = CCTK_FullName (i);
          if (! msg)
          {
            Util_asprintf (&msg, "Periodic 3D output requested for "
                           "'%s'", fullname);
          }
          else
          {
            char *tmp = msg;
            Util_asprintf (&msg, "%s, '%s'", msg, fullname);
            free (tmp);
          }
          free (fullname);
        }
      }
      if (msg)
      {
        CCTK_INFO (msg);
        free (msg);
      }
    }

    /* save the last setting of 'out3D_vars' parameter */
    free (myGH->out3D_vars);
    myGH->out3D_vars = strdup (out3D_vars);
  }
}


/**************************** local functions ******************************/
/* check if this variable can be output (static conditions) */
static int CheckOutputVar (int vindex)
{
  int groupindex;
  cGroup groupinfo;
  char *fullname;
  const char *errormsg;


  /* get the variable group information */
  groupindex = CCTK_GroupIndexFromVarI (vindex);
  CCTK_GroupData (groupindex, &groupinfo);

  /* check if variable is grid array type with 3 dimensions */
  errormsg = NULL;
  if (groupinfo.dim != 3)
  {
    errormsg = "dim != 3";
  }
  else if (groupinfo.grouptype != CCTK_GF && groupinfo.grouptype != CCTK_ARRAY)
  {
    errormsg = "not a grid function or array";
  }
  else if (strncmp (CCTK_VarTypeName (groupinfo.vartype),
                                      "CCTK_VARIABLE_REAL", 18) ||
           CCTK_VarTypeSize (groupinfo.vartype) != sizeof (double))
  {
    errormsg = "not a double-precision floating-point grid function or array";
  }

  if (errormsg)
  {
    fullname = CCTK_FullName (vindex);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "No IOSDF 3D output for '%s': %s", fullname, errormsg);
    free (fullname);
  }

  return (errormsg != NULL);
}


/* callback for CCTK_TraverseString() to set the output flag
   for the given variable */
static void SetOutputFlag (int vindex, const char *optstring, void *arg)
{
  const ioSDFGH *myGH = arg;


  if (CheckOutputVar (vindex) == 0)
  {
    myGH->out3D_every[vindex] = myGH->out3D_every_default;

    if (optstring)
    {
      IOUtil_ParseOutputFrequency ("3D IOSDF", "IOSDF::out3D_vars", 
                                   myGH->stop_on_parse_errors,
                                   vindex, optstring,
                                   &myGH->out3D_every[vindex], NULL);
    }
  }
}
