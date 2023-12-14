#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _WIN32
#include <time.h>
#include <unistd.h>
#endif

#include <ctype.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "CactusPUGH/PUGH/src/include/pugh.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"
#include "IsoSurfacerInit.h"

static const char *rcsid = "$Id: IsoSurfacerInit.c 79 2004-05-17 12:28:58Z goodale $";
CCTK_FILEVERSION(CactusPUGHIO_IsoSurfacer_IsoSurfacerInit_c)

void IsoSurfaceEnd(cGH *GH);
static int IsoSurfacer_ParseIsoString(const char *isostring,isosurfacerGH *myGH);

/***************************************************************/
void *IsoSurfacer_SetupGH (tFleshConfig *config, 
                           int convergence_level, 
                           cGH *GH){
  isosurfacerGH *myGH=(isosurfacerGH*)malloc(sizeof(isosurfacerGH));

  /* initialize values */
  config = config;
  convergence_level = convergence_level;
  GH = GH;

  myGH->funcName=0;
  myGH->formats=0;
  myGH->outfreq=0;
  myGH->firstIteration=0;
  myGH->ComputeNormals=0;
  myGH->isovalue=0.0;
  myGH->perprocessor.verts  = myGH->totals.verts  = NULL;
  myGH->perprocessor.nverts = myGH->totals.nverts = 0;
  myGH->perprocessor.polys  = myGH->totals.polys  = NULL;
  myGH->perprocessor.npolys = myGH->totals.npolys = 0;
  myGH->perprocessor.norms = myGH->totals.norms = NULL;
  myGH->perprocessor.nnorms = myGH->totals.nnorms = 0;
  return myGH;
}

int IsoSurfacer_InitGH (cGH *GH)
{
  const char *my_out_dir;
  DECLARE_CCTK_PARAMETERS
    /*
      The above string declares the following parameters
      char *out_format
      int out_every
      int output_start
      int outer_boundary_cutoff
      char *outdir
     */
  int n,i;
  isosurfacerGH *myGH;
  int Iso_SetupServer(cGH *, isosurfacerGH *, int , int , int , int );
 

  myGH = (isosurfacerGH *) CCTK_GHExtension (GH, "IsoSurfacer");
 
  /*printf("IsoInit\n"); */
  /* initialize values */
  myGH->funcName=0;
  myGH->formats=0;
  myGH->outfreq=out_every >= 0 ? out_every : io_out_every;
  myGH->firstIteration=out_start;
  myGH->ComputeNormals=compute_normals;
  /* printf("*************  compute Normals = %u *****************\n",
     myGH->ComputeNormals); */
  myGH->isovalue=isovalue;
  myGH->perprocessor.verts  = myGH->totals.verts  = NULL;
  myGH->perprocessor.nverts = myGH->totals.nverts = 0;
  myGH->perprocessor.polys  = myGH->totals.polys  = NULL;
  myGH->perprocessor.npolys = myGH->totals.npolys = 0;
  myGH->perprocessor.norms  = myGH->totals.norms = NULL;
  myGH->perprocessor.nnorms = myGH->totals.nnorms = 0;
  myGH->minval=0.0;
  myGH->maxval=1.0;
  
  for(i=0,n=CCTK_NumVars();i<n;i++){ 
    char *fullname = CCTK_FullName (i);
    if(CCTK_Equals (fullname, out_vars)) 
      myGH->funcName=strdup (fullname);
    /* Maybe even set the GF here ? */
    free(fullname);
  }
  
  if(strstr(out_format,"UCD")) myGH->formats|=UCD;
  if(strstr(out_format,"ASCII")) myGH->formats|=ASCII;
  if(strstr(out_format,"BIN")) myGH->formats|=BIN;
  if(strstr(out_format,"SOCK")) myGH->formats|=SOCK;
  if(strstr(out_format,"HDF5")) myGH->formats|=ISOHDF5;
  if(strstr(out_format,"VRML")) myGH->formats|=VRML;

  /* OK, now we test to see if the 'isosurfacer' string 
     overrides everything set up by the new params */
  IsoSurfacer_ParseIsoString(isosurfacer,myGH);

  if(myGH->funcName==0 || myGH->formats==0)
    myGH->RunIsoSurfacer = 0;
  else 
    myGH->RunIsoSurfacer = 1;

  /* get the name for IsoSurfacer output directory and make sure it exists */
  my_out_dir = out_dir;
  if (*my_out_dir)
  {
    i = IOUtil_CreateDirectory (GH, my_out_dir, 0, 0);
    if (i < 0)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Problem creating IsoSurfacer output directory '%s'",
                  my_out_dir);
    }
  }
  else
  {
    my_out_dir = io_out_dir;
    i = 0;
  }
  if (i >= 0 && CCTK_Equals (verbose, "full"))
  {
    CCTK_VInfo (CCTK_THORNSTRING, "IsoSurfacer: Output to directory '%s'",
                my_out_dir);
  }
  Iso_SetupServer(GH,myGH,dataport,controlport, 5, 1); /* needs to move into InitGH */
  /* otherwise, the outdir need not be created if it is '.' */
  return 1;
}

/************************************************************/

void IsoSurfaceEnd(cGH *GH)
{
  isosurfacerGH *myGH;
  myGH = (isosurfacerGH *) GH->extensions [CCTK_GHExtensionHandle ("IsoSurfacer")];
  if( myGH->RunIsoSurfacer == 0 )
    return;
}

/* Parse string from original isosurfacer.
   Typical string is
   "{(wavetoy::phi) (0.35) (SOCK) (1,1,1,1)}"
   
   Since we are limiting ourselves to single-iso-functionality
   right now, we only peel off the first reference to an isosurface
   at this point.  We can ignore the other params since they are
   redundant (already functionally covered by IOBase params) 
*/
static int IsoSurfacer_ParseIsoString(const char *isostring,isosurfacerGH *myGH){
  char *s,*snext,*si,*free_me;
  int len;

  if(!isostring) return 0;
  if((len=strlen(isostring))<4) return 0; /* nothing to write home about here... */
  s = (char *)malloc(len+1);
  strcpy(s,isostring); /* we are going to do some destructive parsing here */
  free_me = s; /* remember this string for when we free it */

  /* Now we parse */
  snext=strchr(s,'('); /* move to first '(' */
  s=snext+1;
  snext=strchr(s,')'); /* find next ',' or ')' which would terminate the list */
  si=strchr(s,',');
  if(si && si<snext) *si='\0';
  else if(snext) *snext='\0';
  else {free(free_me); return 0; } /* parse failure */
  myGH->funcName=(char*)malloc(strlen(s)+1);
  strcpy(myGH->funcName,s); /* got the varname for the output var */
  /* printf("****************IsoSurf[%s]\n",myGH->funcName); */
  /* OK, now we go find the isoval */
  s = strchr(snext+1,'(');
  if(s) s++; else return 1;
  snext=strchr(s,')'); /* find next ',' or ')' which would terminate the list */
  si=strchr(s,',');
  if(si && si<snext) *si='\0';
  else if(snext) *snext='\0';
  else return 1; /* parse failure, but we at least have the right varname.*/
  myGH->isovalue = atof(s);
  /* ******************  */
  s = strchr(snext+1,'(');
  if(s) s++; else return 1;
  snext=strchr(s,')'); /* find next ',' or ')' which would terminate the list */
  if(snext) *snext='\0';
  else return 1; /* parse failure, but we at least have the right varname.*/
   
  if(strstr(s,"UCD")) myGH->formats|=UCD;
  if(strstr(s,"ASCII")) myGH->formats|=ASCII;
  if(strstr(s,"BIN")) myGH->formats|=BIN;
  if(strstr(s,"SOCK")) myGH->formats|=SOCK;
  if(strstr(s,"HDF5")) myGH->formats|=ISOHDF5;
  if(strstr(s,"VRML")) myGH->formats|=VRML;
  /* we are all done now.  The rest is ignored because it contains
     redundant or obsolete information which doesnt really fit into
     the new parser model */
  free(free_me); /* free our temporary storage */
  return 1; /* parse was successful */
}

