 /*@@
   @header    ioSDFGH.h
   @date      Sat 12 June 2004
   @author    Thomas Radke
   @desc 
              The extensions to the GH structure from IOSDF.
   @version   $Header$
 @@*/

#ifndef _IOSDF_IOSDFGH_H_
#define _IOSDF_IOSDFGH_H_ 1

#include "StoreNamedData.h"
#include "bbhutil.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct IOSDFGH
{
  /* default number of times to output */
  int out1D_every_default, out2D_every_default, out3D_every_default;

  /* number of times to output every variable */
  CCTK_INT *out1D_every, *out2D_every, *out3D_every;

  /* lists of variables to output */
  char *out1D_vars, *out2D_vars, *out3D_vars;

  /* directories in which to output */
  char *out1D_dir, *out2D_dir, *out3D_dir;

  /* the last iteration output for var [i] */
  int *out1D_last, *out2D_last, *out3D_last;

  /* database for names of output files that were already created */
  pNamedData *fileList_1D, *fileList_2D, *fileList_3D;

  /* for 1D lines, we define the index where to start the line:
     spxyz[maxdim][XYZdirection][xyzstart_index] */
  int ***spxyz;

  /* for 2D planes, we define the index where to locate the plane
     sp2xyz[maxdim][perpendicular direction] */
  int **sp2xyz;

  /* stop on I/O parameter parsing errors ? */
  int stop_on_parse_errors;

} ioSDFGH;


/* prototypes of functions to be registered as I/O methods */
int IOSDF_Output1DGH (const cGH *GH);
int IOSDF_TriggerOutput1D (const cGH *GH, int);
int IOSDF_TimeFor1D (const cGH *GH, int);
int IOSDF_Output1DVarAs (const cGH *GH, const char *var, const char *alias);
int IOSDF_Output2DGH (const cGH *GH);
int IOSDF_TriggerOutput2D (const cGH *GH, int);
int IOSDF_TimeFor2D (const cGH *GH, int);
int IOSDF_Output2DVarAs (const cGH *GH, const char *var, const char *alias);
int IOSDF_Output3DGH (const cGH *GH);
int IOSDF_TriggerOutput3D (const cGH *GH, int);
int IOSDF_TimeFor3D (const cGH *GH, int);
int IOSDF_Output3DVarAs (const cGH *GH, const char *var, const char *alias);

/* other function prototypes */
void IOSDF_Startup (void);
void IOSDF_Terminate (void);

int IOSDF_Write1D (const cGH *GH, int vindex, const char *alias);
int IOSDF_Write2D (const cGH *GH, int vindex, const char *alias);
int IOSDF_Write3D (const cGH *GH, int vindex, const char *alias);

void IOSDF_Choose1D (const cGH *GH);
void IOSDF_Choose2D (const cGH *GH);

void IOSDF_CheckSteerableParameters1D (ioSDFGH *myGH);
void IOSDF_CheckSteerableParameters2D (ioSDFGH *myGH);
void IOSDF_CheckSteerableParameters3D (ioSDFGH *myGH);

#ifdef __cplusplus
} // extern "C"
#endif

#endif  /* _IOSDF_IOSDFGH_H_ */
