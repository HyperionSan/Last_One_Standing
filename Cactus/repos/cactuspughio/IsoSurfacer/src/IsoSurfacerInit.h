#ifndef _ISOSURFACERINIT_H
#define _ISOSURFACERINIT_H

#include <sys/types.h>

#ifndef _WIN32
#include <sys/times.h>
#else
#include <time.h>
#endif

#include <limits.h>
#include "IsoSurfacerGH.h"

/* prototypes of functions to be registered */
void *IsoSurfacer_SetupGH (tFleshConfig *config, int convergence_level,cGH *GH);
int IsoSurfacer_InitGH (cGH *GH);
int IsoSurfacer (const cGH *GH);
int IsoSurfacer_TriggerOutput (const cGH *GH, int);
int IsoSurfacer_TimeForOutput (const cGH *GH, int);

#endif

