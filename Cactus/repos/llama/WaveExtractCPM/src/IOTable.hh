
#ifndef __iotable_h
#define __iotable_h

#include "cctk.h"

#ifdef __cplusplus
extern "C" {
#endif

void WaveExtractCPM_OutputTableRowASCII(CCTK_ARGUMENTS, const char *name,
                                        int npoints, CCTK_REAL data[]);

void WaveExtractCPM_OutputTableRowHDF5(CCTK_ARGUMENTS, const char *basename,
                                       const char *datasetname,
                                       unsigned int npoints, CCTK_REAL data[]);

#ifdef __cplusplus
} //end extern "C"
#endif
  
#endif
