/*   The extensions to the GH structure from IsoSurfacer. */
#ifndef __IsoSurface_GH_h_
#define __IsoSurface_GH_h_

#include "cctk.h"

#define BIN      1
#define ASCII    2
#define UCD      4
#define VRML     8
#define SOCK     16
#define ISOHDF5  32
#define NONE     64 /* sounds kooky, but its for forcing the 
                       surfacer to run without output for
                       benchmarking purposes */

typedef struct CmdSet {
  char object[64];
  char target[64];
  char value[64];
} CmdSet;

typedef union IsoCommand {
  char buffer[64+64+64];
  CmdSet cmd;
} IsoCommand;

/* Should change name to IsoGeometry */
typedef struct polypatch
{
  CCTK_REAL4    *verts;
  CCTK_INT4      nverts;
  CCTK_INT4     *polys;
  CCTK_INT4      npolys;
  CCTK_REAL4    *norms; /* these are vertex normals */
  CCTK_INT4      nnorms;
} polypatch;

typedef struct isosurfacerGH
{
  char       *funcName; /* will be an array of names & vals later */
  CCTK_REAL   isovalue; /* must be same as GF data for efficiency */
  CCTK_REAL   minval,maxval; /* range of this grid function */
  short int   formats;
  short int   outfreq;
  short int   firstIteration;
  short int  RunIsoSurfacer;
  short int  ComputeNormals;
  polypatch   perprocessor,totals; /* un-collected geometry for this node */
} isosurfacerGH;

/* typedef isosurfacerGH isoparms_st; */

typedef enum {Byte=0,Int8=0,Int16=1,Int32=2,Int64=3,
              Float32=4,Float64=5,
              uInt8=6,uChar=6,uInt16=7,uInt32=8,uInt64=9,
              Char=10,Char8=10,String=10,Unicode=11,
              Char16=11, Error=-1} IsoType;
 
#endif /* __IsoSurface_GH_h_ */
