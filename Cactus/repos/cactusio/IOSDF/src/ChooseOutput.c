/*@@
  @file      ChooseOutput.c
  @author    Thomas Radke
  @date      12 June 2004
  @desc
             Choose what 1D slices and 2D planes to output by IOSDF.
  @enddesc

  @version   $Id$
 @@*/

#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"
#include "ioSDFGH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusIO_IOSDF_ChooseOutput_c)


/********************************************************************
 ********************    Macro Definitions   ************************
 ********************************************************************/
/* macro to choose origin according actual parameter settings:
     1. Indices from IOSDF
     2. Indices from IOUtil
     3. Coords from IOSDF
     4. Coords from IOUtil
 */
#define GET_SLICE(IOSDF_param, IOUtil_param, index, coord)                    \
        {                                                                     \
          index = IOSDF_param##i >= 0 ? IOSDF_param##i : IOUtil_param##i;     \
          coord = IOSDF_param != -424242 ? IOSDF_param : IOUtil_param;        \
        }


/*@@
   @routine   IOSDF_Choose1D
   @author    Thomas Radke
   @date      12 June 2004
   @desc
              Use parameters to choose the 1D slices through the output data.
   @enddesc

   @calls     IOUtil_1DLines

   @var       GH
   @vdesc     pointer to CCTK grid hierarchy
   @vtype     const cGH *
   @vio       in
   @endvar
@@*/
void IOSDF_Choose1D (const cGH *GH)
{
  int i, j, maxdim;
  ioSDFGH *myGH;
  int *origin_index[3];
  CCTK_REAL *origin_phys[3];
  DECLARE_CCTK_PARAMETERS


  /* allocate arrays for origins */
  origin_phys[0] = malloc (3 * 3 * sizeof (CCTK_REAL));
  origin_phys[1] = origin_phys[0] + 3;
  origin_phys[2] = origin_phys[1] + 3;
  origin_index[0] = malloc (3 * 3 * sizeof (int));
  origin_index[1] = origin_index[0] + 3;
  origin_index[2] = origin_index[1] + 3;

  /* get slice points */
  GET_SLICE (out1D_xline_y, out_xline_y, origin_index[0][1], origin_phys[0][1]);
  GET_SLICE (out1D_xline_z, out_xline_z, origin_index[0][2], origin_phys[0][2]);
  GET_SLICE (out1D_yline_x, out_yline_x, origin_index[1][0], origin_phys[1][0]);
  GET_SLICE (out1D_yline_z, out_yline_z, origin_index[1][2], origin_phys[1][2]);
  GET_SLICE (out1D_zline_x, out_zline_x, origin_index[2][0], origin_phys[2][0]);
  GET_SLICE (out1D_zline_y, out_zline_y, origin_index[2][1], origin_phys[2][1]);

  maxdim = CCTK_MaxDim ();
  myGH = CCTK_GHExtension (GH, "IOSDF");
  myGH->spxyz = malloc (maxdim * sizeof (int **));

  for (i = 0; i < maxdim; i++)
  {
    myGH->spxyz[i] = malloc ((i + 1) * sizeof (int *));

    for (j = 0; j <= i; j++)
    {
      myGH->spxyz[i][j] = calloc (i + 1, sizeof (int));
    }

    if (i < 3)
    {
      IOUtil_1DLines (GH, i + 1, origin_index, origin_phys, myGH->spxyz[i]);
    }
  }

  /* free allocated resources */
  free (origin_phys[0]);
  free (origin_index[0]);
}


/*@@
   @routine   IOSDF_Choose2D
   @author    Thomas Radke
   @date      12 June 2004
   @desc
              Use parameters to choose the 2D slices through the output data.
   @enddesc

   @calls     IOUtil_2DPlanes

   @var       GH
   @vdesc     Pointer to CCTK grid hierarchy
   @vtype     const cGH *
   @vio       in
   @endvar
 @@*/
void IOSDF_Choose2D (const cGH *GH)
{
  int i, maxdim;
  ioSDFGH *myGH;
  int origin_index[3];
  CCTK_REAL origin_phys[3];
  DECLARE_CCTK_PARAMETERS


  GET_SLICE (out2D_xyplane_z, out_xyplane_z, origin_index[0], origin_phys[0]);
  GET_SLICE (out2D_xzplane_y, out_xzplane_y, origin_index[1], origin_phys[1]);
  GET_SLICE (out2D_yzplane_x, out_yzplane_x, origin_index[2], origin_phys[2]);

  maxdim = CCTK_MaxDim ();
  myGH = CCTK_GHExtension (GH, "IOSDF");
  myGH->sp2xyz = malloc (3 * sizeof (int *));

  for (i = 0; i < maxdim; i++)
  {
    myGH->sp2xyz[i] = calloc (i + 1, sizeof (int));

    if (i == 2)
    {
      IOUtil_2DPlanes (GH, i + 1, origin_index, origin_phys, myGH->sp2xyz[i]);
    }
  }
}
