/*@@
   @file      Write3D.c
   @date      Sat 12 June 2004
   @author    Thomas Radke
   @desc
              Three-dimensional output of variables in SDF file format.
   @enddesc
   @version   $Id$
 @@*/


#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "util_Table.h"
#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_AdvertisedFiles.h"
#include "ioSDFGH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusIO_IOSDF_Write3D_c)


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static char *OpenFile (const cGH *GH, const char *fullname, const char *alias);

/*@@
   @routine   IOSDF_Write3D
   @date      Sat 12 June 2004
   @author    Thomas Radke
   @desc
              Writes the 3D volume of a variable into a gnuplot SDF file.
   @enddesc
   @calls     Hyperslab_GlobalMappingByIndex
              Hyperslab_FreeMapping
              Hyperslab_Get
              OpenFile

   @var       GH
   @vdesc     Pointer to CCTK GH
   @vtype     const cGH *
   @vio       in
   @endvar
   @var       vindex
   @vdesc     index of variable to output
   @vtype     int
   @vio       in
   @endvar
   @var       alias
   @vdesc     alias name of variable to output
   @vtype     const char *
   @vio       in
   @endvar

   @returntype int
   @returndesc
                0 for success, or<BR>
               -1 if variable has no storage assigned<BR>
               -2 if output file couldn't be opened<BR>
               -3 if hyperslab for coordinates and/or variable couldn't be
                  extracted
   @endreturndesc
@@*/
int IOSDF_Write3D (const cGH *GH, int vindex, const char *alias)
{
  int i, total_hsize;
  int myproc, gindex;
  cGroup gdata;
  CCTK_INT coord_system_handle, coord_handles[3];
  char *fullname, *groupname, *filename;
  double *hdata;
  int extent_int[3];
  int mapping;
  CCTK_INT extent[3], downsample[3], hsize[3];
  double bbox_dbl[2*3];
  CCTK_REAL bbox[2*3], delta[3];
  const double dtime = GH->cctk_time;
  const CCTK_INT origin[] = {0, 0, 0},
                 direction[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  DECLARE_CCTK_PARAMETERS


  /* get the variable group information */
  fullname = CCTK_FullName (vindex);
  gindex = CCTK_GroupIndexFromVarI (vindex);
  CCTK_GroupData (gindex, &gdata);

  /* check if variable has storage assigned */
  if (! CCTK_QueryGroupStorageI (GH, gindex))
  {
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "No IOSDF 3D output for '%s' (no storage)", fullname);
    free (fullname);
    return (-1);
  }

  /* get the coordinate system associated with this grid variable */
  coord_system_handle = -1;
  if (CCTK_IsFunctionAliased ("Coord_GroupSystem"))
  {
    groupname = CCTK_GroupName (gindex);
    coord_system_handle = Coord_GroupSystem (GH, groupname);
    free (groupname);
  }

  /* get the total number of grid points to check for zero-sized variables */
  /* set the extent vector (copy from 'int' to 'CCTK_INT') */
  CCTK_GroupgshVI (GH, 3, extent_int, vindex);
  for (i = 0, total_hsize = 1; i < 3; i++)
  {
    total_hsize *= extent_int[i];
    extent[i] = extent_int[i];
  }
  if (total_hsize <= 0)
  {
    free (fullname);
    return (0);
  }

  downsample[0] = out_downsample_x;
  downsample[1] = out_downsample_y;
  downsample[2] = out_downsample_z;

  /* What processor are we on? */
  myproc = CCTK_MyProc (GH);

  /* Open the file on processor 0 */
  filename = myproc == 0 ? OpenFile (GH, fullname, alias) : NULL;

  /* get the hyperslab mapping */
  mapping = Hyperslab_GlobalMappingByIndex (GH, vindex, 3,
                                            direction, origin, extent,
                                            downsample,
                                            -1,   /* table handle */
                                            NULL  /* conversion fn */,
                                            hsize);
  if (mapping < 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "IOSDF_Write3D: Failed to define hyperslab mapping for "
                "variable '%s'", fullname);
    free (fullname);
    return (-1);
  }
  for (i = 0, total_hsize = 1; i < 3; i++)
  {
    extent_int[i] = hsize[i];
    total_hsize *= hsize[i];
  }
  if (total_hsize <= 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "IOSDF_Write3D: selected hyperslab has zero size for "
                "variable '%s'", fullname);
    Hyperslab_FreeMapping (mapping);
    free (fullname);
    return (-1);
  }

  /* get the bounding box information */
  for (i = 0; i < 3; i++)
  {
    bbox_dbl[2*i + 0] = -1;
    bbox_dbl[2*i + 1] = +1;
  }
  if (coord_system_handle >= 0 &&
      Util_TableGetIntArray (coord_system_handle, 3, coord_handles,
                             "COORDINATES") >= 0)
  {
    for (i = 0; i < 3; i++)
    {
      Util_TableGetReal (coord_handles[i], &bbox[2*i + 0], "COMPMIN");
      Util_TableGetReal (coord_handles[i], &delta[i], "DELTA");

      delta[i] *= downsample[i];
      bbox_dbl[2*i + 0] = bbox[2*i + 0];
      bbox_dbl[2*i + 1] = bbox_dbl[2*i + 0] + delta[i]*(hsize[i]-1);
    }
  }
  gft_out_set_bbox (bbox_dbl, 3);

  /* allocate hyperslab buffer */
  hdata = myproc == 0 ? malloc (total_hsize * sizeof (double)) : NULL;
 
  /* get the hyperslab */
  i = Hyperslab_Get (GH, mapping, 0, vindex, 0, gdata.vartype, hdata);

  /* And dump the data to file */
  if (myproc == 0)
  {
    if (i == 0)
    {
      if (gft_out_brief (filename, dtime, extent_int, 3, hdata) <= 0)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Error writing to 3D IOSDF output file '%s'", filename);
      }
    }
    else
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "IOSDF_Write3D: Failed to extract hyperslab for "
                  "variable '%s'", fullname);
    }

    /* clean up */
    free (hdata);
  } /* end of outputting the data by processor 0 */

  free (fullname);

  return (0);
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
/*@@
   @routine    OpenFile
   @date       Sat 12 June 2004
   @author     Thomas Radke
   @desc
               Opens an SDF file for a given alias name.
               If this is the first time through, it will advertise it
               to IOUtil.
   @enddesc

   @returntype char *
   @returndesc
               the full filename of the SDF output file
   @endreturndesc
 @@*/
static char *OpenFile (const cGH *GH, const char *fullname, const char *alias)
{
  char *filename;
  ioSDFGH *myGH;
  ioAdvertisedFileDesc advertised_file;
  DECLARE_CCTK_PARAMETERS


  /* get handle for IOSDF GH extensions */
  myGH = CCTK_GHExtension (GH, "IOSDF");

  /* see if we are the first time through */
  filename = GetNamedData (myGH->fileList_3D, alias);
  if (! filename)
  {
    filename = malloc (strlen (myGH->out3D_dir) + strlen (alias) + 9);

    /* open/create the file */
    sprintf (filename, "%s%s_3D.sdf", myGH->out3D_dir, alias);

    /* advertise the file for downloading and write file info */
    advertised_file.slice = "";
    advertised_file.thorn = CCTK_THORNSTRING;
    advertised_file.varname = fullname;
    advertised_file.description = "Full-dimensional variable contents";
    advertised_file.mimetype = "application/sdf";

    IOUtil_AdvertiseFile (GH, filename, &advertised_file);

    /* store filename in database */
    StoreNamedData (&myGH->fileList_3D, alias, filename);
  }

  return (filename);
}
