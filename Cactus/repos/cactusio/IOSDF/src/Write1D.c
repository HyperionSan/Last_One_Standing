/*@@
   @file      Write1D.c
   @date      Sat 12 June 2004
   @author    Thomas Radke
   @desc
              Output one-dimensional lines in SDF file format.
   @enddesc
   @version   $Id$
@@*/

#include <math.h>      /* sqrt(3) */
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
CCTK_FILEVERSION(CactusIO_IOSDF_Write1D_c)


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static char *OpenFile (const cGH *GH,
                       const char *fullname,
                       const char *alias,
                       const cGroup *gdata,
                       int dir);


/*@@
   @routine IOSDF_Write1D
   @date    Sat 12 June 2004
   @author  Thomas Radke
   @desc
            This routine does 1D line output along the orthogonals
            and the diagonal (in case of a cubed grid).
            <p>
            It writes to SDF files suitable for gnuplot and xgraph.
            A header telling the physical time prefixes the output data.
   @enddesc
   @calls   Hyperslab_GlobalMappingByIndex
            Hyperslab_FreeMapping
            Hyperslab_GetList
            OpenFile
            WriteData

   @var     GH
   @vdesc   Pointer to CCTK GH
   @vtype   const cGH *
   @vio     in
   @endvar
   @var     vindex
   @vdesc   global index of variable to output
   @vtype   int
   @vio     in
   @endvar
   @var     alias
   @vdesc   alias name (used for creating the output filename)
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
                0 for success, or<BR>
               -1 if variable has no storage assigned
   @endreturndesc
@@*/
int IOSDF_Write1D (const cGH *GH, int vindex, const char *alias)
{
  ioSDFGH *myGH;
  int do_dir[4];
  int i, dir, myproc, gindex, have_coords, mapping;
  int *extent_int;
  cGroup gdata;
  char *fullname, *groupname, *filename;
  CCTK_INT coord_system_handle, coord_handles[3];
  CCTK_REAL minext[3], delta[3];
  CCTK_INT downsample[4];
  CCTK_INT *origin, *direction;
  CCTK_INT hsize, extent;
  double *hdata;
  double bbox[2];
  const double dtime = GH->cctk_time;
  DECLARE_CCTK_PARAMETERS


  /* get the variable's group index and its full name */
  gindex = CCTK_GroupIndexFromVarI (vindex);
  fullname = CCTK_FullName (vindex);

  /* check if variable has storage assigned */
  if (! CCTK_QueryGroupStorageI (GH, gindex))
  {
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "IOSDF_Write1D: No IOSDF_1D output for '%s' (no storage)",
                fullname);
    free (fullname);
    return (-1);
  }

  /* get the handle for IOSDF extensions */
  myGH = CCTK_GHExtension (GH, "IOSDF");

  /* get the variable's group information */
  CCTK_GroupData (gindex, &gdata);

  /* see what slices should be output */
  do_dir[0] = out1D_x && gdata.dim >= 1;
  do_dir[1] = out1D_y && gdata.dim >= 2;
  do_dir[2] = out1D_z && gdata.dim >= 3;
  /* diagonal slice is done only if variable is non-staggered and 3D */
  do_dir[3] = out1D_d && gdata.dim == 3 && gdata.stagtype == 0;
  if (out1D_d && ! do_dir[3] && myGH->out1D_last[vindex] < 0)
  {
    CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,
                "IOSDF_Write1D: No IOSDF_1D diagonal output for '%s' "
                "(only implemented for non-staggered 3D variables)",
                fullname);
  }

  /* return if nothing to do */
  if (! (do_dir[0] || do_dir[1] || do_dir[2] || do_dir[3]))
  {
    free (fullname);
    return (0);
  }

  /* get the coordinate system associated with this grid variable */
  coord_system_handle = -1;
  if (CCTK_IsFunctionAliased ("Coord_GroupSystem"))
  {
    groupname = CCTK_GroupName (gindex);
    coord_system_handle = Coord_GroupSystem (GH, groupname);
    free (groupname);
  }

  myproc = CCTK_MyProc (GH);

  origin     = calloc (2*gdata.dim, sizeof (CCTK_INT));
  direction  = origin + gdata.dim;
  extent_int = malloc ((gdata.dim + 1) * sizeof (int));

  /* set downsampling vector from I/O parameters */
  downsample[0] = out_downsample_x;
  downsample[1] = out_downsample_y;
  downsample[2] = out_downsample_z;
  downsample[3] = 1;

  dir = gdata.dim < 3 ? gdata.dim : 3;

  have_coords = coord_system_handle >= 0 &&
                Util_TableGetIntArray (coord_system_handle, dir,
                                       coord_handles, "COORDINATES") >= 0;
  if (have_coords)
  {
    /* get the coordinates ranges */
    for (i = 0; i < dir; i++)
    {
      Util_TableGetReal (coord_handles[i], &minext[i], "COMPMIN");
      Util_TableGetReal (coord_handles[i], &delta[i], "DELTA");
      delta[i] *= downsample[i];
    }
  }

  /* get the variable's extents, compute the extent for 3D-diagonals as the
     minimum of grid points in each direction */
  CCTK_GroupgshVI (GH, gdata.dim, extent_int, vindex);
  if (gdata.dim == 3)
  {
    extent_int[3] = extent_int[0] < extent_int[1] ?
                    extent_int[0] : extent_int[1];
    if (extent_int[2] < extent_int[3])
    {
      extent_int[3] = extent_int[2];
    }
  }
  /* get the total number of grid points to check for zero-sized variables */
  for (dir = 0, hsize = 1; dir < gdata.dim; dir++)
  {
    hsize *= extent_int[dir];
  }

  /* now do the actual I/O looping over all directions */
  for (dir = 0; dir < 4; dir++)
  {
    /* skip empty slices */
    if (hsize <= 0 || ! do_dir[dir])
    {
      continue;
    }

    /* processor 0 opens the files with the appropriate name */
    filename = NULL;
    if (myproc == 0)
    {
      filename = OpenFile (GH, fullname, alias, &gdata, dir);
    }

    /* set the direction vector */
    for (i = 0; i < gdata.dim; i++)
    {
      direction[i] = (dir == i || dir == 3) ? 1 : 0;
    }

    /* set the extent */
    extent = extent_int[dir];

    /* set the origin of the line */
    if (gdata.grouptype == CCTK_GF && dir < 3)
    {
      for (i = 0; i < gdata.dim; i++)
      {
        origin[i] = myGH->spxyz[gdata.dim-1][dir][i];
      }
      extent -= origin[dir];

      /* correct extent in the case of staggered grids */
      if (CCTK_StaggerDirIndex (dir, gdata.stagtype) == 1)
      {
        extent--;
      }
    }
    else    /* origin for CCTK_ARRAYS is always (0, 0, 0) */
    {
      memset (origin, 0, gdata.dim * sizeof (CCTK_INT));
    }

    mapping = Hyperslab_GlobalMappingByIndex (GH, vindex, 1,
                                              direction, origin, &extent,
                                              &downsample[dir],
                                              -1,   /* table handle */
                                              NULL  /* conversion fn */,
                                              &hsize);
    if (mapping < 0)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "IOSDF_Write1D: Failed to define hyperslab mapping for "
                  "variable '%s'", fullname);
      continue;
    }
    if (hsize <= 0)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "IOSDF_Write1D: selected hyperslab has zero size for "
                  "variable '%s' direction %d", fullname, dir);
      Hyperslab_FreeMapping (mapping);
      continue;
    }

    /* allocate hyperslab buffer on I/O processor */
    hdata = myproc == 0 ? malloc (hsize * sizeof (double)) : NULL;

    /* get the hyperslab */
    i = Hyperslab_Get (GH, mapping, 0, vindex, 0, gdata.vartype, hdata);

    /* release the mapping structure */
    Hyperslab_FreeMapping (mapping);

    /* And dump the data to file */
    if (i != 0)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "IOSDF_Write1D: Failed to extract hyperslab for "
                  "variable '%s'", fullname);
    }
    else if (filename)
    {
      bbox[0] = -1; bbox[1] = +1;
      if (have_coords)
      {
        bbox[0] = minext[dir];
        bbox[1] = dir < 3 ? bbox[0] + delta[dir]*(hsize-1) :
                  (minext[0] + hsize*GH->cctk_delta_space[0]) * sqrt (3);
      }
      gft_out_set_bbox (bbox, 1);

      extent_int[0] = hsize;
      if (gft_out_brief (filename, dtime, extent_int, 1, hdata) <= 0)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Error writing to 1D IOSDF output file '%s'", filename);
      }
    }

    /* clean up */
    free (hdata);

  } /* end of loop through all directions */

  /* free allocated resources */
  free (origin);
  free (fullname);
  free (extent_int);

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
               Opens a set of SDF files for a given alias name.
               If this is the first time through, it will advertise them
               to IOUtil.
   @enddesc
 @@*/
static char *OpenFile (const cGH *GH,
                       const char *fullname,
                       const char *alias,
                       const cGroup *gdata,
                       int dir)
{
  ioSDFGH *myGH;
  int upper, lower;
  char *filename, *nameddata;
  char slicename[40];
  ioAdvertisedFileDesc advertised_file;
  DECLARE_CCTK_PARAMETERS


  /* get handle for and IOSDF GH extensions */
  myGH = CCTK_GHExtension (GH, "IOSDF");

  nameddata = malloc (strlen (alias) + 3);
  sprintf (nameddata, "%s_%d", alias, dir);
  filename = GetNamedData (myGH->fileList_1D, nameddata);
  if (filename)
  {
    free (nameddata);
    return (filename);
  }

  /* get the indices into spxyz[] */
  lower = (dir + 1) % 3;
  upper = (dir + 2) % 3;
  if (upper < lower)
  {
    upper = lower;
    lower = 0;
  }

  if (dir < 3)
  {
    if (gdata->dim == 1)
    {
      strcpy (slicename, "1D");
    }
    else if (gdata->dim == 2)
    {
      /* give the slice origin as range [1 .. n] */
      sprintf (slicename, "%c_%d", 'x' + dir,
               gdata->grouptype == CCTK_GF ?
               myGH->spxyz[gdata->dim-1][dir][lower] : 0);
    }
    else
    {
      /* give the slice origin as range [1 .. n] */
      sprintf (slicename, "%c_%d_%d", 'x' + dir,
               gdata->grouptype == CCTK_GF ?
               myGH->spxyz[gdata->dim-1][dir][lower] : 0,
               gdata->grouptype == CCTK_GF ?
               myGH->spxyz[gdata->dim-1][dir][upper] : 0);
    }
  }
  else
  {
    sprintf (slicename, "%dD_diagonal", gdata->dim);
  }

  filename = malloc (strlen (myGH->out1D_dir) + strlen (alias) +
                     sizeof (slicename) + 6);
  sprintf (filename, "%s%s_%s.sdf", myGH->out1D_dir, alias, slicename);
  StoreNamedData (&myGH->fileList_1D, nameddata, filename);
  free (nameddata);

  /* advertise the output file */
  advertised_file.slice = slicename;
  advertised_file.thorn = CCTK_THORNSTRING;
  advertised_file.varname = fullname;
  advertised_file.description = "One-dimensional line plots";
  advertised_file.mimetype = "application/sdf";

  IOUtil_AdvertiseFile (GH, filename, &advertised_file);

  return (filename);
}
