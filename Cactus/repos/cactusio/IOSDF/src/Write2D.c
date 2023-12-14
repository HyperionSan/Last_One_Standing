/*@@
   @file      Write2D.c
   @date      Sat 12 June 2004
   @author    Thomas Radke
   @desc
              Output two-dimensional slices in SDF file format.
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
static const char *rcsid = "$Id$";
CCTK_FILEVERSION(CactusIO_IOSDF_Write2D_c)


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static char **OpenFile (const cGH *GH,
                        const char *fullname,
                        const char *alias,
                        int dim,
                        int maxdir);


/*@@
   @routine   IOSDF_Write2D
   @date      Sat 12 June 2004
   @author    Thomas Radke
   @desc
              Writes the 2D slices of a variable into separate SDF files.
   @enddesc
   @calls     Hyperslab_GlobalMappingByIndex
              Hyperslab_FreeMapping
              Hyperslab_Get
              OpenFile
              WriteData

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
int IOSDF_Write2D (const cGH *GH, int vindex, const char *alias)
{
  ioSDFGH *myGH;
  int i, total_hsize;
  int dir, dir_i, dir_j, maxdir, myproc, gindex, have_coords;
  int mapping;
  cGroup gdata;
  CCTK_INT coord_system_handle, coord_handles[3];
  char *fullname, *groupname;
  int full_extent[3], slice_extent[2];
  double bbox[2*2];
  CCTK_INT origin[3], extent[2], direction[6], downsample[2], hsize[2];
  CCTK_REAL minext[2], delta[2];
  double *hdata;
  char **filenames;
  const double dtime = GH->cctk_time;
  DECLARE_CCTK_PARAMETERS


  /* get the variable name and group information */
  fullname = CCTK_FullName (vindex);
  gindex = CCTK_GroupIndexFromVarI (vindex);
  CCTK_GroupData (gindex, &gdata);

  /* check if variable has storage assigned */
  if (! CCTK_QueryGroupStorageI (GH, gindex))
  {
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "No IOSDF 2D output for '%s' (no storage)", fullname);
    free (fullname);
    return (-1);
  }

  /* get the handle for IOSDF extensions */
  myGH = CCTK_GHExtension (GH, "IOSDF");

  /* get the number of slices to output */
  /* in general: maxdir = gdata.dim * (gdata.dim - 1) / 2; */
  maxdir = gdata.dim == 2 ? 1 : 3;

  /* get the coordinate system associated with this grid variable */
  coord_system_handle = -1;
  if (CCTK_IsFunctionAliased ("Coord_GroupSystem"))
  {
    groupname = CCTK_GroupName (gindex);
    coord_system_handle = Coord_GroupSystem (GH, groupname);
    free (groupname);
  }

  dir = gdata.dim < 3 ? gdata.dim : 3;

  have_coords = coord_system_handle >= 0 &&
                Util_TableGetIntArray (coord_system_handle, dir,
                                       coord_handles, "COORDINATES") >= 0;

  /* processor 0 opens the files on the first trip through */
  filenames = NULL;
  myproc = CCTK_MyProc (GH);
  if (myproc == 0)
  {
    filenames = OpenFile (GH, fullname, alias, gdata.dim, maxdir);
  }

  /* get the extents of the variable */
  CCTK_GroupgshVI (GH, gdata.dim, full_extent, vindex);

  /* get the total number of grid points to check for zero-sized variables */
  for (dir = 0, hsize[0] = 1; dir < gdata.dim; dir++)
  {
    hsize[0] *= full_extent[dir];
  }

  /* now do the actual I/O looping over all directions */
  for (dir = 0; dir < maxdir; dir++)
  {
    if (hsize[0] <= 0)
    {
      continue;
    }

    /* get the directions to span the hyperslab */
    if (dir == 0)
    {
      dir_i = 0; dir_j = 1;   /* xy */
      downsample[0] = out_downsample_x; downsample[1] = out_downsample_y;
    }
    else if (dir == 1)
    {
      dir_i = 0; dir_j = 2;   /* xz */
      downsample[0] = out_downsample_x; downsample[1] = out_downsample_z;
    }
    else
    {
      dir_i = 1; dir_j = 2;   /* yz */
      downsample[0] = out_downsample_y; downsample[1] = out_downsample_z;
    }

    /* set the extent vector */
    extent[0] = full_extent[dir_i];
    extent[1] = full_extent[dir_j];

    /* set the origin using the slice center from IOUtil */
    memset (origin, 0, sizeof (origin));
    if (have_coords)
    {
      origin[maxdir-dir-1] = myGH->sp2xyz[gdata.dim-1][dir];
    }

    /* set the direction vector */
    memset (direction, 0, sizeof (direction));
    direction[dir_i] = direction[gdata.dim + dir_j] = 1;

    mapping = Hyperslab_GlobalMappingByIndex (GH, vindex, 2,
                                              direction, origin, extent,
                                              downsample,
                                              -1,   /* table handle */
                                              NULL  /* conversion fn */,
                                              hsize);
    if (mapping < 0)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "IOSDF_Write2D: Failed to define hyperslab mapping for "
                  "variable '%s'", fullname);
      continue;
    }
    total_hsize = hsize[0] * hsize[1];
    if (total_hsize <= 0)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "IOSDF_Write2D: selected hyperslab has zero size for "
                  "variable '%s' direction %d", fullname, dir);
      Hyperslab_FreeMapping (mapping);
      continue;
    }

    /* get the bounding box information */
    bbox[0] = bbox[2] = -1;
    bbox[1] = bbox[3] = +1;
    if (have_coords)
    {
      Util_TableGetReal (coord_handles[dir_i], &minext[0], "COMPMIN");
      Util_TableGetReal (coord_handles[dir_i], &delta[0], "DELTA");
      Util_TableGetReal (coord_handles[dir_j], &minext[1], "COMPMIN");
      Util_TableGetReal (coord_handles[dir_j], &delta[1], "DELTA");

      delta[0] *= downsample[0];
      delta[1] *= downsample[1];
      bbox[0] = minext[0];
      bbox[1] = bbox[0] + delta[0]*(hsize[0]-1);
      bbox[2] = minext[1];
      bbox[3] = bbox[2] + delta[1]*(hsize[1]-1);
    }
    gft_out_set_bbox (bbox, 2);

    /* allocate hyperslab buffer */
    hdata = myproc == 0 ? malloc (total_hsize * sizeof (double)) : NULL;

    /* get the hyperslab */
    i = Hyperslab_Get (GH, mapping, 0, vindex, 0, gdata.vartype, hdata);

    /* release the mapping structure */
    Hyperslab_FreeMapping (mapping);

    /* and dump the data to file */
    if (filenames)
    {
      if (i == 0)
      {
        slice_extent[0] = hsize[0];
        slice_extent[1] = hsize[1];

        if (gft_out_brief (filenames[dir], dtime, slice_extent, 2, hdata) <= 0)
        {
          CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Error writing to 2D IOSDF output file '%s'",
                      filenames[dir]);
        }
      }
      else
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "IOSDF_Write2D: Failed to extract hyperslab for "
                    "variable '%s'", fullname);
      }

      /* clean up */
      free (hdata);

    } /* end of outputting the data by processor 0 */

  } /* end of looping through xyz directions */

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
               Opens a set of SDF files for a given alias name.
               If this is the first time through, it will advertise them
               to IOUtil.
   @enddesc

   @returntype char **
   @returndesc
               the set of full filenames
   @endreturndesc
 @@*/
static char **OpenFile (const cGH *GH,
                        const char *fullname,
                        const char *alias,
                        int dim,
                        int maxdir)
{
  int dir;
  char **filenames;
  ioSDFGH *myGH;
  ioAdvertisedFileDesc advertised_file;
  char slicename[30];
  const char *extensions[] = {"xy", "xz", "yz"};
  DECLARE_CCTK_PARAMETERS


  /* get handle for IOSDF GH extensions */
  myGH = CCTK_GHExtension (GH, "IOSDF");

  /* see if we are the first time through */
  filenames = GetNamedData (myGH->fileList_2D, alias);
  if (filenames)
  {
    return (filenames);
  }

  filenames = malloc (maxdir * sizeof (char **));

  /* open/create files for each slice */
  for (dir = 0; dir < maxdir; dir++)
  {
    filenames[dir] = malloc (strlen (myGH->out2D_dir) + strlen (alias) +
                             sizeof (slicename) + 20);

    if (dim == 2)
    {
      strcpy (slicename, "2D");
    }
    else
    {
      /* give the slice origin as range [1 .. n] */
      sprintf (slicename, "%s_%d", extensions[dir],
               myGH->sp2xyz[dim-1][dir]);
    }

    sprintf (filenames[dir], "%s%s_%s.sdf", myGH->out2D_dir, alias, slicename);

    /* advertise the file for downloading and write file info */
    advertised_file.slice = slicename;
    advertised_file.thorn = CCTK_THORNSTRING;
    advertised_file.varname = fullname;
    advertised_file.description = "Two-dimensional slice plots";
    advertised_file.mimetype = "application/sdf";

    IOUtil_AdvertiseFile (GH, filenames[dir], &advertised_file);
  }

  /* store file desriptors in database */
  StoreNamedData (&myGH->fileList_2D, alias, filenames);

  return (filenames);
}
