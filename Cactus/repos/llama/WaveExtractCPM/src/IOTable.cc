#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string>
#include <map>
#include <sys/stat.h>
#include <errno.h>
#include <iostream>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

// Provide support for outputting tabular data
// Copied from Multipole
// Eventually should be moved into its own thorn, e.g. IOTable
// ASCII and HDF5 both supported
// ASCII format should be simpler and easier to parse than CarpetIOASCII


#ifdef HAVE_CAPABILITY_HDF5
// We currently support the HDF5 1.6 API (and when using 1.8 the
// compatibility mode introduced by H5_USE_16_API).  Several machines
// in SimFactory use HDF5 1.6, so we cannot drop support for it.  It
// seems it is hard to support both the 1.6 and 1.8 API
// simultaneously; for example H5Fopen takes a different number of
// arguments in the two versions.
#define H5_USE_16_API
#include <hdf5.h>
#endif

#include "IOTable.hh"

// check return code of HDF5 call abort with an error message if there was an error.
// adapted from CarpetIOHDF5/src/CarpetIOHDF5.hh.
#define HDF5_ERROR(fn_call)                                                   \
        do {                                                                  \
          hid_t _error_code = fn_call;                                        \
                                                                              \
                                                                              \
          if (_error_code < 0)                                                \
          {                                                                   \
            CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,              \
                        "HDF5 call '%s' returned error code %d",              \
                        #fn_call, (int)_error_code);                          \
          }                                                                   \
        } while (0)


using namespace std;

////////////////////////////////////////////////////////////////////////
// File manipulation
////////////////////////////////////////////////////////////////////////

FILE *WaveExtractCPM_OpenOutputFile(CCTK_ARGUMENTS, const char *name_p)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  string name(name_p);
  
  bool first_time = cctk_iteration == 0;
  const char *mode = first_time ? "w" : "a";

  string output_name(string(io_out_dir) + "/" + name);

  FILE *fp = fopen(output_name.c_str(), mode);

  if (fp == 0)
  {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, "%s", (string("Could not open output file ") + output_name).c_str());
  }

  return fp;
}

////////////////////////////////////////////////////////////////////////
// ASCII IO
////////////////////////////////////////////////////////////////////////

void WaveExtractCPM_OutputTableRowASCII(CCTK_ARGUMENTS, const char *name_p, int npoints, CCTK_REAL data[])
{
  DECLARE_CCTK_ARGUMENTS;
  const string name(name_p);

  if (FILE *fp = WaveExtractCPM_OpenOutputFile(CCTK_PASS_CTOC, name.c_str()))
  {
    for (int i=0; i < npoints; i++)
    {
      fprintf(fp, "%.19g", data[i]);
      if (i < npoints-1)
      {
        fprintf(fp, "\t");
      }
    }
    fprintf(fp, "\n");
    fclose(fp);
  }
}

////////////////////////////////////////////////////////////////////////
// HDF5 table output
////////////////////////////////////////////////////////////////////////

#ifdef HAVE_CAPABILITY_HDF5

static bool file_exists(const string &name)
{
  struct stat sts;
  return !(stat(name.c_str(), &sts) == -1 && errno == ENOENT);
}

static bool dataset_exists(hid_t file, const string &dataset_name)
{
  // To test whether a dataset exists, the recommended way in API 1.6
  // is to use H5Gget_objinfo, but this prints an error to stderr if
  // the dataset does not exist.  We explicitly avoid this by wrapping
  // the call in H5E_BEGIN_TRY/H5E_END_TRY statements.  In 1.8,
  // H5Gget_objinfo is deprecated, and H5Lexists does the job.  See
  // http://www.mail-archive.com/hdf-forum@hdfgroup.org/msg00125.html

  #if 1
  bool exists;
  H5E_BEGIN_TRY
  {
    exists = H5Gget_objinfo(file, dataset_name.c_str(), 1, NULL) >= 0;
  }
  H5E_END_TRY;
  return exists;
  #else
  return H5Lexists(file, dataset_name.c_str(), H5P_DEFAULT);
  #endif
}

extern "C"
void WaveExtractCPM_OutputTableRowHDF5(CCTK_ARGUMENTS, const char *basename_p,
                                       const char *datasetname_p,
                                       unsigned int npoints, CCTK_REAL data[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  string basename(basename_p);
  string datasetname(datasetname_p);

  if (strlen(io_out_dir) == 0)
  {
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "IO::out_dir is empty; cannot write HDF5 file"
                );
  }
  
  string output_name = string(io_out_dir) + string("/") + basename;
  static map<string,bool> checked; // Has the given file been checked
                                   // for truncation? map<*,bool>
                                   // defaults to false

  hid_t file;

  if (!file_exists(output_name) || (!checked[output_name] && IO_TruncateOutputFiles(cctkGH)))
  {
    file = H5Fcreate(output_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }
  else
  {
    file = H5Fopen(output_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  }

  checked[output_name] = true;

  hid_t dataset = -1;

  if (dataset_exists(file, datasetname))
  {
    dataset = H5Dopen(file, datasetname.c_str());
  }
  else
  {
    hsize_t dims[2] = {0,npoints};
    hsize_t maxdims[2] = {H5S_UNLIMITED, npoints};
    hid_t dataspace = H5Screate_simple(2, dims, maxdims);

    hid_t cparms = -1;
    hsize_t chunk_dims[2] = {(hsize_t)hdf5_chunk_size,npoints};
    cparms = H5Pcreate (H5P_DATASET_CREATE);
    HDF5_ERROR(H5Pset_chunk(cparms, 2, chunk_dims));

    dataset = H5Dcreate(file, datasetname.c_str(), H5T_NATIVE_DOUBLE, dataspace, cparms);
    H5Pclose(cparms);
  }

  hid_t filespace = H5Dget_space (dataset);
  hsize_t filedims[2];
  hsize_t maxdims[2];
  HDF5_ERROR(H5Sget_simple_extent_dims(filespace, filedims, maxdims));

  if (filedims[1] != npoints)
  {
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot change the number of columns in dataset %s in file %s from %d to %d",
                datasetname.c_str(), output_name.c_str(), (int) filedims[1], npoints);
  }
  
  filedims[0] += 1;
  hsize_t size[2] = {filedims[0], filedims[1]};
  HDF5_ERROR(H5Dextend (dataset, size));
  HDF5_ERROR(H5Sclose(filespace));

  /* Select a hyperslab  */
  hsize_t offset[2] = {filedims[0]-1, 0};
  hsize_t dims2[2] = {1,npoints};
  filespace = H5Dget_space (dataset);
  HDF5_ERROR(H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL,
                                  dims2, NULL));

  hid_t memdataspace = H5Screate_simple(2, dims2, NULL);

  /* Write the data to the hyperslab  */
  HDF5_ERROR(H5Dwrite (dataset, H5T_NATIVE_DOUBLE, memdataspace, filespace,
                       H5P_DEFAULT, data));

  HDF5_ERROR(H5Dclose(dataset));
  HDF5_ERROR(H5Sclose(filespace));
  HDF5_ERROR(H5Sclose(memdataspace));

  HDF5_ERROR(H5Fclose(file));
}

#else

void WaveExtractCPM_OutputArrayToH5File(CCTK_ARGUMENTS, const string &basename, const string &datasetname,
                                        int npoints, CCTK_REAL data[])
{
  CCTK_WARN(0,"HDF5 output has been requested but Cactus has been compiled without HDF5 support");
}

#endif

