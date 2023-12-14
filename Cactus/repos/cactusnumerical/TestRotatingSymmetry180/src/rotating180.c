#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include <assert.h>
#include <string.h>

#define DIM(v) ((sizeof(v))/sizeof(v[0]))

enum action {init, diff};

static void do_work(int init_or_diff, CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
  
  const CCTK_REAL *xyz[3] = {x,y,z};
  const CCTK_REAL *rxyz[4] = {r,x,y,z};
  // these are the normal tensor grid functions
  const char *tensor_groups[] = {
    "TestRotatingSymmetry180::gf_u_single", "TestRotatingSymmetry180::gf_u_vector", "TestRotatingSymmetry180::gf_d_single", "TestRotatingSymmetry180::gf_d_vector", "TestRotatingSymmetry180::gf_4u_single", "TestRotatingSymmetry180::gf_4u_vector", "TestRotatingSymmetry180::gf_4d_single", "TestRotatingSymmetry180::gf_4d_vector", "TestRotatingSymmetry180::gf_uu_sym_single", "TestRotatingSymmetry180::gf_uu_sym_vector", "TestRotatingSymmetry180::gf_dd_sym_single", "TestRotatingSymmetry180::gf_dd_sym_vector", "TestRotatingSymmetry180::gf_uu_single", "TestRotatingSymmetry180::gf_uu_vector", "TestRotatingSymmetry180::gf_dd_single", "TestRotatingSymmetry180::gf_dd_vector", "TestRotatingSymmetry180::gf_du_single", "TestRotatingSymmetry180::gf_du_vector", "TestRotatingSymmetry180::gf_ud_single", "TestRotatingSymmetry180::gf_ud_vector", "TestRotatingSymmetry180::gf_ddd_sym_single", "TestRotatingSymmetry180::gf_ddd_sym_vector", "TestRotatingSymmetry180::gf_4uu_sym_single", "TestRotatingSymmetry180::gf_4uu_sym_vector", "TestRotatingSymmetry180::gf_4dd_sym_single", "TestRotatingSymmetry180::gf_4dd_sym_vector",
  };
  // scalars of various sorts
  const char *scalar_groups[] = {
    "TestRotatingSymmetry180::gf_none_single", "TestRotatingSymmetry180::gf_none_vector", "TestRotatingSymmetry180::gf_scalar_single", "TestRotatingSymmetry180::gf_scalar_vector", "TestRotatingSymmetry180::gf_4scalar_single", "TestRotatingSymmetry180::gf_4scalar_vector",
  };
  // weyl scalars
  const char *weyl_scalar_groups[] = {
    "TestRotatingSymmetry180::gf_weylscalars_real_single", "TestRotatingSymmetry180::gf_weylscalars_real_vector",
  };
  int const weylparities[10][3] =
    {{+1,+1,+1},
     {-1,-1,-1},
     {+1,+1,+1},
     {-1,-1,-1},
     {+1,+1,+1},
     {-1,-1,-1},
     {+1,+1,+1},
     {-1,-1,-1},
     {+1,+1,+1},
     {-1,-1,-1}};
  // the special case of the velocity
  const char *vel_groups[] = {
    "TestRotatingSymmetry180::gf_veld","TestRotatingSymmetry180::gf_velu",
  };
  const char *vel4_groups[] = {
    "TestRotatingSymmetry180::gf_vel4d","TestRotatingSymmetry180::gf_vel4u",
  };

  assert(cctk_nghostzones[0] == 1 && cctk_nghostzones[1] == 1 && cctk_nghostzones[2] == 1); 
  
  CCTK_INT symbnd[6];
  int ierr = GetSymmetryBoundaries(cctkGH, DIM(symbnd), symbnd);
  assert(!ierr);
  const int imin[3] = {init_or_diff == init ? 1 - (cctk_bbox[0] && !symbnd[0]) : 0, 
                       init_or_diff == init ? 1 - (cctk_bbox[2] && !symbnd[2]) : 0, 
                       init_or_diff == init ? 1 - (cctk_bbox[4] && !symbnd[4]) : 0}; 
  const int imax[3] = {cctk_lsh[0] - (init_or_diff == init ? 1 - (cctk_bbox[1] && !symbnd[1]): 0), 
                       cctk_lsh[1] - (init_or_diff == init ? 1 - (cctk_bbox[3] && !symbnd[3]): 0), 
                       cctk_lsh[2] - (init_or_diff == init ? 1 - (cctk_bbox[5] && !symbnd[5]): 0)}; 

  ptrdiff_t npoints = cctk_ash[0]*cctk_ash[1]*cctk_ash[2];

  if(init_or_diff == init)
    *num_diffs = 0;

  for(size_t g = 0 ; g < DIM(tensor_groups) ; g++) {
    const int firstvar = CCTK_FirstVarIndex(tensor_groups[g]);
    assert(firstvar >= 0);
    const int nvars = CCTK_NumVarsInGroup(tensor_groups[g]);
    assert(nvars > 0);
    cGroup groupdata;
    int ierr = CCTK_GroupData(CCTK_GroupIndex(tensor_groups[g]), &groupdata);
    assert(!ierr);
    const int vectorlength = groupdata.vectorlength;
    for(int vi = 0 ; vi < nvars ; vi++) {
      const int vectorindex = vi % vectorlength;
      const char *varname = CCTK_VarName(firstvar+vi);
      assert(varname);
      // make use of all variables ending in the tensor index they use
      const char *type = strrchr(varname, '_')+1;
      assert(type);
      const char *bracket = strchr(varname, '[');
      const size_t typelen = bracket ? bracket - type : strlen(type);
      assert(typelen <= strlen(type));
      CCTK_REAL *vardataptr = CCTK_VarDataPtrI(cctkGH, 0, firstvar+vi);
      assert(vardataptr);
      if(init_or_diff == init)
        memset(vardataptr, -1, sizeof(*vardataptr)*npoints); // initialize (to NaN)
      for(int k = imin[2] ; k < imax[2] ; k++)
      for(int j = imin[1] ; j < imax[1] ; j++)
      for(int i = imin[0] ; i < imax[0] ; i++) {
        int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
        CCTK_REAL data = vectorindex+1.; // initialize each vector of vars with different values
        for(size_t t = 0 ; t < typelen ; t++) {
          if(type[t] == 't') {
            data *= r[idx];
          } else {
            // yes this can fail on machines where xyz is not contiguous...
            assert(type[t] == 'x' || type[t] == 'y' || type[t] == 'z');
            data *= xyz[type[t]-'x'][idx];
          }
        } // for typelen
        vardataptr[idx] = init_or_diff == init ? data : vardataptr[idx] - data;
        if(init_or_diff == diff && vardataptr[idx] != 0)
          *num_diffs += 1;
      } // for npoints
    } // for nvars
  } // for tensor_groups

  for(size_t g = 0 ; g < DIM(scalar_groups) ; g++) {
    const int firstvar = CCTK_FirstVarIndex(scalar_groups[g]);
    assert(firstvar >= 0);
    const int nvars = CCTK_NumVarsInGroup(scalar_groups[g]);
    assert(nvars > 0);
    cGroup groupdata;
    int ierr = CCTK_GroupData(CCTK_GroupIndex(scalar_groups[g]), &groupdata);
    assert(!ierr);
    const int vectorlength = groupdata.vectorlength;
    for(int vi = 0 ; vi < nvars ; vi++) {
      const int vectorindex = vi % vectorlength;
      CCTK_REAL *vardataptr = CCTK_VarDataPtrI(cctkGH, 0, firstvar+vi);
      assert(vardataptr);
      if(init_or_diff == init)
        memset(vardataptr, -1, sizeof(*vardataptr)*npoints); // initialize (to NaN)
      for(int k = imin[2] ; k < imax[2] ; k++)
      for(int j = imin[1] ; j < imax[1] ; j++)
      for(int i = imin[0] ; i < imax[0] ; i++) {
        int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
        CCTK_REAL data = (vectorindex+1) * r[idx]; // initialize each vector of vars with different values
        vardataptr[idx] = init_or_diff == init ? data : vardataptr[idx] - data;
        if(init_or_diff == diff && vardataptr[idx] != 0)
          *num_diffs += 1;
      } // for npoints
    } // for nvars
  } // for scalars

  for(size_t g = 0 ; g < DIM(weyl_scalar_groups) ; g++) {
    const int firstvar = CCTK_FirstVarIndex(weyl_scalar_groups[g]);
    assert(firstvar >= 0);
    const int nvars = CCTK_NumVarsInGroup(weyl_scalar_groups[g]);
    assert(nvars > 0);
    cGroup groupdata;
    int ierr = CCTK_GroupData(CCTK_GroupIndex(weyl_scalar_groups[g]), &groupdata);
    assert(!ierr);
    const int vectorlength = groupdata.vectorlength;
    for(int vi = 0 ; vi < nvars ; vi++) {
      const int vectorindex = vi % vectorlength;
      const int index = vi / vectorlength;
      CCTK_REAL *vardataptr = CCTK_VarDataPtrI(cctkGH, 0, firstvar+vi);
      assert(vardataptr);
      if(init_or_diff == init)
        memset(vardataptr, -1, sizeof(*vardataptr)*npoints); // initialize (to NaN)
      for(int k = imin[2] ; k < imax[2] ; k++)
      for(int j = imin[1] ; j < imax[1] ; j++)
      for(int i = imin[0] ; i < imax[0] ; i++) {
        int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
        CCTK_REAL data = (vectorindex+1); // initialize each vector of vars with different values
        for(int d = 0 ; d < 3 ; d++) {
          data *= weylparities[index][d] == +1 ? r[idx] : xyz[d][idx];
      } // for npoints
        vardataptr[idx] = init_or_diff == init ? data : vardataptr[idx] - data;
        if(init_or_diff == diff && vardataptr[idx] != 0)
          *num_diffs += 1;
      } // for npoints
    } // for nvars
  } // for weyl_scalar_groups

  for(size_t g = 0 ; g < DIM(vel_groups) ; g++) {
    const int firstvar = CCTK_FirstVarIndex(vel_groups[g]);
    assert(firstvar >= 0);
    const int nvars = CCTK_NumVarsInGroup(vel_groups[g]);
    assert(nvars > 0);
    assert(nvars == 3);
    cGroup groupdata;
    int ierr = CCTK_GroupData(CCTK_GroupIndex(vel_groups[g]), &groupdata);
    assert(!ierr);
    const int vectorlength = groupdata.vectorlength;
    assert(vectorlength == 3);
    for(int vi = 0 ; vi < nvars ; vi++) {
      CCTK_REAL *vardataptr = CCTK_VarDataPtrI(cctkGH, 0, firstvar+vi);
      assert(vardataptr);
      if(init_or_diff == init)
        memset(vardataptr, -1, sizeof(*vardataptr)*npoints); // initialize (to NaN)
      for(int k = imin[2] ; k < imax[2] ; k++)
      for(int j = imin[1] ; j < imax[1] ; j++)
      for(int i = imin[0] ; i < imax[0] ; i++) {
        int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
        CCTK_REAL data = xyz[vi][idx];
        vardataptr[idx] = init_or_diff == init ? data : vardataptr[idx] - data;
        if(init_or_diff == diff && vardataptr[idx] != 0)
          *num_diffs += 1;
      } // for npoints
    } // for nvars
  } // for vel_groups

  for(size_t g = 0 ; g < DIM(vel4_groups) ; g++) {
    const int firstvar = CCTK_FirstVarIndex(vel4_groups[g]);
    assert(firstvar >= 0);
    const int nvars = CCTK_NumVarsInGroup(vel4_groups[g]);
    assert(nvars > 0);
    assert(nvars == 4);
    cGroup groupdata;
    int ierr = CCTK_GroupData(CCTK_GroupIndex(vel4_groups[g]), &groupdata);
    assert(!ierr);
    const int vectorlength = groupdata.vectorlength;
    assert(vectorlength == 4);
    for(int vi = 0 ; vi < nvars ; vi++) {
      CCTK_REAL *vardataptr = CCTK_VarDataPtrI(cctkGH, 0, firstvar+vi);
      assert(vardataptr);
      if(init_or_diff == init)
        memset(vardataptr, -1, sizeof(*vardataptr)*npoints); // initialize (to NaN)
      for(int k = imin[2] ; k < imax[2] ; k++)
      for(int j = imin[1] ; j < imax[1] ; j++)
      for(int i = imin[0] ; i < imax[0] ; i++) {
        int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
        CCTK_REAL data = rxyz[vi][idx];
        vardataptr[idx] = init_or_diff == init ? data : vardataptr[idx] - data;
        if(init_or_diff == diff && vardataptr[idx] != 0)
          *num_diffs += 1;
      } // for npoints
    } // for nvars
  } // for vel4_groups

}

void TestRotatingSymmetry180_Initialize(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  do_work(init, CCTK_PASS_CTOC);
}

void TestRotatingSymmetry180_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  const char *groups[] = {
    "TestRotatingSymmetry180::gf_none_single", "TestRotatingSymmetry180::gf_none_vector", "TestRotatingSymmetry180::gf_scalar_single", "TestRotatingSymmetry180::gf_scalar_vector", "TestRotatingSymmetry180::gf_4scalar_single", "TestRotatingSymmetry180::gf_4scalar_vector", "TestRotatingSymmetry180::gf_u_single", "TestRotatingSymmetry180::gf_u_vector", "TestRotatingSymmetry180::gf_d_single", "TestRotatingSymmetry180::gf_d_vector", "TestRotatingSymmetry180::gf_4u_single", "TestRotatingSymmetry180::gf_4u_vector", "TestRotatingSymmetry180::gf_4d_single", "TestRotatingSymmetry180::gf_4d_vector", "TestRotatingSymmetry180::gf_uu_sym_single", "TestRotatingSymmetry180::gf_uu_sym_vector", "TestRotatingSymmetry180::gf_dd_sym_single", "TestRotatingSymmetry180::gf_dd_sym_vector", "TestRotatingSymmetry180::gf_uu_single", "TestRotatingSymmetry180::gf_uu_vector", "TestRotatingSymmetry180::gf_dd_single", "TestRotatingSymmetry180::gf_dd_vector", "TestRotatingSymmetry180::gf_du_single", "TestRotatingSymmetry180::gf_du_vector", "TestRotatingSymmetry180::gf_ud_single", "TestRotatingSymmetry180::gf_ud_vector", "TestRotatingSymmetry180::gf_ddd_sym_single", "TestRotatingSymmetry180::gf_ddd_sym_vector", "TestRotatingSymmetry180::gf_4uu_sym_single", "TestRotatingSymmetry180::gf_4uu_sym_vector", "TestRotatingSymmetry180::gf_4dd_sym_single", "TestRotatingSymmetry180::gf_4dd_sym_vector", "TestRotatingSymmetry180::gf_weylscalars_real_single", "TestRotatingSymmetry180::gf_weylscalars_real_vector", "TestRotatingSymmetry180::gf_veld", "TestRotatingSymmetry180::gf_velu", "TestRotatingSymmetry180::gf_vel4d", "TestRotatingSymmetry180::gf_vel4u",
  };

  for(size_t g = 0 ; g < DIM(groups) ; g++) {
    int ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, groups[g], "none");
    assert(!ierr);
  }
}

void TestRotatingSymmetry180_Compare(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  do_work(diff, CCTK_PASS_CTOC);
}
