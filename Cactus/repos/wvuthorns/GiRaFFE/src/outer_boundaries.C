/*******************************************************
 * Outer boundaries are handled as follows:
 * (-1) Update RHS quantities, leave RHS quantities zero on all outer ghostzones (including outer AMR refinement, processor, and outer boundaries)
 * ( 0) Let MoL update all evolution variables
 * ( 1) Apply outer boundary conditions (BCs) on A_{\mu}
 * ( 2) Compute B^i from A_i everywhere, synchronize B^i
 * ( 3) Call con2prim to get primitives on interior pts
 * ( 4) Apply outer BCs on {vx,vy,vz}.
 * ( 5) (optional) set conservatives on outer boundary.
 *******************************************************/

#include "cctk.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "GiRaFFE_headers.h"
#include "inlined_functions.C"

#define IDX(i,j,k) CCTK_GFINDEX3D(cctkGH,(i),(j),(k))

#define XMAX_OB_LINEAR_EXTRAP(FUNC,imax) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) FUNC[IDX(imax,j,k)] = 2.0 * FUNC[IDX(imax-1,j,k)] - FUNC[IDX(imax-2,j,k)];
#define YMAX_OB_LINEAR_EXTRAP(FUNC,jmax) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,jmax,k)] = 2.0 * FUNC[IDX(i,jmax-1,k)] - FUNC[IDX(i,jmax-2,k)];
#define ZMAX_OB_LINEAR_EXTRAP(FUNC,kmax) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,j,kmax)] = 2.0 * FUNC[IDX(i,j,kmax-1)] - FUNC[IDX(i,j,kmax-2)];

#define XMIN_OB_LINEAR_EXTRAP(FUNC,imin) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) FUNC[IDX(imin,j,k)] = 2.0 * FUNC[IDX(imin+1,j,k)] - FUNC[IDX(imin+2,j,k)];
#define YMIN_OB_LINEAR_EXTRAP(FUNC,jmin) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,jmin,k)] = 2.0 * FUNC[IDX(i,jmin+1,k)] - FUNC[IDX(i,jmin+2,k)];
#define ZMIN_OB_LINEAR_EXTRAP(FUNC,kmin) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,j,kmin)] = 2.0 * FUNC[IDX(i,j,kmin+1)] - FUNC[IDX(i,j,kmin+2)];

/*********************************************
 * Apply outer boundary conditions on A_{\mu}
 ********************************************/
extern "C" void GiRaFFE_outer_boundaries_on_A_mu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_EQUALS(EM_BC,"frozen")) {
    CCTK_VInfo(CCTK_THORNSTRING,"Warning: EM_BC=frozen WILL NOT NECESSARILY result in frozen B-fields at the boundaries, since A-fields may evolve even with B fixed in time!\n");
    return;
  }

  bool Symmetry_none=false; if(CCTK_EQUALS(Symmetry,"none")) Symmetry_none=true;

  const int levelnumber = GetRefinementLevel(cctkGH);

  GiRaFFE_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij(cctkGH,cctk_lsh,  gxx,gxy,gxz,gyy,gyz,gzz,alp,
                                                                        gtxx,gtxy,gtxz,gtyy,gtyz,gtzz,
                                                                        gtupxx,gtupxy,gtupxz,gtupyy,gtupyz,gtupzz,
                                                                        phi_bssn,psi_bssn,lapm1);

  // Don't apply approximate outer boundary conditions on initial data, which should be defined everywhere, or on levels != [coarsest level].
  if(cctk_iteration==0 || levelnumber!=0) return;

  if(cctk_nghostzones[0]!=cctk_nghostzones[1] || cctk_nghostzones[0]!=cctk_nghostzones[2])
    CCTK_VError(VERR_DEF_PARAMS,"ERROR: GiRaFFE outer BC driver does not support unequal number of ghostzones in different directions!");
  for(int which_bdry_pt=0;which_bdry_pt<cctk_nghostzones[0];which_bdry_pt++) {
    const int imax=cctk_lsh[0]-cctk_nghostzones[0]+which_bdry_pt; // for cctk_nghostzones==3, this goes {cctk_lsh-3,cctk_lsh-2,cctk_lsh-1}; outer bdry pt is at cctk_lsh-1
    const int jmax=cctk_lsh[1]-cctk_nghostzones[1]+which_bdry_pt;
    const int kmax=cctk_lsh[2]-cctk_nghostzones[2]+which_bdry_pt;

    const int imin=cctk_nghostzones[0]-which_bdry_pt-1; // for cctk_nghostzones==3, this goes {2,1,0}
    const int jmin=cctk_nghostzones[1]-which_bdry_pt-1;
    const int kmin=cctk_nghostzones[2]-which_bdry_pt-1;

    if(cctk_bbox[1]) { XMAX_OB_LINEAR_EXTRAP(Ax,imax); XMAX_OB_LINEAR_EXTRAP(Ay,imax); XMAX_OB_LINEAR_EXTRAP(Az,imax); XMAX_OB_LINEAR_EXTRAP(psi6phi,imax); }
    if(cctk_bbox[3]) { YMAX_OB_LINEAR_EXTRAP(Ax,jmax); YMAX_OB_LINEAR_EXTRAP(Ay,jmax); YMAX_OB_LINEAR_EXTRAP(Az,jmax); YMAX_OB_LINEAR_EXTRAP(psi6phi,jmax); }
    if(cctk_bbox[5]) { ZMAX_OB_LINEAR_EXTRAP(Ax,kmax); ZMAX_OB_LINEAR_EXTRAP(Ay,kmax); ZMAX_OB_LINEAR_EXTRAP(Az,kmax); ZMAX_OB_LINEAR_EXTRAP(psi6phi,kmax); }

    if(cctk_bbox[0]) {                    XMIN_OB_LINEAR_EXTRAP(Ax,imin); XMIN_OB_LINEAR_EXTRAP(Ay,imin); XMIN_OB_LINEAR_EXTRAP(Az,imin); XMIN_OB_LINEAR_EXTRAP(psi6phi,imin); }
    if(cctk_bbox[2]) {                    YMIN_OB_LINEAR_EXTRAP(Ax,jmin); YMIN_OB_LINEAR_EXTRAP(Ay,jmin); YMIN_OB_LINEAR_EXTRAP(Az,jmin); YMIN_OB_LINEAR_EXTRAP(psi6phi,jmin); }
    if((cctk_bbox[4]) && Symmetry_none) { ZMIN_OB_LINEAR_EXTRAP(Ax,kmin); ZMIN_OB_LINEAR_EXTRAP(Ay,kmin); ZMIN_OB_LINEAR_EXTRAP(Az,kmin); ZMIN_OB_LINEAR_EXTRAP(psi6phi,kmin); }
  }
}

#define XMAX_OB_SIMPLE_COPY(FUNC,imax) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) FUNC[IDX(imax,j,k)] = FUNC[IDX(imax-1,j,k)];
#define YMAX_OB_SIMPLE_COPY(FUNC,jmax) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,jmax,k)] = FUNC[IDX(i,jmax-1,k)];
#define ZMAX_OB_SIMPLE_COPY(FUNC,kmax) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,j,kmax)] = FUNC[IDX(i,j,kmax-1)];

#define XMIN_OB_SIMPLE_COPY(FUNC,imin) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) FUNC[IDX(imin,j,k)] = FUNC[IDX(imin+1,j,k)];
#define YMIN_OB_SIMPLE_COPY(FUNC,jmin) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,jmin,k)] = FUNC[IDX(i,jmin+1,k)];
#define ZMIN_OB_SIMPLE_COPY(FUNC,kmin) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,j,kmin)] = FUNC[IDX(i,j,kmin+1)];


#define XMAX_INFLOW_CHECK(vx,imax) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) if(vx[IDX(imax,j,k)]<0.) vx[IDX(imax,j,k)]=0.;
#define YMAX_INFLOW_CHECK(vy,jmax) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) if(vy[IDX(i,jmax,k)]<0.) vy[IDX(i,jmax,k)]=0.;
#define ZMAX_INFLOW_CHECK(vz,kmax) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) if(vz[IDX(i,j,kmax)]<0.) vz[IDX(i,j,kmax)]=0.;

#define XMIN_INFLOW_CHECK(vx,imin) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) if(vx[IDX(imin,j,k)]>0.) vx[IDX(imin,j,k)]=0.;
#define YMIN_INFLOW_CHECK(vy,jmin) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) if(vy[IDX(i,jmin,k)]>0.) vy[IDX(i,jmin,k)]=0.;
#define ZMIN_INFLOW_CHECK(vz,kmin) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) if(vz[IDX(i,j,kmin)]>0.) vz[IDX(i,j,kmin)]=0.;


/*******************************************************
 * Apply outer boundary conditions on {vx,vy,vz}
 * It is better to apply BCs on primitives than conservs,
 * because small errors in conservs can be greatly
 * amplified in con2prim, sometimes leading to unphysical
 * primitives & unnecessary fixes.
 *******************************************************/
extern "C" void GiRaFFE_outer_boundaries_on_vx_vy_vz(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_EQUALS(Velocity_BC,"frozen")) return;

  bool Symmetry_none=false; if(CCTK_EQUALS(Symmetry,"none")) Symmetry_none=true;

  const int levelnumber = GetRefinementLevel(cctkGH);

  // Don't apply approximate outer boundary conditions on initial data, which should be defined everywhere, or on levels != [coarsest level].
  if(cctk_iteration==0 || levelnumber!=0) return;

  int OUTFLOW=1;
  if(CCTK_EQUALS(Velocity_BC,"copy")) OUTFLOW=0;

  GiRaFFE_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij(cctkGH,cctk_lsh,  gxx,gxy,gxz,gyy,gyz,gzz,alp,
                                                                        gtxx,gtxy,gtxz,gtyy,gtyz,gtzz,
                                                                        gtupxx,gtupxy,gtupxz,gtupyy,gtupyz,gtupzz,
                                                                        phi_bssn,psi_bssn,lapm1);

  if(cctk_nghostzones[0]!=cctk_nghostzones[1] || cctk_nghostzones[0]!=cctk_nghostzones[2])
    CCTK_VError(VERR_DEF_PARAMS,"ERROR: GiRaFFE outer BC driver does not support unequal number of ghostzones in different directions!");
  for(int which_bdry_pt=0;which_bdry_pt<cctk_nghostzones[0];which_bdry_pt++) {
    const int imax=cctk_lsh[0]-cctk_nghostzones[0]+which_bdry_pt; // for cctk_nghostzones==3, this goes {cctk_lsh-3,cctk_lsh-2,cctk_lsh-1}; outer bdry pt is at cctk_lsh-1
    const int jmax=cctk_lsh[1]-cctk_nghostzones[1]+which_bdry_pt;
    const int kmax=cctk_lsh[2]-cctk_nghostzones[2]+which_bdry_pt;

    const int imin=cctk_nghostzones[0]-which_bdry_pt-1; // for cctk_nghostzones==3, this goes {2,1,0}
    const int jmin=cctk_nghostzones[1]-which_bdry_pt-1;
    const int kmin=cctk_nghostzones[2]-which_bdry_pt-1;

    // Order here is for compatibility with old version of this code.
    /* XMIN & XMAX */
    // i=imax=outer boundary
    if(cctk_bbox[1]) { XMAX_OB_SIMPLE_COPY(vx,imax); XMAX_OB_SIMPLE_COPY(vy,imax); XMAX_OB_SIMPLE_COPY(vz,imax); if(OUTFLOW) XMAX_INFLOW_CHECK(vx,imax); }
    // i=imin=outer boundary
    if(cctk_bbox[0]) {
      XMIN_OB_SIMPLE_COPY(vx,imin); XMIN_OB_SIMPLE_COPY(vy,imin); XMIN_OB_SIMPLE_COPY(vz,imin); if(OUTFLOW) XMIN_INFLOW_CHECK(vx,imin); }

    /* YMIN & YMAX */
    // j=jmax=outer boundary
    if(cctk_bbox[3]) { YMAX_OB_SIMPLE_COPY(vx,jmax); YMAX_OB_SIMPLE_COPY(vy,jmax); YMAX_OB_SIMPLE_COPY(vz,jmax); if(OUTFLOW) YMAX_INFLOW_CHECK(vy,jmax); }
    // j=jmin=outer boundary
    if(cctk_bbox[2]) {
      YMIN_OB_SIMPLE_COPY(vx,jmin); YMIN_OB_SIMPLE_COPY(vy,jmin); YMIN_OB_SIMPLE_COPY(vz,jmin); if(OUTFLOW) YMIN_INFLOW_CHECK(vy,jmin); }

    /* ZMIN & ZMAX */
    // k=kmax=outer boundary
    if(cctk_bbox[5]) { ZMAX_OB_SIMPLE_COPY(vx,kmax); ZMAX_OB_SIMPLE_COPY(vy,kmax); ZMAX_OB_SIMPLE_COPY(vz,kmax); if(OUTFLOW) ZMAX_INFLOW_CHECK(vz,kmax); }
    // k=kmin=outer boundary
    if((cctk_bbox[4]) && Symmetry_none) {
      ZMIN_OB_SIMPLE_COPY(vx,kmin); ZMIN_OB_SIMPLE_COPY(vy,kmin); ZMIN_OB_SIMPLE_COPY(vz,kmin); if(OUTFLOW) ZMIN_INFLOW_CHECK(vz,kmin); }
  }

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        if(((cctk_bbox[0]) && i<cctk_nghostzones[0]) ||
           ((cctk_bbox[1]) && i>=cctk_lsh[0]-cctk_nghostzones[0]) ||
           ((cctk_bbox[2]) && j<cctk_nghostzones[1]) ||
           ((cctk_bbox[3]) && j>=cctk_lsh[1]-cctk_nghostzones[1]) ||
           ((cctk_bbox[4]) && k<cctk_nghostzones[2] && CCTK_EQUALS(Symmetry,"none")) ||
           ((cctk_bbox[5]) && k>=cctk_lsh[2]-cctk_nghostzones[2])) {
          const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
          int ww;

          CCTK_REAL METRIC[NUMVARS_FOR_METRIC],dummy=-1e300;
          ww=0;
          //psi[index] = exp(phi[index]);
          METRIC[ww] = phi_bssn[index];ww++;
          METRIC[ww] = dummy;          ww++; // Don't need to set psi.
          METRIC[ww] = gtxx[index];    ww++;
          METRIC[ww] = gtxy[index];    ww++;
          METRIC[ww] = gtxz[index];    ww++;
          METRIC[ww] = gtyy[index];    ww++;
          METRIC[ww] = gtyz[index];    ww++;
          METRIC[ww] = gtzz[index];    ww++;
          METRIC[ww] = lapm1[index];   ww++;
          METRIC[ww] = betax[index];   ww++;
          METRIC[ww] = betay[index];   ww++;
          METRIC[ww] = betaz[index];   ww++;
          METRIC[ww] = gtupxx[index];  ww++;
          METRIC[ww] = gtupyy[index];  ww++;
          METRIC[ww] = gtupzz[index];  ww++;
          METRIC[ww] = gtupxy[index];  ww++;
          METRIC[ww] = gtupxz[index];  ww++;
          METRIC[ww] = gtupyz[index];  ww++;

          CCTK_REAL PRIMS[MAXNUMVARS];
          ww=0;
          PRIMS[ww] = vx[index];    ww++;
          PRIMS[ww] = vy[index];    ww++;
          PRIMS[ww] = vz[index];    ww++;
          PRIMS[ww] = Bx[index];    ww++;
          PRIMS[ww] = By[index];    ww++;
          PRIMS[ww] = Bz[index];    ww++;

          CCTK_REAL CONSERVS[NUM_CONSERVS];

          GiRaFFE_compute_conservatives(PRIMS,METRIC, CONSERVS);

          vx[index]    = PRIMS[VX];
          vy[index]    = PRIMS[VY];
          vz[index]    = PRIMS[VZ];

          mhd_st_x[index]=CONSERVS[STILDEX];
          mhd_st_y[index]=CONSERVS[STILDEY];
          mhd_st_z[index]=CONSERVS[STILDEZ];
        }
      }
}

#include "compute_conservatives_FFE.C"
