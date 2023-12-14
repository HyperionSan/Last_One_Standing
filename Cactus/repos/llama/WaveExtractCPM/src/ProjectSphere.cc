#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include <vector>
#include "carpetinterp2.hh"

using namespace std;
using namespace CarpetInterp2;

extern "C" void WavExtrCPM_MetricToSphericalBasis(CCTK_ARGUMENTS);

struct interp_setup_t_CPM {
  fasterp_setup_t * fasterp_setup;
  int npoints;

  interp_setup_t_CPM (int npoints_) : fasterp_setup (NULL), npoints (npoints_) {}
  ~interp_setup_t_CPM () { if (fasterp_setup) delete fasterp_setup; }
};
vector<interp_setup_t_CPM *> interp_setups_CPM;   // has to be the same number as there are detectors!


#define SINDEX(lsh, i, j) j*lsh[0]+i

extern "C" void WavExtrCPM_ProjectSphere(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   
   int ierr = 0;
   int di;
   int lsh[2], lbnd[2];
   CCTK_REAL dtheta, dphi, dthetainv, dphiinv;
   
   if (verbose>2) CCTK_INFO("Interpolating metric and derivatives onto sphere");

   if (*do_nothing == 1) return;

   if (cctk_iteration != 0)
   {
       if (cctk_iteration % my_out_every_det[*current_detector-1] != 0)
       {
         if (verbose>2) CCTK_INFO("No time for this detector");
         return;
      }
   }

   *current_detector_radius = detector_radius[*current_detector-1];
   di=*current_detector-1;

   // Sanity Check
   if (*current_detector_radius < 1.e-10)
   {
      *do_nothing = 1;
      CCTK_WARN(1,"This should never happen: The detector radius is 0!");
      return;
   }

   if (calc_when_necessary == 1)
   {
      if (cctk_time < *current_detector_radius-50)
      {
         if (verbose>2) CCTK_INFO("No time for this detector");
         return;
      }
      int istat = CCTK_IsFunctionAliased("MergerHandler_WeHaveMerger");
      if (istat == 1)
      {
         if (MergerHandler_WeHaveMerger() == 1)
         {
            if (cctk_time > MergerHandler_MergerTime()+*current_detector_radius+ringdown_margin)
            {
               if (verbose>2) CCTK_INFO("No time for this detector");
               return;
            }
         }
      }
   }
    

   if (verbose>1)
      cout << "Analysing Detector No.: " << *current_detector << " Radius " << *current_detector_radius << endl;



   // local shape of the 2D grid arrays
   ierr = CCTK_GrouplshGN(cctkGH, 2, lsh, "WaveExtractCPM::surface_arrays");
   if ( ierr < 0 )
      CCTK_WARN(0, "cannot get local size for surface arrays");

   //if (size(ctheta,1)<2) call CCTK_WARN (0, "internal error")
   //if (size(cphi,2)<2) call CCTK_WARN (0, "internal error")
   dtheta = ctheta[SINDEX(lsh,1,0)] - ctheta[SINDEX(lsh,0,0)];
   dphi = cphi[SINDEX(lsh,0,1)] - cphi[SINDEX(lsh,0,0)];

   
//   int number_of_vars = 12;
// SHH CHANGE
    int number_of_vars = 27;
    
    vector<string> varnames(number_of_vars);
   
    varnames[0] = "admbase::gxx";
    varnames[1] = "admbase::gxy";
    varnames[2] = "admbase::gxz";
    varnames[3] = "admbase::gyy";
    varnames[4] = "admbase::gyz";
    varnames[5] = "admbase::gzz";
    varnames[6] = "admderivatives::gxx_dr";
    varnames[7] = "admderivatives::gxy_dr";
    varnames[8] = "admderivatives::gxz_dr";
    varnames[9] = "admderivatives::gyy_dr";
    varnames[10] = "admderivatives::gyz_dr";
    varnames[11] = "admderivatives::gzz_dr";
// SHH CHANGE
    varnames[12] = "admbase::betax";
    varnames[13] = "admbase::betay";
    varnames[14] = "admbase::betaz";
    varnames[15] = "admbase::dtbetax";
    varnames[16] = "admbase::dtbetay";
    varnames[17] = "admbase::dtbetaz";

    varnames[18] = "admderivatives::gxx_dt";
    varnames[19] = "admderivatives::gxy_dt";
    varnames[20] = "admderivatives::gxz_dt";
    varnames[21] = "admderivatives::gyy_dt";
    varnames[22] = "admderivatives::gyz_dt";
    varnames[23] = "admderivatives::gzz_dt";

    varnames[24] = "admderivatives::betax_dr";
    varnames[25] = "admderivatives::betay_dr";
    varnames[26] = "admderivatives::betaz_dr";
    

   // get the input variable's index
   vector<int> varindices (number_of_vars, -1);
   
   for (int i=0; i < number_of_vars; ++i)
   {
      varindices[i] = CCTK_VarIndex(varnames[i].c_str());
      if (varindices[i] < 0)
         CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "couldn't get index of slice variable '%s'", varnames[i].c_str());
   }
   
   vector<CCTK_REAL*> values (number_of_vars, static_cast<CCTK_REAL*>(NULL));
   
   // set output pointers
   values[0] = gxxi;
   values[1] = gxyi;
   values[2] = gxzi;
   values[3] = gyyi;
   values[4] = gyzi;
   values[5] = gzzi;
   values[6] = dr_gxxi;
   values[7] = dr_gxyi;
   values[8] = dr_gxzi;
   values[9] = dr_gyyi;
   values[10] = dr_gyzi;
   values[11] = dr_gzzi;
    
    // SHH CHANGE
    values[12] = gtxi;
    values[13] = gtyi;
    values[14] = gtzi;
    values[15] = dt_gtxi;
    values[16] = dt_gtyi;
    values[17] = dt_gtzi;
    
    values[18] = dt_gxxi;
    values[19] = dt_gxyi;
    values[20] = dt_gxzi;
    values[21] = dt_gyyi;
    values[22] = dt_gyzi;
    values[23] = dt_gzzi;
    
    values[24] = dr_gtxi;
    values[25] = dr_gtyi;
    values[26] = dr_gtzi;
   
   int interpolator_order = interpolation_stencil;
   
   // set up the interpolation
   interp_setups_CPM.resize(maximum_detector_number, NULL);
   assert(interp_setups_CPM.size() > di);
   interp_setup_t_CPM* &interp_setup = interp_setups_CPM[di];
   if (not interp_setup or interp_setup->fasterp_setup->outofdate()) {
      if (interp_setup)
        delete interp_setup;
      interp_setup = new interp_setup_t_CPM(lsh[0]*lsh[1]);

      // allocate storage for coordinates
      fasterp_glocs_t locations (interp_setup->npoints);

      // get Cartesian coordinate values of gridpoints on spherical surface
      CCTK_REAL rad = *current_detector_radius;
      for (int ip=0; ip < lsh[1]; ++ip) {
         for (int it=0; it < lsh[0]; ++it) {
            const int l = SINDEX(lsh, it, ip);
            locations.coords[0][l] = origin_x + rad * sintheta[l] * cosphi[l];
            locations.coords[1][l] = origin_y + rad * sintheta[l] * sinphi[l];
            locations.coords[2][l] = origin_z + rad * costheta[l];
         }
      }

      interp_setup->fasterp_setup =
         new fasterp_setup_t(cctkGH, locations, interpolator_order);
   }

   // do the interpolation
   assert(interp_setup->fasterp_setup);
   assert(interp_setup->npoints == lsh[0]*lsh[1]);
   interp_setup->fasterp_setup->interpolate(cctkGH, varindices, values);
}

extern "C" void WavExtrCPM_MetricToSphericalBasis(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS

   // Convert to spherical coord. system
   // FIXME : WARN, THESE EQUATIONS PROBABLY CANT WORK FOR PHYSICAL METRIC
   CCTK_REAL rad = *current_detector_radius;
   CCTK_REAL r2=rad*rad;
   int lsh[2], lbnd[2];
   int ierr = 0;

   ierr = CCTK_GrouplshGN(cctkGH, 2, lsh, "WaveExtractCPM::surface_arrays");
   if ( ierr < 0 )
      CCTK_WARN(0, "cannot get local size for surface arrays");

   // Find out the lower bounds of the distributed integration grid arrays.
   ierr = CCTK_GrouplbndGN(cctkGH,2,lbnd,"WaveExtractCPM::surface_integrands");
   if ( ierr < 0 )
      CCTK_WARN(0, "cannot get lower bounds for surface integrands");

   
   for (int ip = 0; ip < lsh[1]; ++ip)
   {
      for (int it = 0; it < lsh[0]; ++it)
      {
         const int l = SINDEX(lsh, it, ip);
         
         CCTK_REAL costh  = costheta[l];
         CCTK_REAL costh2 = costh*costh;
         CCTK_REAL sinth  = sintheta[l];
         CCTK_REAL sinth2 = sinth*sinth;
         CCTK_REAL cosph  = cosphi[l];
         CCTK_REAL cosph2 = cosph*cosph;
         CCTK_REAL sinph  = sinphi[l];
         CCTK_REAL sinph2 = sinph*sinph;
         
          CCTK_REAL one = 1.0;
          CCTK_REAL two = 2.0;
          CCTK_REAL four = 4.0;
         
          CCTK_REAL g_tx    = gtxi[l];
          CCTK_REAL g_ty    = gtyi[l];
          CCTK_REAL g_tz    = gtzi[l];
          CCTK_REAL g_xx    = gxxi[l];
          CCTK_REAL g_xy    = gxyi[l];
          CCTK_REAL g_xz    = gxzi[l];
          CCTK_REAL g_yy    = gyyi[l];
          CCTK_REAL g_yz    = gyzi[l];
          CCTK_REAL g_zz    = gzzi[l];
          
          CCTK_REAL dr_g_tx = dr_gtxi[l];
          CCTK_REAL dr_g_ty = dr_gtyi[l];
          CCTK_REAL dr_g_tz = dr_gtzi[l];

          CCTK_REAL dr_g_xx = dr_gxxi[l];
          CCTK_REAL dr_g_xy = dr_gxyi[l];
          CCTK_REAL dr_g_xz = dr_gxzi[l];
          CCTK_REAL dr_g_yy = dr_gyyi[l];
          CCTK_REAL dr_g_yz = dr_gyzi[l];
          CCTK_REAL dr_g_zz = dr_gzzi[l];

          CCTK_REAL dt_g_tx = dt_gtxi[l];
          CCTK_REAL dt_g_ty = dt_gtyi[l];
          CCTK_REAL dt_g_tz = dt_gtzi[l];
          CCTK_REAL dt_g_xx = dt_gxxi[l];
          CCTK_REAL dt_g_xy = dt_gxyi[l];
          CCTK_REAL dt_g_xz = dt_gxzi[l];
          CCTK_REAL dt_g_yy = dt_gyyi[l];
          CCTK_REAL dt_g_yz = dt_gyzi[l];
          CCTK_REAL dt_g_zz = dt_gzzi[l];


          gtr[l] = g_tz*costh + g_tx*cosph*sinth + g_ty*sinth*sinph;
          
          gtth[l] = rad*(g_tx*costh*cosph - g_tz*sinth + g_ty*costh*sinph);
          
          gtphi[l] = rad*(g_ty*cosph*sinth - g_tx*sinth*sinph);
          
          grr[l] = g_zz*costh2 + two*g_xz*costh*cosph*sinth + g_xx*cosph2*sinth2 + two*g_yz*costh*sinth*sinph + two*g_xy*cosph*sinth2*sinph + g_yy*sinth2*sinph2;
          

          grth[l] = rad*(-(g_zz*costh*sinth) + g_xx*costh*cosph2*sinth + g_xz*cosph*(costh2 - sinth2) + two*g_xy*costh*cosph*sinth*sinph + g_yz*(costh2 - sinth2)*sinph + g_yy*costh*sinth*sinph2);
          

          grphi[l] = rad*(g_yz*costh*cosph*sinth - g_xz*costh*sinth*sinph - g_xx*cosph*sinth2*sinph + g_yy*cosph*sinth2*sinph + g_xy*sinth2*(cosph2 - sinph2));
          
          gthth[l] = r2*(g_xx*costh2*cosph2 - two*g_xz*costh*cosph*sinth + g_zz*sinth2 + two*g_xy*costh2*cosph*sinph - two*g_yz*costh*sinth*sinph + g_yy*costh2*sinph2);
          
          gthphi[l] = r2*(-(g_yz*cosph*sinth2) - g_xx*costh*cosph*sinth*sinph + g_yy*costh*cosph*sinth*sinph + g_xz*sinth2*sinph + g_xy*costh*sinth*(cosph2 - sinph2));
          
          gphiphi[l] = r2*(g_yy*cosph2*sinth2 - two*g_xy*cosph*sinth2*sinph + g_xx*sinth2*sinph2);
          
          dr_gtth[l] = (one/rad) * gtth[l]
                      +
                      rad*(dr_g_tx*costh*cosph - dr_g_tz*sinth + dr_g_ty*costh*sinph);
          
          dr_gtphi[l] = (one/rad) * gtphi[l]
                        +
                        rad*(dr_g_ty*cosph*sinth - dr_g_tx*sinth*sinph);
                   
          dr_gthth[l] = (two/rad) * gthth[l]
                        +
                        r2*(dr_g_xx*costh2*cosph2
                            - two*dr_g_xz*costh*cosph*sinth
                            + dr_g_zz*sinth2
                            + two*dr_g_xy*costh2*cosph*sinph
                            - two*dr_g_yz*costh*sinth*sinph
                            + dr_g_yy*costh2*sinph2);
          
          dr_gthphi[l] =  (two/rad) * gthphi[l]
                        +
                        r2*(-(dr_g_yz*cosph*sinth2)
                            - dr_g_xx*costh*cosph*sinth*sinph
                            + dr_g_yy*costh*cosph*sinth*sinph
                            + dr_g_xz*sinth2*sinph
                            + dr_g_xy*costh*sinth*(cosph2 - sinph2));
          
          dr_gphiphi[l] = (two/rad) * gphiphi[l]
                        +
                        r2*(dr_g_yy*cosph2*sinth2
                            - two*dr_g_xy*cosph*sinth2*sinph
                            + dr_g_xx*sinth2*sinph2);
          
                    
          dt_grth[l] = rad*(-(dt_g_zz*costh*sinth) + dt_g_xx*costh*cosph2*sinth + dt_g_xz*cosph*(costh2 - sinth2) + two*dt_g_xy*costh*cosph*sinth*sinph + dt_g_yz*(costh2 - sinth2)*sinph + dt_g_yy*costh*sinth*sinph2);
          
          dt_grphi[l] = rad*(dt_g_yz*costh*cosph*sinth - dt_g_xz*costh*sinth*sinph - dt_g_xx*cosph*sinth2*sinph + dt_g_yy*cosph*sinth2*sinph + dt_g_xy*sinth2*(cosph2 - sinph2));
          
          dt_gthth[l] = r2*(dt_g_xx*costh2*cosph2 - two*dt_g_xz*costh*cosph*sinth + dt_g_zz*sinth2 + two*dt_g_xy*costh2*cosph*sinph - two*dt_g_yz*costh*sinth*sinph + dt_g_yy*costh2*sinph2);
          
          dt_gthphi[l] = r2*(-(dt_g_yz*cosph*sinth2) - dt_g_xx*costh*cosph*sinth*sinph + dt_g_yy*costh*cosph*sinth*sinph + dt_g_xz*sinth2*sinph + dt_g_xy*costh*sinth*(cosph2 - sinph2));
          
          dt_gphiphi[l] = r2*(dt_g_yy*cosph2*sinth2 - two*dt_g_xy*cosph*sinth2*sinph + dt_g_xx*sinth2*sinph2);
          

         if (it+lbnd[0]>=int_ntheta[*current_detector-1] ||
            ip+lbnd[1]>=int_nphi[*current_detector-1] )
         {
             
             gtr[l] = 0;
             gtth[l] = 0;
             gtphi[l] = 0;

            grr[l]=0;
            grth[l]=0;
            grphi[l]=0;
            gthth[l]=0;
            gthphi[l]=0;
            gphiphi[l]=0;
             
            dr_gthth[l]=0;
            dr_gthphi[l]=0;
            dr_gphiphi[l]=0;
             
             dr_gtth[l] = 0;
             dr_gtphi[l] = 0;
             
             dt_grth[l] = 0;
             dt_grphi[l] = 0;
             dt_gthth[l] = 0;
             dt_gthphi[l] = 0;
             dt_gphiphi[l] = 0;
         }
      }
   }
}
