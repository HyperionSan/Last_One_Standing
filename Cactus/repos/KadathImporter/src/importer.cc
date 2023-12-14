#include <cmath>
#include <functional>
#include <string>
#include <chrono>
#include <vector>
#include <array>
#include <type_traits>

#ifdef _STANDALONE_
#include "importer.h"
#include <iostream>
#else
#include <importer.h>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

//template function that handles moving the space-time values to their cctk couterparts
template<typename T>
void import_data(T& exported_vals, CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INFO("Moving Vacuum Data to output vars");
  int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  #pragma omp parallel for
  for (int i = 0; i < npoints; ++i) {
    // FIXME this expects GXX to be conformally flat
    alp[i] = (puncture_lapse) ? 1. / std::sqrt(exported_vals[GXX][i]) : exported_vals[ALPHA][i];

    // we don't set beta
    betax[i] = 0;
    betay[i] = 0;
    betaz[i] = 0;

    gxx[i] = exported_vals[GXX][i];
    gxy[i] = exported_vals[GXY][i];
    gxz[i] = exported_vals[GXZ][i];
    gyy[i] = exported_vals[GYY][i];
    gyz[i] = exported_vals[GYZ][i];
    gzz[i] = exported_vals[GZZ][i];

    kxx[i] = exported_vals[KXX][i];
    kxy[i] = exported_vals[KXY][i];
    kxz[i] = exported_vals[KXZ][i];
    kyy[i] = exported_vals[KYY][i];
    kyz[i] = exported_vals[KYZ][i];
    kzz[i] = exported_vals[KZZ][i];
  } // for i
}

//template function that handles moving the space-time values to their cctk couterparts
template<typename T>
void import_data_wmatter(T& exported_vals, CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INFO("Moving Matter Data to output vars");
  int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  #pragma omp parallel for
  for (int i = 0; i < npoints; ++i) {
    // FIXME this expects GXX to be conformally flat
    alp[i] = (puncture_lapse) ? 1. / std::sqrt(exported_vals[GXX][i]) : exported_vals[ALPHA][i];

    // we don't set beta
    betax[i] = 0;
    betay[i] = 0;
    betaz[i] = 0;

    gxx[i] = exported_vals[GXX][i];
    gxy[i] = exported_vals[GXY][i];
    gxz[i] = exported_vals[GXZ][i];
    gyy[i] = exported_vals[GYY][i];
    gyz[i] = exported_vals[GYZ][i];
    gzz[i] = exported_vals[GZZ][i];

    kxx[i] = exported_vals[KXX][i];
    kxy[i] = exported_vals[KXY][i];
    kxz[i] = exported_vals[KXZ][i];
    kyy[i] = exported_vals[KYY][i];
    kyz[i] = exported_vals[KYZ][i];
    kzz[i] = exported_vals[KZZ][i];
    
    rho[i] = exported_vals[RHO][i];
    eps[i] = exported_vals[EPS][i];
    press[i] = exported_vals[PRESS][i];

    vel[i] = exported_vals[VELX][i];
    vel[i + npoints] = exported_vals[VELY][i];
    vel[i + 2 * npoints] = exported_vals[VELZ][i];
  } // for i
}
#endif

//Generalized importer to reduce duplicate code
#ifdef _STANDALONE_
template<typename T>
void print_var(T ary, const size_t var) {
  for(const auto& e :  ary[var])
    std::cout << e << std::endl;
}

extern "C" void KadathImporter(const int npoints, std::string id_type, std::vector<double> xx, 
  std::vector<double> yy, std::vector<double> zz, const int interp_order, const double interpolation_offset, const double delta_r_rel) {
    const char* filename{"id.info"};
#else
extern "C" void KadathImporter(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO("Setting up KADATH initial data");

  CCTK_INFO("Setting up coordinates");

  int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  std::vector<double> xx(npoints), yy(npoints), zz(npoints);

  #pragma omp parallel for
  for (int i = 0; i < npoints; ++i) {
    xx[i] = x[i];
    yy[i] = y[i];
    zz[i] = z[i];
  }
  std::string id_type{type};
  CCTK_INFO("Starting data import");
#endif
  auto start = std::chrono::system_clock::now();
  if(id_type == "BH") {
    auto exported_vals = std::move(KadathExportBH(npoints, xx.data(), 
        yy.data(), zz.data(), filename, interpolation_offset, interp_order, delta_r_rel));
    
    #ifndef _STANDALONE_
    import_data(exported_vals, CCTK_PASS_CTOC);
    #endif
  } else if(id_type == "BBH") {
    auto exported_vals = std::move(KadathExportBBH(npoints, xx.data(), 
        yy.data(), zz.data(), filename, interpolation_offset, interp_order, delta_r_rel));
    
    #ifndef _STANDALONE_
    import_data(exported_vals, CCTK_PASS_CTOC);
    #endif
  } else if(id_type == "NS") {
    auto exported_vals = std::move(KadathExportNS(npoints, xx.data(), yy.data(), zz.data(), filename));
    
    #ifndef _STANDALONE_
    import_data_wmatter(exported_vals, CCTK_PASS_CTOC);
    #endif
  } else if(id_type == "BNS") {
    auto exported_vals = std::move(KadathExportBNS(npoints, xx.data(), yy.data(), zz.data(), filename));
    
    #ifndef _STANDALONE_
    import_data_wmatter(exported_vals, CCTK_PASS_CTOC);
    #else
    print_var(exported_vals, RHO);
    #endif
  } else if(id_type == "BHNS") {
    auto exported_vals = std::move(KadathExportBHNS(npoints, xx.data(), yy.data(), zz.data(), filename, 
          interpolation_offset, interp_order, delta_r_rel));
    
    #ifndef _STANDALONE_
    import_data_wmatter(exported_vals, CCTK_PASS_CTOC);
    #endif
  }
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  
  #ifndef _STANDALONE_
  CCTK_VInfo(CCTK_THORNSTRING, "Filling took %g sec", elapsed_seconds.count());
  #else
  printf("Filling took %g sec\n", elapsed_seconds.count());
  #endif
}
