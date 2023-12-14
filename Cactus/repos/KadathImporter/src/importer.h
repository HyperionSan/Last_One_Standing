/**
 * file: importer.hh
 * author(s)/maintainer(s): 
 *  Samuel David Tootle <tootle@itp.uni-frankfurt.de>
 *  Ludwig Jens Papenfort <papenfort@th.physik.uni-frankfurt.de>
 * License: GPL
 *
 * The following is a minimal header file to include the necessary enumerators and function defintions
 * that reflect those used by the exporters included in the FUKA branch of Kadath.  This makes it
 * possible to only need to compile the Einstein Toolkit with the static libary stored in ${HOME_KADATH}/lib
 * without needing to include or compile with Kadath header files.  This is especially important
 * since features in FUKA require C++17, whereas, the toolkit is primarily based on C++11/14.
 * Therefore, this circumvents the need to compile with C++17 when compiling the Einstein Toolkit.
 */

#include <array>
#include <vector>

// enumeration for quantities to export to evolution codes
enum sim_vac_quants {
  ALPHA,
  BETAX,
  BETAY,
  BETAZ,
  GXX,
  GXY,
  GXZ,
  GYY,
  GYZ,
  GZZ,
  KXX,
  KXY,
  KXZ,
  KYY,
  KYZ,
  KZZ,
  NUM_VOUT
};

// chained enums together to include matter quantities
enum sim_matter_quants {
  RHO=NUM_VOUT,
  EPS,
  PRESS,
  VELX,
  VELY,
  VELZ,
  NUM_OUT
};


constexpr size_t VAC = 1;
constexpr size_t MATTER = 2;
std::array<std::vector<double>,NUM_VOUT> KadathExportBH(int const npoints,
                                                        double const * xx, double const * yy, double const * zz,
                                                        char const * fn,
                                                        double const interpolation_offset, int const interp_order,
                                                        double const delta_r_rel);

std::array<std::vector<double>,NUM_OUT> KadathExportNS(int const npoints,
                                                        double const * xx, double const * yy, double const * zz,
                                                        char const * fn);

std::array<std::vector<double>,NUM_OUT> KadathExportBNS(int const npoints,
                                                        double const * xx, double const * yy, double const * zz,
                                                        char const * fn);

std::array<std::vector<double>,NUM_VOUT> KadathExportBBH(int const npoints,
                                                        double const * xx, double const * yy, double const * zz,
                                                        char const * fn,
                                                        double const interpolation_offset, int const interp_order,
                                                        double const delta_r_rel);

std::array<std::vector<double>,NUM_OUT> KadathExportBHNS(int const npoints,
                                                        double const * xx, double const * yy, double const * zz,
                                                        char const * fn,
                                                        double const interpolation_offset, int const interp_order,
                                                        double const delta_r_rel);

#ifdef _STANDALONE_
extern "C" void KadathImporter(const int npoints, std::string id_type, std::vector<double> xx, 
  std::vector<double> yy, std::vector<double> zz, const int interp_order = 8, const double interpolation_offset=0., const double delta_r_rel = 0.3); 
#endif
