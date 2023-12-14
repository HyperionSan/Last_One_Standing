#include <cstddef>
#include "gsl/gsl_sf_legendre.h"

namespace gsl_interface {

  extern "C" {
    void gsl_sf_legendre_sphPlm_deriv ( const int lmax, const double x, double res[], double res_d[] ) {

      size_t mylmax{lmax};
      gsl_sf_legendre_deriv_alt_array_e ( GSL_SF_LEGENDRE_SPHARM, mylmax, x, -1.0, res, res_d );
    }

    void gsl_sf_legendre_sphPlm_deriv2 ( const int lmax, const double x, double res[], double res_d[], double res_d2[] ) {

      size_t mylmax{lmax};
      gsl_sf_legendre_deriv2_alt_array_e ( GSL_SF_LEGENDRE_SPHARM, mylmax, x, -1.0, res, res_d, res_d2 );
    }

    int gsl_sf_legendre_index ( const int l, const int m ) {

      size_t myl{l};
      size_t mym{m};
      size_t myindex{gsl_sf_legendre_array_index(myl,mym)};
      int index{myindex+1};
      return index;
    }

    int gsl_sf_legendre_n ( const int lmax ) {
   
      size_t mylmax{lmax};
      size_t myn{gsl_sf_legendre_array_n(mylmax)};
      int n{myn};
      return n;
    }
  }

}
