#include <complex>
#include <cassert>
#include <iostream>
#include <vector>
#include "EffectiveSource.h"

namespace rwz_source_interface {
  static std::vector<RWZ_EffectiveSource*> effsource;

  extern "C" {
    void init_rwz_source ( const int* nmodes, const int l[],
		       const int m[], const double* M , const int mparity[]) {
      assert(effsource.empty());
      for (int i=0; i<*nmodes; i++) {
        RWZ_EffectiveSource::parity_t parity{(mparity[i]%2==0)?RWZ_EffectiveSource::EvenParity : RWZ_EffectiveSource::OddParity};
        effsource.push_back(new RWZ_EffectiveSource(l[i], m[i], *M, parity));
      }
    }

/*    void init_source_q ( const int* nmodes, const int l[],
		       const int m[], const long double* M ) {
      assert(effsource.empty());
      double Md = *M;
      for (int i=0; i<*nmodes; i++) {
        effsource.push_back(new EffectiveSource(l[i], m[i], Md));
      }
    } */

/*    void set_window_q ( const long double* r1, const long double* w1,
		        const long double* q1, const long double* s1,
			const long double* r2, const long double* w2,
		        const long double* q2, const long double* s2,
			const int* nmodes ) {
      assert(effsource.size()==*nmodes);
      double r1d = *r1;
      double w1d = *w1;
      double q1d = *q1;
      double s1d = *s1;
      double r2d = *r2;
      double w2d = *w2;
      double q2d = *q2;
      double s2d = *s2;
      for (auto& x : effsource) 
	(*x).set_window ( r1d, w1d, q1d, s1d, r2d, w2d, q2d, s2d );
    } */

/*    void calc_window_q ( const int* n, const long double r[],
		       double long Win[], double long dWin[],
		       double long d2Win[] ) {
      assert(!effsource.empty());
      double rd[*n];
      double Wind[*n];
      double dWind[*n];
      double d2Wind[*n];
      for (int i=0; i++; i<*n) {
	rd[i] = r[i];
      }
      (*effsource.at(0)).calc_window ( *n, rd, Wind, dWind, d2Wind );
      for (int i=0; i++; i<*n) {
        Win[i] = Wind[i];
        dWin[i] = dWind[i];
        d2Win[i] = d2Wind[i];
      }
    } */
     
    void set_rwz_time_window ( const double* T, const double* dT_dt,
		           const double* d2T_dt2, 
			   const int* do_smooth_after_lmax, 
			   const int* nmodes ) {
      assert(static_cast<int>(effsource.size())==*nmodes);
//std::cout << "nmodes= " << *nmodes << ", do smooth= " << *do_smooth_after_lmax << std::endl;
      for (auto& x : effsource)
        if ( (*x).get_l() <= *do_smooth_after_lmax ) {
//std::cout << "we're in" << std::endl;
          (*x).set_rwz_time_window ( 1.0, 0.0, 0.0 );
	} else {
//std::cout << "we're out" << std::endl;
          (*x).set_rwz_time_window ( *T, *dT_dt, *d2T_dt2 );
        }
    }

/*    void set_time_window_q ( const long double* T, const long double* dT_dt,
		           const long double* d2T_dt2, const int* nmodes ) {
      assert(effsource.size()==*nmodes);
      double Td = *T;
      double dT_dtd = *dT_dt;
      double d2T_dt2d = *d2T_dt2;
      for (auto& x : effsource)
        (*x).set_time_window ( Td, dT_dtd, d2T_dt2d );
    } */

    void set_rwz_particle ( const double* r, const double* psi,
		        const double* ur, const double* En,
			const double* Lz, const double* ar,
			const double* apsi, const double* dardt,
			const double* dapsidt, const double* d2ardt2,
			const double* d2apsidt2, const int* nmodes ) {
      assert(static_cast<int>(effsource.size())==*nmodes);
#pragma omp parallel for shared(r,psi,ur,En,Lz,ar,apsi,dardt,dapsidt,d2ardt2,d2apsidt2,effsource)
/*      for (auto& x : effsource) 
	(*x).set_particle ( *r, *psi, *ur, *En, *Lz, *ar, *apsi, *dardt,
                            *dapsidt, *d2ardt2, *d2apsidt2); */

        for (int i=0; i<*nmodes;i++) {
          (*effsource.at(i)).set_rwz_particle ( *r, *psi, *ur, *En, *Lz, *ar, *apsi,
			                 *dardt, *dapsidt, *d2ardt2, *d2apsidt2);
        }
    }

/*    void set_particle_q ( const long double* r, const long double* psi,
		        const long double* ur, const long double* En,
		        const long double* Lz, const long double* ar,
			const long double* apsi, const long double* dardt,
                        const long double* dapsidt, const long double* d2ardt2,
                        const long double* d2apsidt2, const int* nmodes ) {
      assert(effsource.size()==*nmodes);
      double rd = *r;
      double psid = *psi;
      double urd = *ur;
      double End = *En;
      double Lzd = *Lz;
      double ard = *ar;
      double apsid = *apsi;
      double dardtd = *dardt;
      double dapsidtd = *dapsidt;
      double d2ardt2d = *d2ardt2;
      double d2apsidt2d = *d2apsidt2;
      for (auto& x : effsource) 
	(*x).set_particle ( rd, psid, urd, End, Lzd, ard, apsid, dardtd,
                            dapsidtd, d2ardt2d, d2apsidt2d );
    } */

/*    void set_particle ( const double* p, const double* e,
		        const double* chi, const double* psi,
			const double* r, const double* ur,
			const int* use_osc, const double* ar,
			const double* apsi, const int* nmodes ) {
      assert(effsource.size()==*nmodes);
      for (auto& x : effsource) 
	(*x).set_particle ( *p, *e, *chi, *psi, *r, *ur, *use_osc, *ar, *apsi);
    }

    void set_particle_q ( const long double* p, const long double* e,
		        const long double* chi, const long double* psi,
		        const long double* r, const long double* ur,
                        const int* use_osc, const long double* ar,
			const long double* apsi, const int* nmodes ) {
      assert(effsource.size()==*nmodes);
      double pd = *p;
      double ed = *e;
      double chid = *chi;
      double psid = *psi;
      double rd = *r;
      double urd = *ur;
      double ard = *ar;
      double apsid = *apsi;
      for (auto& x : effsource) 
	(*x).set_particle ( pd, ed, chid, psid, rd, urd, *use_osc, ard, apsid );
    } */

/*    void set_particle ( const double* p, const double* e,
		        const double* chi, const double* psi,
			const double* r, const double* ur,
			const int* use_osc, const int* nmodes ) {
      assert(effsource.size()==*nmodes);
      for (auto& x : effsource) 
	(*x).set_particle ( *p, *e, *chi, *psi, *r, *ur, *use_osc );
    }

    void set_particle_q ( const long double* p, const long double* e,
		        const long double* chi, const long double* psi,
		        const long double* r, const long double* ur,
                        const int* use_osc, const int* nmodes ) {
      assert(effsource.size()==*nmodes);
      double pd = *p;
      double ed = *e;
      double chid = *chi;
      double psid = *psi;
      double rd = *r;
      double urd = *ur;
      for (auto& x : effsource) 
	(*x).set_particle ( pd, ed, chid, psid, rd, urd, *use_osc );
    } */

    void eval_rwz_source (const int* mode,  const double* r,
		      double* sre, double* sim ) {
      /* I first check that the mode number is zero or positive */
      assert((*mode)>=0);
      /* I then convert it to size_t and compare it with the size of the
         vector of effective sources. Hopefully, in the event of a compiler
         where the size of a vector is not of type size_t, the compiler will
         complain. */
      assert(effsource.size()>=static_cast<size_t>(*mode));
      std::complex<double> src = (*(effsource.at(*mode)))(*r);
      *sre = std::real(src);
      *sim = std::imag(src);
    }

/*    void eval_source_q (const int* mode,  const long double* r,
		      long double* sre, long double* sim ) {
      assert(effsource.size()>=*mode);
      double rd = *r;
      std::complex<double> src = (*(effsource.at(*mode)))(rd);
      *sre = std::real(src);
      *sim = std::imag(src);
    } */

    void eval_rwz_source_all (const int* mode,  const int* n, const double r[],
                          const double Win[], const double dWin[],
                          const double d2Win[], double sre[], double sim[] ) {
/*      printf("mode = %d\n", *mode);
      printf("n = %d\n", *n); */
      assert((*mode)>=0);
      assert(effsource.size()>=static_cast<size_t>(*mode));

      (*(effsource.at(*mode)))(*n, r, Win, dWin, d2Win, sre, sim);
    }

/*    void eval_source_all_q (const int* mode,  const int* n, 
		            const long double r[], const long double Win[],
			    const long double dWin[], const long double d2Win[],
			    long double sre[], long double sim[] ) {
      assert(effsource.size()>=*mode);
      double rd[*n];
      double Wind[*n];
      double dWind[*n];
      double d2Wind[*n];
      double sred[*n];
      double simd[*n];
      for (int i=0; i++; i<*n) {
	rd[i] = r[i];
	Wind[i] = Win[i];
	dWind[i] = dWin[i];
	d2Wind[i] = d2Win[i];
      }
      (*(effsource.at(*mode)))(*n, rd, Wind, dWind, d2Wind, sred, simd);
      for (int i=0; i++; i<*n) {
	sre[i] = sred[i];
	sim[i] = simd[i];
      }
    } */

    void clean_rwz_source () {
      /* Use integers of type size_t in order to avoid type conversions */
      size_t s{effsource.size()};
      for (size_t i=0; i<s; i++) {
        assert(effsource.at(s-1-i)!=nullptr);
        delete effsource.at(s-1-i);
        effsource.pop_back(); 
      }
      assert(effsource.size()==0);
    }

    void Psi ( const int* mode, const double* r, const int* sign,
               double* psire, double* psiim ) {
      assert((*mode)>=0);
      assert(effsource.size()>=static_cast<size_t>(*mode));
      std::complex<double> psi{(*(effsource.at(*mode))).Psi(*r, *sign)};
      *psire = std::real(psi);
      *psiim = std::imag(psi);
    }

/*    void Psi_q ( const int* mode, const long double* r,
               long double* psire, long double* psiim ) {
      assert(effsource.size()>=*mode);
      double rd = *r;
      std::complex<double> psi{(*(effsource.at(*mode))).Psi(rd)};
      *psire = std::real(psi);
      *psiim = std::imag(psi);
    } */

    void dPsi_dr ( const int* mode, const double* r, const int* sign,
                   double* dpsidrre, double* dpsidrim ) {
      assert((*mode)>=0);
      assert(effsource.size()>=static_cast<size_t>(*mode));
      std::complex<double> dpsidr{(*(effsource.at(*mode))).dPsi_dr(*r, *sign)};
      *dpsidrre = std::real(dpsidr);
      *dpsidrim = std::imag(dpsidr);
    }

/*    void dPsi_dr_q ( const int* mode, const long double* r,
                   long double* dpsidrre, long double* dpsidrim ) {
      assert(effsource.size()>=*mode);
      double rd = *r;
      std::complex<double> dpsidr{(*(effsource.at(*mode))).dPsi_dr(rd)};
      *dpsidrre = std::real(dpsidr);
      *dpsidrim = std::imag(dpsidr);
    } */

    void dPsi_dt ( const int* mode, const double* r, const int* sign,
                   double* dpsidtre, double* dpsidtim ) {
      assert((*mode)>=0);
      assert(effsource.size()>=static_cast<size_t>(*mode));
      std::complex<double> dpsidt{(*(effsource.at(*mode))).dPsi_dt(*r, *sign)};
      *dpsidtre = std::real(dpsidt);
      *dpsidtim = std::imag(dpsidt);
    }

/*    void dPsi_dt_q ( const int* mode, const long double* r,
                   long double* dpsidtre, long double* dpsidtim ) {
      assert(effsource.size()>=*mode);
      double rd = *r;
      std::complex<double> dpsidt{(*(effsource.at(*mode))).dPsi_dt(rd)};
      *dpsidtre = std::real(dpsidt);
      *dpsidtim = std::imag(dpsidt);
    } */
  }
}
