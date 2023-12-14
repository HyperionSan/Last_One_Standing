/*
 * Copyright 2016 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#ifndef EFFECTIVESOURCE_H
#define EFFECTIVESOURCE_H

#include <complex>
#include <array>

class RWZ_EffectiveSource
{
public:
  enum parity_t { EvenParity, OddParity };
  RWZ_EffectiveSource(const int l, const int m, const double M, const parity_t parity);

  void operator()(const int n, const double r[],
    const double W[], const double dW_dr[], const double d2W_dr2[],
    double sre[], double sim[]);
  std::complex<double> operator()(const double r);
  std::complex<double> Psi(const double r, const int sign);
  std::complex<double> dPsi_dr(const double r, const int sign);
  std::complex<double> dPsi_dt(const double r, const int sign);

  void set_rwz_time_window(const double T, const double dT_dt, const double d2T_dt2);
  void set_rwz_particle(const double r, const double phi,
                    const double ur, const double E, const double L,
                    const double ar, const double aphi,
                    const double dardt, const double daphidt,
                    const double d2ardt2, const double d2aphidt2);
  int get_l();

private:
  const int l, m;
  const double M;
  const parity_t parity;

  double r0, phi0;
  double T, dT_dt, d2T_dt2;
  double rt, rtt, phit, phitt, c, ct, ctt;
  std::complex<double> w0, w1p, w1m, w2p, w2m, w3p, w3m;
  std::array<std::array<double,4>,4> coeffs, dcoeffs_dt, d2coeffs_dt2;
  std::array<std::array<double,4>,4> abscoeffs, dabscoeffs_dt, d2abscoeffs_dt2;
};

#endif /* EFFECTIVESOURCE_H */
