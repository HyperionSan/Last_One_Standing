#include "math.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
/*
 * Converting from ADM to BSSN quantities is required in the Einstein Toolkit,
 * as initial data are given in terms of ADM quantities, and Baikal evolves the BSSN quantities.
 */
void ADM_to_BSSN_Baikal(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS_ADM_to_BSSN_Baikal;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL *restrict alphaSphorCartGF = alp;
  const CCTK_REAL *restrict BSphorCartU0GF = dtbetax;
  const CCTK_REAL *restrict BSphorCartU1GF = dtbetay;
  const CCTK_REAL *restrict BSphorCartU2GF = dtbetaz;
  const CCTK_REAL *restrict KSphorCartDD00GF = kxx;
  const CCTK_REAL *restrict KSphorCartDD01GF = kxy;
  const CCTK_REAL *restrict KSphorCartDD02GF = kxz;
  const CCTK_REAL *restrict KSphorCartDD11GF = kyy;
  const CCTK_REAL *restrict KSphorCartDD12GF = kyz;
  const CCTK_REAL *restrict KSphorCartDD22GF = kzz;
  const CCTK_REAL *restrict betaSphorCartU0GF = betax;
  const CCTK_REAL *restrict betaSphorCartU1GF = betay;
  const CCTK_REAL *restrict betaSphorCartU2GF = betaz;
  const CCTK_REAL *restrict gammaSphorCartDD00GF = gxx;
  const CCTK_REAL *restrict gammaSphorCartDD01GF = gxy;
  const CCTK_REAL *restrict gammaSphorCartDD02GF = gxz;
  const CCTK_REAL *restrict gammaSphorCartDD11GF = gyy;
  const CCTK_REAL *restrict gammaSphorCartDD12GF = gyz;
  const CCTK_REAL *restrict gammaSphorCartDD22GF = gzz;
  #pragma omp parallel for
  for (int i2 = 0; i2 < cctk_lsh[2]; i2++) {
    for (int i1 = 0; i1 < cctk_lsh[1]; i1++) {
      for (int i0 = 0; i0 < cctk_lsh[0]; i0++) {
        /*
         * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
         */
        const double alphaSphorCart = alphaSphorCartGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double betaSphorCartU0 = betaSphorCartU0GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double betaSphorCartU1 = betaSphorCartU1GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double betaSphorCartU2 = betaSphorCartU2GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double BSphorCartU0 = BSphorCartU0GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double BSphorCartU1 = BSphorCartU1GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double BSphorCartU2 = BSphorCartU2GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double gammaSphorCartDD00 = gammaSphorCartDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double gammaSphorCartDD01 = gammaSphorCartDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double gammaSphorCartDD02 = gammaSphorCartDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double gammaSphorCartDD11 = gammaSphorCartDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double gammaSphorCartDD12 = gammaSphorCartDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double gammaSphorCartDD22 = gammaSphorCartDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double KSphorCartDD00 = KSphorCartDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double KSphorCartDD01 = KSphorCartDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double KSphorCartDD02 = KSphorCartDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double KSphorCartDD11 = KSphorCartDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double KSphorCartDD12 = KSphorCartDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double KSphorCartDD22 = KSphorCartDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        /*
         * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
         */
        const double FDPart3_1 = gammaSphorCartDD00*((gammaSphorCartDD12)*(gammaSphorCartDD12));
        const double FDPart3_3 = ((gammaSphorCartDD01)*(gammaSphorCartDD01))*gammaSphorCartDD22;
        const double FDPart3_5 = ((gammaSphorCartDD02)*(gammaSphorCartDD02))*gammaSphorCartDD11;
        const double FDPart3_6 = -FDPart3_1 - FDPart3_3 - FDPart3_5 + gammaSphorCartDD00*gammaSphorCartDD11*gammaSphorCartDD22 + 2*gammaSphorCartDD01*gammaSphorCartDD02*gammaSphorCartDD12;
        const double FDPart3_7 = (1.0/(FDPart3_6));
        const double FDPart3_8 = cbrt(FDPart3_7);
        const double FDPart3_9 = 2*FDPart3_7;
        const double FDPart3_10 = FDPart3_7*KSphorCartDD00*(gammaSphorCartDD11*gammaSphorCartDD22 - ((gammaSphorCartDD12)*(gammaSphorCartDD12))) + FDPart3_7*KSphorCartDD11*(gammaSphorCartDD00*gammaSphorCartDD22 - ((gammaSphorCartDD02)*(gammaSphorCartDD02))) + FDPart3_7*KSphorCartDD22*(gammaSphorCartDD00*gammaSphorCartDD11 - ((gammaSphorCartDD01)*(gammaSphorCartDD01))) + FDPart3_9*KSphorCartDD01*(-gammaSphorCartDD01*gammaSphorCartDD22 + gammaSphorCartDD02*gammaSphorCartDD12) + FDPart3_9*KSphorCartDD02*(gammaSphorCartDD01*gammaSphorCartDD12 - gammaSphorCartDD02*gammaSphorCartDD11) + FDPart3_9*KSphorCartDD12*(-gammaSphorCartDD00*gammaSphorCartDD12 + gammaSphorCartDD01*gammaSphorCartDD02);
        const double FDPart3_11 = (1.0/3.0)*FDPart3_10;
        hDD00GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_8*gammaSphorCartDD00 - 1;
        hDD01GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_8*gammaSphorCartDD01;
        hDD02GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_8*gammaSphorCartDD02;
        hDD11GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_8*gammaSphorCartDD11 - 1;
        hDD12GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_8*gammaSphorCartDD12;
        hDD22GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_8*gammaSphorCartDD22 - 1;
        aDD00GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_8*(-FDPart3_11*gammaSphorCartDD00 + KSphorCartDD00);
        aDD01GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_8*(-FDPart3_11*gammaSphorCartDD01 + KSphorCartDD01);
        aDD02GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_8*(-FDPart3_11*gammaSphorCartDD02 + KSphorCartDD02);
        aDD11GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_8*(-FDPart3_11*gammaSphorCartDD11 + KSphorCartDD11);
        aDD12GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_8*(-FDPart3_11*gammaSphorCartDD12 + KSphorCartDD12);
        aDD22GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_8*(-FDPart3_11*gammaSphorCartDD22 + KSphorCartDD22);
        trKGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_10;
        vetU0GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = betaSphorCartU0;
        vetU1GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = betaSphorCartU1;
        vetU2GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = betaSphorCartU2;
        betU0GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = BSphorCartU0;
        betU1GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = BSphorCartU1;
        betU2GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = BSphorCartU2;
        alphaGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = alphaSphorCart;
        cfGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = pow(FDPart3_6/(-FDPart3_1*FDPart3_7 - FDPart3_3*FDPart3_7 - FDPart3_5*FDPart3_7 + FDPart3_7*gammaSphorCartDD00*gammaSphorCartDD11*gammaSphorCartDD22 + 2*FDPart3_7*gammaSphorCartDD01*gammaSphorCartDD02*gammaSphorCartDD12), -1.0/6.0);
        
      } // END LOOP: for (int i0 = 0; i0 < cctk_lsh[0]; i0++)
    } // END LOOP: for (int i1 = 0; i1 < cctk_lsh[1]; i1++)
  } // END LOOP: for (int i2 = 0; i2 < cctk_lsh[2]; i2++)

  const CCTK_REAL invdx0 = 1.0/CCTK_DELTA_SPACE(0);
  const CCTK_REAL invdx1 = 1.0/CCTK_DELTA_SPACE(1);
  const CCTK_REAL invdx2 = 1.0/CCTK_DELTA_SPACE(2);
  if(FD_order == 2) {
#pragma omp parallel for
    for (int i2 = cctk_nghostzones[2]; i2 < cctk_lsh[2]-cctk_nghostzones[2]; i2++) {
      for (int i1 = cctk_nghostzones[1]; i1 < cctk_lsh[1]-cctk_nghostzones[1]; i1++) {
        for (int i0 = cctk_nghostzones[0]; i0 < cctk_lsh[0]-cctk_nghostzones[0]; i0++) {
            /*
             * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
             */
            const double hDD00_i0_i1_i2m1 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
            const double hDD00_i0_i1m1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
            const double hDD00_i0m1_i1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
            const double hDD00 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD00_i0p1_i1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
            const double hDD00_i0_i1p1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
            const double hDD00_i0_i1_i2p1 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
            const double hDD01_i0_i1_i2m1 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
            const double hDD01_i0_i1m1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
            const double hDD01_i0m1_i1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
            const double hDD01 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD01_i0p1_i1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
            const double hDD01_i0_i1p1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
            const double hDD01_i0_i1_i2p1 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
            const double hDD02_i0_i1_i2m1 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
            const double hDD02_i0_i1m1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
            const double hDD02_i0m1_i1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
            const double hDD02 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD02_i0p1_i1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
            const double hDD02_i0_i1p1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
            const double hDD02_i0_i1_i2p1 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
            const double hDD11_i0_i1_i2m1 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
            const double hDD11_i0_i1m1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
            const double hDD11_i0m1_i1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
            const double hDD11 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD11_i0p1_i1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
            const double hDD11_i0_i1p1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
            const double hDD11_i0_i1_i2p1 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
            const double hDD12_i0_i1_i2m1 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
            const double hDD12_i0_i1m1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
            const double hDD12_i0m1_i1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
            const double hDD12 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD12_i0p1_i1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
            const double hDD12_i0_i1p1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
            const double hDD12_i0_i1_i2p1 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
            const double hDD22_i0_i1_i2m1 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
            const double hDD22_i0_i1m1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
            const double hDD22_i0m1_i1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
            const double hDD22 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD22_i0p1_i1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
            const double hDD22_i0_i1p1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
            const double hDD22_i0_i1_i2p1 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
            const double FDPart1_Rational_1_2 = 1.0/2.0;
            const double FDPart1_0 = FDPart1_Rational_1_2*invdx0;
            const double FDPart1_1 = FDPart1_Rational_1_2*invdx1;
            const double FDPart1_2 = FDPart1_Rational_1_2*invdx2;
            const double hDD_dD000 = FDPart1_0*(-hDD00_i0m1_i1_i2 + hDD00_i0p1_i1_i2);
            const double hDD_dD001 = FDPart1_1*(-hDD00_i0_i1m1_i2 + hDD00_i0_i1p1_i2);
            const double hDD_dD002 = FDPart1_2*(-hDD00_i0_i1_i2m1 + hDD00_i0_i1_i2p1);
            const double hDD_dD010 = FDPart1_0*(-hDD01_i0m1_i1_i2 + hDD01_i0p1_i1_i2);
            const double hDD_dD011 = FDPart1_1*(-hDD01_i0_i1m1_i2 + hDD01_i0_i1p1_i2);
            const double hDD_dD012 = FDPart1_2*(-hDD01_i0_i1_i2m1 + hDD01_i0_i1_i2p1);
            const double hDD_dD020 = FDPart1_0*(-hDD02_i0m1_i1_i2 + hDD02_i0p1_i1_i2);
            const double hDD_dD021 = FDPart1_1*(-hDD02_i0_i1m1_i2 + hDD02_i0_i1p1_i2);
            const double hDD_dD022 = FDPart1_2*(-hDD02_i0_i1_i2m1 + hDD02_i0_i1_i2p1);
            const double hDD_dD110 = FDPart1_0*(-hDD11_i0m1_i1_i2 + hDD11_i0p1_i1_i2);
            const double hDD_dD111 = FDPart1_1*(-hDD11_i0_i1m1_i2 + hDD11_i0_i1p1_i2);
            const double hDD_dD112 = FDPart1_2*(-hDD11_i0_i1_i2m1 + hDD11_i0_i1_i2p1);
            const double hDD_dD120 = FDPart1_0*(-hDD12_i0m1_i1_i2 + hDD12_i0p1_i1_i2);
            const double hDD_dD121 = FDPart1_1*(-hDD12_i0_i1m1_i2 + hDD12_i0_i1p1_i2);
            const double hDD_dD122 = FDPart1_2*(-hDD12_i0_i1_i2m1 + hDD12_i0_i1_i2p1);
            const double hDD_dD220 = FDPart1_0*(-hDD22_i0m1_i1_i2 + hDD22_i0p1_i1_i2);
            const double hDD_dD221 = FDPart1_1*(-hDD22_i0_i1m1_i2 + hDD22_i0_i1p1_i2);
            const double hDD_dD222 = FDPart1_2*(-hDD22_i0_i1_i2m1 + hDD22_i0_i1_i2p1);
            /*
             * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
             */
            const double FDPart3_0 = hDD22 + 1;
            const double FDPart3_1 = -FDPart3_0*hDD01 + hDD02*hDD12;
            const double FDPart3_4 = hDD11 + 1;
            const double FDPart3_6 = hDD00 + 1;
            const double FDPart3_7 = (1.0/(FDPart3_0*FDPart3_4*FDPart3_6 - FDPart3_0*((hDD01)*(hDD01)) - FDPart3_4*((hDD02)*(hDD02)) - FDPart3_6*((hDD12)*(hDD12)) + 2*hDD01*hDD02*hDD12));
            const double FDPart3_8 = (1.0/2.0)*FDPart3_7;
            const double FDPart3_9 = FDPart3_1*FDPart3_8;
            const double FDPart3_10 = -FDPart3_4*hDD02 + hDD01*hDD12;
            const double FDPart3_11 = FDPart3_10*FDPart3_8;
            const double FDPart3_12 = hDD_dD012 + hDD_dD021 - hDD_dD120;
            const double FDPart3_13 = FDPart3_7*(FDPart3_0*FDPart3_4 - ((hDD12)*(hDD12)));
            const double FDPart3_14 = (1.0/2.0)*FDPart3_13;
            const double FDPart3_15 = -FDPart3_6*hDD12 + hDD01*hDD02;
            const double FDPart3_16 = 2*FDPart3_7;
            const double FDPart3_17 = FDPart3_15*FDPart3_16;
            const double FDPart3_18 = hDD_dD012 - hDD_dD021 + hDD_dD120;
            const double FDPart3_19 = FDPart3_10*FDPart3_16;
            const double FDPart3_20 = -hDD_dD012 + hDD_dD021 + hDD_dD120;
            const double FDPart3_21 = FDPart3_1*FDPart3_16;
            const double FDPart3_22 = 2*hDD_dD122 - hDD_dD221;
            const double FDPart3_23 = 2*hDD_dD022 - hDD_dD220;
            const double FDPart3_24 = FDPart3_7*(FDPart3_4*FDPart3_6 - ((hDD01)*(hDD01)));
            const double FDPart3_25 = -hDD_dD112 + 2*hDD_dD121;
            const double FDPart3_26 = 2*hDD_dD011 - hDD_dD110;
            const double FDPart3_27 = FDPart3_7*(FDPart3_0*FDPart3_6 - ((hDD02)*(hDD02)));
            const double FDPart3_28 = -hDD_dD001 + 2*hDD_dD010;
            const double FDPart3_29 = -hDD_dD002 + 2*hDD_dD020;
            const double FDPart3_30 = FDPart3_15*FDPart3_8;
            const double FDPart3_31 = (1.0/2.0)*FDPart3_27;
            const double FDPart3_32 = (1.0/2.0)*FDPart3_24;
            lambdaU0GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_13*(FDPart3_11*FDPart3_29 + FDPart3_14*hDD_dD000 + FDPart3_28*FDPart3_9) + FDPart3_17*(FDPart3_11*hDD_dD221 + FDPart3_12*FDPart3_14 + FDPart3_9*hDD_dD112) + FDPart3_19*(FDPart3_11*hDD_dD220 + FDPart3_14*hDD_dD002 + FDPart3_18*FDPart3_9) + FDPart3_21*(FDPart3_11*FDPart3_20 + FDPart3_14*hDD_dD001 + FDPart3_9*hDD_dD110) + FDPart3_24*(FDPart3_11*hDD_dD222 + FDPart3_14*FDPart3_23 + FDPart3_22*FDPart3_9) + FDPart3_27*(FDPart3_11*FDPart3_25 + FDPart3_14*FDPart3_26 + FDPart3_9*hDD_dD111);
            lambdaU1GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_13*(FDPart3_28*FDPart3_31 + FDPart3_29*FDPart3_30 + FDPart3_9*hDD_dD000) + FDPart3_17*(FDPart3_12*FDPart3_9 + FDPart3_30*hDD_dD221 + FDPart3_31*hDD_dD112) + FDPart3_19*(FDPart3_18*FDPart3_31 + FDPart3_30*hDD_dD220 + FDPart3_9*hDD_dD002) + FDPart3_21*(FDPart3_20*FDPart3_30 + FDPart3_31*hDD_dD110 + FDPart3_9*hDD_dD001) + FDPart3_24*(FDPart3_22*FDPart3_31 + FDPart3_23*FDPart3_9 + FDPart3_30*hDD_dD222) + FDPart3_27*(FDPart3_25*FDPart3_30 + FDPart3_26*FDPart3_9 + FDPart3_31*hDD_dD111);
            lambdaU2GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_13*(FDPart3_11*hDD_dD000 + FDPart3_28*FDPart3_30 + FDPart3_29*FDPart3_32) + FDPart3_17*(FDPart3_11*FDPart3_12 + FDPart3_30*hDD_dD112 + FDPart3_32*hDD_dD221) + FDPart3_19*(FDPart3_11*hDD_dD002 + FDPart3_18*FDPart3_30 + FDPart3_32*hDD_dD220) + FDPart3_21*(FDPart3_11*hDD_dD001 + FDPart3_20*FDPart3_32 + FDPart3_30*hDD_dD110) + FDPart3_24*(FDPart3_11*FDPart3_23 + FDPart3_22*FDPart3_30 + FDPart3_32*hDD_dD222) + FDPart3_27*(FDPart3_11*FDPart3_26 + FDPart3_25*FDPart3_32 + FDPart3_30*hDD_dD111);
          
        } // END LOOP: for (int i0 = cctk_nghostzones[0]; i0 < cctk_lsh[0]-cctk_nghostzones[0]; i0++)
      } // END LOOP: for (int i1 = cctk_nghostzones[1]; i1 < cctk_lsh[1]-cctk_nghostzones[1]; i1++)
    } // END LOOP: for (int i2 = cctk_nghostzones[2]; i2 < cctk_lsh[2]-cctk_nghostzones[2]; i2++)
  }
  if(FD_order == 4) {
#pragma omp parallel for
    for (int i2 = cctk_nghostzones[2]; i2 < cctk_lsh[2]-cctk_nghostzones[2]; i2++) {
      for (int i1 = cctk_nghostzones[1]; i1 < cctk_lsh[1]-cctk_nghostzones[1]; i1++) {
        for (int i0 = cctk_nghostzones[0]; i0 < cctk_lsh[0]-cctk_nghostzones[0]; i0++) {
            /*
             * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
             */
            const double hDD00_i0_i1_i2m2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-2)];
            const double hDD00_i0_i1_i2m1 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
            const double hDD00_i0_i1m2_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1-2,i2)];
            const double hDD00_i0_i1m1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
            const double hDD00_i0m2_i1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0-2,i1,i2)];
            const double hDD00_i0m1_i1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
            const double hDD00 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD00_i0p1_i1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
            const double hDD00_i0p2_i1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0+2,i1,i2)];
            const double hDD00_i0_i1p1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
            const double hDD00_i0_i1p2_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1+2,i2)];
            const double hDD00_i0_i1_i2p1 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
            const double hDD00_i0_i1_i2p2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+2)];
            const double hDD01_i0_i1_i2m2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-2)];
            const double hDD01_i0_i1_i2m1 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
            const double hDD01_i0_i1m2_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1-2,i2)];
            const double hDD01_i0_i1m1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
            const double hDD01_i0m2_i1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0-2,i1,i2)];
            const double hDD01_i0m1_i1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
            const double hDD01 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD01_i0p1_i1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
            const double hDD01_i0p2_i1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0+2,i1,i2)];
            const double hDD01_i0_i1p1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
            const double hDD01_i0_i1p2_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1+2,i2)];
            const double hDD01_i0_i1_i2p1 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
            const double hDD01_i0_i1_i2p2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+2)];
            const double hDD02_i0_i1_i2m2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-2)];
            const double hDD02_i0_i1_i2m1 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
            const double hDD02_i0_i1m2_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1-2,i2)];
            const double hDD02_i0_i1m1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
            const double hDD02_i0m2_i1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0-2,i1,i2)];
            const double hDD02_i0m1_i1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
            const double hDD02 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD02_i0p1_i1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
            const double hDD02_i0p2_i1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0+2,i1,i2)];
            const double hDD02_i0_i1p1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
            const double hDD02_i0_i1p2_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1+2,i2)];
            const double hDD02_i0_i1_i2p1 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
            const double hDD02_i0_i1_i2p2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+2)];
            const double hDD11_i0_i1_i2m2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-2)];
            const double hDD11_i0_i1_i2m1 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
            const double hDD11_i0_i1m2_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1-2,i2)];
            const double hDD11_i0_i1m1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
            const double hDD11_i0m2_i1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0-2,i1,i2)];
            const double hDD11_i0m1_i1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
            const double hDD11 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD11_i0p1_i1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
            const double hDD11_i0p2_i1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0+2,i1,i2)];
            const double hDD11_i0_i1p1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
            const double hDD11_i0_i1p2_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1+2,i2)];
            const double hDD11_i0_i1_i2p1 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
            const double hDD11_i0_i1_i2p2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+2)];
            const double hDD12_i0_i1_i2m2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-2)];
            const double hDD12_i0_i1_i2m1 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
            const double hDD12_i0_i1m2_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1-2,i2)];
            const double hDD12_i0_i1m1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
            const double hDD12_i0m2_i1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0-2,i1,i2)];
            const double hDD12_i0m1_i1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
            const double hDD12 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD12_i0p1_i1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
            const double hDD12_i0p2_i1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0+2,i1,i2)];
            const double hDD12_i0_i1p1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
            const double hDD12_i0_i1p2_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1+2,i2)];
            const double hDD12_i0_i1_i2p1 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
            const double hDD12_i0_i1_i2p2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+2)];
            const double hDD22_i0_i1_i2m2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-2)];
            const double hDD22_i0_i1_i2m1 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
            const double hDD22_i0_i1m2_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1-2,i2)];
            const double hDD22_i0_i1m1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
            const double hDD22_i0m2_i1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0-2,i1,i2)];
            const double hDD22_i0m1_i1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
            const double hDD22 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD22_i0p1_i1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
            const double hDD22_i0p2_i1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0+2,i1,i2)];
            const double hDD22_i0_i1p1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
            const double hDD22_i0_i1p2_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1+2,i2)];
            const double hDD22_i0_i1_i2p1 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
            const double hDD22_i0_i1_i2p2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+2)];
            const double FDPart1_Rational_2_3 = 2.0/3.0;
            const double FDPart1_Rational_1_12 = 1.0/12.0;
            const double hDD_dD000 = invdx0*(FDPart1_Rational_1_12*(hDD00_i0m2_i1_i2 - hDD00_i0p2_i1_i2) + FDPart1_Rational_2_3*(-hDD00_i0m1_i1_i2 + hDD00_i0p1_i1_i2));
            const double hDD_dD001 = invdx1*(FDPart1_Rational_1_12*(hDD00_i0_i1m2_i2 - hDD00_i0_i1p2_i2) + FDPart1_Rational_2_3*(-hDD00_i0_i1m1_i2 + hDD00_i0_i1p1_i2));
            const double hDD_dD002 = invdx2*(FDPart1_Rational_1_12*(hDD00_i0_i1_i2m2 - hDD00_i0_i1_i2p2) + FDPart1_Rational_2_3*(-hDD00_i0_i1_i2m1 + hDD00_i0_i1_i2p1));
            const double hDD_dD010 = invdx0*(FDPart1_Rational_1_12*(hDD01_i0m2_i1_i2 - hDD01_i0p2_i1_i2) + FDPart1_Rational_2_3*(-hDD01_i0m1_i1_i2 + hDD01_i0p1_i1_i2));
            const double hDD_dD011 = invdx1*(FDPart1_Rational_1_12*(hDD01_i0_i1m2_i2 - hDD01_i0_i1p2_i2) + FDPart1_Rational_2_3*(-hDD01_i0_i1m1_i2 + hDD01_i0_i1p1_i2));
            const double hDD_dD012 = invdx2*(FDPart1_Rational_1_12*(hDD01_i0_i1_i2m2 - hDD01_i0_i1_i2p2) + FDPart1_Rational_2_3*(-hDD01_i0_i1_i2m1 + hDD01_i0_i1_i2p1));
            const double hDD_dD020 = invdx0*(FDPart1_Rational_1_12*(hDD02_i0m2_i1_i2 - hDD02_i0p2_i1_i2) + FDPart1_Rational_2_3*(-hDD02_i0m1_i1_i2 + hDD02_i0p1_i1_i2));
            const double hDD_dD021 = invdx1*(FDPart1_Rational_1_12*(hDD02_i0_i1m2_i2 - hDD02_i0_i1p2_i2) + FDPart1_Rational_2_3*(-hDD02_i0_i1m1_i2 + hDD02_i0_i1p1_i2));
            const double hDD_dD022 = invdx2*(FDPart1_Rational_1_12*(hDD02_i0_i1_i2m2 - hDD02_i0_i1_i2p2) + FDPart1_Rational_2_3*(-hDD02_i0_i1_i2m1 + hDD02_i0_i1_i2p1));
            const double hDD_dD110 = invdx0*(FDPart1_Rational_1_12*(hDD11_i0m2_i1_i2 - hDD11_i0p2_i1_i2) + FDPart1_Rational_2_3*(-hDD11_i0m1_i1_i2 + hDD11_i0p1_i1_i2));
            const double hDD_dD111 = invdx1*(FDPart1_Rational_1_12*(hDD11_i0_i1m2_i2 - hDD11_i0_i1p2_i2) + FDPart1_Rational_2_3*(-hDD11_i0_i1m1_i2 + hDD11_i0_i1p1_i2));
            const double hDD_dD112 = invdx2*(FDPart1_Rational_1_12*(hDD11_i0_i1_i2m2 - hDD11_i0_i1_i2p2) + FDPart1_Rational_2_3*(-hDD11_i0_i1_i2m1 + hDD11_i0_i1_i2p1));
            const double hDD_dD120 = invdx0*(FDPart1_Rational_1_12*(hDD12_i0m2_i1_i2 - hDD12_i0p2_i1_i2) + FDPart1_Rational_2_3*(-hDD12_i0m1_i1_i2 + hDD12_i0p1_i1_i2));
            const double hDD_dD121 = invdx1*(FDPart1_Rational_1_12*(hDD12_i0_i1m2_i2 - hDD12_i0_i1p2_i2) + FDPart1_Rational_2_3*(-hDD12_i0_i1m1_i2 + hDD12_i0_i1p1_i2));
            const double hDD_dD122 = invdx2*(FDPart1_Rational_1_12*(hDD12_i0_i1_i2m2 - hDD12_i0_i1_i2p2) + FDPart1_Rational_2_3*(-hDD12_i0_i1_i2m1 + hDD12_i0_i1_i2p1));
            const double hDD_dD220 = invdx0*(FDPart1_Rational_1_12*(hDD22_i0m2_i1_i2 - hDD22_i0p2_i1_i2) + FDPart1_Rational_2_3*(-hDD22_i0m1_i1_i2 + hDD22_i0p1_i1_i2));
            const double hDD_dD221 = invdx1*(FDPart1_Rational_1_12*(hDD22_i0_i1m2_i2 - hDD22_i0_i1p2_i2) + FDPart1_Rational_2_3*(-hDD22_i0_i1m1_i2 + hDD22_i0_i1p1_i2));
            const double hDD_dD222 = invdx2*(FDPart1_Rational_1_12*(hDD22_i0_i1_i2m2 - hDD22_i0_i1_i2p2) + FDPart1_Rational_2_3*(-hDD22_i0_i1_i2m1 + hDD22_i0_i1_i2p1));
            /*
             * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
             */
            const double FDPart3_0 = hDD22 + 1;
            const double FDPart3_1 = -FDPart3_0*hDD01 + hDD02*hDD12;
            const double FDPart3_4 = hDD11 + 1;
            const double FDPart3_6 = hDD00 + 1;
            const double FDPart3_7 = (1.0/(FDPart3_0*FDPart3_4*FDPart3_6 - FDPart3_0*((hDD01)*(hDD01)) - FDPart3_4*((hDD02)*(hDD02)) - FDPart3_6*((hDD12)*(hDD12)) + 2*hDD01*hDD02*hDD12));
            const double FDPart3_8 = (1.0/2.0)*FDPart3_7;
            const double FDPart3_9 = FDPart3_1*FDPart3_8;
            const double FDPart3_10 = -FDPart3_4*hDD02 + hDD01*hDD12;
            const double FDPart3_11 = FDPart3_10*FDPart3_8;
            const double FDPart3_12 = hDD_dD012 + hDD_dD021 - hDD_dD120;
            const double FDPart3_13 = FDPart3_7*(FDPart3_0*FDPart3_4 - ((hDD12)*(hDD12)));
            const double FDPart3_14 = (1.0/2.0)*FDPart3_13;
            const double FDPart3_15 = -FDPart3_6*hDD12 + hDD01*hDD02;
            const double FDPart3_16 = 2*FDPart3_7;
            const double FDPart3_17 = FDPart3_15*FDPart3_16;
            const double FDPart3_18 = hDD_dD012 - hDD_dD021 + hDD_dD120;
            const double FDPart3_19 = FDPart3_10*FDPart3_16;
            const double FDPart3_20 = -hDD_dD012 + hDD_dD021 + hDD_dD120;
            const double FDPart3_21 = FDPart3_1*FDPart3_16;
            const double FDPart3_22 = 2*hDD_dD122 - hDD_dD221;
            const double FDPart3_23 = 2*hDD_dD022 - hDD_dD220;
            const double FDPart3_24 = FDPart3_7*(FDPart3_4*FDPart3_6 - ((hDD01)*(hDD01)));
            const double FDPart3_25 = -hDD_dD112 + 2*hDD_dD121;
            const double FDPart3_26 = 2*hDD_dD011 - hDD_dD110;
            const double FDPart3_27 = FDPart3_7*(FDPart3_0*FDPart3_6 - ((hDD02)*(hDD02)));
            const double FDPart3_28 = -hDD_dD001 + 2*hDD_dD010;
            const double FDPart3_29 = -hDD_dD002 + 2*hDD_dD020;
            const double FDPart3_30 = FDPart3_15*FDPart3_8;
            const double FDPart3_31 = (1.0/2.0)*FDPart3_27;
            const double FDPart3_32 = (1.0/2.0)*FDPart3_24;
            lambdaU0GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_13*(FDPart3_11*FDPart3_29 + FDPart3_14*hDD_dD000 + FDPart3_28*FDPart3_9) + FDPart3_17*(FDPart3_11*hDD_dD221 + FDPart3_12*FDPart3_14 + FDPart3_9*hDD_dD112) + FDPart3_19*(FDPart3_11*hDD_dD220 + FDPart3_14*hDD_dD002 + FDPart3_18*FDPart3_9) + FDPart3_21*(FDPart3_11*FDPart3_20 + FDPart3_14*hDD_dD001 + FDPart3_9*hDD_dD110) + FDPart3_24*(FDPart3_11*hDD_dD222 + FDPart3_14*FDPart3_23 + FDPart3_22*FDPart3_9) + FDPart3_27*(FDPart3_11*FDPart3_25 + FDPart3_14*FDPart3_26 + FDPart3_9*hDD_dD111);
            lambdaU1GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_13*(FDPart3_28*FDPart3_31 + FDPart3_29*FDPart3_30 + FDPart3_9*hDD_dD000) + FDPart3_17*(FDPart3_12*FDPart3_9 + FDPart3_30*hDD_dD221 + FDPart3_31*hDD_dD112) + FDPart3_19*(FDPart3_18*FDPart3_31 + FDPart3_30*hDD_dD220 + FDPart3_9*hDD_dD002) + FDPart3_21*(FDPart3_20*FDPart3_30 + FDPart3_31*hDD_dD110 + FDPart3_9*hDD_dD001) + FDPart3_24*(FDPart3_22*FDPart3_31 + FDPart3_23*FDPart3_9 + FDPart3_30*hDD_dD222) + FDPart3_27*(FDPart3_25*FDPart3_30 + FDPart3_26*FDPart3_9 + FDPart3_31*hDD_dD111);
            lambdaU2GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_13*(FDPart3_11*hDD_dD000 + FDPart3_28*FDPart3_30 + FDPart3_29*FDPart3_32) + FDPart3_17*(FDPart3_11*FDPart3_12 + FDPart3_30*hDD_dD112 + FDPart3_32*hDD_dD221) + FDPart3_19*(FDPart3_11*hDD_dD002 + FDPart3_18*FDPart3_30 + FDPart3_32*hDD_dD220) + FDPart3_21*(FDPart3_11*hDD_dD001 + FDPart3_20*FDPart3_32 + FDPart3_30*hDD_dD110) + FDPart3_24*(FDPart3_11*FDPart3_23 + FDPart3_22*FDPart3_30 + FDPart3_32*hDD_dD222) + FDPart3_27*(FDPart3_11*FDPart3_26 + FDPart3_25*FDPart3_32 + FDPart3_30*hDD_dD111);
          
        } // END LOOP: for (int i0 = cctk_nghostzones[0]; i0 < cctk_lsh[0]-cctk_nghostzones[0]; i0++)
      } // END LOOP: for (int i1 = cctk_nghostzones[1]; i1 < cctk_lsh[1]-cctk_nghostzones[1]; i1++)
    } // END LOOP: for (int i2 = cctk_nghostzones[2]; i2 < cctk_lsh[2]-cctk_nghostzones[2]; i2++)
  }

  ExtrapolateGammas(cctkGH, lambdaU0GF);
  ExtrapolateGammas(cctkGH, lambdaU1GF);
  ExtrapolateGammas(cctkGH, lambdaU2GF);
}
