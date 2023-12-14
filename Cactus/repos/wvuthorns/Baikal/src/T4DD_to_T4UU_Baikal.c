#include "math.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
/*
 * Compute T4UU from T4DD (provided by TmunuBase),
 * using BSSN quantities as inputs for the 4D raising operation
 * 
 * WARNING: Do not enable SIMD here, as it is not guaranteed that
 *          cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2] is a multiple of
 *          SIMD_width!
 */
void T4DD_to_T4UU_Baikal(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS_T4DD_to_T4UU_Baikal;
  DECLARE_CCTK_PARAMETERS;

  // First read in TmunuBase gridfunctions
  const CCTK_REAL *restrict T4DD00GF = eTtt;
  const CCTK_REAL *restrict T4DD01GF = eTtx;
  const CCTK_REAL *restrict T4DD02GF = eTty;
  const CCTK_REAL *restrict T4DD03GF = eTtz;
  const CCTK_REAL *restrict T4DD11GF = eTxx;
  const CCTK_REAL *restrict T4DD12GF = eTxy;
  const CCTK_REAL *restrict T4DD13GF = eTxz;
  const CCTK_REAL *restrict T4DD22GF = eTyy;
  const CCTK_REAL *restrict T4DD23GF = eTyz;
  const CCTK_REAL *restrict T4DD33GF = eTzz;

  // TmunuBase provides T_{\alpha \beta}, but BSSN needs T^{\alpha \beta}.
  //  Here we perform the raising operation using BSSN quantities as input:
  #pragma omp parallel for
  for (int i2 = 0; i2 < cctk_lsh[2]; i2++) {
    for (int i1 = 0; i1 < cctk_lsh[1]; i1++) {
      for (int i0 = 0; i0 < cctk_lsh[0]; i0++) {
        /*
         * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
         */
        const double hDD00 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double hDD01 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double hDD02 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double hDD11 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double hDD12 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double hDD22 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double vetU0 = vetU0GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double vetU1 = vetU1GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double vetU2 = vetU2GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double cf = cfGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double alpha = alphaGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double T4DD00 = T4DD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double T4DD01 = T4DD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double T4DD02 = T4DD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double T4DD03 = T4DD03GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double T4DD11 = T4DD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double T4DD12 = T4DD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double T4DD13 = T4DD13GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double T4DD22 = T4DD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double T4DD23 = T4DD23GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double T4DD33 = T4DD33GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        /*
         * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
         */
        const double FDPart3_0 = (1.0/((alpha)*(alpha)*(alpha)*(alpha)));
        const double FDPart3_1 = FDPart3_0*T4DD00;
        const double FDPart3_2 = 2*FDPart3_0;
        const double FDPart3_3 = T4DD01*vetU0;
        const double FDPart3_4 = T4DD02*vetU1;
        const double FDPart3_5 = T4DD03*vetU2;
        const double FDPart3_6 = ((vetU0)*(vetU0));
        const double FDPart3_8 = ((vetU1)*(vetU1));
        const double FDPart3_10 = ((vetU2)*(vetU2));
        const double FDPart3_13 = FDPart3_1*vetU0;
        const double FDPart3_15 = (1.0/((alpha)*(alpha)));
        const double FDPart3_16 = FDPart3_15*vetU0;
        const double FDPart3_17 = (1.0/((cf)*(cf)*(cf)*(cf)));
        const double FDPart3_18 = hDD22 + 1;
        const double FDPart3_19 = pow(cf, -6);
        const double FDPart3_22 = hDD11 + 1;
        const double FDPart3_24 = hDD00 + 1;
        const double FDPart3_25 = (1.0/(FDPart3_18*FDPart3_19*FDPart3_22*FDPart3_24 - FDPart3_18*FDPart3_19*((hDD01)*(hDD01)) - FDPart3_19*FDPart3_22*((hDD02)*(hDD02)) - FDPart3_19*FDPart3_24*((hDD12)*(hDD12)) + 2*FDPart3_19*hDD01*hDD02*hDD12));
        const double FDPart3_26 = -FDPart3_16*vetU1 + FDPart3_25*(-FDPart3_17*FDPart3_18*hDD01 + FDPart3_17*hDD02*hDD12);
        const double FDPart3_27 = FDPart3_15*T4DD02;
        const double FDPart3_28 = FDPart3_26*FDPart3_27;
        const double FDPart3_31 = -FDPart3_16*vetU2 + FDPart3_25*(-FDPart3_17*FDPart3_22*hDD02 + FDPart3_17*hDD01*hDD12);
        const double FDPart3_32 = FDPart3_15*T4DD03;
        const double FDPart3_33 = FDPart3_31*FDPart3_32;
        const double FDPart3_34 = FDPart3_26*T4DD12;
        const double FDPart3_35 = FDPart3_31*T4DD13;
        const double FDPart3_36 = FDPart3_15*vetU1;
        const double FDPart3_37 = FDPart3_26*T4DD22;
        const double FDPart3_38 = FDPart3_31*T4DD23;
        const double FDPart3_39 = FDPart3_26*T4DD23;
        const double FDPart3_40 = FDPart3_15*vetU2;
        const double FDPart3_41 = FDPart3_31*T4DD33;
        const double FDPart3_42 = -FDPart3_15*FDPart3_6 + FDPart3_25*(FDPart3_17*FDPart3_18*FDPart3_22 - FDPart3_17*((hDD12)*(hDD12)));
        const double FDPart3_43 = FDPart3_15*T4DD01;
        const double FDPart3_44 = FDPart3_42*FDPart3_43;
        const double FDPart3_45 = FDPart3_42*T4DD11;
        const double FDPart3_46 = FDPart3_42*T4DD12;
        const double FDPart3_47 = FDPart3_42*T4DD13;
        const double FDPart3_50 = FDPart3_26*FDPart3_43;
        const double FDPart3_51 = FDPart3_25*(-FDPart3_17*FDPart3_24*hDD12 + FDPart3_17*hDD01*hDD02) - FDPart3_36*vetU2;
        const double FDPart3_52 = FDPart3_32*FDPart3_51;
        const double FDPart3_54 = FDPart3_16*FDPart3_51;
        const double FDPart3_55 = FDPart3_26*T4DD13;
        const double FDPart3_58 = -FDPart3_15*FDPart3_8 + FDPart3_25*(FDPart3_17*FDPart3_18*FDPart3_24 - FDPart3_17*((hDD02)*(hDD02)));
        const double FDPart3_59 = FDPart3_27*FDPart3_58;
        const double FDPart3_60 = FDPart3_58*T4DD12;
        const double FDPart3_62 = FDPart3_58*T4DD23;
        const double FDPart3_64 = FDPart3_31*FDPart3_43;
        const double FDPart3_65 = FDPart3_27*FDPart3_51;
        const double FDPart3_67 = -FDPart3_10*FDPart3_15 + FDPart3_25*(FDPart3_17*FDPart3_22*FDPart3_24 - FDPart3_17*((hDD01)*(hDD01)));
        const double FDPart3_68 = FDPart3_32*FDPart3_67;
        const double FDPart3_71 = ((FDPart3_26)*(FDPart3_26));
        const double FDPart3_72 = ((FDPart3_31)*(FDPart3_31));
        const double FDPart3_73 = 2*FDPart3_26;
        const double FDPart3_76 = FDPart3_26*FDPart3_31;
        const double FDPart3_77 = ((FDPart3_51)*(FDPart3_51));
        const double FDPart3_78 = 2*vetU1;
        const double FDPart3_79 = 2*FDPart3_51;
        const double FDPart3_80 = 2*vetU2;
        T4UU00GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_0*FDPart3_10*T4DD33 + FDPart3_0*FDPart3_6*T4DD11 + FDPart3_0*FDPart3_8*T4DD22 + FDPart3_1 - FDPart3_2*FDPart3_3 - FDPart3_2*FDPart3_4 - FDPart3_2*FDPart3_5 + FDPart3_2*T4DD12*vetU0*vetU1 + FDPart3_2*T4DD13*vetU0*vetU2 + FDPart3_2*T4DD23*vetU1*vetU2;
        T4UU01GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_0*FDPart3_4*vetU0 + FDPart3_0*FDPart3_5*vetU0 + FDPart3_0*FDPart3_6*T4DD01 - FDPart3_13 + FDPart3_16*FDPart3_34 + FDPart3_16*FDPart3_35 + FDPart3_16*FDPart3_45 - FDPart3_28 - FDPart3_33 + FDPart3_36*FDPart3_37 + FDPart3_36*FDPart3_38 + FDPart3_36*FDPart3_46 + FDPart3_39*FDPart3_40 + FDPart3_40*FDPart3_41 + FDPart3_40*FDPart3_47 - FDPart3_44;
        T4UU02GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_0*FDPart3_3*vetU1 + FDPart3_0*FDPart3_5*vetU1 + FDPart3_0*FDPart3_8*T4DD02 - FDPart3_1*vetU1 + FDPart3_16*FDPart3_26*T4DD11 + FDPart3_16*FDPart3_60 + FDPart3_34*FDPart3_36 + FDPart3_36*FDPart3_51*T4DD23 + FDPart3_36*FDPart3_58*T4DD22 + FDPart3_40*FDPart3_51*T4DD33 + FDPart3_40*FDPart3_55 + FDPart3_40*FDPart3_62 - FDPart3_50 - FDPart3_52 + FDPart3_54*T4DD13 - FDPart3_59;
        T4UU03GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_0*FDPart3_10*T4DD03 + FDPart3_0*FDPart3_3*vetU2 + FDPart3_0*FDPart3_4*vetU2 - FDPart3_1*vetU2 + FDPart3_16*FDPart3_31*T4DD11 + FDPart3_16*FDPart3_67*T4DD13 + FDPart3_31*FDPart3_36*T4DD12 + FDPart3_35*FDPart3_40 + FDPart3_36*FDPart3_51*T4DD22 + FDPart3_36*FDPart3_67*T4DD23 + FDPart3_40*FDPart3_51*T4DD23 + FDPart3_40*FDPart3_67*T4DD33 + FDPart3_54*T4DD12 - FDPart3_64 - FDPart3_65 - FDPart3_68;
        T4UU11GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_1*FDPart3_6 + 2*FDPart3_16*FDPart3_31*T4DD03 + FDPart3_16*FDPart3_73*T4DD02 + 2*FDPart3_31*FDPart3_47 + ((FDPart3_42)*(FDPart3_42))*T4DD11 + 2*FDPart3_44*vetU0 + FDPart3_46*FDPart3_73 + FDPart3_71*T4DD22 + FDPart3_72*T4DD33 + 2*FDPart3_76*T4DD23;
        T4UU12GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_13*vetU1 + FDPart3_16*FDPart3_58*T4DD02 + FDPart3_26*FDPart3_45 + FDPart3_28*vetU1 + FDPart3_33*vetU1 + FDPart3_37*FDPart3_58 + FDPart3_38*FDPart3_58 + FDPart3_39*FDPart3_51 + FDPart3_41*FDPart3_51 + FDPart3_44*vetU1 + FDPart3_46*FDPart3_58 + FDPart3_47*FDPart3_51 + FDPart3_50*vetU0 + FDPart3_54*T4DD03 + FDPart3_71*T4DD12 + FDPart3_76*T4DD13;
        T4UU13GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_13*vetU2 + FDPart3_16*FDPart3_67*T4DD03 + FDPart3_28*vetU2 + FDPart3_31*FDPart3_45 + FDPart3_33*vetU2 + FDPart3_37*FDPart3_51 + FDPart3_38*FDPart3_51 + FDPart3_39*FDPart3_67 + FDPart3_41*FDPart3_67 + FDPart3_44*vetU2 + FDPart3_46*FDPart3_51 + FDPart3_47*FDPart3_67 + FDPart3_54*T4DD02 + FDPart3_64*vetU0 + FDPart3_72*T4DD13 + FDPart3_76*T4DD12;
        T4UU22GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_1*FDPart3_8 + FDPart3_50*FDPart3_78 + FDPart3_52*FDPart3_78 + FDPart3_55*FDPart3_79 + ((FDPart3_58)*(FDPart3_58))*T4DD22 + FDPart3_59*FDPart3_78 + FDPart3_60*FDPart3_73 + FDPart3_62*FDPart3_79 + FDPart3_71*T4DD11 + FDPart3_77*T4DD33;
        T4UU23GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_1*vetU1*vetU2 + FDPart3_31*FDPart3_60 + FDPart3_34*FDPart3_51 + FDPart3_35*FDPart3_51 + FDPart3_50*vetU2 + FDPart3_51*FDPart3_58*T4DD22 + FDPart3_51*FDPart3_67*T4DD33 + FDPart3_52*vetU2 + FDPart3_55*FDPart3_67 + FDPart3_59*vetU2 + FDPart3_62*FDPart3_67 + FDPart3_64*vetU1 + FDPart3_65*vetU1 + FDPart3_68*vetU1 + FDPart3_76*T4DD11 + FDPart3_77*T4DD23;
        T4UU33GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_1*FDPart3_10 + FDPart3_31*FDPart3_79*T4DD12 + 2*FDPart3_35*FDPart3_67 + FDPart3_64*FDPart3_80 + FDPart3_65*FDPart3_80 + ((FDPart3_67)*(FDPart3_67))*T4DD33 + FDPart3_67*FDPart3_79*T4DD23 + FDPart3_68*FDPart3_80 + FDPart3_72*T4DD11 + FDPart3_77*T4DD22;
        
      } // END LOOP: for (int i0 = 0; i0 < cctk_lsh[0]; i0++)
    } // END LOOP: for (int i1 = 0; i1 < cctk_lsh[1]; i1++)
  } // END LOOP: for (int i2 = 0; i2 < cctk_lsh[2]; i2++)
}
