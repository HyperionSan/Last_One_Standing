fp t110;
fp t122;
fp t121;
fp t113;
fp t108;
fp t120;
fp t112;
fp t106;
fp t119;
fp t105;
fp t118;
fp t117;
fp t116;
fp t107;
fp t115;
fp t109;
      t110 = RATIONAL(1.0,9.0)*y;
      t122 = RATIONAL(-1.0,18.0)+t110;
      t121 = t110+RATIONAL(1.0,18.0);
      t113 = RATIONAL(1.0,12.0);
      t108 = t113*x;
      t120 = t108+t122;
      t112 = RATIONAL(-1.0,12.0);
      t106 = t112*x;
      t119 = t106+t122;
      t105 = t112*z;
      t118 = t105+t122;
      t117 = t106+t121;
      t116 = t108+t121;
      t107 = t113*z;
      t115 = t107+t121;
      t109 = RATIONAL(-2.0,9.0)*y;
      coeffs_dy->coeff_m1_m1_m1 = t107+t120;
      coeffs_dy->coeff_0_m1_m1 = t107+t122;
      coeffs_dy->coeff_p1_m1_m1 = t107+t119;
      coeffs_dy->coeff_m1_0_m1 = t109;
      coeffs_dy->coeff_0_0_m1 = t109;
      coeffs_dy->coeff_p1_0_m1 = t109;
      coeffs_dy->coeff_m1_p1_m1 = t105+t117;
      coeffs_dy->coeff_0_p1_m1 = t105+t121;
      coeffs_dy->coeff_p1_p1_m1 = t105+t116;
      coeffs_dy->coeff_m1_m1_0 = t120;
      coeffs_dy->coeff_0_m1_0 = t122;
      coeffs_dy->coeff_p1_m1_0 = t119;
      coeffs_dy->coeff_m1_0_0 = t109;
      coeffs_dy->coeff_0_0_0 = t109;
      coeffs_dy->coeff_p1_0_0 = t109;
      coeffs_dy->coeff_m1_p1_0 = t117;
      coeffs_dy->coeff_0_p1_0 = t121;
      coeffs_dy->coeff_p1_p1_0 = t116;
      coeffs_dy->coeff_m1_m1_p1 = t108+t118;
      coeffs_dy->coeff_0_m1_p1 = t118;
      coeffs_dy->coeff_p1_m1_p1 = t106+t118;
      coeffs_dy->coeff_m1_0_p1 = t109;
      coeffs_dy->coeff_0_0_p1 = t109;
      coeffs_dy->coeff_p1_0_p1 = t109;
      coeffs_dy->coeff_m1_p1_p1 = t106+t115;
      coeffs_dy->coeff_0_p1_p1 = t115;
      coeffs_dy->coeff_p1_p1_p1 = t108+t115;