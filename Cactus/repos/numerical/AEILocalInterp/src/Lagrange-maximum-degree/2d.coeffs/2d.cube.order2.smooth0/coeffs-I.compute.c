fp t27;
fp t25;
fp t18;
fp t28;
fp t21;
fp t34;
fp t33;
fp t22;
fp t17;
fp t23;
fp t32;
fp t15;
fp t31;
fp t14;
fp t30;
fp t24;
fp t20;
fp t29;
fp t19;
fp t16;
fp t13;
fp t12;
      t27 = x*x;
      t25 = RATIONAL(1.0,6.0);
      t18 = t25*t27;
      t28 = y*y;
      t21 = t25*t28;
      t34 = t18+t21+RATIONAL(-1.0,9.0);
      t33 = x*y;
      t22 = RATIONAL(-1.0,3.0);
      t17 = t22*t27;
      t23 = RATIONAL(2.0,9.0);
      t32 = t17+t21+t23;
      t15 = t22*t28;
      t31 = t15+t18+t23;
      t14 = t25*y;
      t30 = t14+t34;
      t24 = RATIONAL(-1.0,6.0);
      t20 = t24*y;
      t29 = t20+t34;
      t19 = t24*x;
      t16 = t25*x;
      t13 = RATIONAL(-1.0,4.0)*t33;
      t12 = RATIONAL(1.0,4.0)*t33;
      coeffs_I->coeff_m1_m1 = t19+t12+t29;
      coeffs_I->coeff_0_m1 = t20+t32;
      coeffs_I->coeff_p1_m1 = t13+t16+t29;
      coeffs_I->coeff_m1_0 = t19+t31;
      coeffs_I->coeff_0_0 = t15+t17+RATIONAL(5.0,9.0);
      coeffs_I->coeff_p1_0 = t16+t31;
      coeffs_I->coeff_m1_p1 = t19+t13+t30;
      coeffs_I->coeff_0_p1 = t14+t32;
      coeffs_I->coeff_p1_p1 = t16+t12+t30;