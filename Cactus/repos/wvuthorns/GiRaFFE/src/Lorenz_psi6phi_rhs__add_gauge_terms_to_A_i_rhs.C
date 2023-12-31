static inline CCTK_REAL avg(const CCTK_REAL f[PLUS2+1][PLUS2+1][PLUS2+1],const int imin,const int imax, const int jmin,const int jmax, const int kmin,const int kmax);

static void Lorenz_psi6phi_rhs__add_gauge_terms_to_A_i_rhs(const cGH *cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,const CCTK_REAL *dX,CCTK_REAL **in_vars,const CCTK_REAL *psi6phi,
                                                           CCTK_REAL *shiftx_iphjphkph,CCTK_REAL *shifty_iphjphkph,CCTK_REAL *shiftz_iphjphkph,
                                                           CCTK_REAL *alpha_iphjphkph,CCTK_REAL *alpha_Phi_minus_betaj_A_j_iphjphkph,CCTK_REAL *alpha_sqrtg_Ax_interp,
                                                           CCTK_REAL *alpha_sqrtg_Ay_interp,CCTK_REAL *alpha_sqrtg_Az_interp,
                                                           CCTK_REAL *psi6phi_rhs,CCTK_REAL *Ax_rhs,CCTK_REAL *Ay_rhs,CCTK_REAL *Az_rhs) {
  DECLARE_CCTK_PARAMETERS;
  /* Compute
   * \partial_t psi6phi = -\partial_j ( \alpha \sqrt{\gamma} A^j - \beta^j psi6phi)
   *    (Eq 13 of http://arxiv.org/pdf/1110.4633.pdf), using Lorenz gauge.
   * Note that the RHS consists of a shift advection term on psi6phi and
   *    a term depending on the vector potential.
   * psi6phi is defined at (i+1/2,j+1/2,k+1/2), but instead of reconstructing
   *    to compute the RHS of \partial_t psi6phi, we instead use standard
   *    interpolations.
   */

  const CCTK_REAL dXm1=1.0/dX[0];
  const CCTK_REAL dYm1=1.0/dX[1];
  const CCTK_REAL dZm1=1.0/dX[2];

  // The stencil here is {-1,1},{-1,1},{-1,1} for x,y,z directions, respectively.
  //     Note that ALL input variables are defined at ALL gridpoints, so no
  //     worries about ghostzones.
#pragma omp parallel for
  for(int k=1;k<cctk_lsh[2]-1;k++) for(int j=1;j<cctk_lsh[1]-1;j++) for(int i=1;i<cctk_lsh[0]-1;i++) {
        const int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        CCTK_REAL INTERP_VARS[MAXNUMINTERP][PLUS2+1][PLUS2+1][PLUS2+1];

        // First compute \partial_j \alpha \sqrt{\gamma} A^j (RHS of \partial_i psi6phi)
        // FIXME: Would be much cheaper & easier to unstagger A_i, raise, then interpolate A^i.
        //        However, we keep it this way to be completely compatible with the original
        //        Illinois GRMHD thorn, called mhd_evolve.
        //
        //Step 1) j=x: Need to raise A_i, but to do that, we must have all variables at the same gridpoints:
        // The goal is to compute \partial_j (\alpha \sqrt{\gamma} A^j) at (i+1/2,j+1/2,k+1/2)
        //    We do this by first interpolating (RHS1x) = (\alpha \sqrt{\gamma} A^x) at
        //    (i,j+1/2,k+1/2)and (i+1,j+1/2,k+1/2), then taking \partial_x (RHS1x) =
        //    [ RHS1x(i+1,j+1/2,k+1/2) - RHS1x(i,j+1/2,k+1/2) ]/dX.
        // First bring gup's, psi, and alpha to (i,j+1/2,k+1/2):
        int num_vars_to_interp;
        int vars_to_interpolate[MAXNUMINTERP] = {GUPXXI,GUPXYI,GUPXZI,GUPYYI,GUPYZI,GUPZZI,LAPM1I,PSII,SHIFTXI,SHIFTYI,SHIFTZI};
        num_vars_to_interp = 11;
        // We may set interp_limits to be more general than we need.
        int interp_limits[6] = {-1,1,-1,1,-1,1}; SET_INDEX_ARRAYS_3DBLOCK(interp_limits);
        //SET_INDEX_ARRAYS_3DBLOCK(interp_limits);
        for(int ww=0;ww<num_vars_to_interp;ww++) {
          int whichvar=vars_to_interpolate[ww];
          // Read in variable at interp. stencil points from main memory, store in INTERP_VARS.
          for(int kk=PLUS0;kk<=PLUS1;kk++) for(int jj=PLUS0;jj<=PLUS1;jj++) for(int ii=PLUS0;ii<=PLUS1;ii++) {
                INTERP_VARS[whichvar][kk][jj][ii] = in_vars[whichvar][index_arr_3DB[kk][jj][ii]]; }
        }

        // Next set \alpha at (i+1/2,j+1/2,k+1/2). Will come in handy when computing damping term later.
        alpha_iphjphkph[index] = avg(INTERP_VARS[LAPM1I] , PLUS0,PLUS1, PLUS0,PLUS1, PLUS0,PLUS1)+1.0;

        //A_x needs a stencil s.t. interp_limits={0,1,-1,1,-1,1}:
        for(int kk=MINUS1;kk<=PLUS1;kk++) for(int jj=MINUS1;jj<=PLUS1;jj++) for(int ii=PLUS0;ii<=PLUS1;ii++) {
              INTERP_VARS[A_XI][kk][jj][ii] = in_vars[A_XI][index_arr_3DB[kk][jj][ii]]; }
        //A_y needs a stencil s.t. interp_limits={-1,1,0,1,-1,1}:
        for(int kk=MINUS1;kk<=PLUS1;kk++) for(int jj=PLUS0;jj<=PLUS1;jj++) for(int ii=MINUS1;ii<=PLUS1;ii++) {
              INTERP_VARS[A_YI][kk][jj][ii] = in_vars[A_YI][index_arr_3DB[kk][jj][ii]]; }
        //A_z needs a stencil s.t. interp_limits={-1,1,-1,1,0,1}:
        for(int kk=PLUS0;kk<=PLUS1;kk++) for(int jj=MINUS1;jj<=PLUS1;jj++) for(int ii=MINUS1;ii<=PLUS1;ii++) {
              INTERP_VARS[A_ZI][kk][jj][ii] = in_vars[A_ZI][index_arr_3DB[kk][jj][ii]]; }

        // FIRST DO A^X TERM (interpolate to (i,j+1/2,k+1/2) )
        // \alpha \sqrt{\gamma} A^x = \alpha psi^6 A^x (RHS of \partial_i psi6phi)
        // Note that gupij is \tilde{\gamma}^{ij}, so we need to multiply by \psi^{-4}.
        const CCTK_REAL gupxx_jphkph = avg(INTERP_VARS[GUPXXI], PLUS0,PLUS0, PLUS0,PLUS1, PLUS0,PLUS1);
        const CCTK_REAL gupxy_jphkph = avg(INTERP_VARS[GUPXYI], PLUS0,PLUS0, PLUS0,PLUS1, PLUS0,PLUS1);
        const CCTK_REAL gupxz_jphkph = avg(INTERP_VARS[GUPXZI], PLUS0,PLUS0, PLUS0,PLUS1, PLUS0,PLUS1);

        for(int kk=PLUS0;kk<=PLUS1;kk++) for(int jj=PLUS0;jj<=PLUS1;jj++) for(int ii=PLUS0;ii<=PLUS1;ii++) {
              const CCTK_REAL Psi2 = INTERP_VARS[PSII][kk][jj][ii]*INTERP_VARS[PSII][kk][jj][ii];
              const CCTK_REAL alpha = INTERP_VARS[LAPM1I][kk][jj][ii]+1.0;
              INTERP_VARS[LAPSE_PSI2I][kk][jj][ii]=alpha*Psi2;
              INTERP_VARS[LAPSE_OVER_PSI6I][kk][jj][ii]=alpha/(Psi2*Psi2*Psi2);
            }

        const CCTK_REAL lapse_Psi2_jphkph = avg(INTERP_VARS[LAPSE_PSI2I], PLUS0,PLUS0, PLUS0,PLUS1, PLUS0,PLUS1);

        const CCTK_REAL A_x_jphkph   = avg(INTERP_VARS[A_XI], PLUS0,PLUS0, PLUS0,PLUS0, PLUS0,PLUS0); // @ (i,j+1/2,k+1/2)
        const CCTK_REAL A_y_jphkph   = avg(INTERP_VARS[A_YI],MINUS1,PLUS0, PLUS0,PLUS1, PLUS0,PLUS0); // @ (i+1/2,j,k+1/2)
        const CCTK_REAL A_z_jphkph   = avg(INTERP_VARS[A_ZI],MINUS1,PLUS0, PLUS0,PLUS0, PLUS0,PLUS1); // @ (i+1/2,j+1/2,k)

        alpha_sqrtg_Ax_interp[index] = lapse_Psi2_jphkph*
          ( gupxx_jphkph*A_x_jphkph + gupxy_jphkph*A_y_jphkph + gupxz_jphkph*A_z_jphkph );


        // DO A^Y TERM (interpolate to (i+1/2,j,k+1/2) )
        // \alpha \sqrt{\gamma} A^y = \alpha psi^6 A^y (RHS of \partial_i psi6phi)
        // Note that gupij is \tilde{\gamma}^{ij}, so we need to multiply by \psi^{-4}.
        const CCTK_REAL gupxy_iphkph = avg(INTERP_VARS[GUPXYI], PLUS0,PLUS1, PLUS0,PLUS0, PLUS0,PLUS1);
        const CCTK_REAL gupyy_iphkph = avg(INTERP_VARS[GUPYYI], PLUS0,PLUS1, PLUS0,PLUS0, PLUS0,PLUS1);
        const CCTK_REAL gupyz_iphkph = avg(INTERP_VARS[GUPYZI], PLUS0,PLUS1, PLUS0,PLUS0, PLUS0,PLUS1);

        const CCTK_REAL lapse_Psi2_iphkph = avg(INTERP_VARS[LAPSE_PSI2I], PLUS0,PLUS1, PLUS0,PLUS0, PLUS0,PLUS1);
        //CCTK_REAL lapse_iphkph = avg(INTERP_VARS[LAPM1I], PLUS0,PLUS1, PLUS0,PLUS0, PLUS0,PLUS1)+1.0;
        //CCTK_REAL psi_iphkph   = avg(INTERP_VARS[PSII  ], PLUS0,PLUS1, PLUS0,PLUS0, PLUS0,PLUS1);

        const CCTK_REAL A_x_iphkph   = avg(INTERP_VARS[A_XI], PLUS0,PLUS1,MINUS1,PLUS0, PLUS0,PLUS0); // @ (i,j+1/2,k+1/2)
        const CCTK_REAL A_y_iphkph   = avg(INTERP_VARS[A_YI], PLUS0,PLUS0, PLUS0,PLUS0, PLUS0,PLUS0); // @ (i+1/2,j,k+1/2)
        const CCTK_REAL A_z_iphkph   = avg(INTERP_VARS[A_ZI], PLUS0,PLUS0,MINUS1,PLUS0, PLUS0,PLUS1); // @ (i+1/2,j+1/2,k)

        alpha_sqrtg_Ay_interp[index] = lapse_Psi2_iphkph*
          ( gupxy_iphkph*A_x_iphkph + gupyy_iphkph*A_y_iphkph + gupyz_iphkph*A_z_iphkph );

        // DO A^Z TERM (interpolate to (i+1/2,j+1/2,k) )
        // \alpha \sqrt{\gamma} A^z = \alpha psi^6 A^z (RHS of \partial_i psi6phi)
        // Note that gupij is \tilde{\gamma}^{ij}, so we need to multiply by \psi^{-4}.
        const CCTK_REAL gupxz_iphjph = avg(INTERP_VARS[GUPXZI], PLUS0,PLUS1, PLUS0,PLUS1, PLUS0,PLUS0);
        const CCTK_REAL gupyz_iphjph = avg(INTERP_VARS[GUPYZI], PLUS0,PLUS1, PLUS0,PLUS1, PLUS0,PLUS0);
        const CCTK_REAL gupzz_iphjph = avg(INTERP_VARS[GUPZZI], PLUS0,PLUS1, PLUS0,PLUS1, PLUS0,PLUS0);
        //CCTK_REAL lapse_iphjph = avg(INTERP_VARS[LAPM1I], PLUS0,PLUS1, PLUS0,PLUS1, PLUS0,PLUS0)+1.0;
        //CCTK_REAL psi_iphjph   = avg(INTERP_VARS[PSII  ], PLUS0,PLUS1, PLUS0,PLUS1, PLUS0,PLUS0);

        const CCTK_REAL lapse_Psi2_iphjph = avg(INTERP_VARS[LAPSE_PSI2I], PLUS0,PLUS1, PLUS0,PLUS1, PLUS0,PLUS0);

        const CCTK_REAL A_x_iphjph   = avg(INTERP_VARS[A_XI], PLUS0,PLUS1, PLUS0,PLUS0,MINUS1,PLUS0); // @ (i,j+1/2,k+1/2)
        const CCTK_REAL A_y_iphjph   = avg(INTERP_VARS[A_YI], PLUS0,PLUS0, PLUS0,PLUS1,MINUS1,PLUS0); // @ (i+1/2,j,k+1/2)
        const CCTK_REAL A_z_iphjph   = avg(INTERP_VARS[A_ZI], PLUS0,PLUS0, PLUS0,PLUS0, PLUS0,PLUS0); // @ (i+1/2,j+1/2,k)

        alpha_sqrtg_Az_interp[index] = lapse_Psi2_iphjph*
          ( gupxz_iphjph*A_x_iphjph + gupyz_iphjph*A_y_iphjph + gupzz_iphjph*A_z_iphjph );


        // Next set \alpha \Phi - \beta^j A_j at (i+1/2,j+1/2,k+1/2):
        //   We add a "L" suffix to shifti_iphjphkph to denote "LOCAL", as we set
        //      shifti_iphjphkph[] gridfunction below.
        const CCTK_REAL shiftx_iphjphkphL = avg(INTERP_VARS[SHIFTXI], PLUS0,PLUS1, PLUS0,PLUS1, PLUS0,PLUS1);
        const CCTK_REAL shifty_iphjphkphL = avg(INTERP_VARS[SHIFTYI], PLUS0,PLUS1, PLUS0,PLUS1, PLUS0,PLUS1);
        const CCTK_REAL shiftz_iphjphkphL = avg(INTERP_VARS[SHIFTZI], PLUS0,PLUS1, PLUS0,PLUS1, PLUS0,PLUS1);
        const CCTK_REAL lapse_over_Psi6_iphjphkphL = avg(INTERP_VARS[LAPSE_OVER_PSI6I], PLUS0,PLUS1, PLUS0,PLUS1, PLUS0,PLUS1);
        //CCTK_REAL psi_iphjphkph = avg(INTERP_VARS[PSII  ], PLUS0,PLUS1, PLUS0,PLUS1, PLUS0,PLUS1);
        //CCTK_REAL psi2_iphjphkph= psi_iphjphkph*psi_iphjphkph;
        //CCTK_REAL psi6_iphjphkph= psi2_iphjphkph*psi2_iphjphkph*psi2_iphjphkph;
        const CCTK_REAL A_x_iphjphkph = avg(INTERP_VARS[A_XI], PLUS0,PLUS1, PLUS0,PLUS0, PLUS0,PLUS0); // @ (i,j+1/2,k+1/2)
        const CCTK_REAL A_y_iphjphkph = avg(INTERP_VARS[A_YI], PLUS0,PLUS0, PLUS0,PLUS1, PLUS0,PLUS0); // @ (i+1/2,j,k+1/2)
        const CCTK_REAL A_z_iphjphkph = avg(INTERP_VARS[A_ZI], PLUS0,PLUS0, PLUS0,PLUS0, PLUS0,PLUS1); // @ (i+1/2,j+1/2,k)

        alpha_Phi_minus_betaj_A_j_iphjphkph[index] = psi6phi[index]*lapse_over_Psi6_iphjphkphL
          - (shiftx_iphjphkphL*A_x_iphjphkph + shifty_iphjphkphL*A_y_iphjphkph + shiftz_iphjphkphL*A_z_iphjphkph);

        // Finally, save shifti_iphjphkph, for \partial_j \beta^j psi6phi
        shiftx_iphjphkph[index]=shiftx_iphjphkphL;
        shifty_iphjphkph[index]=shifty_iphjphkphL;
        shiftz_iphjphkph[index]=shiftz_iphjphkphL;
      }

  // This loop requires two additional ghostzones in every direction. Hence the following loop definition:
#pragma omp parallel for
  for(int k=cctk_nghostzones[2];k<cctk_lsh[2]-cctk_nghostzones[2];k++) for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-cctk_nghostzones[1];j++) for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-cctk_nghostzones[0];i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // \partial_t A_i = [reconstructed stuff] + [gauge stuff],
        //    where [gauge stuff] = -\partial_i (\alpha \Phi - \beta^j A_j)
        const CCTK_REAL alpha_Phi_minus_betaj_A_j_iphjphkphL = alpha_Phi_minus_betaj_A_j_iphjphkph[index];
        // - partial_i -> - (A_{i} - A_{i-1})/dX = (A_{i-1} - A_{i})/dX, for Ax
        Ax_rhs[index] += dXm1*(alpha_Phi_minus_betaj_A_j_iphjphkph[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - alpha_Phi_minus_betaj_A_j_iphjphkphL);
        Ay_rhs[index] += dYm1*(alpha_Phi_minus_betaj_A_j_iphjphkph[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - alpha_Phi_minus_betaj_A_j_iphjphkphL);
        Az_rhs[index] += dZm1*(alpha_Phi_minus_betaj_A_j_iphjphkph[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - alpha_Phi_minus_betaj_A_j_iphjphkphL);

        // \partial_t psi6phi = [shift advection term] + \partial_j (\alpha \sqrt{\gamma} A^j)
        // Here we compute [shift advection term] = \partial_j (\beta^j psi6phi)
        // Cache misses are likely more expensive than branch mispredictions here,
        //       which is why we use if() statements and array lookups inside the if()'s.
        CCTK_REAL psi6phi_rhsL=0.0;
        const CCTK_REAL psi6phiL=psi6phi[index];
        const CCTK_REAL shiftx_iphjphkphL=shiftx_iphjphkph[index];
        const CCTK_REAL shifty_iphjphkphL=shifty_iphjphkph[index];
        const CCTK_REAL shiftz_iphjphkphL=shiftz_iphjphkph[index];

        // \partial_x (\beta^x psi6phi) :
        if(shiftx_iphjphkphL < 0.0) {
          psi6phi_rhsL+=0.5*dXm1*(+    shiftx_iphjphkph[CCTK_GFINDEX3D(cctkGH,i-2,j,k)]*psi6phi[CCTK_GFINDEX3D(cctkGH,i-2,j,k)]
                                  -4.0*shiftx_iphjphkph[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]*psi6phi[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]
                                  +3.0*shiftx_iphjphkphL*                               psi6phiL);
        } else {
          psi6phi_rhsL+=0.5*dXm1*(-    shiftx_iphjphkph[CCTK_GFINDEX3D(cctkGH,i+2,j,k)]*psi6phi[CCTK_GFINDEX3D(cctkGH,i+2,j,k)]
                                  +4.0*shiftx_iphjphkph[CCTK_GFINDEX3D(cctkGH,i+1,j,k)]*psi6phi[CCTK_GFINDEX3D(cctkGH,i+1,j,k)]
                                  -3.0*shiftx_iphjphkphL*                               psi6phiL);
        }

        // \partial_y (\beta^y psi6phi) :
        if(shifty_iphjphkphL < 0.0) {
          psi6phi_rhsL+=0.5*dYm1*(+    shifty_iphjphkph[CCTK_GFINDEX3D(cctkGH,i,j-2,k)]*psi6phi[CCTK_GFINDEX3D(cctkGH,i,j-2,k)]
                                  -4.0*shifty_iphjphkph[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]*psi6phi[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]
                                  +3.0*shifty_iphjphkphL*                               psi6phiL);
        } else {
          psi6phi_rhsL+=0.5*dYm1*(-    shifty_iphjphkph[CCTK_GFINDEX3D(cctkGH,i,j+2,k)]*psi6phi[CCTK_GFINDEX3D(cctkGH,i,j+2,k)]
                                  +4.0*shifty_iphjphkph[CCTK_GFINDEX3D(cctkGH,i,j+1,k)]*psi6phi[CCTK_GFINDEX3D(cctkGH,i,j+1,k)]
                                  -3.0*shifty_iphjphkphL*                               psi6phiL);
        }

        // \partial_z (\beta^z psi6phi) :
        if(shiftz_iphjphkphL < 0.0) {
          psi6phi_rhsL+=0.5*dZm1*(+    shiftz_iphjphkph[CCTK_GFINDEX3D(cctkGH,i,j,k-2)]*psi6phi[CCTK_GFINDEX3D(cctkGH,i,j,k-2)]
                                  -4.0*shiftz_iphjphkph[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]*psi6phi[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]
                                  +3.0*shiftz_iphjphkphL*                               psi6phiL);
        } else {
          psi6phi_rhsL+=0.5*dZm1*(-    shiftz_iphjphkph[CCTK_GFINDEX3D(cctkGH,i,j,k+2)]*psi6phi[CCTK_GFINDEX3D(cctkGH,i,j,k+2)]
                                  +4.0*shiftz_iphjphkph[CCTK_GFINDEX3D(cctkGH,i,j,k+1)]*psi6phi[CCTK_GFINDEX3D(cctkGH,i,j,k+1)]
                                  -3.0*shiftz_iphjphkphL*                               psi6phiL);
        }

        // Next we add \partial_j (\alpha \sqrt{\gamma} A^j) to \partial_t psi6phi:
        psi6phi_rhsL+=dXm1*(alpha_sqrtg_Ax_interp[index] - alpha_sqrtg_Ax_interp[CCTK_GFINDEX3D(cctkGH,i+1,j,k)])
          +           dYm1*(alpha_sqrtg_Ay_interp[index] - alpha_sqrtg_Ay_interp[CCTK_GFINDEX3D(cctkGH,i,j+1,k)])
          +           dZm1*(alpha_sqrtg_Az_interp[index] - alpha_sqrtg_Az_interp[CCTK_GFINDEX3D(cctkGH,i,j,k+1)]);

        // *GENERALIZED* LORENZ GAUGE:
        // Finally, add damping factor to \partial_t psi6phi
        //subtract lambda * alpha psi^6 Phi
        psi6phi_rhsL+=-damp_lorenz*alpha_iphjphkph[index]*psi6phiL;

        psi6phi_rhs[index] = psi6phi_rhsL;
      }
}

static inline CCTK_REAL avg(const CCTK_REAL f[PLUS2+1][PLUS2+1][PLUS2+1],const int imin,const int imax, const int jmin,const int jmax, const int kmin,const int kmax) {
  CCTK_REAL retval=0.0,num_in_sum=0.0;
  for(int kk=kmin;kk<=kmax;kk++) for(int jj=jmin;jj<=jmax;jj++) for(int ii=imin;ii<=imax;ii++) {
        retval+=f[kk][jj][ii]; num_in_sum++;
      }
  return retval/num_in_sum;
}
