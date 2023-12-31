/*****************************************
 * PPM Reconstruction Interface.
 * Zachariah B. Etienne (2013)
 *
 * This version of PPM implements the standard
 * Colella & Woodward PPM, but in the GRFFE
 * limit, where P=rho=0. Thus, e.g., ftilde=0.
 *****************************************/

#define MINUS2 0
#define MINUS1 1
#define PLUS0  2
#define PLUS1  3
#define PLUS2  4
#define MAXNUMINDICES 5
//    ^^^^^^^^^^^^^ Be _sure_ to define MAXNUMINDICES appropriately!

// You'll find the #define's for LOOP_DEFINE and SET_INDEX_ARRAYS inside:
#include "loop_defines_reconstruction.h"

static inline CCTK_REAL slope_limit(const CCTK_REAL dU,const CCTK_REAL dUp1);
static inline void monotonize(const CCTK_REAL U,CCTK_REAL &Ur,CCTK_REAL &Ul);

static void reconstruct_set_of_prims_PPM_GRFFE(const cGH *cctkGH,const int *cctk_lsh,const int flux_dirn,const int num_prims_to_reconstruct,const int *which_prims_to_reconstruct,
                                               const gf_and_gz_struct *in_prims,gf_and_gz_struct *out_prims_r,gf_and_gz_struct *out_prims_l, CCTK_REAL *temporary) {

  CCTK_REAL U[MAXNUMVARS][MAXNUMINDICES],dU[MAXNUMVARS][MAXNUMINDICES],slope_lim_dU[MAXNUMVARS][MAXNUMINDICES],
    Ur[MAXNUMVARS][MAXNUMINDICES],Ul[MAXNUMVARS][MAXNUMINDICES];
  int ijkgz_lo_hi[4][2];

  for(int ww=0;ww<num_prims_to_reconstruct;ww++) {
    const int whichvar=which_prims_to_reconstruct[ww];

    if(in_prims[whichvar].gz_lo[flux_dirn]!=0 || in_prims[whichvar].gz_hi[flux_dirn]!=0) {
      CCTK_VError(VERR_DEF_PARAMS,"TOO MANY GZ'S! WHICHVAR=%d: %d %d %d : %d %d %d DIRECTION %d",whichvar,
		  in_prims[whichvar].gz_lo[1],in_prims[whichvar].gz_lo[2],in_prims[whichvar].gz_lo[3],
		  in_prims[whichvar].gz_hi[1],in_prims[whichvar].gz_hi[2],in_prims[whichvar].gz_hi[3],flux_dirn);
    }

    // *** LOOP 1: Interpolate to Ur and Ul, which are face values ***
    //  You will find that Ur depends on U at MINUS1,PLUS0, PLUS1,PLUS2, and
    //                     Ul depends on U at MINUS2,MINUS1,PLUS0,PLUS1.
    //  However, we define the below loop from MINUS2 to PLUS2. Why not split
    //     this up and get additional points? Maybe we should. In GRMHD, the
    //     reason is that later on, Ur and Ul depend on ftilde, which is
    //     defined from MINUS2 to PLUS2, so we would lose those points anyway.
    //     But in GRFFE, ftilde is set to zero, so there may be a potential
    //     for boosting performance here.
    LOOP_DEFINE(2,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
      SET_INDEX_ARRAYS(-2,2,flux_dirn);
      /* *** LOOP 1a: READ INPUT *** */
      // Read in a primitive at all gridpoints between m = MINUS2 & PLUS2, where m's direction is given by flux_dirn. Store to U.
      for(int ii=MINUS2;ii<=PLUS2;ii++) U[whichvar][ii] = in_prims[whichvar].gf[index_arr[flux_dirn][ii]];

      /* *** LOOP 1b: DO COMPUTATION *** */
      /* First, compute simple dU = U(i) - U(i-1), where direction of i
       *         is given by flux_dirn, and U is a primitive variable:
       *         {vx,vy,vz,Bx,By,Bz}. */
      // Note that for Ur and Ul at i, we must compute dU(i-1),dU(i),dU(i+1),
      //         and dU(i+2)
      dU[whichvar][MINUS1] = U[whichvar][MINUS1]- U[whichvar][MINUS2];
      dU[whichvar][PLUS0]  = U[whichvar][PLUS0] - U[whichvar][MINUS1];
      dU[whichvar][PLUS1]  = U[whichvar][PLUS1] - U[whichvar][PLUS0];
      dU[whichvar][PLUS2]  = U[whichvar][PLUS2] - U[whichvar][PLUS1];

      // Then, compute slope-limited dU, using MC slope limiter:
      slope_lim_dU[whichvar][MINUS1]=slope_limit(dU[whichvar][MINUS1],dU[whichvar][PLUS0]);
      slope_lim_dU[whichvar][PLUS0] =slope_limit(dU[whichvar][PLUS0], dU[whichvar][PLUS1]);
      slope_lim_dU[whichvar][PLUS1] =slope_limit(dU[whichvar][PLUS1], dU[whichvar][PLUS2]);

      // Finally, compute face values Ur and Ul based on the PPM prescription
      //   (Eq. A1 in http://arxiv.org/pdf/astro-ph/0503420.pdf, but using standard 1/6=(1.0/6.0) coefficient)
      // Ur[PLUS0] represents U(i+1/2)
      // We applied a simplification to the following line: Ur=U+0.5*(U(i+1)-U) + ... = 0.5*(U(i+1)+U) + ...
      Ur[whichvar][PLUS0] = 0.5*(U[whichvar][PLUS1] + U[whichvar][PLUS0] ) + (1.0/6.0)*(slope_lim_dU[whichvar][PLUS0]  - slope_lim_dU[whichvar][PLUS1]);
      // Ul[PLUS0] represents U(i-1/2)
      // We applied a simplification to the following line: Ul=U(i-1)+0.5*(U-U(i-1)) + ... = 0.5*(U+U(i-1)) + ...
      Ul[whichvar][PLUS0] = 0.5*(U[whichvar][PLUS0] + U[whichvar][MINUS1]) + (1.0/6.0)*(slope_lim_dU[whichvar][MINUS1] - slope_lim_dU[whichvar][PLUS0]);

      /* *** LOOP 1c: WRITE OUTPUT *** */
      // Store right face values to {vxr,vyr,vzr,Bxr,Byr,Bzr},
      //    and left face values to {vxl,vyl,vzl,Bxl,Byl,Bzl}
      out_prims_r[whichvar].gf[index_arr[flux_dirn][PLUS0]] = Ur[whichvar][PLUS0];
      out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]] = Ul[whichvar][PLUS0];
    }

    // *** LOOP 2 (REMOVED): STEEPEN RHOB. RHOB DOES NOT EXIST IN GRFFE EQUATIONS ***
  }

  // *** LOOP 3: FLATTEN BASED ON FTILDE AND MONOTONIZE ***
  for(int ww=0;ww<num_prims_to_reconstruct;ww++) {
    const int whichvar=which_prims_to_reconstruct[ww];
    // ftilde() depends on P(MINUS2,MINUS1,PLUS1,PLUS2), THUS IS SET TO ZERO IN GRFFE
    LOOP_DEFINE(2,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
      SET_INDEX_ARRAYS(0,0,flux_dirn);

      U[whichvar][PLUS0]  = in_prims[whichvar].gf[index_arr[flux_dirn][PLUS0]];
      Ur[whichvar][PLUS0] = out_prims_r[whichvar].gf[index_arr[flux_dirn][PLUS0]];
      Ul[whichvar][PLUS0] = out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]];

      // ftilde_gf was computed in the function compute_ftilde_gf(), called before this routine
      //CCTK_REAL ftilde = ftilde_gf[index_arr[flux_dirn][PLUS0]];
      // ...and then flatten (local operation)
      Ur[whichvar][PLUS0]   = Ur[whichvar][PLUS0];
      Ul[whichvar][PLUS0]   = Ul[whichvar][PLUS0];

      // Then monotonize
      monotonize(U[whichvar][PLUS0],Ur[whichvar][PLUS0],Ul[whichvar][PLUS0]);

      out_prims_r[whichvar].gf[index_arr[flux_dirn][PLUS0]] = Ur[whichvar][PLUS0];
      out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]] = Ul[whichvar][PLUS0];
    }
    // Note: ftilde=0 in GRFFE. Ur depends on ftilde, which depends on points of U between MINUS2 and PLUS2
    out_prims_r[whichvar].gz_lo[flux_dirn]+=2;
    out_prims_r[whichvar].gz_hi[flux_dirn]+=2;
    // Note: ftilde=0 in GRFFE. Ul depends on ftilde, which depends on points of U between MINUS2 and PLUS2
    out_prims_l[whichvar].gz_lo[flux_dirn]+=2;
    out_prims_l[whichvar].gz_hi[flux_dirn]+=2;
  }

  // *** LOOP 4: SHIFT Ur AND Ul ***
  /* Currently face values are set so that
   *      a) Ur(i) represents U(i+1/2), and
   *      b) Ul(i) represents U(i-1/2)
   *    Here, we shift so that the indices are consistent:
   *      a) U(i-1/2+epsilon) = oldUl(i)   = newUr(i)
   *      b) U(i-1/2-epsilon) = oldUr(i-1) = newUl(i)
   *    Note that this step is not strictly necessary if you keep
   *      track of indices when computing the flux. */
  for(int ww=0;ww<num_prims_to_reconstruct;ww++) {
    const int whichvar=which_prims_to_reconstruct[ww];
    LOOP_DEFINE(3,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
      SET_INDEX_ARRAYS(-1,0,flux_dirn);
      temporary[index_arr[flux_dirn][PLUS0]] = out_prims_r[whichvar].gf[index_arr[flux_dirn][MINUS1]];
    }

    LOOP_DEFINE(3,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
      SET_INDEX_ARRAYS(0,0,flux_dirn);
      // Then shift so that Ur represents the gridpoint at i-1/2+epsilon,
      //                and Ul represents the gridpoint at i-1/2-epsilon.
      // Ur(i-1/2) = Ul(i-1/2)     = U(i-1/2+epsilon)
      // Ul(i-1/2) = Ur(i+1/2 - 1) = U(i-1/2-epsilon)
      out_prims_r[whichvar].gf[index_arr[flux_dirn][PLUS0]] = out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]];
      out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]] = temporary[index_arr[flux_dirn][PLUS0]];
    }
    // Ul was just shifted, so we lost another ghostzone.
    out_prims_l[whichvar].gz_lo[flux_dirn]+=1;
    out_prims_l[whichvar].gz_hi[flux_dirn]+=0;
    // As for Ur, we didn't need to get rid of another ghostzone,
    //    but we did ... seems wasteful!
    out_prims_r[whichvar].gz_lo[flux_dirn]+=1;
    out_prims_r[whichvar].gz_hi[flux_dirn]+=0;

  }
}

// Set SLOPE_LIMITER_COEFF = 2.0 for MC, 1 for minmod
#define SLOPE_LIMITER_COEFF 2.0

//Eq. 60 in JOURNAL OF COMPUTATIONAL PHYSICS 123, 1-14 (1996)
//   [note the factor of 2 missing in the |a_{j+1} - a_{j}| term].
//   Recall that dU = U_{i} - U_{i-1}.
static inline CCTK_REAL slope_limit(const CCTK_REAL dU,const CCTK_REAL dUp1) {
  if(dU*dUp1 > 0.0) {
    //delta_m_U=0.5 * [ (u_(i+1)-u_i) + (u_i-u_(i-1)) ] = (u_(i+1) - u_(i-1))/2  <-- first derivative, second-order; this should happen most of the time (smooth flows)
    const CCTK_REAL delta_m_U = 0.5*(dU + dUp1);
    // EXPLANATION OF BELOW LINE OF CODE.
    // In short, sign_delta_a_j = sign(delta_m_U) = (0.0 < delta_m_U) - (delta_m_U < 0.0).
    //    If delta_m_U>0, then (0.0 < delta_m_U)==1, and (delta_m_U < 0.0)==0, so sign_delta_a_j=+1
    //    If delta_m_U<0, then (0.0 < delta_m_U)==0, and (delta_m_U < 0.0)==1, so sign_delta_a_j=-1
    //    If delta_m_U==0,then (0.0 < delta_m_U)==0, and (delta_m_U < 0.0)==0, so sign_delta_a_j=0
    const int sign_delta_m_U = (0.0 < delta_m_U) - (delta_m_U < 0.0);
    //Decide whether to use 2nd order derivative or first-order derivative, limiting slope.
    return sign_delta_m_U*MIN(fabs(delta_m_U),MIN(SLOPE_LIMITER_COEFF*fabs(dUp1),SLOPE_LIMITER_COEFF*fabs(dU)));
  }
  return 0.0;
}

static inline void monotonize(const CCTK_REAL U,CCTK_REAL &Ur,CCTK_REAL &Ul) {
  const CCTK_REAL dU = Ur - Ul;
  const CCTK_REAL mU = 0.5*(Ur+Ul);

  if ( (Ur-U)*(U-Ul) <= 0.0) {
    Ur = U;
    Ul = U;
    return;
  }
  if ( dU*(U-mU) > (1.0/6.0)*SQR(dU)) {
    Ul = 3.0*U - 2.0*Ur;
    return;
  }
  if ( dU*(U-mU) < -(1.0/6.0)*SQR(dU)) {
    Ur = 3.0*U - 2.0*Ul;
    return;
  }
}
