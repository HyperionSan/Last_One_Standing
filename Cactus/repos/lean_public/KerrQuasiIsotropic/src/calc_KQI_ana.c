/* intial data thorn: KQI_analytic */
/*======================================================*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "KQI_ana_utils.h"

/* -------------------------------------------------------------------*/
void KQI_analytic(CCTK_ARGUMENTS);
void
KQI_analytic (CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

    //Structure of this function:
    //---------------------------------------------------------
    // 1) Define parameters
    // 2) Define coordinates
    // 3) Define metric functions 
    // 4) Initialize gauge functions
    // 4) Define metric, curvature etc in Cartesian coordinates (x,y,z) 


//  CCTK_INFO("=== Begin KQI initial data ===");

  /*=== define BH parameters ===*/
  /*----------------------------*/
  CCTK_REAL mass, mass2, spin, spin2, rBLp, rBLm;
  mass  = m_plus;
  mass2 = mass * mass;
  spin  = spin_plus;
  spin2 = spin  * spin;
  rBLp  = mass + sqrt( mass2 - spin2 );
  rBLm  = mass - sqrt( mass2 - spin2 );

  /*----------------------------*/

  /*=== define grid length ===*/
  /*--------------------------*/
  CCTK_INT imin[3], imax[3];
  for (int d = 0; d < 3; ++ d)
  {
    imin[d] = 0;
    imax[d] = cctk_lsh[d];
  }

  /*--------------------------*/

/*=== loops over full grid ===*/
/*----------------------------*/
//#pragma omp parallel for
  for (int k = imin[2]; k < imax[2]; ++k)
  {
   for (int j = imin[1]; j < imax[1]; ++j)
   {
    for (int i = imin[0]; i < imax[0]; ++i)
    {

     const int ind = CCTK_GFINDEX3D (cctkGH, i, j, k);


    /*=== define position parameters ===*/
    /*----------------------------------*/
    // Define coordinates
    CCTK_REAL xx, yy, zz;
    xx = x[ind] - pos_plus[0];
    yy = y[ind] - pos_plus[1];
    zz = z[ind] - pos_plus[2];

    CCTK_REAL RR, RR2;
    RR2 = xx * xx + yy * yy + zz * zz;
    if( RR2 < pow( eps_r, 2 ) ) 
        RR2 = pow( eps_r, 2 );
    RR  = sqrt( RR2 );

    CCTK_REAL rho, rho2, rho3;
    rho2 = xx * xx + yy * yy; 
    if( rho2 < pow( eps_r, 2 ) )
        rho2 = pow( eps_r, 2 );
    rho  = sqrt( rho2 );
    rho3 = rho2 * rho;

    CCTK_REAL ctheta, ctheta2;
    ctheta  = zz / RR;
    ctheta2 = ctheta * ctheta;

    CCTK_REAL stheta, stheta2;
    stheta  = rho / RR;
    stheta2 = stheta  * stheta;
    /*----------------------------------*/

    /*=== define metric functions ======*/ 
    /*----------------------------------*/
    CCTK_REAL rBL, rBL2;
    rBL  = RR * ( 1.0 + 0.25 * rBLp / RR ) * ( 1.0 + 0.25 * rBLp / RR );
    rBL2 = rBL * rBL;

    CCTK_REAL drBLdR;
    drBLdR = 1.0 - rBLp*rBLp / ( 16.0 * RR2 );

    CCTK_REAL Delt, Sigm, Sigm2, fctFF;
    Delt  = rBL2 + spin2 - 2 * mass * rBL;
    Sigm  = rBL2 + spin2 * ctheta2;
    Sigm2 = Sigm * Sigm;
    fctFF = ( rBL2 + spin2 ) * ( rBL2 + spin2 ) - Delt * spin2 * stheta2;

    CCTK_REAL psi04;
    psi04 = Sigm / RR2;

    CCTK_REAL fctGG, fctHH;
    fctGG = rBLm / ( RR2 * ( rBL - rBLm ) );
    fctHH = ( 2.0 * mass * rBL + Sigm ) / ( RR2 * Sigm2 );

    CCTK_REAL detgij;
    detgij = pow( psi04, 3) * ( 1.0 + RR2 * fctGG ) * ( 1.0 + spin2 * rho2 * fctHH );

    /*----------------------------------*/

    /*=== initialize gauge functions ===*/
    /*----------------------------------*/
    if( CCTK_Equals(initial_lapse, "KQI_ana" ))
    // analytic continuation
    {
      alp[ind] = ( 4.0*RR - rBLp) * sqrt( rBL -rBLm )  / sqrt( 16.0*RR * ( rBL2 + spin2 * ( 1.0 + 2.0*mass*rBL * stheta2 / Sigm ) ) );
    }
    else if( CCTK_Equals(initial_lapse, "KQI_abs" ))
    // absolute value of lapse
    // Discouraged. Better to use the analytic continuation of KQI_ana.
    {
      alp[ind] = sqrt( Delt * Sigm / fctFF );
    }
    else if( CCTK_Equals(initial_lapse, "KQI_prec" ))
    // precollapsed lapse
    {
      alp[ind] = pow( detgij, -1./6 ); 
    }
    else
    {
      alp[ind] = 1.0;
    }


    if (CCTK_Equals(initial_shift, "KQI_ana"))
    {
      CCTK_REAL bphi;
      bphi = 2.0 * spin * mass * rBL / fctFF;

      betax[ind] =   yy * bphi;
      betay[ind] = - xx * bphi;
      betaz[ind] = 0.0;

    }
    else
    {
      betax[ind] = 0.0;
      betay[ind] = 0.0;
      betaz[ind] = 0.0;
    }

    /*----------------------------------*/

    /*=== conformal metric in Cartesian coords ===*/
    /*--------------------------------------------*/

    gxx[ind] = psi04 * ( 1.0 + xx*xx * fctGG + spin2 * yy*yy * fctHH );
    gxy[ind] = psi04 * (       xx*yy * fctGG - spin2 * xx*yy * fctHH );
    gxz[ind] = psi04 * (       xx*zz * fctGG );
    gyy[ind] = psi04 * ( 1.0 + yy*yy * fctGG + spin2 * xx*xx * fctHH );
    gyz[ind] = psi04 * (       yy*zz * fctGG );
    gzz[ind] = psi04 * ( 1.0 + zz*zz * fctGG );

    /*--------------------------------------------*/

    /*=== Aij in Cartesian coords ===*/
    /*-------------------------------*/
//HW: TODO: take off
    // analytic continuation
////    alpKQI = ( 4.0*RR - rBLp) * sqrt( rBL -rBLm )  / sqrt( 16.0*RR * ( rBL2 + spin2 * ( 1.0 + 2.0*mass*rBL * stheta2 / Sigm ) ) );
//    CCTK_REAL alpKQI = alp[ind];

    CCTK_REAL auxKij, facKij, facKijRho, facKijZ;
    auxKij    = 2.0 * rBL2 * ( rBL2 + spin2 ) + Sigm * ( rBL2 - spin2 );
    facKij    = alp[ind] * spin * mass * stheta2 / ( RR2 * rho3 * Delt * Sigm2 );
    facKijRho = 2.0 * zz  * spin2 * rBL * Delt * ctheta * stheta - rho * RR * drBLdR * auxKij;
    facKijZ   = 2.0 * rho * spin2 * rBL * Delt * ctheta * stheta + zz  * RR * drBLdR * auxKij;

    kxx[ind] =   2.0 * xx * yy   * facKij * facKijRho;
    kxy[ind] = ( yy*yy - xx*xx ) * facKij * facKijRho;
    kxz[ind] = - yy * rho        * facKij * facKijZ;
    kyy[ind] = - 2.0 * xx * yy   * facKij * facKijRho;
    kyz[ind] =   xx * rho        * facKij * facKijZ;
    kzz[ind] =   0.0;

    /*---------------------------------------*/

    } /* for i */
   }  /* for j */
  }   /* for k */
/*=== end of loops over grid ===*/
/*------------------------------*/

//  CCTK_INFO("=== End KQI initial data ===");

}
/* -------------------------------------------------------------------*/
