/*
 * The following sets up a vector potential of the form
 * A_\phi = \varpi^2 A_b max[(EE-EE_cut),0],
 * where \varpi is the cylindrical radius: sqrt(x^2+y^2),
 * and EE \in {\rho, P} is the variable P or rho, specifying
 *   whether the vector potential is proportional to P or rho
 *   in the region greater than the cutoff.
 *
 * This is inspired by the choice of, e.g., the Illinois
 * NR group (see, e.g., 2nd full paragraph of pg 2
 * in https://arxiv.org/pdf/astro-ph/0510653.pdf)
 *
 * You may disable the \varpi^2 multiplication by simply
 * setting enable_varpi_squared_multiplication = false .
 *
 * This formulation assumes that A_r and A_\theta = 0;
 * only A_\phi can be nonzero. Thus the coordinate
 * transformations are as follows:
 *
 * A_x = dphi/dx A_phi
 *     = d[atan(y/x)]/dx A_phi
 *     = -y/(x^2+y^2) A_phi
 *     = -y/varpi^2 A_phi
 * A_y = dphi/dy A_phi
 *     = d[atan(y/x)]/dy A_phi
 *     =  x/(x^2+y^2) A_phi
 *     =  x/varpi^2 A_phi
 *
 * Note that both staggered (IllinoisGRMHD-compatible) and unstaggered
 * A-field options are supported.
 */
#include <cstdio>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )

extern "C"
void Seed_Magnetic_Fields (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO ("Seeding magnetic fields");

  CCTK_REAL dX = CCTK_DELTA_SPACE(0);
  CCTK_REAL dY = CCTK_DELTA_SPACE(1);



  if(CCTK_EQUALS(Afield_type,"Pressure_prescription")) {
    // A_phi proportional to (pressure - p_cut)^n_s
    if(enable_IllinoisGRMHD_staggered_A_fields) {
#pragma omp parallel for
      for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
            int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
            int indexip1=CCTK_GFINDEX3D(cctkGH,MIN(i+1,cctk_lsh[0]-1),j,k);
            int indexjp1=CCTK_GFINDEX3D(cctkGH,i,MIN(j+1,cctk_lsh[1]-1),k);
            int indexkp1=CCTK_GFINDEX3D(cctkGH,i,j,MIN(k+1,cctk_lsh[2]-1));
            int indexip1kp1=CCTK_GFINDEX3D(cctkGH,MIN(i+1,cctk_lsh[0]-1),j,MIN(k+1,cctk_lsh[2]-1));
            int indexjp1kp1=CCTK_GFINDEX3D(cctkGH,i,MIN(j+1,cctk_lsh[1]-1),MIN(k+1,cctk_lsh[2]-1));

            CCTK_REAL xL = x[index];
            CCTK_REAL yL = y[index];
            // zL unused

            CCTK_REAL PL = press[index]; // Assumes HydroBase pressure is set!

            CCTK_REAL PLip1 = press[indexip1];
            CCTK_REAL PLjp1 = press[indexjp1];
            CCTK_REAL PLkp1 = press[indexkp1];
            CCTK_REAL PLip1kp1 = press[indexip1kp1];
            CCTK_REAL PLjp1kp1 = press[indexjp1kp1];


            CCTK_REAL Pressure_at_Ax_stagger = 0.25*(PL + PLjp1 + PLkp1 + PLjp1kp1);
            Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] = -(yL + 0.5*dY)*A_b*pow(MAX(Pressure_at_Ax_stagger-P_cut,0.0),n_s);
            if(enable_varpi_squared_multiplication == false) {
              CCTK_REAL x2_plus_y2 = xL*xL + (yL + 0.5*dY)*(yL + 0.5*dY);
              Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] /= x2_plus_y2;
            }

            CCTK_REAL Pressure_at_Ay_stagger = 0.25*(PL + PLip1 + PLkp1 + PLip1kp1);
            Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] =  (xL + 0.5*dX)*A_b*pow(MAX(Pressure_at_Ay_stagger-P_cut,0.0),n_s);
            if(enable_varpi_squared_multiplication == false) {
              CCTK_REAL x2_plus_y2 = (xL + 0.5*dX)*(xL + 0.5*dX) + yL*yL;
              Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] /= x2_plus_y2;
            }
            Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)] =  0.0;
            Aphi[index]  = 0.0;

          }
    } else { /* NOTE: We assume vertex-centered grids if enable_IllinoisGRMHD_staggered_A_fields==false */
#pragma omp parallel for
      for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
            int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
            CCTK_REAL PL = press[index]; // Assumes HydroBase pressure is set!
            Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] = -y[index]*A_b*pow(MAX(PL-P_cut,0.0),n_s);
            Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] =  x[index]*A_b*pow(MAX(PL-P_cut,0.0),n_s);
            if(enable_varpi_squared_multiplication == false) {
              CCTK_REAL xL = x[index];
              CCTK_REAL yL = y[index];
              CCTK_REAL x2_plus_y2 = xL*xL + yL*yL;
              Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] /= x2_plus_y2;
              Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] /= x2_plus_y2;
            }
            Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)] =  0.0;
            Aphi[index]  = 0.0;
          }
    }
  } else if(CCTK_EQUALS(Afield_type,"Density_prescription")) {
    // A_phi proportional to (rho - rho_cut)
    if(enable_IllinoisGRMHD_staggered_A_fields) {
#pragma omp parallel for
      for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
            int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
            int indexip1=CCTK_GFINDEX3D(cctkGH,MIN(i+1,cctk_lsh[0]-1),j,k);
            int indexjp1=CCTK_GFINDEX3D(cctkGH,i,MIN(j+1,cctk_lsh[1]-1),k);
            int indexkp1=CCTK_GFINDEX3D(cctkGH,i,j,MIN(k+1,cctk_lsh[2]-1));
            int indexip1kp1=CCTK_GFINDEX3D(cctkGH,MIN(i+1,cctk_lsh[0]-1),j,MIN(k+1,cctk_lsh[2]-1));
            int indexjp1kp1=CCTK_GFINDEX3D(cctkGH,i,MIN(j+1,cctk_lsh[1]-1),MIN(k+1,cctk_lsh[2]-1));

            CCTK_REAL xL = x[index];
            CCTK_REAL yL = y[index];
            // zL unused

            CCTK_REAL rhoL = rho[index]; // Assumes HydroBase rho is set!

            CCTK_REAL rhoLip1 = rho[indexip1];
            CCTK_REAL rhoLjp1 = rho[indexjp1];
            CCTK_REAL rhoLkp1 = rho[indexkp1];
            CCTK_REAL rhoLip1kp1 = rho[indexip1kp1];
            CCTK_REAL rhoLjp1kp1 = rho[indexjp1kp1];

            CCTK_REAL rho_at_Ax_stagger = 0.25*(rhoL + rhoLjp1 + rhoLkp1 + rhoLjp1kp1);
            Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] = -(yL + 0.5*dY)*A_b*MAX(rho_at_Ax_stagger-rho_cut, 0.0);
            if(enable_varpi_squared_multiplication == false) {
              CCTK_REAL x2_plus_y2 = xL*xL + (yL + 0.5*dY)*(yL + 0.5*dY);
              Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] /= x2_plus_y2;
            }

            CCTK_REAL rho_at_Ay_stagger = 0.25*(rhoL + rhoLip1 + rhoLkp1 + rhoLip1kp1);
            Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] =  (xL + 0.5*dX)*A_b*MAX(rho_at_Ay_stagger-rho_cut, 0.0);
            if(enable_varpi_squared_multiplication == false) {
              CCTK_REAL x2_plus_y2 = (xL + 0.5*dX)*(xL + 0.5*dX) + yL*yL;
              Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] /= x2_plus_y2;
            }

            Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)] =  0.0;
            Aphi[index]  = 0.0;
          }
    } else { /* NOTE: We assume vertex-centered grids if enable_IllinoisGRMHD_staggered_A_fields==false */
#pragma omp parallel for
      for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
            int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
            CCTK_REAL rhoL = rho[index]; // Assumes HydroBase rho is set!
            Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] = -y[index]*A_b*MAX(rhoL-rho_cut, 0.0);
            Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] =  x[index]*A_b*MAX(rhoL-rho_cut, 0.0);
            if(enable_varpi_squared_multiplication == false) {
              CCTK_REAL xL = x[index];
              CCTK_REAL yL = y[index];
              CCTK_REAL x2_plus_y2 = xL*xL + yL*yL;
              Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] /= x2_plus_y2;
              Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] /= x2_plus_y2;
            }

            Avec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)] =  0.0;
            Aphi[index]  = 0.0;
          }
    }
  }
}
