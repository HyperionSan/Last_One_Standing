#include <cstdio>
#include <cassert>
#include <vector>
#include <ios>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <algorithm>

using std::min;
using std::max;

extern "C"
void Seed_Magnetic_Fields_Privt (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO ("Seeding magnetic fields");

  double dX = CCTK_DELTA_SPACE(0);
  double dY = CCTK_DELTA_SPACE(1);
  double Ax_yshift_staggering = 0.;
  double Ay_xshift_staggering = 0.;

  if(enable_IllinoisGRMHD_staggered_A_fields) {
#pragma omp parallel for
    for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
          int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
          int indexip1=CCTK_GFINDEX3D(cctkGH,min(i+1,cctk_lsh[0]-1),j,k);
          int indexjp1=CCTK_GFINDEX3D(cctkGH,i,min(j+1,cctk_lsh[1]-1),k);
          int indexkp1=CCTK_GFINDEX3D(cctkGH,i,j,min(k+1,cctk_lsh[2]-1));
          int indexip1kp1=CCTK_GFINDEX3D(cctkGH,min(i+1,cctk_lsh[0]-1),j,min(k+1,cctk_lsh[2]-1));
          int indexjp1kp1=CCTK_GFINDEX3D(cctkGH,i,min(j+1,cctk_lsh[1]-1),min(k+1,cctk_lsh[2]-1));

          double xL = x[index];
          double yL = y[index];
          double zL = z[index];
          double x1 = xL - x_c1;
          double x2 = xL - x_c2;
          
          double r1 = sqrt(x1*x1 + yL*yL + zL*zL);
          double r2 = sqrt(x2*x2 + yL*yL + zL*zL);

	  if(CCTK_EQUALS(A_field_type, "poloidal_A_interior"))
	    {
	      double PL = press[index]; // Assumes HydroBase pressure is set!

	      double PLip1 = press[indexip1];
	      double PLjp1 = press[indexjp1];
	      double PLkp1 = press[indexkp1];
	      double PLip1kp1 = press[indexip1kp1];
	      double PLjp1kp1 = press[indexjp1kp1];
          
	      double Pressure_at_Ax_stagger = 0.25*(PL + PLjp1 + PLkp1 + PLjp1kp1);
	      Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = -(y[index] + 0.5*dY)*A_b*pow(max(Pressure_at_Ax_stagger-P_cut,0.0),n_s);

	      double Pressure_at_Ay_stagger = 0.25*(PL + PLip1 + PLkp1 + PLip1kp1);

	      if(!have_two_NSs_along_x_axis) Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] =  (xL + 0.5*dX)*A_b*pow(max(Pressure_at_Ay_stagger-P_cut,0.0),n_s);
	      else {
		if(r1<=r_NS1) Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] =  (x1 + 0.5*dX)*A_b*pow(max(Pressure_at_Ay_stagger-P_cut,0.0),n_s);
		if(r2<=r_NS2) Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] =  (x2 + 0.5*dX)*A_b*pow(max(Pressure_at_Ay_stagger-P_cut,0.0),n_s);
		if(r1>r_NS1 && r2>r_NS2) Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = 0.0; // No external B-field.
	      }
	    }
	  // Modified to include poloidal vector potential option
	  else if(CCTK_EQUALS(A_field_type, "dipolar_A_everywhere"))
	    {
	      if(!enable_IllinoisGRMHD_staggered_A_fields)
		{
		  CCTK_WARN(1, "dipolar_A_everywhere requires a staggered grid");
		}

	      double pi = 3.141592653589793;

	      // \varpi^2
	      double varpi1_2 = x1 * x1 + yL * yL;
	      double varpi2_2 = x2 * x2 + yL * yL;

	      // Current loop radius
	      double r01 = r_zero_NS1 * r_NS1;
	      double r02 = r_zero_NS2 * r_NS2;

	      // phi-component of vector potential in spherical basis
	      // Eq (2) in Paschalidis et al PRD 88 021504(R) (2013)
	      double Ap1 = pi * r01 * r01 * I_zero_NS1 * pow(r01 * r01 + r1 * r1, -1.5) * (1.0 + 1.875 * r01 * r01 * (r01 * r01 + varpi1_2) * pow(r01 * r01 + r1 * r1, -2.0));
	      double Ap2 = pi * r02 * r02 * I_zero_NS2 * pow(r02 * r02 + r2 * r2, -1.5) * (1.0 + 1.875 * r02 * r02 * (r02 * r02 + varpi2_2) * pow(r02 * r02 + r2 * r2, -2.0));

	      // x-component of vector potential in Cartesian basis
	      Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = - (yL + 0.5 * dY) * (Ap1 + Ap2);

	      // y-component of vector potential in Cartesian basis
	      Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = (x1 + 0.5 * dX) * Ap1 + (x2 + 0.5 * dX) * Ap2;
	    }

	  // z-component of vector potential in Cartesian basis
          Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] =  0.0;
          Aphi[index]  = 0.0;
          
        }
  } else {
#pragma omp parallel for
    for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
          int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
          double PL = press[index]; // Assumes HydroBase pressure is set!
          Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = -y[index]*A_b*pow(max(PL-P_cut,0.0),n_s);
          Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] =  x[index]*A_b*pow(max(PL-P_cut,0.0),n_s);
          Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] =  0.0;
          Aphi[index]  = 0.0;
        }
  }
}

