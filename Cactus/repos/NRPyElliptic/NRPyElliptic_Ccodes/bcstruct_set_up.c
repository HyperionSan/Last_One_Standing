#include "././NRPy_basic_defines.h"
#include "././NRPy_function_prototypes.h"
/*
 * EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt():
 *   A coordinate system's "eigencoordinate" is the simplest member
 *   of its family; all spherical-like coordinate systems have
 *   Spherical as their eigencoordinate. The same is true for
 *   cylindrical-like (Cylindrical is eigencoordinate),
 *   Cartesian-like (Cartesian is the eigencoordinate), and
 *   SymTP-like (SymTP is the eigencoordinate) coordinates.
 * 
 *   For a given gridpoint (i0,i1,i2) and corresponding coordinate
 *   (x0,x1,x2), this function performs the dual mapping
 *   (x0,x1,x2) -> (Cartx,Carty,Cartz) -> (x0,x1,x2)'
 *   Note that (x0,x1,x2) IS NOT ALWAYS equal to (x0,x1,x2)';
 *   For example consider in Spherical coordinates
 *   (x0,x1,x2)=(r,theta,phi)=(-0.1,pi/4,pi/4).
 *   This point will map to (x0,x1,x2)', in which x0>0,
 *   because the inversion r=sqrt(Cartx^2+Carty^2+Cartz^2)
 *   is always positive. In this case, (x0,x1,x2) is considered
 *   an *inner* boundary point, and on a cell-centered grid
 *   is guaranteed to map to a grid point in the grid interior;
 *   filling in this point requires copying data, and possibly
 *   multiplying by a +/- 1 if the data is from a gridfunction
 *   storing tensors/vectors.
 */
void EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(const paramstruct *restrict params, REAL *restrict xx[3],
                                                               const int i0, const int i1, const int i2,
                                                               REAL x0x1x2_inbounds[3], int i0i1i2_inbounds[3]) {
#include "./set_Cparameters.h"


  // This is a 3-step algorithm:
  // Step 1: (x0,x1,x2) -> (Cartx,Carty,Cartz)
  //         Find the Cartesian coordinate that (x0,x1,x2)
  //         maps to, assuming (x0,x1,x2) is the eigen-
  //         coordinate. Note that we assume (x0,x1,x2)
  //         has the same grid boundaries in both the
  //         original coordinate and the eigencoordinate.
  // Step 2: (Cartx,Carty,Cartz) -> (x0,x1,x2)'
  //         Find the interior eigencoordinate point
  //         (x0,x1,x2)' to which (Cartx,Carty,Cartz)
  //         maps, as well as the corresponding
  //         gridpoint integer index (i0,i1,i2). For
  //         cell-centered grids, (x0,x1,x2) will always
  //         overlap exactly (to roundoff error) a point
  //         on the numerical grid.
  // Step 3: Sanity check
  //         Convert x0(i0_inbounds),x1(i1_inbounds),x2(i2_inbounds) -> (Cartx,Carty,Cartz),
  //         and check that
  //         (Cartx,Carty,Cartz) == (Cartx(x0(i0)),Cartx(x1(i1)),Cartx(x2(i2)))
  //         If not, error out!

  // Step 1: Convert the (curvilinear) coordinate (x0,x1,x2) to Cartesian coordinates
  REAL xCart[3];  // where (x,y,z) is output
  {
    // xx_to_Cart for EigenCoordinate SymTP (orig coord = SinhSymTP):
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
    /*
     *  Original SymPy expressions:
     *  "[xCart[0] = xx0*sin(xx1)*cos(xx2),
     *    xCart[1] = xx0*sin(xx1)*sin(xx2),
     *    xCart[2] = sqrt(bScale**2 + xx0**2)*cos(xx1)]"
     */
    {
      const double tmp_0 = xx0*sin(xx1);
      xCart[0] = tmp_0*cos(xx2);
      xCart[1] = tmp_0*sin(xx2);
      xCart[2] = sqrt(((bScale)*(bScale)) + ((xx0)*(xx0)))*cos(xx1);
    }
  }

  REAL Cartx = xCart[0];
  REAL Carty = xCart[1];
  REAL Cartz = xCart[2];
  // Step 2: Find the (i0_inbounds,i1_inbounds,i2_inbounds) corresponding to the above Cartesian coordinate.
  //   If (i0_inbounds,i1_inbounds,i2_inbounds) is in a ghost zone, then it must equal (i0,i1,i2), and
  //      the point is an outer boundary point.
  //   Otherwise (i0_inbounds,i1_inbounds,i2_inbounds) is in the grid interior, and data at (i0,i1,i2)
  //      must be replaced with data at (i0_inbounds,i1_inbounds,i2_inbounds), but multiplied by the
  //      appropriate parity condition (+/- 1).
  REAL Cart_to_xx0_inbounds,Cart_to_xx1_inbounds,Cart_to_xx2_inbounds;
  // Cart_to_xx for EigenCoordinate SymTP (orig coord = SinhSymTP);
  /*
   *  Original SymPy expressions:
   *  "[Cart_to_xx0_inbounds = M_SQRT1_2*sqrt(Cartx**2 + Carty**2 + Cartz**2 - bScale**2 + sqrt(-4*Cartz**2*bScale**2 + bScale**4 + 2*bScale**2*(Cartx**2 + Carty**2 + Cartz**2) + (Cartx**2 + Carty**2 + Cartz**2)**2)),
   *    Cart_to_xx1_inbounds = acos(M_SQRT1_2*sqrt(1 + (Cartx**2 + Carty**2 + Cartz**2)/bScale**2 - sqrt(-4*Cartz**2*bScale**2 + bScale**4 + 2*bScale**2*(Cartx**2 + Carty**2 + Cartz**2) + (Cartx**2 + Carty**2 + Cartz**2)**2)/bScale**2)*sign(Cartz)),
   *    Cart_to_xx2_inbounds = atan2(Carty, Cartx)]"
   */
  {
    const double tmp_1 = ((bScale)*(bScale));
    const double tmp_2 = ((Cartx)*(Cartx)) + ((Carty)*(Carty)) + ((Cartz)*(Cartz));
    const double tmp_3 = sqrt(-4*((Cartz)*(Cartz))*tmp_1 + ((bScale)*(bScale)*(bScale)*(bScale)) + 2*tmp_1*tmp_2 + ((tmp_2)*(tmp_2)));
    const double tmp_4 = (1.0/(tmp_1));
    Cart_to_xx0_inbounds = M_SQRT1_2*sqrt(-tmp_1 + tmp_2 + tmp_3);
    Cart_to_xx1_inbounds = acos(M_SQRT1_2*sqrt(tmp_2*tmp_4 - tmp_3*tmp_4 + 1)*(((Cartz) > 0) - ((Cartz) < 0)));
    Cart_to_xx2_inbounds = atan2(Carty, Cartx);
  }

  // Next compute xxmin[i]. By definition,
  //    xx[i][j] = xxmin[i] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*dxxi;
  // -> xxmin[i] = xx[i][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxxi
  const REAL xxmin[3] = {
    xx[0][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx0,
    xx[1][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx1,
    xx[2][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx2 };

  // Finally compute i{0,1,2}_inbounds (add 0.5 to account for rounding down)
  const int i0_inbounds = (int)( (Cart_to_xx0_inbounds - xxmin[0] - (1.0/2.0)*dxx0 + ((REAL)NGHOSTS)*dxx0)/dxx0 + 0.5 );
  const int i1_inbounds = (int)( (Cart_to_xx1_inbounds - xxmin[1] - (1.0/2.0)*dxx1 + ((REAL)NGHOSTS)*dxx1)/dxx1 + 0.5 );
  const int i2_inbounds = (int)( (Cart_to_xx2_inbounds - xxmin[2] - (1.0/2.0)*dxx2 + ((REAL)NGHOSTS)*dxx2)/dxx2 + 0.5 );

  // Step 3: Convert x0(i0_inbounds),x1(i1_inbounds),x2(i2_inbounds) -> (Cartx,Carty,Cartz),
  //         and check that
  //         (Cartx,Carty,Cartz) == (Cartx(x0(i0)),Cartx(x1(i1)),Cartx(x2(i2)))
  //         If not, error out!

  // Step 3.a: Compute x0(i0_inbounds),x1(i1_inbounds),x2(i2_inbounds):
  const REAL x0_inbounds = xx[0][i0_inbounds];
  const REAL x1_inbounds = xx[1][i1_inbounds];
  const REAL x2_inbounds = xx[2][i2_inbounds];

  // Step 3.b: Compute {x,y,z}Cart_from_xx, as a
  //           function of i0,i1,i2
  REAL xCart_from_xx, yCart_from_xx, zCart_from_xx;
  {
    // xx_to_Cart for Coordinate SinhSymTP):
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
    /*
     *  Original SymPy expressions:
     *  "[xCart_from_xx = AMAX*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*sin(xx1)*cos(xx2)/(exp(1/SINHWAA) - exp(-1/SINHWAA)),
     *    yCart_from_xx = AMAX*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*sin(xx1)*sin(xx2)/(exp(1/SINHWAA) - exp(-1/SINHWAA)),
     *    zCart_from_xx = sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)]"
     */
    const double tmp_0 = (1.0/(SINHWAA));
    const double tmp_1 = exp(tmp_0) - exp(-tmp_0);
    const double tmp_3 = exp(tmp_0*xx0) - exp(-tmp_0*xx0);
    const double tmp_4 = AMAX*tmp_3*sin(xx1)/tmp_1;
    xCart_from_xx = tmp_4*cos(xx2);
    yCart_from_xx = tmp_4*sin(xx2);
    zCart_from_xx = sqrt(((AMAX)*(AMAX))*((tmp_3)*(tmp_3))/((tmp_1)*(tmp_1)) + ((bScale)*(bScale)))*cos(xx1);
  }

  // Step 3.c: Compute {x,y,z}Cart_from_xx_inbounds, as a
  //           function of i0_inbounds,i1_inbounds,i2_inbounds
  REAL xCart_from_xx_inbounds, yCart_from_xx_inbounds, zCart_from_xx_inbounds;
  {
    // xx_to_Cart_inbounds for Coordinate SinhSymTP):
    REAL xx0 = xx[0][i0_inbounds];
    REAL xx1 = xx[1][i1_inbounds];
    REAL xx2 = xx[2][i2_inbounds];
    /*
     *  Original SymPy expressions:
     *  "[xCart_from_xx_inbounds = AMAX*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*sin(xx1)*cos(xx2)/(exp(1/SINHWAA) - exp(-1/SINHWAA)),
     *    yCart_from_xx_inbounds = AMAX*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*sin(xx1)*sin(xx2)/(exp(1/SINHWAA) - exp(-1/SINHWAA)),
     *    zCart_from_xx_inbounds = sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)]"
     */
    const double tmp_0 = (1.0/(SINHWAA));
    const double tmp_1 = exp(tmp_0) - exp(-tmp_0);
    const double tmp_3 = exp(tmp_0*xx0) - exp(-tmp_0*xx0);
    const double tmp_4 = AMAX*tmp_3*sin(xx1)/tmp_1;
    xCart_from_xx_inbounds = tmp_4*cos(xx2);
    yCart_from_xx_inbounds = tmp_4*sin(xx2);
    zCart_from_xx_inbounds = sqrt(((AMAX)*(AMAX))*((tmp_3)*(tmp_3))/((tmp_1)*(tmp_1)) + ((bScale)*(bScale)))*cos(xx1);
  }

  // Step 3.d: Compare xCart_from_xx to xCart_from_xx_inbounds;
  //           they should be identical!!!
#define EPS_REL 1e-8
  const REAL norm_factor = sqrt(xCart_from_xx*xCart_from_xx + yCart_from_xx*yCart_from_xx + zCart_from_xx*zCart_from_xx) + 1e-15;
  if(fabs( (double)(xCart_from_xx - xCart_from_xx_inbounds) ) > EPS_REL * norm_factor ||
     fabs( (double)(yCart_from_xx - yCart_from_xx_inbounds) ) > EPS_REL * norm_factor ||
     fabs( (double)(zCart_from_xx - zCart_from_xx_inbounds) ) > EPS_REL * norm_factor) {
    fprintf(stderr,"Error in SinhSymTP coordinate system: Cartesian disagreement: ( %.15e %.15e %.15e ) != ( %.15e %.15e %.15e ) | xx: %e %e %e -> %e %e %e | %d %d %d\n",
            (double)xCart_from_xx,(double)yCart_from_xx,(double)zCart_from_xx,
            (double)xCart_from_xx_inbounds,(double)yCart_from_xx_inbounds,(double)zCart_from_xx_inbounds,
            xx[0][i0],xx[1][i1],xx[2][i2],
            xx[0][i0_inbounds],xx[1][i1_inbounds],xx[2][i2_inbounds],
            Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2);
    exit(1);
  }

  // Step 4: Set output arrays.
  x0x1x2_inbounds[0] = xx[0][i0_inbounds];
  x0x1x2_inbounds[1] = xx[1][i1_inbounds];
  x0x1x2_inbounds[2] = xx[2][i2_inbounds];
  i0i1i2_inbounds[0] = i0_inbounds;
  i0i1i2_inbounds[1] = i1_inbounds;
  i0i1i2_inbounds[2] = i2_inbounds;
}
/*
 * set_parity_for_inner_boundary_single_pt():
 *   Given (x0,x1,x2)=(xx0,xx1,xx2) and
 *   (x0,x1,x2)'=(x0x1x2_inbounds[0],x0x1x2_inbounds[1],x0x1x2_inbounds[2])
 *   (see description of
 *   EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt()
 *   above for more details), here we compute the parity conditions
 *   for all 10 tensor types supported by NRPy+.
 */
void set_parity_for_inner_boundary_single_pt(const paramstruct *restrict params, const REAL xx0,const REAL xx1,const REAL xx2,
                                             const REAL x0x1x2_inbounds[3], const int idx,
                                             innerpt_bc_struct *restrict innerpt_bc_arr) {
#include "./set_Cparameters.h"


  const REAL xx0_inbounds = x0x1x2_inbounds[0];
  const REAL xx1_inbounds = x0x1x2_inbounds[1];
  const REAL xx2_inbounds = x0x1x2_inbounds[2];

  REAL REAL_parity_array[10];
  {
    // Evaluate dot products needed for setting parity
    //     conditions at a given point (xx0,xx1,xx2),
    //     using C code generated by NRPy+

    // NRPy+ Curvilinear Boundary Conditions: Unit vector dot products for all
    //      ten parity conditions, in given coordinate system.
    //      Needed for automatically determining sign of tensor across coordinate boundary.
    // Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb
    /*
     *  Original SymPy expressions:
     *  "[REAL_parity_array[0] = 1,
     *    REAL_parity_array[1] = AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)),
     *    REAL_parity_array[2] = AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)),
     *    REAL_parity_array[3] = sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds),
     *    REAL_parity_array[4] = (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)))**2,
     *    REAL_parity_array[5] = (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)))*(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2))),
     *    REAL_parity_array[6] = (sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds))*(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2))),
     *    REAL_parity_array[7] = (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)))**2,
     *    REAL_parity_array[8] = (sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds))*(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))*cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**2) + sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*sin(xx1)*sin(xx1_inbounds)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)*sqrt(AMAX**2*(exp(xx0_inbounds/SINHWAA) - exp(-xx0_inbounds/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1_inbounds)**2))),
     *    REAL_parity_array[9] = (sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds))**2]"
     */
    {
      const double tmp_0 = (1.0/(SINHWAA));
      const double tmp_2 = exp(tmp_0*xx0) - exp(-tmp_0*xx0);
      const double tmp_4 = exp(tmp_0*xx0_inbounds) - exp(-tmp_0*xx0_inbounds);
      const double tmp_5 = ((AMAX)*(AMAX))/((exp(tmp_0) - exp(-tmp_0))*(exp(tmp_0) - exp(-tmp_0)));
      const double tmp_6 = ((bScale)*(bScale));
      const double tmp_7 = sin(xx1);
      const double tmp_8 = ((tmp_2)*(tmp_2))*tmp_5;
      const double tmp_9 = sin(xx1_inbounds);
      const double tmp_10 = ((tmp_4)*(tmp_4))*tmp_5;
      const double tmp_11 = 1/(sqrt(tmp_10 + tmp_6*((tmp_9)*(tmp_9)))*sqrt(tmp_6*((tmp_7)*(tmp_7)) + tmp_8));
      const double tmp_12 = tmp_11*tmp_2*tmp_4*tmp_5*cos(xx1)*cos(xx1_inbounds);
      const double tmp_13 = sin(xx2)*sin(xx2_inbounds);
      const double tmp_14 = tmp_11*tmp_7*tmp_9*sqrt(tmp_10 + tmp_6)*sqrt(tmp_6 + tmp_8);
      const double tmp_15 = cos(xx2)*cos(xx2_inbounds);
      const double tmp_16 = tmp_12 + tmp_13*tmp_14 + tmp_14*tmp_15;
      const double tmp_17 = tmp_12*tmp_13 + tmp_12*tmp_15 + tmp_14;
      const double tmp_18 = tmp_13 + tmp_15;
      REAL_parity_array[0] = 1;
      REAL_parity_array[1] = tmp_16;
      REAL_parity_array[2] = tmp_17;
      REAL_parity_array[3] = tmp_18;
      REAL_parity_array[4] = ((tmp_16)*(tmp_16));
      REAL_parity_array[5] = tmp_16*tmp_17;
      REAL_parity_array[6] = tmp_16*tmp_18;
      REAL_parity_array[7] = ((tmp_17)*(tmp_17));
      REAL_parity_array[8] = tmp_17*tmp_18;
      REAL_parity_array[9] = ((tmp_18)*(tmp_18));
    }

  }

  // Next perform sanity check on parity array output: should be +1 or -1 to within 8 significant digits:
  for(int whichparity=0;whichparity<10;whichparity++) {
    //printf("Good? Parity %d evaluated to %e\n",whichparity,(double)REAL_parity_array[whichparity]);
    if( fabs(REAL_parity_array[whichparity]) < 1 - 1e-8 || fabs(REAL_parity_array[whichparity]) > 1 + 1e-8 ) {
      fprintf(stderr,"Error at point (%e %e %e), which maps to (%e %e %e).\n",
              xx0,xx1,xx2, xx0_inbounds,xx1_inbounds,xx2_inbounds);
      fprintf(stderr,"Parity evaluated to %e , which is not within 8 significant digits of +1 or -1.\n",
              REAL_parity_array[whichparity]);
      exit(1);
    }
    // The typecast (int8_t)REAL_parity_array[whichparity] *does not work*.
    //  Thankfully, we've already checked whether REAL_parity_array[whichparity]
    //  is within 1e-8 of +/- 1, so here we just check the sign of the
    //  REAL_parity_array to find the correct value of innerpt_bc_arr[idx].parity[parity].
    for(int parity=0;parity<10;parity++) {
      innerpt_bc_arr[idx].parity[parity] = 1;
      if(REAL_parity_array[parity] < 0) innerpt_bc_arr[idx].parity[parity] = -1;
    }
  } // END for(int whichparity=0;whichparity<10;whichparity++)
}

/*
 * At each coordinate point (x0,x1,x2) situated at grid index (i0,i1,i2):
 * 
 * Step 1: Set up inner boundary structs bcstruct->inner_bc_array[].
 *   Recall that at each inner boundary point we must set innerpt_bc_struct:
 *     typedef struct __innerpt_bc_struct__ {
 *       int dstpt;  // dstpt is the 3D grid index IDX3S(i0,i1,i2) of the inner boundary point (i0,i1,i2)
 *       int srcpt;  // srcpt is the 3D grid index (a la IDX3S) to which the inner boundary point maps
 *       int8_t parity[10];  // parity[10] is a calculation of dot products for the 10 independent parity types
 *     } innerpt_bc_struct;
 *   At each ghostzone (i.e., each point within NGHOSTS points from grid boundary):
 *     Call EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt().
 *         This function converts the curvilinear coordinate (x0,x1,x2) to the corresponding
 *         Cartesian coordinate (x,y,z), then finds the grid point
 *         (i0_inbounds,i1_inbounds,i2_inbounds) in the grid interior or outer boundary
 *         corresponding to this Cartesian coordinate (x,y,z).
 *     If (i0,i1,i2) *is not* the same as (i0_inbounds,i1_inbounds,i2_inbounds),
 *         then we are at an inner boundary point. We must set
 *         Set bcstruct->inner_bc_array for this point, which requires we specify
 *         both (i0_inbounds,i1_inbounds,i2_inbounds) [just found!] and parity
 *         conditions for this gridpoint. The latter is found & specified within the
 *         function set_parity_for_inner_boundary_single_pt().
 *     If (i0,i1,i2) *is* the same as (i0_inbounds,i1_inbounds,i2_inbounds),
 *         then we are at an outer boundary point. Take care of outer BCs in Step 2.
 * Step 2: Set up outer boundary structs bcstruct->outer_bc_array[which_gz][face][idx2d]:
 *   Recall that at each inner boundary point we must set outerpt_bc_struct:
 *     typedef struct __outerpt_bc_struct__ {
 *       short i0,i1,i2;  // the outer boundary point grid index (i0,i1,i2), on the 3D grid
 *       int8_t FACEX0,FACEX1,FACEX2;  // 1-byte integers that store
 *       //                               FACEX0,FACEX1,FACEX2 = +1, 0, 0 if on the i0=i0min face,
 *       //                               FACEX0,FACEX1,FACEX2 = -1, 0, 0 if on the i0=i0max face,
 *       //                               FACEX0,FACEX1,FACEX2 =  0,+1, 0 if on the i1=i2min face,
 *       //                               FACEX0,FACEX1,FACEX2 =  0,-1, 0 if on the i1=i1max face,
 *       //                               FACEX0,FACEX1,FACEX2 =  0, 0,+1 if on the i2=i2min face, or
 *       //                               FACEX0,FACEX1,FACEX2 =  0, 0,-1 if on the i2=i2max face,
 *     } outerpt_bc_struct;
 *   Outer boundary points are filled from the inside out, two faces at a time.
 *     E.g., consider a Cartesian coordinate grid that has 14 points in each direction,
 *     including the ghostzones, with NGHOSTS=2.
 *     We first fill in the lower x0 face with (i0=1,i1={2,11},i2={2,11}). We fill these
 *     points in first, since they will in general (at least in the case of extrapolation
 *     outer BCs) depend on e.g., i0=2 and i0=3 points.
 *     Simultaneously we can fill in the upper x0 face with (i0=12,i1={2,11},i2={2,11}),
 *     since these points depend only on e.g., i0=11 and i0=10 (again assuming extrap. BCs).
 *     Next we can fill in the lower x1 face: (i0={1,12},i1=2,i2={2,11}). Notice these
 *     depend on i0 min and max faces being filled. The remaining pattern goes like this:
 *     Upper x1 face: (i0={1,12},i1=12,i2={2,11})
 *     Lower x2 face: (i0={1,12},i1={1,12},i2=1)
 *     Upper x2 face: (i0={1,12},i1={1,12},i2=12)
 *     Lower x0 face: (i0=0,i1={1,12},i2={1,12})
 *     Upper x0 face: (i0=13,i1={1,12},i2={1,12})
 *     Lower x1 face: (i0={0,13},i1=0,i2={2,11})
 *     Upper x1 face: (i0={0,13},i1=13,i2={2,11})
 *     Lower x2 face: (i0={0,13},i1={0,13},i2=0)
 *     Upper x2 face: (i0={0,13},i1={0,13},i2=13)
 *   Note that we allocate a outerpt_bc_struct at *all* boundary points,
 *     regardless of whether the point is an outer or inner point. However
 *     the struct is set only at outer boundary points. This is slightly
 *     wasteful, but only in memory, not in CPU.
 */
void bcstruct_set_up(const paramstruct *restrict params, REAL *restrict xx[3], bc_struct *restrict bcstruct) {
#include "./set_Cparameters.h"


  ////////////////////////////////////////
  // STEP 1: SET UP INNER BOUNDARY STRUCTS
  {
    // First count the number of inner points.
    int num_inner = 0;
    LOOP_OMP("omp parallel for reduction(+:num_inner)",
             i0,0,Nxx_plus_2NGHOSTS0,  i1,0,Nxx_plus_2NGHOSTS1,  i2,0,Nxx_plus_2NGHOSTS2) {
      const int i0i1i2[3] = { i0,i1,i2 };
      if(!IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2, NGHOSTS)) {
        REAL x0x1x2_inbounds[3];
        int i0i1i2_inbounds[3];
        EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
        if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
          // this is a pure outer boundary point.
        } else {
          // this is an inner boundary point, which maps either
          //  to the grid interior or to an outer boundary point
          num_inner++;
        }
      }
    }
    // Store num_inner to bc_info:
    bcstruct->bc_info.num_inner_boundary_points = num_inner;

    // Next allocate memory for inner_boundary_points:
    bcstruct->inner_bc_array = (innerpt_bc_struct *restrict)malloc( sizeof(innerpt_bc_struct)*num_inner );
  }

  // Then set inner_bc_array:
  {
    int which_inner = 0;
    LOOP_NOOMP(i0,0,Nxx_plus_2NGHOSTS0,  i1,0,Nxx_plus_2NGHOSTS1,  i2,0,Nxx_plus_2NGHOSTS2) {
      const int i0i1i2[3] = { i0,i1,i2 };
      if(!IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2, NGHOSTS)) {
        REAL x0x1x2_inbounds[3];
        int i0i1i2_inbounds[3];
        EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
        if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
          // this is a pure outer boundary point.
        } else {
          bcstruct->inner_bc_array[which_inner].dstpt = IDX3S(i0,i1,i2);
          bcstruct->inner_bc_array[which_inner].srcpt = IDX3S(i0i1i2_inbounds[0],i0i1i2_inbounds[1],i0i1i2_inbounds[2]);
          //printf("%d / %d\n",which_inner, bc_info->num_inner_boundary_points);
          set_parity_for_inner_boundary_single_pt(params, xx[0][i0],xx[1][i1],xx[2][i2],
                                                  x0x1x2_inbounds, which_inner, bcstruct->inner_bc_array);

          which_inner++;
        }
      }
    }
  }

  ////////////////////////////////////////
  // STEP 2: SET UP OUTER BOUNDARY STRUCTS
  // First set up loop bounds for outer boundary condition updates,
  //   store to bc_info->bc_loop_bounds[which_gz][face][]. Also
  //   allocate memory for outer_bc_array[which_gz][face][]:
  int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
  int imax[3] = { Nxx_plus_2NGHOSTS0-NGHOSTS, Nxx_plus_2NGHOSTS1-NGHOSTS, Nxx_plus_2NGHOSTS2-NGHOSTS };
  for(int which_gz=0;which_gz<NGHOSTS;which_gz++) {
    const int x0min_face_range[6] = { imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2] };  imin[0]--;
    const int x0max_face_range[6] = { imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2] };  imax[0]++;
    const int x1min_face_range[6] = { imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2] };  imin[1]--;
    const int x1max_face_range[6] = { imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2] };  imax[1]++;
    const int x2min_face_range[6] = { imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2] };  imin[2]--;
    const int x2max_face_range[6] = { imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1 };  imax[2]++;

    int face=0;
    ////////////////////////
    // x0min and x0max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x0min and x0max faces have exactly the same size.
    bcstruct->pure_outer_bc_array[which_gz + NGHOSTS*(face/2)] = (outerpt_bc_struct *restrict)malloc(sizeof(outerpt_bc_struct) * 2 *
                                                                                                     ((x0min_face_range[1]-x0min_face_range[0]) *
                                                                                                      (x0min_face_range[3]-x0min_face_range[2]) *
                                                                                                      (x0min_face_range[5]-x0min_face_range[4])));
    // x0min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x0min_face_range[i]; }
    face++;
    // x0max face: Set loop bounds & allocate memory for outer_bc_array:
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x0max_face_range[i]; }
    face++;
    ////////////////////////

    ////////////////////////
    // x1min and x1max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x1min and x1max faces have exactly the same size.
    bcstruct->pure_outer_bc_array[which_gz + NGHOSTS*(face/2)] = (outerpt_bc_struct *restrict)malloc(sizeof(outerpt_bc_struct) * 2 *
                                                                                                     ((x1min_face_range[1]-x1min_face_range[0]) *
                                                                                                      (x1min_face_range[3]-x1min_face_range[2]) *
                                                                                                      (x1min_face_range[5]-x1min_face_range[4])));
    // x1min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x1min_face_range[i]; }
    face++;
    // x1max face: Set loop bounds & allocate memory for outer_bc_array:
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x1max_face_range[i]; }
    face++;
    ////////////////////////


    ////////////////////////
    // x2min and x2max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x2min and x2max faces have exactly the same size.
    bcstruct->pure_outer_bc_array[which_gz + NGHOSTS*(face/2)] = (outerpt_bc_struct *restrict)malloc(sizeof(outerpt_bc_struct) * 2 *
                                                                                                     ((x2min_face_range[1]-x2min_face_range[0]) *
                                                                                                      (x2min_face_range[3]-x2min_face_range[2]) *
                                                                                                      (x2min_face_range[5]-x2min_face_range[4])));
    // x2min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x2min_face_range[i]; }
    face++;
    // x2max face: Set loop bounds & allocate memory for outer_bc_array:
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x2max_face_range[i]; }
    face++;
    ////////////////////////
  }

  for(int which_gz=0;which_gz<NGHOSTS;which_gz++) for(int dirn=0;dirn<3;dirn++) {
      int idx2d = 0;
      // LOWER FACE: dirn=0 -> x0min; dirn=1 -> x1min; dirn=2 -> x2min
      {
        const int face = dirn*2;
#define IDX2D_BCS(i0,i0min,i0max, i1,i1min,i1max ,i2,i2min,i2max)       \
        ( ((i0)-(i0min)) + ((i0max)-(i0min)) * ( ((i1)-(i1min)) + ((i1max)-(i1min)) * ((i2)-(i2min)) ) )
        const int FACEX0=(face==0) - (face==1); // +1 if face==0 (x0min) ; -1 if face==1 (x0max). Otherwise 0.
        const int FACEX1=(face==2) - (face==3); // +1 if face==2 (x1min) ; -1 if face==3 (x1max). Otherwise 0.
        const int FACEX2=(face==4) - (face==5); // +1 if face==4 (x2min) ; -1 if face==5 (x2max). Otherwise 0.
        LOOP_NOOMP(i0,bcstruct->bc_info.bc_loop_bounds[which_gz][face][0],bcstruct->bc_info.bc_loop_bounds[which_gz][face][1],
                   i1,bcstruct->bc_info.bc_loop_bounds[which_gz][face][2],bcstruct->bc_info.bc_loop_bounds[which_gz][face][3],
                   i2,bcstruct->bc_info.bc_loop_bounds[which_gz][face][4],bcstruct->bc_info.bc_loop_bounds[which_gz][face][5]) {
          REAL x0x1x2_inbounds[3];
          int i0i1i2_inbounds[3];
          EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
          if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0 = i0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1 = i1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2 = i2;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0 = FACEX0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1 = FACEX1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2 = FACEX2;
            idx2d++;
          }
        }
      }
      // UPPER FACE: dirn=0 -> x0max; dirn=1 -> x1max; dirn=2 -> x2max
      {
        const int face = dirn*2+1;
        const int FACEX0=(face==0) - (face==1); // +1 if face==0 ; -1 if face==1. Otherwise 0.
        const int FACEX1=(face==2) - (face==3); // +1 if face==2 ; -1 if face==3. Otherwise 0.
        const int FACEX2=(face==4) - (face==5); // +1 if face==4 ; -1 if face==5. Otherwise 0.
        LOOP_NOOMP(i0,bcstruct->bc_info.bc_loop_bounds[which_gz][face][0],bcstruct->bc_info.bc_loop_bounds[which_gz][face][1],
                   i1,bcstruct->bc_info.bc_loop_bounds[which_gz][face][2],bcstruct->bc_info.bc_loop_bounds[which_gz][face][3],
                   i2,bcstruct->bc_info.bc_loop_bounds[which_gz][face][4],bcstruct->bc_info.bc_loop_bounds[which_gz][face][5]) {
          REAL x0x1x2_inbounds[3];
          int i0i1i2_inbounds[3];
          EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
          if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0 = i0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1 = i1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2 = i2;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0 = FACEX0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1 = FACEX1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2 = FACEX2;
            idx2d++;
          }
        }
      }
      bcstruct->bc_info.num_pure_outer_boundary_points[which_gz][dirn] = idx2d;
    }
}
