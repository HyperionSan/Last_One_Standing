/*
    Copyright 2017 Philippe Grandclement

    This file is part of Kadath.

    Kadath is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Kadath is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Kadath.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "headcpp.hpp"

#include "utilities.hpp"
#include "spheric.hpp"
#include "point.hpp"
#include "val_domain.hpp"
namespace Kadath {
// Standard constructor
Domain_shell::Domain_shell (int num, int ttype, double rint, double rext, const Point& cr, const Dim_array& nbr) : Domain(num, ttype, nbr),
							 alpha((rext-rint)/2.), beta((rext+rint)/2.), center(cr) {
    
     assert (nbr.get_ndim()==3) ;
     assert (cr.get_ndim()==3) ;    // Affectation de type_point :
     do_coloc() ;

}


// Constructor by copy
Domain_shell::Domain_shell (const Domain_shell& so) : Domain(so), alpha(so.alpha), 
						beta(so.beta), center(so.center){}



Domain_shell::Domain_shell (int num, FILE* fd) : Domain(num, fd), center(fd) {
	fread_be (&alpha, sizeof(double), 1, fd) ;
	fread_be (&beta, sizeof(double), 1, fd) ;
	do_coloc() ;
}

// Destructor
Domain_shell::~Domain_shell() {}

void Domain_shell::save (FILE* fd) const {
	nbr_points.save(fd) ;
	nbr_coefs.save(fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;
	fwrite_be (&type_base, sizeof(int), 1, fd) ;	
	center.save(fd) ;
	fwrite_be (&alpha, sizeof(double), 1, fd) ;	
	fwrite_be (&beta, sizeof(double), 1, fd) ;
}

ostream& Domain_shell::print (ostream& o) const {
  o << "Shell" << endl ;
  o <<  beta-alpha << " < r < " << beta+alpha << endl ;
  o << "Center  = " << center << endl ;
  o << "Nbr pts = " << nbr_points << endl ;
  o << endl ;
  return o ;
}



Val_domain Domain_shell::der_normal (const Val_domain& so, int bound) const {

	Val_domain res (so.der_var(1)) ;
	switch (bound) {
		case OUTER_BC :
			res /= alpha ;
			break ;
		case INNER_BC :
			res /= alpha ;
			break ;
		default:
			cerr << "Unknown boundary case in Domain_shell::der_normal" << endl ;
			abort() ;
		}
return res ;
}

// Computes the Cartesian coordinates
void Domain_shell::do_absol () const  {
	for (int i=0 ; i<3 ; i++)
	   assert (coloc[i] != 0x0) ;
	for (int i=0 ; i<3 ; i++)
	   assert (absol[i] == 0x0) ;
	for (int i=0 ; i<3 ; i++) {
	   absol[i] = new Val_domain(this) ;
	   absol[i]->allocate_conf() ;
	   }
	Index index (nbr_points) ;

	do  {
		absol[0]->set(index) = (alpha* ((*coloc[0])(index(0))) +beta)*
		 sin((*coloc[1])(index(1)))*cos((*coloc[2])(index(2))) + center(1) ;
		absol[1]->set(index) = (alpha* ((*coloc[0])(index(0))) +beta) *
		 sin((*coloc[1])(index(1)))*sin((*coloc[2])(index(2))) + center(2) ;
		absol[2]->set(index) = (alpha* ((*coloc[0])(index(0))) + beta) * cos((*coloc[1])(index(1))) + center(3) ;
			}
	while (index.inc())  ;
}

// Computes the radius
void Domain_shell::do_radius () const  {
	for (int i=0 ; i<3 ; i++)
	   assert (coloc[i] != 0x0) ;
	assert (radius == 0x0) ;
	radius = new Val_domain(this) ;
	radius->allocate_conf() ;
	Index index (nbr_points) ;
	do 
		radius->set(index) = alpha* ((*coloc[0])(index(0))) + beta ;
	while (index.inc())  ;
}

// Computes the Cartesian coordinates
void Domain_shell::do_cart () const  {
	for (int i=0 ; i<3 ; i++)
	   assert (coloc[i] != 0x0) ;
	for (int i=0 ; i<3 ; i++)
	   assert (cart[i] == 0x0) ;
	for (int i=0 ; i<3 ; i++) {
	   cart[i] = new Val_domain(this) ;
	   cart[i]->allocate_conf() ;
	   }
	Index index (nbr_points) ;

	do  {
		cart[0]->set(index) = (alpha* ((*coloc[0])(index(0))) +beta)*
		 sin((*coloc[1])(index(1)))*cos((*coloc[2])(index(2))) + center(1) ;
		cart[1]->set(index) = (alpha* ((*coloc[0])(index(0))) +beta) *
		 sin((*coloc[1])(index(1)))*sin((*coloc[2])(index(2))) + center(2) ;
		cart[2]->set(index) = (alpha* ((*coloc[0])(index(0))) + beta) * cos((*coloc[1])(index(1))) + center(3) ;
			}
	while (index.inc())  ;
}

// Computes the Cartesian coordinates
void Domain_shell::do_cart_surr () const  {
	for (int i=0 ; i<3 ; i++)
	   assert (coloc[i] != 0x0) ;
	for (int i=0 ; i<3 ; i++)
	   assert (cart_surr[i] == 0x0) ;
	for (int i=0 ; i<3 ; i++) {
	   cart_surr[i] = new Val_domain(this) ;
	   cart_surr[i]->allocate_conf() ;
	   }
	Index index (nbr_points) ;

	do  {
		cart_surr[0]->set(index) = sin((*coloc[1])(index(1)))*cos((*coloc[2])(index(2))) ;
		cart_surr[1]->set(index) = sin((*coloc[1])(index(1)))*sin((*coloc[2])(index(2))) ;
		cart_surr[2]->set(index) =  cos((*coloc[1])(index(1))) ;
			}
	while (index.inc())  ;
}

// Check if a point is inside this domain
bool Domain_shell::is_in (const Point& xx, double prec) const {

	assert (xx.get_ndim()==3) ;
	
	double x_loc = xx(1) - center(1) ;
	double y_loc = xx(2) - center(2) ;
	double z_loc = xx(3) - center(3) ;
	double air_loc = sqrt (x_loc*x_loc + y_loc*y_loc + z_loc*z_loc) ;

	bool res = ((air_loc / (alpha+beta) - 1  <= prec) && (air_loc/ (beta-alpha) - 1 >= -prec)) ? true : false ;
	return res ;
}

//Computes the numerical coordinates, as a function of the absolute ones.
const Point Domain_shell::absol_to_num(const Point& abs) const {

	assert (is_in(abs)) ;
	Point num(3) ;
	
	double x_loc = abs(1) - center(1) ;
	double y_loc = abs(2) - center(2) ;
	double z_loc = abs(3) - center(3) ;
	double air = sqrt(x_loc*x_loc+y_loc*y_loc+z_loc*z_loc) ;
	num.set(1) = (air-beta)/alpha ;
	double rho = sqrt(x_loc*x_loc+y_loc*y_loc) ;
	
	if (rho==0) {
	    // On the axis ?
	    num.set(2) = (z_loc>=0) ? 0 : M_PI ;
	    num.set(3) = 0 ;
	}
 	else {
	    num.set(2) = atan(rho/z_loc) ;
	    num.set(3) = atan2 (y_loc, x_loc) ;
        }	
	
	if (num(2) <0)
	    num.set(2) = M_PI + num(2) ;

	return num ;
}

// Convert absolute coordinates to numerical ones
const Point Domain_shell::absol_to_num_bound(const Point& abs, int bound) const {

	assert ((bound==OUTER_BC) || (bound==INNER_BC)) ;
	assert (is_in(abs, 1e-3)) ;
	Point num(3) ;
	
	double x_loc = abs(1) - center(1) ;
	double y_loc = abs(2) - center(2) ;
	double z_loc = abs(3) - center(3) ;
	
	switch (bound) {
	  case INNER_BC:
	    num.set(1) = -1 ;
	    break ;
	  case OUTER_BC:
	    num.set(1) = 1 ;
	    break ;
	  default:
	      cerr << "unknown boundary in Domain_shell::absol_to_num" << endl ;
	      abort() ;
	}
	
	double rho = sqrt(x_loc*x_loc+y_loc*y_loc) ;
	
	if (rho==0) {
	    // Sur l'axe
	    num.set(2) = (z_loc>=0) ? 0 : M_PI ;
	    num.set(3) = 0 ;
	}
 	else {
	    num.set(2) = atan(rho/z_loc) ;
	    num.set(3) = atan2 (y_loc, x_loc) ;
        }	
	
	if (num(2) <0)
	    num.set(2) = M_PI + num(2) ;
	
	return num ;
}
double coloc_leg (int, int) ;
void Domain_shell::do_coloc () {

	switch (type_base) {
		case CHEB_TYPE:
			nbr_coefs = nbr_points ;
			nbr_coefs.set(2) += 2 ;
			del_deriv() ;
			for (int i=0 ; i<ndim ; i++)
				coloc[i] = new Array<double> (nbr_points(i)) ;
			for (int i=0 ; i<nbr_points(0) ; i++)
				coloc[0]->set(i) = -cos(M_PI*i/(nbr_points(0)-1)) ;
			for (int j=0 ; j<nbr_points(1) ; j++)
				coloc[1]->set(j) = M_PI/2.*j/(nbr_points(1)-1) ;
			for (int k=0 ; k<nbr_points(2) ; k++)
				coloc[2]->set(k) = M_PI*2.*k/nbr_points(2) ;
			break ;
		case LEG_TYPE:
			nbr_coefs = nbr_points ;
			nbr_coefs.set(2) += 2 ;
			del_deriv() ;
			for (int i=0 ; i<ndim ; i++) 
				coloc[i] = new Array<double> (nbr_points(i)) ;
			for (int i=0 ; i<nbr_points(0) ; i++)
				coloc[0]->set(i) = coloc_leg(i, nbr_points(0)) ;
			for (int j=0 ; j<nbr_points(1) ; j++)
				coloc[1]->set(j) = M_PI/2.*j/(nbr_points(1)-1) ;
			for (int k=0 ; k<nbr_points(2) ; k++)
				coloc[2]->set(k) = M_PI*2.*k/nbr_points(2) ;
			break ;
		default :
			cerr << "Unknown type of basis in Domain_shell::do_coloc" << endl ;
			abort() ;
	}
}

// standard base for a symetric function in z, using Chebyshev
void Domain_shell::set_cheb_base(Base_spectral& base) const  {

	int m ;

	assert (type_base == CHEB_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? COS_EVEN : SIN_ODD ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = CHEB ;
		 }
	}	
}

void Domain_shell::set_cheb_r_base(Base_spectral& base) const  {
  set_cheb_base(base) ;
}

void Domain_shell::set_legendre_r_base(Base_spectral& base) const  {
  set_legendre_base(base) ;
}
  
// standard base for a anti-symetric function in z, using Chebyshev
void Domain_shell::set_anti_cheb_base(Base_spectral& base) const  {

	int m ;

	assert (type_base == CHEB_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? COS_ODD : SIN_EVEN ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = CHEB ;
		 }
	}	
}

void Domain_shell::set_cheb_base_r_spher(Base_spectral& base) const  {

	int m ;

	assert (type_base == CHEB_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? COS_EVEN : SIN_ODD ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = CHEB ;
		 }
	}	
}

void Domain_shell::set_cheb_base_t_spher(Base_spectral& base) const  {

	int m ;

	assert (type_base == CHEB_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? SIN_EVEN : COS_ODD ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = CHEB ;
		 }
	}	
}

void Domain_shell::set_cheb_base_p_spher(Base_spectral& base) const  {

	int m ;

	assert (type_base == CHEB_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? SIN_ODD : COS_EVEN ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = CHEB ;
		 }
	}	
}
void Domain_shell::set_cheb_base_r_mtz(Base_spectral& base) const  {

	int m ;

	assert (type_base == CHEB_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? COS_EVEN : SIN_ODD ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = CHEB ;
		 }
	}	
}

void Domain_shell::set_cheb_base_t_mtz(Base_spectral& base) const  {

	int m ;

	assert (type_base == CHEB_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? SIN_ODD : COS_EVEN ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = CHEB ;
		 }
	}	
}

void Domain_shell::set_cheb_base_p_mtz(Base_spectral& base) const  {

	int m ;

	assert (type_base == CHEB_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? SIN_EVEN: COS_ODD ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = CHEB ;
		 }
	}	
}

// standard base for a symetric function in z, using Legendre
void Domain_shell::set_legendre_base(Base_spectral& base) const {

	int m ;

	assert (type_base == LEG_TYPE) ;
	base.allocate(nbr_coefs) ;
	
	base.def = true ;    
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? COS_EVEN : SIN_ODD ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = LEG ;
		 }
	}
}

// standard base for a anti-symetric function in z, using Legendre
void Domain_shell::set_anti_legendre_base(Base_spectral& base) const {

	int m ;
	assert (type_base == LEG_TYPE) ;

	base.allocate(nbr_coefs) ;
	
	base.def = true ;    
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? COS_ODD : SIN_EVEN ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = LEG ;
		 }
	}
}

void Domain_shell::set_legendre_base_r_spher(Base_spectral& base) const  {

	int m ;

	assert (type_base == LEG_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? COS_EVEN : SIN_ODD ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = LEG ;
		 }
	}	
}

void Domain_shell::set_legendre_base_t_spher(Base_spectral& base) const  {

	int m ;

	assert (type_base == LEG_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? SIN_EVEN : COS_ODD ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = LEG ;
		 }
	}	
}

void Domain_shell::set_legendre_base_p_spher(Base_spectral& base) const  {

	int m ;

	assert (type_base == LEG_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? SIN_ODD : COS_EVEN ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = LEG ;
		 }
	}	
}

void Domain_shell::set_legendre_base_r_mtz(Base_spectral& base) const  {

	int m ;

	assert (type_base == LEG_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? COS_EVEN : SIN_ODD ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = LEG ;
		 }
	}	
}

void Domain_shell::set_legendre_base_t_mtz(Base_spectral& base) const  {

	int m ;

	assert (type_base == LEG_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? SIN_ODD : COS_EVEN ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = LEG ;
		 }
	}	
}

void Domain_shell::set_legendre_base_p_mtz(Base_spectral& base) const  {

	int m ;

	assert (type_base == LEG_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	base.bases_1d[2]->set(0) = COSSIN ;
	
	Index index (base.bases_1d[0]->get_dimensions()) ;
	
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
	        m = (k%2==0) ? k/2 : (k-1)/2 ;
		base.bases_1d[1]->set(k) = (m%2==0) ? SIN_EVEN : COS_ODD ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = LEG ;
		 }
	}	
}

// Computes the derivatives with respect to XYZ as a function of the numerical ones.
void Domain_shell::do_der_abs_from_der_var(Val_domain** der_var, Val_domain** der_abs) const {

	// d/dx :
	Val_domain sintdr (der_var[0]->mult_sin_theta()/alpha) ;	
	Val_domain dtsr (*der_var[1]/get_radius()) ;
	dtsr.set_base() = der_var[1]->get_base() ;
	Val_domain dpsr (*der_var[2]/get_radius()) ;	
	dpsr.set_base() = der_var[2]->get_base() ;
	Val_domain costdtsr (dtsr.mult_cos_theta()) ;

	Val_domain dpsrssint (dpsr.div_sin_theta()) ;
	der_abs[0] = new Val_domain ((sintdr+costdtsr).mult_cos_phi() - dpsrssint.mult_sin_phi()) ;

	// d/dy :
	der_abs[1] = new Val_domain ((sintdr+costdtsr).mult_sin_phi() + dpsrssint.mult_cos_phi()) ;
	// d/dz :
	der_abs[2] = new Val_domain (der_var[0]->mult_cos_theta()/alpha - dtsr.mult_sin_theta()) ;
}

// Rules for the multiplication of two basis.
Base_spectral Domain_shell::mult (const Base_spectral& a, const Base_spectral& b) const {

	assert (a.ndim==3) ;
	assert (b.ndim==3) ;
	
	Base_spectral res(3) ;
	bool res_def = true ;

	if (!a.def)
		res_def=false ;
	if (!b.def)
		res_def=false ;
		
	if (res_def) {

	// Base in phi :
	res.bases_1d[2] = new Array<int> (a.bases_1d[2]->get_dimensions()) ;
	switch ((*a.bases_1d[2])(0)) {
		case COSSIN:
			switch ((*b.bases_1d[2])(0)) {
				case COSSIN:
					res.bases_1d[2]->set(0) = COSSIN ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		default:
			res_def = false ;
			break ;
	}

	// Bases in theta :
	// On check l'alternance :
	Index index_1 (a.bases_1d[1]->get_dimensions()) ;
	res.bases_1d[1] = new Array<int> (a.bases_1d[1]->get_dimensions()) ;
	do {
	switch ((*a.bases_1d[1])(index_1)) {
		case COS_EVEN:
			switch ((*b.bases_1d[1])(index_1)) {
				case COS_EVEN:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? COS_EVEN : SIN_ODD ;
					break ;
				case COS_ODD:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? COS_ODD : SIN_EVEN ;
					break ;
				case SIN_EVEN:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? SIN_EVEN : COS_ODD ;
					break ;
				case SIN_ODD:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? SIN_ODD : COS_EVEN ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		case COS_ODD:
			switch ((*b.bases_1d[1])(index_1)) {
				case COS_EVEN:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? COS_ODD : SIN_EVEN ;
					break ;
				case COS_ODD:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? COS_EVEN : SIN_ODD ;
					break ;
				case SIN_EVEN:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? SIN_ODD : COS_EVEN ;
					break ;
				case SIN_ODD:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? SIN_EVEN : COS_ODD ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		case SIN_EVEN:
			switch ((*b.bases_1d[1])(index_1)) {
				case COS_EVEN:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? SIN_EVEN : COS_ODD ;
					break ;
				case COS_ODD:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? SIN_ODD : COS_EVEN ;
					break ;
				case SIN_EVEN:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? COS_EVEN : SIN_ODD ;
					break ;
				case SIN_ODD:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? COS_ODD : SIN_EVEN ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		case SIN_ODD:
			switch ((*b.bases_1d[1])(index_1)) {
				case COS_EVEN:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? SIN_ODD : COS_EVEN ;
					break ;
				case COS_ODD:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? SIN_EVEN : COS_ODD ;
					break ;
				case SIN_EVEN:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? COS_ODD : SIN_EVEN ;
					break ;
				case SIN_ODD:
					res.bases_1d[1]->set(index_1) = (index_1(0)%4<2) ? COS_EVEN : SIN_ODD ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		default:
			res_def = false ;
			break ;
		}
	}
	while (index_1.inc()) ;


	// Base in r :
	Index index_0 (a.bases_1d[0]->get_dimensions()) ;
	res.bases_1d[0] = new Array<int> (a.bases_1d[0]->get_dimensions()) ;
	do {
	switch ((*a.bases_1d[0])(index_0)) {
		case CHEB:
			switch ((*b.bases_1d[0])(index_0)) {
				case CHEB:
					res.bases_1d[0]->set(index_0) = CHEB ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		case LEG:
			switch ((*b.bases_1d[0])(index_0)) {
				case LEG:
					res.bases_1d[0]->set(index_0) = LEG ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		default:
			res_def = false ;
			break ;
		}
	}
	while (index_0.inc()) ;
	}
		
	if (!res_def) 
		for (int dim=0 ; dim<a.ndim ; dim++)
			if (res.bases_1d[dim]!= 0x0) {
				delete res.bases_1d[dim] ;
				res.bases_1d[dim] = 0x0 ;
				}
	res.def = res_def ;
	return res ;
}

int Domain_shell::give_place_var (char* p) const {
    int res = -1 ;
    if (strcmp(p,"R ")==0)
	res = 0 ;
    if (strcmp(p,"T ")==0)
	res = 1 ;
    if (strcmp(p,"P ")==0)
	res = 2 ;
    return res ;
}

double Domain_shell::integ(const Val_domain& so, int bound) const
{
   Val_domain rrso (mult_r(mult_r(mult_sin_theta(so)))) ;
	
	double res = 0 ;
	if (!so.check_if_zero())
	{

	  int baset = (*rrso.get_base().bases_1d[1]) (0) ;
	  Index pcf (nbr_coefs) ;
	switch (baset) {
	  case COS_ODD :
	break ;
	case SIN_EVEN :
	  break ;
	case COS_EVEN : {
	res += M_PI*val_boundary(bound, rrso, pcf) ;
	break ;
    }
    case SIN_ODD : {
	  for (int j=0 ; j<nbr_coefs(1) ; j++) {
	    pcf.set(1) = j ;
	    res += 2./(2*double(j)+1) * val_boundary(bound, rrso, pcf) ;
	  }
      break ;
    }
    
    default : 
      cerr << "Case not yet implemented in Domain_shell::integ" << endl ;
      abort() ;
  }
    res *= 2*M_PI ;
	}
    return res ;
}
}
