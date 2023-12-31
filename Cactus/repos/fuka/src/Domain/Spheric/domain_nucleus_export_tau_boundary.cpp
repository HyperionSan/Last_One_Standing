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

#include "spheric.hpp"
#include "scalar.hpp"
#include "tensor.hpp"

namespace Kadath {
void Domain_nucleus::export_tau_val_domain_boundary_vr (const Val_domain& so, int bound, Array<double>& sec, int& pos_sec, int ncond) const {
	
	if (so.check_if_zero())
		pos_sec += ncond ;
	else {
	so.coef() ;
	Index pos_cf (nbr_coefs) ;
	Index pos_galerkin (nbr_coefs) ;

	//k==0 ;
	{
	  pos_cf.set(2) = 0 ; 
	int baset = (*so.get_base().bases_1d[1]) (0) ;
	assert (baset==COS_EVEN) ;
	for (int j=0 ; j<nbr_coefs(1) ; j++) {
		pos_cf.set(1) = j ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf) ;
		pos_sec ++ ;
	}
	}
	
	// k==2
	if (nbr_coefs(2)-1>2)
	{
	  pos_cf.set(2) = 2 ; 
	int baset = (*so.get_base().bases_1d[1]) (2) ;
	assert (baset==SIN_ODD) ;
	for (int j=0 ; j<nbr_coefs(1)-1 ; j++) {
		pos_cf.set(1) = j ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf) ;
		pos_sec ++ ;
	}
	}
	
	// k==3
	if (nbr_coefs(2)-1>3)
	{
	  pos_cf.set(2) = 3 ; 
	int baset = (*so.get_base().bases_1d[1]) (3) ;
	assert (baset==SIN_ODD) ;
	for (int j=0 ; j<nbr_coefs(1)-1 ; j++) {
		pos_cf.set(1) = j ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf) ;
		pos_sec ++ ;
	}
	}
	
	// others
	for (int k=4 ; k<nbr_coefs(2)-1 ; k++) {
	  
	  pos_cf.set(2) = k ;
	  int baset = (*so.get_base().bases_1d[1]) (k) ;
	  int mquant = (k%2==0) ? int(k/2) : int((k-1)/2) ;
	  
	  if (mquant%2==0) {
	    assert (baset==COS_EVEN) ;
	    for (int j=1 ; j<nbr_coefs(1) ; j++) {
	      	pos_cf.set(1) = j ;
		pos_galerkin = pos_cf ;
		pos_galerkin.set(1) = 0 ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf) -2.*val_boundary(bound, so, pos_galerkin) ;
		pos_sec ++ ;
	  }
	  }
	  
	  if (mquant%2==1) {
	    assert (baset==SIN_ODD) ;
	    for (int j=0 ; j<nbr_coefs(1)-1 ; j++) {
	      	pos_cf.set(1) = j ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf)  ;
		pos_sec ++ ;
	  }
	  }
	}
	}
}

void Domain_nucleus::export_tau_val_domain_boundary_vt (const Val_domain& so, int bound, Array<double>& sec, int& pos_sec, int ncond) const {
	
	if (so.check_if_zero())
		pos_sec += ncond ;
	else {
	  
	so.coef() ;
	Index pos_cf (nbr_coefs) ;
	Index pos_galerkin (nbr_coefs) ;

	//k==0 ;
	{
	  pos_cf.set(2) = 0 ; 
	int baset = (*so.get_base().bases_1d[1]) (0) ;
	assert (baset==SIN_EVEN) ;
	for (int j=1 ; j<nbr_coefs(1)-1 ; j++) {
		pos_cf.set(1) = j ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf) ; 
		pos_sec ++ ;
	}
	}
	
	// k==2
	if (nbr_coefs(2)-1>2)
	{
	  pos_cf.set(2) = 2 ; 
	int baset = (*so.get_base().bases_1d[1]) (2) ;
	assert (baset==COS_ODD) ;
	for (int j=0 ; j<nbr_coefs(1)-1 ; j++) {
		pos_cf.set(1) = j ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf) ;
		pos_sec ++ ;
	}
	}
	
	// k==3
	if (nbr_coefs(2)-1>3)
	{
	  pos_cf.set(2) = 3 ; 
	int baset = (*so.get_base().bases_1d[1]) (3) ;
	assert (baset==COS_ODD) ;
	for (int j=0 ; j<nbr_coefs(1)-1 ; j++) {
		pos_cf.set(1) = j ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf) ;
		pos_sec ++ ;
	}
	}
	
	// others
	for (int k=4 ; k<nbr_coefs(2)-1 ; k++) {
	  
	  pos_cf.set(2) = k ;
	  int baset = (*so.get_base().bases_1d[1]) (k) ;
	  int mquant = (k%2==0) ? int(k/2) : int((k-1)/2) ;
	  
	  if (mquant%2==0) {
	    assert (baset==SIN_EVEN) ;
	    for (int j= 1 ; j<nbr_coefs(1)-1 ; j++) {
	      	pos_cf.set(1) = j ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf)  ;
		pos_sec ++ ;
	  }
	  }
	  
	  if (mquant%2==1) {
	    assert (baset==COS_ODD) ;
	      for (int j=1 ; j<nbr_coefs(1)-1 ; j++) {
	      	pos_cf.set(1) = j ;
		pos_galerkin = pos_cf ;
		pos_galerkin.set(1) = 0 ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf) -val_boundary(bound, so, pos_galerkin) ;
		pos_sec ++ ;
	  }
	  }
	  }
	}
	
	}



void Domain_nucleus::export_tau_val_domain_boundary_vp (const Val_domain& so, int bound, Array<double>& sec, int& pos_sec, int ncond) const {
	
	if (so.check_if_zero())
		pos_sec += ncond ;
	else {
	so.coef() ;
	Index pos_cf (nbr_coefs) ;
	Index pos_galerkin (nbr_coefs) ;

	// k==0
	{
	  pos_cf.set(2) = 0 ; 
	int baset = (*so.get_base().bases_1d[1]) (0) ;
	assert (baset==SIN_ODD) ;
	for (int j=0 ; j<nbr_coefs(1)-1 ; j++) {
		pos_cf.set(1) = j ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf) ;
		pos_sec ++ ;
	}
	}
	
	// k==2
	if (nbr_coefs(2)-1>2)
	{
	  pos_cf.set(2) = 2 ; 
	int baset = (*so.get_base().bases_1d[1]) (2) ;
	assert (baset==COS_EVEN) ;
	for (int j=0 ; j<nbr_coefs(1) ; j++) {
		pos_cf.set(1) = j ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf) ;
		pos_sec ++ ;
	}
	}
	
	// k==3
	if (nbr_coefs(2)-1>3)
	{
	  pos_cf.set(2) = 3 ; 
	int baset = (*so.get_base().bases_1d[1]) (3) ;
	assert (baset==COS_EVEN) ;
	for (int j=0 ; j<nbr_coefs(1) ; j++) {
		pos_cf.set(1) = j ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf) ;
		pos_sec ++ ;
	}
	}
	
	// others
	for (int k=4 ; k<nbr_coefs(2)-1 ; k++) {
	  
	  pos_cf.set(2) = k ;
	  int baset = (*so.get_base().bases_1d[1]) (k) ;
	  int mquant = (k%2==0) ? int(k/2) : int((k-1)/2) ;
	  
	  if (mquant%2==0) {
	    assert (baset==SIN_ODD) ;
	    for (int j= 0 ; j<nbr_coefs(1)-1 ; j++) {
	      	pos_cf.set(1) = j ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf)  ;
		pos_sec ++ ;
	  }
	  }
	  
	  if (mquant%2==1) {
	    assert (baset==COS_EVEN) ;
	   for (int j=1 ; j<nbr_coefs(1) ; j++) {
	      	pos_cf.set(1) = j ;
		pos_galerkin = pos_cf ;
		pos_galerkin.set(1) = 0 ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf) -2.*val_boundary(bound, so, pos_galerkin) ;
		pos_sec ++ ;
	  }
	  }
	}
	}
}

void Domain_nucleus::export_tau_val_domain_boundary (const Val_domain& so, int mlim, int bound, Array<double>& sec, int& pos_sec, int ncond) const {
	
	if (so.check_if_zero())
		pos_sec += ncond ;
	else {
	so.coef() ;
	int kmin = 2*mlim + 2 ;
	Index pos_cf (nbr_coefs) ;
	Index pos_galerkin (nbr_coefs) ;

	// Loop on phi :
	for (int k=0 ; k<nbr_coefs(2)-1 ; k++)
		if (k!=1) {
			pos_cf.set(2) = k ; 
			// Loop on theta
			int baset = (*so.get_base().bases_1d[1]) (k) ;
			for (int j=0 ; j<nbr_coefs(1) ; j++) {
				pos_cf.set(1) = j ;
				switch (baset) {
					case COS_EVEN:
						if (k<kmin) {
								sec.set(pos_sec) = val_boundary(bound, so, pos_cf) ;
								pos_sec ++ ;
							}
							else if (j!=0) {
								// Galerkin base
								pos_galerkin = pos_cf ;
								pos_galerkin.set(1) = 0 ;
								sec.set(pos_sec) = val_boundary(bound, so, pos_cf) 
									-2.*val_boundary(bound, so, pos_galerkin) ;
								pos_sec ++ ;
							}
						break ;
					case COS_ODD:
						if (j!=nbr_coefs(1)-1) {
							if (k<kmin) {
								sec.set(pos_sec) = val_boundary(bound, so, pos_cf) ;
								pos_sec ++ ;
							}
							else if (j!=0) {
								// Galerkin base
								pos_galerkin = pos_cf ;
								pos_galerkin.set(1) = 0 ;
								sec.set(pos_sec) = val_boundary(bound, so, pos_cf) 
									-val_boundary(bound, so, pos_galerkin) ;
								pos_sec ++ ;
							}}
						break ;
					case SIN_EVEN:
						if ((j!=0) && (j!=nbr_coefs(1)-1)) {
							if (k<kmin+2) {
							sec.set(pos_sec) = val_boundary(bound, so, pos_cf) ;
							pos_sec ++ ;
							}
							else if (j!=1) {
								// Galerkin base
								// Galerkin base
								pos_galerkin = pos_cf ;
								pos_galerkin.set(1) = 1 ;
								sec.set(pos_sec) = val_boundary(bound, so, pos_cf) 
									-j*val_boundary(bound, so, pos_galerkin) ;
								pos_sec ++ ;
							}
							}
						break ;
					case SIN_ODD:
						if (j!=nbr_coefs(1)-1) {
							if (k<kmin+2) {
								sec.set(pos_sec) = val_boundary(bound, so, pos_cf) ;
								pos_sec ++ ;
							}
							else if (j!=0) {
								// Galerkin base
								pos_galerkin = pos_cf ;
								pos_galerkin.set(1) = 0 ;
								sec.set(pos_sec) = val_boundary(bound, so, pos_cf) 
									-(2*j+1)*val_boundary(bound, so, pos_galerkin) ;
								pos_sec ++ ;
							}}
						break ;
					default:
						cerr << "Unknow theta basis in Domain_nucleus::export_tau_val_domain_boundary" << endl ;
						abort() ;
					}
				}
		}
	}
}

void Domain_nucleus::export_tau_boundary (const Tensor& tt, int dom, int bound, Array<double>& res, int& pos_res, const Array<int>& ncond,
										int n_cmp, Array<int>** p_cmp) const {

	// Check boundary
	if (bound!=OUTER_BC) {
		cerr << "Unknown boundary in Domain_nucleus::export_tau_boundary" << endl ;
		abort() ;
	}

	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			if (!tt.is_m_order_affected())
			  export_tau_val_domain_boundary (tt()(dom), 0, bound, res, pos_res, ncond(0)) ;
			else 
			    export_tau_val_domain_boundary (tt()(dom), tt.get_parameters()->get_m_order(), bound, res, pos_res, ncond(0)) ;
			break ;
		case 1 : {
			bool found = false ;
			// Cartesian basis
			if (tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) {
				if (n_cmp==-1) {
					export_tau_val_domain_boundary (tt(1)(dom), 0,  bound, res, pos_res, ncond(0)) ;
					export_tau_val_domain_boundary (tt(2)(dom), 0, bound,  res, pos_res, ncond(1)) ;
					export_tau_val_domain_boundary (tt(3)(dom), 0,  bound, res, pos_res, ncond(2)) ;
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if ((*p_cmp[i])(0)==1)
						export_tau_val_domain_boundary (tt(1)(dom), 0,  bound, res, pos_res, ncond(i)) ;
					if ((*p_cmp[i])(0)==2)
						export_tau_val_domain_boundary (tt(2)(dom), 0,  bound, res, pos_res, ncond(i)) ;
					if ((*p_cmp[i])(0)==3)
						export_tau_val_domain_boundary (tt(3)(dom), 0,  bound, res, pos_res, ncond(i)) ;
				}
				found = true ;
			}
			// Spherical coordinates
			if (tt.get_basis().get_basis(dom)==SPHERICAL_BASIS) {
				if (n_cmp==-1) {
					export_tau_val_domain_boundary_vr (tt(1)(dom),  bound, res, pos_res, ncond(0)) ;
					export_tau_val_domain_boundary_vt (tt(2)(dom),  bound, res, pos_res, ncond(1)) ;
					export_tau_val_domain_boundary_vp (tt(3)(dom),  bound, res, pos_res, ncond(2)) ;
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if ((*p_cmp[i])(0)==1)
						export_tau_val_domain_boundary_vr (tt(1)(dom),  bound, res, pos_res, ncond(i)) ;
					if ((*p_cmp[i])(0)==2)
						export_tau_val_domain_boundary_vt (tt(2)(dom),  bound, res, pos_res, ncond(i)) ;
					if ((*p_cmp[i])(0)==3)
						export_tau_val_domain_boundary_vp (tt(3)(dom),  bound, res, pos_res, ncond(i)) ;
				}
				found = true ;
			}
			if (!found) {
				cerr << "Unknown type of vector Domain_nucleus::export_tau_boundary" << endl ;
				abort() ;
			}
		}
			break ;
		case 2 : {
			bool found = false ;
			// Cartesian basis and symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==6)) {
				if (n_cmp==-1) {
					export_tau_val_domain_boundary (tt(1,1)(dom), 0,  bound, res, pos_res, ncond(0)) ;
					export_tau_val_domain_boundary (tt(1,2)(dom), 0,  bound, res, pos_res, ncond(1)) ;
					export_tau_val_domain_boundary (tt(1,3)(dom), 0,  bound, res, pos_res, ncond(2)) ;
					export_tau_val_domain_boundary (tt(2,2)(dom), 0,  bound, res, pos_res, ncond(3)) ;
					export_tau_val_domain_boundary (tt(2,3)(dom), 0,  bound, res, pos_res, ncond(4)) ;
					export_tau_val_domain_boundary (tt(3,3)(dom), 0, bound,  res, pos_res, ncond(5)) ;
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==1))
						export_tau_val_domain_boundary (tt(1, 1)(dom), 0, bound,  res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==2))
						export_tau_val_domain_boundary (tt(1, 2)(dom), 0, bound,  res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain_boundary (tt(1, 3)(dom), 0,  bound, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==2))
						export_tau_val_domain_boundary (tt(2, 2)(dom), 0,  bound, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain_boundary (tt(2, 3)(dom), 0,  bound, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain_boundary (tt(3, 3)(dom), 0,  bound, res, pos_res, ncond(i)) ;
				}
				found = true ;
			}
			// Cartesian basis and not symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==9)) {
				if (n_cmp==-1) {
					export_tau_val_domain_boundary (tt(1,1)(dom), 0,  bound, res, pos_res, ncond(0)) ;
					export_tau_val_domain_boundary (tt(1,2)(dom), 0,  bound, res, pos_res, ncond(1)) ;
					export_tau_val_domain_boundary (tt(1,3)(dom), 0,  bound, res, pos_res, ncond(2)) ;
					export_tau_val_domain_boundary (tt(2,1)(dom), 0,  bound, res, pos_res, ncond(3)) ;
					export_tau_val_domain_boundary (tt(2,2)(dom), 0,  bound, res, pos_res, ncond(4)) ;
					export_tau_val_domain_boundary (tt(2,3)(dom), 0, bound,  res, pos_res, ncond(5)) ;
					export_tau_val_domain_boundary (tt(3,1)(dom), 0,  bound, res, pos_res, ncond(6)) ;
					export_tau_val_domain_boundary (tt(3,2)(dom), 0,  bound, res, pos_res, ncond(7)) ;
					export_tau_val_domain_boundary (tt(3,3)(dom), 0, bound,  res, pos_res, ncond(8)) ;
					
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==1))
						export_tau_val_domain_boundary (tt(1, 1)(dom), 0, bound,  res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==2))
						export_tau_val_domain_boundary (tt(1, 2)(dom), 0, bound,  res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain_boundary (tt(1, 3)(dom), 0,  bound, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==1))
						export_tau_val_domain_boundary (tt(2, 1)(dom), 0, bound,  res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==2))
						export_tau_val_domain_boundary (tt(2, 2)(dom), 0,  bound, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain_boundary (tt(2, 3)(dom), 0,  bound, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==1))
						export_tau_val_domain_boundary (tt(3, 1)(dom), 0, bound,  res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==2))
						export_tau_val_domain_boundary (tt(3, 2)(dom), 0,  bound, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain_boundary (tt(3, 3)(dom), 0,  bound, res, pos_res, ncond(i)) ;
				}
				found = true ;
			}
			if (!found) {
				cerr << "Unknown type of 2-tensor Domain_nucleus::export_tau_boundary" << endl ;
				abort() ;
			}
		}
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_nucleus::export_tau_boundary" << endl ;
			break ;
	}
}}
