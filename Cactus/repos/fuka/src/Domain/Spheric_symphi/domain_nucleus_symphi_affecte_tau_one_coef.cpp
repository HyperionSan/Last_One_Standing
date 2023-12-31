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
#include "spheric_symphi.hpp"
#include "array.hpp"
#include "scalar.hpp"
#include "tensor.hpp"

namespace Kadath {
void Domain_nucleus_symphi::affecte_tau_one_coef_val_domain (Val_domain& so, int cc, int& conte) const {

	so.is_zero = false ;
	so.allocate_coef() ;
	*so.cf=0. ;
	Index pos_cf(nbr_coefs) ;

	bool found = false ;

	// Positions of the Galerkin basis
	Index pos_gal_t (nbr_coefs) ;
	Index pos_gal_r (nbr_coefs) ;
	Index pos_gal_rt (nbr_coefs) ;
	double fact_t, fact_r, fact_rt ;

	int kmin, kmax ;
		// Base in phi 
		int basep = (*so.get_base().bases_1d[2]) (0) ;
		switch (basep) {
					case COS_EVEN:
						kmin = 0 ;
						kmax = nbr_coefs(2)-1 ;
						break ;
					case COS_ODD:
						kmin = 0 ;
						kmax = nbr_coefs(2)-2 ;
						break ;
					case SIN_EVEN:
						kmin = 1 ;
						kmax = nbr_coefs(2)-2 ;
						break ;
					case SIN_ODD:
						kmin = 0 ;
						kmax = nbr_coefs(2)-2 ;
						break ;
					default:
						cerr << "Unknow phi basis in Domain_nucleus_symphi::affecte_tau_one_coef_val_domain" << endl ;
						abort() ;
		}



	int lquant ;
	// Loop on phi :
	for (int k=kmin ; k<=kmax ; k++) {

		int mquant ;
		switch (basep) {
					case COS_EVEN:
						mquant = 2*k ;
						break ;
					case COS_ODD:
						mquant = 2*k+1 ;
						break ;
					case SIN_EVEN:
						mquant = 2*k ;
						break ;
					case SIN_ODD:
						mquant = 2*k+1 ;
						break ;
					default:
						cerr << "Unknow phi basis in Domain_nucleus_symphi::affecte_tau_one_coef_val_domain" << endl ;
						abort() ;
		}



		pos_cf.set(2) = k ;
			// Loop on theta
			int baset = (*so.get_base().bases_1d[1]) (k) ;
			for (int j=0 ; j<nbr_coefs(1) ; j++) {	
				int baser = (*so.get_base().bases_1d[0]) (j, k) ;
				pos_cf.set(1) = j ;
				// Loop on r :
				for (int i=0 ; i<nbr_coefs(0) ; i++) {
					pos_cf.set(0) = i ;
					switch (baset) {
						case COS_EVEN :
							lquant = 2*j ;
							// No galerkin :
							if ((mquant==0) && (lquant==0))  {
								if (conte==cc)  {
								  found = true ;
								  so.cf->set(pos_cf) = 1. ;
								  }
								  conte ++ ;
							}
							else if (mquant==0) {
							  if (i!=0) {
								if (conte==cc) {
								found = true ;
								// Galerkin base in r only
								pos_gal_r = pos_cf ;
								pos_gal_r.set(0) = 0 ;
								switch (baser) {
									case CHEB_EVEN :
									  fact_r = - pow(-1, i) ;
									  break ;
									case LEG_EVEN : {
									  fact_r = -1. ;
									  for (int t=0 ; t<i ; t++)
									    fact_r *= -double(2*t+1)/double(2*t+2) ;
									  }
									  break ;
									default :
									  cerr << "Strange base in Domain_nucleus_symphi::affecte_one_coef_val_domain" << endl ;
									  abort()  ;
								}
								so.cf->set(pos_cf) = 1 ;
								so.cf->set(pos_gal_r) += fact_r ;
								}
								conte ++ ;
							  }
							}
							else if ((j!=0) && (i!=0)) {
							    if (conte==cc) {
							    found = true ;
							    // Need to use two_dimensional Galerkin basis (aouch !)
							    pos_gal_r = pos_cf ;
							    pos_gal_r.set(0) = 0 ;
							    pos_gal_t = pos_cf ;
							    pos_gal_t.set(1) = 0 ;
							    pos_gal_rt = pos_cf ;
							    pos_gal_rt.set(0) = 0 ;
							    pos_gal_rt.set(1) = 0 ;
							    switch (baser) {
							      case CHEB_EVEN :
								fact_r = -pow(-1, i) ;
								fact_t = -1. ;
								fact_rt = pow(-1, i) ;
								break ;
							      case LEG_EVEN : {
								double l0 = 1 ;
								 for (int t=0 ; t<i ; t++)
								  l0 *= -double(2*t+1)/double(2*t+2) ;
								fact_r = - l0 ;
								fact_t = -1. ;
								fact_rt = l0 ;
								}
								break ;
							    default :
									  cerr << "Strange base in Domain_nucleus_symphi::affecte_one_coef_val_domain" << endl ;
									  abort()  ;
								}
							      so.cf->set(pos_cf) = 1. ;
							      so.cf->set(pos_gal_r) = fact_r ;
							      so.cf->set(pos_gal_t) = fact_t ;
							      so.cf->set(pos_gal_rt) = fact_rt ;
							      }
							      conte ++ ;
							    }
						break ;
					case COS_ODD:
						lquant = 2*j+1 ;
						if ((j!=nbr_coefs(1)-1) && (i!=nbr_coefs(0)-1))  {
							      if ((mquant==0) && (lquant<=1)) {
								  if (conte==cc) {
								      found = true ;
								      so.cf->set(pos_cf) = 1. ;
								      }
								  conte++  ;
							      }
							      else {
								if ((mquant==0) && (i!=0)) {
								  if (conte==cc) {
								    found = true ;
								    pos_gal_r = pos_cf ;
								    pos_gal_r.set(0) = 0 ;
								    switch (baser) {
								      case CHEB_ODD :
									fact_r = - (2*i+1) * pow(-1, i) ;
									break ;
								      case LEG_ODD : {
									fact_r = -1. ;
									for (int t=0 ; t<i ; t++)
									  fact_r *= -double(2*t+3)/double(2*t+2) ;
									}
									break ;
								      default :
									cerr << "Strange base in Domain_nucleus_symphi::affecte_one_coef_val_domain" << endl ;
									abort()  ;
								      }

								  so.cf->set(pos_cf) = 1. ;
								  so.cf->set(pos_gal_r) = fact_r ;
								}
							       conte ++ ;
							}
							else if ((j!=0) && (i!=0)) {
							    if (conte==cc) {
								found = true ;
								 // Need to use two_dimensional Galerkin basis (aouch !)
								pos_gal_r = pos_cf ;
								pos_gal_r.set(0) = 0 ;
								pos_gal_t = pos_cf ;
								pos_gal_t.set(1) = 0 ;
								pos_gal_rt = pos_cf ;
								pos_gal_rt.set(0) = 0 ;
								pos_gal_rt.set(1) = 0 ;
								switch (baser) {
								  case CHEB_ODD :
								    fact_r = -pow(-1, i)*(2*i+1) ;
								    fact_t = -1. ;
								    fact_rt = pow(-1, i)*(2*i+1) ;
								    break ;
								  case LEG_ODD : {
								    double l0 = 1 ;
								    for (int t=0 ; t<i ; t++)
								      l0 *= -double(2*t+3)/double(2*t+2) ;
								    fact_r = - l0 ;
								    fact_t = -1. ;
								    fact_rt = l0 ;
								    }
								    break ;
								default :
									  cerr << "Strange base in Domain_nucleus_symphi::affecte_one_coef_val_domain" << endl ;
									  abort()  ;
								}  
							      so.cf->set(pos_cf) = 1. ;
							      so.cf->set(pos_gal_r) = fact_r ;
							      so.cf->set(pos_gal_t) = fact_t ;
							      so.cf->set(pos_gal_rt) = fact_rt ;
							      }
							      conte ++ ;
							}
						    }
						}
						break ;
					case SIN_EVEN:
						lquant = 2*j ;
						if ((j!=0) && (j!=nbr_coefs(1)-1)) { 
						if ((mquant<=1) && (lquant==0)) {
						    if (conte==cc) {
								      found = true ;
								      so.cf->set(pos_cf) = 1. ;
								      }
						    conte ++ ;
						}
						else {
						if ((mquant<=1) && (i!=0)) {
							// Galerkin base in r only
							 if (conte==cc) {
								found = true ;
							pos_gal_r = pos_cf ;
							pos_gal_r.set(0) = 0 ;
							switch (baser) {
								case CHEB_EVEN :
								  fact_r = - pow(-1, i) ;
								  break ;
								case LEG_EVEN : {
								  fact_r = -1. ;
								  for (int t=0 ; t<i ; t++)
								  fact_r *= -double(2*t+1)/double(2*t+2) ;
								  }
								  break ;
								default :
								  cerr << "Strange base in Domain_nucleus_symphi::affecte_tau_val_domain" << endl ;
								  abort()  ;
								}	  
							    so.cf->set(pos_cf) = 1. ;
								so.cf->set(pos_gal_r) = fact_r ;
								}
							    conte ++ ;
							}
						else  {
						//Double Galerkin
						if ((j!=1) && (i!=0)) {

								 if (conte==cc) {
								    found = true ;
								 // Need to use two_dimensional Galerkin basis (aouch !)
							    pos_gal_r = pos_cf ;
							    pos_gal_r.set(0) = 0 ;
							    pos_gal_t = pos_cf ;
							    pos_gal_t.set(1) = 1 ;
							    pos_gal_rt = pos_cf ;
							    pos_gal_rt.set(0) = 0 ;
							    pos_gal_rt.set(1) = 1 ;
							     switch (baser) {
							      case CHEB_EVEN :
								fact_r = -pow(-1, i) ;
								fact_t = -j ;
								fact_rt = pow(-1, i)*j ;
								break ;
							      case LEG_EVEN : {
								double l0 = 1 ;
								 for (int t=0 ; t<i ; t++)
								  l0 *= -double(2*t+1)/double(2*t+2) ;
								fact_r = - l0 ;
								fact_t = -j ;
								fact_rt = l0*j ;
								}
								break ;
							    default :
									  cerr << "Strange base in Domain_nucleus_symphi::affecte_tau_val_domain" << endl ;
									  abort()  ;
								}  
							      so.cf->set(pos_cf) = 1 ;
							      so.cf->set(pos_gal_r) = fact_r ;
							      so.cf->set(pos_gal_t) = fact_t ;
							      so.cf->set(pos_gal_rt) = fact_rt ;
								}
							      conte ++ ;
							}
						      }
						}
						}
						break ;
					case SIN_ODD:
						lquant = 2*j+1 ;
						if ((j!=nbr_coefs(1)-1) && (i!=nbr_coefs(0)-1))  {
							      if ((mquant<=1) && (lquant<=1)) {
								  if (conte==cc) {
								      found = true ;
								      so.cf->set(pos_cf) = 1. ;
								      }
								  conte++  ;
							      }
							      else {
								if ((mquant<=1) && (i!=0)) {
								  if (conte==cc) {
								    found = true ;
								    pos_gal_r = pos_cf ;
								    pos_gal_r.set(0) = 0 ;
								    switch (baser) {
								      case CHEB_ODD :
									fact_r = - (2*i+1) * pow(-1, i) ;
									break ;
								      case LEG_ODD : {
									fact_r = -1. ;
									for (int t=0 ; t<i ; t++)
									  fact_r *= -double(2*t+3)/double(2*t+2) ;
									}
									break ;
								      default :
									cerr << "Strange base in Domain_nucleus_symphi::affecte_one_coef_val_domain" << endl ;
									abort()  ;
								      }

								  so.cf->set(pos_cf) = 1. ;
								  so.cf->set(pos_gal_r) = fact_r ;
								}
							       conte ++ ;
							}
							else if ((j!=0) && (i!=0)) {
							    if (conte==cc) {
								found = true ;
								 // Need to use two_dimensional Galerkin basis (aouch !)
								pos_gal_r = pos_cf ;
								pos_gal_r.set(0) = 0 ;
								pos_gal_t = pos_cf ;
								pos_gal_t.set(1) = 0 ;
								pos_gal_rt = pos_cf ;
								pos_gal_rt.set(0) = 0 ;
								pos_gal_rt.set(1) = 0 ;
								switch (baser) {
								  case CHEB_ODD :
								    fact_r = -pow(-1, i)*(2*i+1) ;
								    fact_t = -(2*j+1) ;
								    fact_rt = pow(-1, i)*(2*i+1)*(2*j+1) ;
								    break ;
								  case LEG_ODD : {
								    double l0 = 1 ;
								    for (int t=0 ; t<i ; t++)
								      l0 *= -double(2*t+3)/double(2*t+2) ;
								    fact_r = - l0 ;
								    fact_t = -(2*j+1) ;
								    fact_rt = l0*(2*j+1) ;
								    }
								    break ;
								default :
									  cerr << "Strange base in Domain_nucleus_symphi::affecte_one_coef_val_domain" << endl ;
									  abort()  ;
								}  
							      so.cf->set(pos_cf) = 1. ;
							      so.cf->set(pos_gal_r) = fact_r ;
							      so.cf->set(pos_gal_t) = fact_t ;
							      so.cf->set(pos_gal_rt) = fact_rt ;
							      }
							      conte ++ ;
							}
						    }
						}
						break ;


					default:
						cerr << "Unknow theta basis in Domain_nucleus_symphi::affecte_coef_val_domain" << endl ;
						abort() ;
					}
				}
			}
	}
	// If not found put to zero :
	if (!found)
		so.set_zero() ;
}

void Domain_nucleus_symphi::affecte_tau_one_coef (Tensor& tt, int dom, int cc, int& pos_cf) const {

	// Check right domain
	assert (tt.get_space().get_domain(dom)==this) ;

	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			affecte_tau_one_coef_val_domain (tt.set().set_domain(dom), cc, pos_cf) ;
			break ;
		case 1 : {
			bool found = false ;
			// Cartesian basis
			if (tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) {
				affecte_tau_one_coef_val_domain (tt.set(1).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(3).set_domain(dom), cc, pos_cf) ;
				found = true ;
			}
			
			if (!found) {
				cerr << "Unknown type of vector Domain_nucleus_symphi::affecte_tau_one_coef" << endl ;
				abort() ;
			}
		}
			break ;
		case 2 : {
			bool found = false ;
			// Cartesian basis and symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==6)) {
				affecte_tau_one_coef_val_domain (tt.set(1,1).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(1,2).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(1,3).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2,2).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2,3).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(3,3).set_domain(dom), cc, pos_cf) ;
				found = true ;
			}
			// Cartesian basis and not symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==9)) {
				affecte_tau_one_coef_val_domain (tt.set(1,1).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(1,2).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(1,3).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2,1).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2,2).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2,3).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(3,1).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(3,2).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(3,3).set_domain(dom), cc, pos_cf) ;
				found = true ;
			}
			if (!found) {
				cerr << "Unknown type of 2-tensor Domain_nucleus_symphi::affecte_tau_one_coef" << endl ;
				abort() ;
			}
		}
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_nucleus_symphi::affecte_tau" << endl ;
			break ;
	}
}}
