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
int Domain_nucleus::nbr_conditions_val_domain_vr (const Val_domain& so, int order) const {
	
	int res = 0 ;
	
	Index pos (nbr_coefs) ;
	do {
		bool indic = true ;
		// True coef in phi ?
		if ((pos(2)==1) || (pos(2)==nbr_coefs(2)-1))
			indic = false ;
		// Get base in theta :
		int baset = (*so.get_base().bases_1d[1]) (pos(2)) ;
		switch (baset) {
					case COS_EVEN:
						if ((pos(1)==0) && (pos(2)!=0))
							indic = false ;
						break ;
					case SIN_ODD:
						if (pos(1)==nbr_coefs(1)-1)
							indic = false ;
						break ;
					default:
						cerr << "Unknow theta basis in Domain_nucleus::nbr_conditions_val_domain_vr" << endl ;
						abort() ;
		}
	
		int max = 0 ;
		if (indic) {
			// Base in r :
			int baser = (*so.get_base().bases_1d[0]) (pos(1), pos(2)) ;
			switch (baser) {
				case CHEB_EVEN :
					if (pos(0)==0)
					  indic = false ;
					max = nbr_coefs(0) ;
					break ;
				case LEG_EVEN :
					if  (pos(0)==0)
					  indic = false ;
					max = nbr_coefs(0) ;
					break ;
				case CHEB_ODD :
					if (pos(0)==nbr_coefs(0)-1)
						indic = false ;
					if ((pos(2)==0) && (pos(1)>1) && (pos(0)==0))
					  indic = false ;
					if ((pos(2)>0) && (pos(0)==0))
					  indic = false ;
					max = nbr_coefs(0)-1 ;
					break ;
				case LEG_ODD :
					if (pos(0)==nbr_coefs(0)-1)
						indic = false ;
					if ((pos(2)==0) && (pos(1)>1) && (pos(0)==0))
					  indic = false ;
					if ((pos(2)>0) && (pos(0)==0))
					  indic = false ;
					max = nbr_coefs(0)-1 ;
					break ;
				default:
					cerr << "Unknow radial basis in Domain_nucleus::nbr_conditions_val_domain_vr" << endl ;
					abort() ;
			}
		}
		// Order with respect to r :
		int lim = 0 ;
		switch (order) {
		  case 2 :
		      lim = max-2 ;
		      break ;
		  case 0 :
		      lim = max-1 ;
		      break ;
		  default :
		      cerr << "Unknown case in Domain_nucleus_nbr_conditions" << endl ;
		      abort() ;
		}
		if (pos(0)>lim)
			indic = false ;
		if (indic)
			res ++ ;
	}
	while (pos.inc()) ;

	return res ;
}

int Domain_nucleus::nbr_conditions_val_domain_vt (const Val_domain& so, int order) const {
	
	int res = 0 ;
	
	Index pos (nbr_coefs) ;
	do {
		bool indic = true ;
		// True coef in phi ?
		if ((pos(2)==1) || (pos(2)==nbr_coefs(2)-1))
			indic = false ;
		// Get base in theta :
		int baset = (*so.get_base().bases_1d[1]) (pos(2)) ;
		switch (baset) {
					case SIN_EVEN:
						if ((pos(1)==nbr_coefs(1)-1) || (pos(1)==0))
						  indic= false ;
						break ;
					case COS_ODD:
						if (pos(1)==nbr_coefs(1)-1)
							indic = false ;
						if ((pos(1)==0) && (pos(2)>3))
							indic = false ;
						break ;
					default:
						cerr << "Unknow theta basis in Domain_nucleus::nbr_conditions_val_domain_vt" << endl ;
						abort() ;
		}
	
		int max = 0 ;
		if (indic) {
			// Base in r :
			int baser = (*so.get_base().bases_1d[0]) (pos(1), pos(2)) ;
			switch (baser) {
				case CHEB_EVEN :
					if (((pos(2)==2) || (pos(2)==3)) && (pos(1)>0) && (pos(0)==0))
					  indic = false ;
					if ((pos(2)>3) && (pos(0)==0))
					  indic = false ;
					max = nbr_coefs(0) ;
					break ;
				case LEG_EVEN :
					if (((pos(2)==2) || (pos(2)==3)) && (pos(1)>0) && (pos(0)==0))
					  indic = false ;
					if ((pos(2)>3) && (pos(0)==0))
					  indic = false ;
					max = nbr_coefs(0) ;
					break ;
				case CHEB_ODD :
					if (pos(0)==nbr_coefs(0)-1)
						indic = false ;
					if (pos(0)==0)
					  indic = false ;
					max = nbr_coefs(0)-1 ;
					break ;
				case LEG_ODD :
					if (pos(0)==nbr_coefs(0)-1)
						indic = false ;
					if (pos(0)==0)
					  indic = false ;
					max = nbr_coefs(0)-1 ;
					break ;
				default:
					cerr << "Unknow radial basis in Domain_nucleus::nbr_conditions_val_domain_vt" << endl ;
					abort() ;
			}
		}
		// Order with respect to r :
		int lim = 0 ;
		switch (order) {
		  case 2 :
		      lim = max-2 ;
		      break ;
		  case 0 :
		      lim = max-1 ;
		      break ;
		  default :
		      cerr << "Unknown case in Domain_nucleus_nbr_conditions" << endl ;
		      abort() ;
		}
		if (pos(0)>lim)
			indic = false ;
		if (indic)
			res ++ ;
	}
	while (pos.inc()) ;

	return res ;
}

int Domain_nucleus::nbr_conditions_val_domain_vp (const Val_domain& so, int order) const {
	
	int res = 0 ;
	
	Index pos (nbr_coefs) ;
	do {
		bool indic = true ;
		// True coef in phi ?
		if ((pos(2)==1) || (pos(2)==nbr_coefs(2)-1))
			indic = false ;
		// Get base in theta :
		int baset = (*so.get_base().bases_1d[1]) (pos(2)) ;
		switch (baset) {
					case COS_EVEN:
					  if ((pos(2)>3) && (pos(1)==0))
						  indic = false ;
						break ;
					case SIN_ODD:
						if (pos(1)==nbr_coefs(1)-1)
							indic = false ;
						break ;
					default:
						cerr << "Unknow theta basis in Domain_nucleus::nbr_conditions_val_domain_vp" << endl ;
						abort() ;
		}
	
		int max = 0 ;
		if (indic) {
			// Base in r :
			int baser = (*so.get_base().bases_1d[0]) (pos(1), pos(2)) ;
			switch (baser) {
				case CHEB_EVEN :
					if (((pos(2)==2) || (pos(2)==3)) && (pos(1)>0) && (pos(0)==0))
					  indic = false ;
					if ((pos(2)>3) && (pos(0)==0))
					  indic = false ;
					max = nbr_coefs(0) ;
					break ;
				case LEG_EVEN :
					if (((pos(2)==2) || (pos(2)==3)) && (pos(1)>0) && (pos(0)==0))
					  indic = false ;
					if ((pos(2)>3) && (pos(0)==0))
					  indic = false ;
					max = nbr_coefs(0) ;
					break ;
				case CHEB_ODD :
					if (pos(0)==nbr_coefs(0)-1)
						indic = false ;
					if ((pos(2)==0) && (pos(1)>0) && (pos(0)==0))
					  indic = false ;
					if (((pos(2)==4) || (pos(2)==5)) && (pos(1)>0) && (pos(0)==0))
					  indic = false ;
					if ((pos(2)>5) && (pos(0)==0))
					  indic = false ;
					max = nbr_coefs(0)-1 ;
					break ;
				case LEG_ODD :
					if (pos(0)==nbr_coefs(0)-1)
						indic = false ;
					if ((pos(2)==0) && (pos(1)>0) && (pos(0)==0))
					  indic = false ;
					if (((pos(2)==4) || (pos(2)==5)) && (pos(1)>0) && (pos(0)==0))
					  indic = false ;
					if ((pos(2)>5) && (pos(0)==0))
					  indic = false ;
					max = nbr_coefs(0)-1 ;
					break ;
				default:
					cerr << "Unknow radial basis in Domain_nucleus::nbr_conditions_val_domain_vp" << endl ;
					abort() ;
			}
		}
		// Order with respect to r :
		int lim = 0 ;
		switch (order) {
		  case 2 :
		      lim = max-2 ;
		      break ;
		  case 0 :
		      lim = max-1 ;
		      break ;
		  default :
		      cerr << "Unknown case in Domain_nucleus_nbr_conditions" << endl ;
		      abort() ;
		}
		if (pos(0)>lim)
			indic = false ;
		if (indic)
			res ++ ;
	}
	while (pos.inc()) ;

	return res ;
}


int Domain_nucleus::nbr_conditions_val_domain (const Val_domain& so, int mlim, int llim, int order) const {
	
	int res = 0 ;
	int kmin = 2*mlim + 2 ;

	Index pos (nbr_coefs) ;
	do {
		bool indic = true ;
		// True coef in phi ?
		if ((pos(2)==1) || (pos(2)==nbr_coefs(2)-1))
			indic = false ;
		// Get base in theta :
		int baset = (*so.get_base().bases_1d[1]) (pos(2)) ;
		int lquant ;
		switch (baset) {
					case COS_EVEN:
						if ((pos(1)==0) && (pos(2)>=kmin))
							indic = false ;
						lquant = 2*pos(1) ;
						break ;
					case COS_ODD:
						if ((pos(1)==nbr_coefs(1)-1) || ((pos(1)==0) && (pos(2)>=kmin)))
							indic = false ;
						lquant = 2*pos(1)+1 ;
						break ;
					case SIN_EVEN:
						if (((pos(1)==1) && (pos(2)>=kmin+2)) || (pos(1)==0) || (pos(1)==nbr_coefs(1)-1)) 
							indic = false  ;
						lquant = 2*pos(1) ;
						break ;
					case SIN_ODD:
						if (((pos(1)==0) && (pos(2)>=kmin+2)) || (pos(1)==nbr_coefs(1)-1))
							indic = false ;
						lquant = 2*pos(1)+1 ;
						break ;
					default:
						cerr << "Unknow theta basis in Domain_nucleus::nbr_unknowns_val_domain" << endl ;
						abort() ;
		}

		int max = 0;
		if (indic) {
			// Base in r :
			int baser = (*so.get_base().bases_1d[0]) (pos(1), pos(2)) ;
			
			switch (baser) {
				case CHEB_EVEN :
					if ((lquant>llim) && (pos(0)==0))
						indic = false ;
					max = nbr_coefs(0) ;
					break ;
				case LEG_EVEN :
					if ((lquant>llim) && (pos(0)==0))
						indic = false ;
					max = nbr_coefs(0) ;
					break ;
				case CHEB_ODD :
					if ((lquant>llim+1) && (pos(0)==0))
						indic = false ;
					max = nbr_coefs(0)-1 ;
					break ;
				case LEG_ODD :
					if ((lquant>llim+1) && (pos(0)==0))
						indic = false ;
					max = nbr_coefs(0)-1 ;
					break ;
				default:
					cerr << "Unknow radial basis in Domain_nucleus::nbr_unknowns_val_domain" << endl ;
					abort() ;
			}
		}

		// Order with respect to r :
		int lim = 0 ;
		switch (order) {
		  case 2 :
		      lim = max-2 ;
		      break ;
		  case 1 :
		  {
		      if ((pos(1)==0) && (pos(2)==0))
			lim = max-2 ;
		      else
			lim = max -1 ;}
		      break ;
		  case 0 :
		      lim = max-1 ;
		      break ;
		  default :
		      cerr << "Unknown case in Domain_nucleus_nbr_conditions" << endl ;
		      abort() ;
		}
		
		
		if (pos(0)>lim)
			indic = false ;

		if (indic)
			res ++ ;
	}
	while (pos.inc()) ;
	return res ;
}

Array<int> Domain_nucleus::nbr_conditions (const Tensor& tt, int dom, int order, int n_cmp, Array<int>** p_cmp) const {

	int size = (n_cmp==-1) ? tt.get_n_comp() : n_cmp ;
	Array<int> res (size) ;
	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			if (!tt.is_m_order_affected())
			  res.set(0) = nbr_conditions_val_domain (tt()(dom), 0, 0, order) ;
			else
			  res.set(0) = nbr_conditions_val_domain (tt()(dom), tt.get_parameters()->get_m_order(), 0, order) ;
			break ;
		case 1 : {
			bool found = false ;
			// Cartesian basis
			if (tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) {
				if (n_cmp==-1) {
					res.set(0) = nbr_conditions_val_domain (tt(1)(dom), 0, 0, order) ;
					res.set(1) = nbr_conditions_val_domain (tt(2)(dom), 0, 0, order) ;
					res.set(2) = nbr_conditions_val_domain (tt(3)(dom), 0, 0, order) ;
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if ((*p_cmp[i])(0)==1)
						res.set(i) = nbr_conditions_val_domain (tt(1)(dom), 0, 0, order) ;
					if ((*p_cmp[i])(0)==2)
						res.set(i) = nbr_conditions_val_domain (tt(2)(dom), 0, 0, order) ;
					if ((*p_cmp[i])(0)==3)
						res.set(i) = nbr_conditions_val_domain (tt(3)(dom), 0, 0, order) ;
				}
				found = true ;
			}
			// Spherical coordinates
			if (tt.get_basis().get_basis(dom)==SPHERICAL_BASIS) {
				if (n_cmp==-1) {
					res.set(0) = nbr_conditions_val_domain_vr (tt(1)(dom), order) ;
					res.set(1) = nbr_conditions_val_domain_vt (tt(2)(dom), order) ;
					res.set(2) = nbr_conditions_val_domain_vp (tt(3)(dom), order) ;
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if ((*p_cmp[i])(0)==1)
						res.set(i) = nbr_conditions_val_domain_vr (tt(1)(dom), order) ;
					if ((*p_cmp[i])(0)==2)
						res.set(i) = nbr_conditions_val_domain_vt (tt(2)(dom), order) ;
					if ((*p_cmp[i])(0)==3)
						res.set(i) = nbr_conditions_val_domain_vp (tt(3)(dom), order) ;
				}
				found = true ;
			}
			if (!found) {
				cerr << "Unknown type of vector Domain_nucleus::nbr_conditions" << endl ;
				abort() ;
			}
		}
		      break ;
		case 2 : {
			bool found = false ;
			// Cartesian basis and symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==6)) {
				if (n_cmp==-1) {
					res.set(0) = nbr_conditions_val_domain (tt(1,1)(dom), 0, 0, order) ;
					res.set(1) = nbr_conditions_val_domain (tt(1,2)(dom), 0, 0, order) ;
					res.set(2) = nbr_conditions_val_domain (tt(1,3)(dom), 0, 0, order) ;
					res.set(3) = nbr_conditions_val_domain (tt(2,2)(dom), 0, 0, order) ;
					res.set(4) = nbr_conditions_val_domain (tt(2,3)(dom), 0, 0, order) ;
					res.set(5) = nbr_conditions_val_domain (tt(3,3)(dom), 0, 0, order) ;
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==1))
						res.set(i) = nbr_conditions_val_domain (tt(1, 1)(dom), 0, 0, order) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==2))
						res.set(i) = nbr_conditions_val_domain (tt(1, 2)(dom), 0, 0, order) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==3))
						res.set(i) = nbr_conditions_val_domain (tt(1, 3)(dom), 0, 0, order) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==2))
						res.set(i) = nbr_conditions_val_domain (tt(2, 2)(dom), 0, 0, order) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==3))
						res.set(i) = nbr_conditions_val_domain (tt(2, 3)(dom), 0, 0, order) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==3))
						res.set(i) = nbr_conditions_val_domain (tt(3, 3)(dom), 0, 0, order) ;
				}
				found = true ;
			}
			// Cartesian basis and not symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==9)) {
				if (n_cmp==-1) {
					res.set(0) = nbr_conditions_val_domain (tt(1,1)(dom), 0, 0, order) ;
					res.set(1) = nbr_conditions_val_domain (tt(1,2)(dom), 0, 0, order) ;
					res.set(2) = nbr_conditions_val_domain (tt(1,3)(dom), 0, 0, order) ;
					res.set(3) = nbr_conditions_val_domain (tt(2,1)(dom), 0, 0, order) ;
					res.set(4) = nbr_conditions_val_domain (tt(2,2)(dom), 0, 0, order) ;
					res.set(5) = nbr_conditions_val_domain (tt(2,3)(dom), 0, 0, order) ;
					res.set(6) = nbr_conditions_val_domain (tt(3,1)(dom), 0, 0, order) ;
					res.set(7) = nbr_conditions_val_domain (tt(3,2)(dom), 0, 0, order) ;
					res.set(8) = nbr_conditions_val_domain (tt(3,3)(dom), 0, 0, order) ;
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==1))
						res.set(i) = nbr_conditions_val_domain (tt(1, 1)(dom), 0, 0, order) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==2))
						res.set(i) = nbr_conditions_val_domain (tt(1, 2)(dom), 0, 0, order) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==3))
						res.set(i) = nbr_conditions_val_domain (tt(1, 3)(dom), 0, 0, order) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==1))
						res.set(i) = nbr_conditions_val_domain (tt(2, 1)(dom), 0, 0, order) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==2))
						res.set(i) = nbr_conditions_val_domain (tt(2, 2)(dom), 0, 0, order) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==3))
						res.set(i) = nbr_conditions_val_domain (tt(2, 3)(dom), 0, 0, order) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==1))
						res.set(i) = nbr_conditions_val_domain (tt(3, 1)(dom), 0, 0, order) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==2))
						res.set(i) = nbr_conditions_val_domain (tt(3, 2)(dom), 0, 0, order) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==3))
						res.set(i) = nbr_conditions_val_domain (tt(3, 3)(dom), 0, 0, order) ;
				}
				found = true ;
			}
			if (!found) {
				cerr << "Unknown type of 2-tensor Domain_nucleus::nbr_conditions" << endl ;
				abort() ;
			}
		}
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_nucleus::nbr_conditions" << endl ;
			break ;
	}
	return res ;
}}
