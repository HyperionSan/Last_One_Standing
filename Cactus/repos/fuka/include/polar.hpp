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

#ifndef __POLAR_HPP_
#define __POLAR_HPP_

#include "space.hpp"
namespace Kadath {


/**
* Class for a 2-dimensional spherical domain containing the origin and a symetry with respect to the plane \f$ z=0 \f$.
* It is intended to deal with systems with a symmetry with respect to \f$\varphi\f$.
* \li centered on the point \c center \f$(X_c, Z_c)\f$
* \li The numerical coordinates are :
*
* \f$ 0 \leq x \leq 1 \f$
*
* \f$ 0 \leq \theta^\star \leq \pi/2 \f$

* \li Standard spherical coordinates :
*
* \f$ \rho = \alpha x \f$
*
* \f$ \theta = \theta^\star \f$
*
* \ingroup domain
*/
class Domain_polar_nucleus : public Domain {

 private:
  double alpha ; ///< Relates the numerical to the physical radii.
  Point center ; ///< Absolute coordinates of the center.
 
 public: 
  /**
  * Standard constructor :
  * @param num : number of the domain (used by the \c Space).
  * @param ttype : Chebyshev or Legendre type of spectral expansion.
  * @param radius : radius of the nucleus.
  * @param cr : center of the spherical coordinates.
  * @param nbr : number of points in each dimension.
  */
  Domain_polar_nucleus (int nim, int ttype, double radius, const Point& cr, const Dim_array& nbr) ;
  Domain_polar_nucleus (const Domain_polar_nucleus& so) ; ///< Copy constructor.
  /**
  * Constructor from a file
  * @param num : number of the domain (used by the \c Space).
  * @param ff: file containd the domain, generated by the save function.
  */
  Domain_polar_nucleus (int num, FILE* ff) ;

  virtual ~Domain_polar_nucleus() ;
  virtual void save (FILE*) const ;

  private:    
    virtual void do_absol ()  const ;
    virtual void do_radius () const ;
    virtual void do_cart ()  const ; 
    virtual void do_cart_surr ()  const ; 

  private:
     virtual void set_cheb_base(Base_spectral&) const ;       
     virtual void set_legendre_base(Base_spectral&) const ;
     virtual void set_anti_cheb_base(Base_spectral&) const ;       
     virtual void set_anti_legendre_base(Base_spectral&) const ;
     virtual void set_cheb_base_with_m(Base_spectral&, int m) const ;       
     virtual void set_legendre_base_with_m(Base_spectral&, int m) const ;
     virtual void set_anti_cheb_base_with_m(Base_spectral&, int m) const ;       
     virtual void set_anti_legendre_base_with_m(Base_spectral&, int m) const ;
     virtual void do_coloc () ;
     virtual int give_place_var (char*) const ;

  public:
     virtual Point get_center () const {return center ;} ;
     virtual bool is_in(const Point&xx, double prec=1e-13) const ;
     virtual const Point absol_to_num(const Point&) const;
     virtual void do_der_abs_from_der_var(Val_domain** der_var, Val_domain** der_abs) const ;
     virtual Base_spectral mult (const Base_spectral&, const Base_spectral&) const ;

  public:
    
     virtual Val_domain mult_cos_theta (const Val_domain&) const ;
     virtual Val_domain mult_sin_theta (const Val_domain&) const ;
     virtual Val_domain div_sin_theta (const Val_domain&) const ;
     virtual Val_domain div_x (const Val_domain&) const ;     
     virtual Val_domain mult_r (const Val_domain&) const ;
     virtual Val_domain div_r (const Val_domain&) const ;
     virtual Val_domain laplacian (const Val_domain&, int) const ;    
     virtual Val_domain laplacian2 (const Val_domain&, int) const ;
     virtual Val_domain der_r (const Val_domain&) const ; 
     virtual Val_domain dt (const Val_domain&) const ; 
     virtual Val_domain srdr (const Val_domain&) const ;
     virtual double integrale (const Val_domain&) const ;
     virtual double integ_volume (const Val_domain&) const ;

   
     virtual double val_boundary (int, const Val_domain&, const Index&) const ;
     virtual void find_other_dom (int, int, int&, int&) const ;
     virtual Val_domain der_normal (const Val_domain&, int) const ;

     virtual int nbr_unknowns (const Tensor&, int) const ;
	/**
	* Computes the number of true unknowns of a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions to determine the precise number of degrees of freedom.
	* @param so : the field.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param llim: limit for the regularity (quantum number wrt \f$\theta\f$).
	* @returns the number of true unknowns.
	*/
     int nbr_unknowns_val_domain (const Val_domain& so, int mquant, int llim) const ;
     virtual Array<int> nbr_conditions (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given tensorial equation in the bulk.
	* It takes into account the various Galerkin basis used.
	* @param so : the residual of the equation.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param llim : limit for the regularity (quantum number wrt \f$\theta\f$).
	* @param order : order of the equation (i.e. 2 for a Laplacian for instance)
	* @returns the number of true unknowns.
	*/
     int nbr_conditions_val_domain (const Val_domain& so, int mquant, int llim, int order) const ;
     virtual Array<int> nbr_conditions_boundary (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given equation on a boundary.
	* It takes into account the various Galerkin basis used.
	* It is used for implementing boundary conditions and matching ones.
	* @param eq : the residual of the equation.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @returns the number of true conditions.
	*/
     int nbr_conditions_val_domain_boundary (const Val_domain& eq, int mquant) const ;
     virtual void export_tau (const Tensor&, int, int, Array<double>&, int&, const Array<int>&,int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports a residual equation in the bulk.
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param llim : limit for the regularity (quantum number wrt \f$\theta\f$).
	* @param order : describes the order of the equation (2 for a Laplacian for instance).
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond :  the corresponding number of equations. It is used when the equation is null.
	*/
     void export_tau_val_domain (const Val_domain& eq, int mquant, int llim, int order, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void export_tau_boundary (const Tensor&, int, int, Array<double>&, int&, const Array<int>&, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports all the residual equations corresponding to a tensorial one on a given boundary
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param bound : the boundary at which the equation is enforced.
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond :  the corresponding number of equations. It is used when the residual is null.
	*/
     void export_tau_val_domain_boundary (const Val_domain& eq, int mquant, int bound, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void affecte_tau (Tensor&, int, const Array<double>&, int&) const ;
	/**
	* Affects some coefficients to a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the field to be affected.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param llim : limit for the regularity (quantum number wrt \f$\theta\f$).
	* @param cf : \c Array of the coefficients used.
	* @param pos_cf : current position in the array of coefficients.
	*/
     void affecte_tau_val_domain (Val_domain& so, int mquant, int llim, const Array<double>& cf, int& pos_cf) const ;
     virtual void affecte_tau_one_coef (Tensor&, int, int, int&) const ;
	/**
	* Sets at most one coefficient of a \c Val_domain to 1.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the \c Val_domain to be affected. It is set to zero if cc does not corresponds to another field.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param llim : limit for the regularity (quantum number wrt \f$\theta\f$).
	* @param cc : location, in the overall system, of the coefficient to be set to 1.
	* @param pos_cf : current position.
	*/
      void affecte_tau_one_coef_val_domain (Val_domain& so, int mquant, int llim, int cc, int& pos_cf) const ;

public:
     virtual ostream& print (ostream& o) const ;
} ;


/**
* Class for a 2-dimensional spherical shell and a symmetry with respect to the plane \f$ z=0 \f$.
* It is intended to deal with systems with a symmetry with respect to \f$\varphi\f$.
* \li 2 dimensions.
* \li centered on the point \c center \f$(X_c, Z_c)\f$
* \li The numerical coordinates are :
*
* \f$ -1 \leq x \leq 1 \f$
*
* \f$ 0 \leq \theta^\star \leq \pi/2 \f$
*
* \li Standard spherical coordinates :
*
* \f$ \rho = \alpha x + \beta \f$
*
* \f$ \theta = \theta^\star \f$
*
* \ingroup domain
*/
class Domain_polar_shell : public Domain {

 private:
  double alpha ; ///< Relates the numerical to the physical radii.
  double beta ; ///< Relates the numerical to the physical radii.
  Point center ; ///< Absolute coordinates of the center.
 
 public:  
 /**
  * Standard constructor :
  * @param num : number of the domain (used by the \c Space).
  * @param ttype : Chebyshev or Legendre type of spectral expansion.
  * @param r_int : inner radius of the shell.
  * @param r_ext : outer radius of the shell.
  * @param cr : center of the spherical coordinates.
  * @param nbr : number of points in each dimension.
  */
  Domain_polar_shell (int num, int ttype, double r_int, double r_ext, const Point& cr, const Dim_array& nbr) ;
  Domain_polar_shell (const Domain_polar_shell& so) ; ///< Copy constructor.
  /**
  * Constructor from a file
  * @param num : number of the domain (used by the \c Space).
  * @param ff: file containd the domain, generated by the save function.
  */
  Domain_polar_shell (int num, FILE* ff) ;

  virtual ~Domain_polar_shell() ;
  virtual void save (FILE*) const ;

  private:    
    virtual void do_absol ()  const ;
    virtual void do_radius () const ;
    virtual void do_cart ()  const ; 
    virtual void do_cart_surr ()  const ; 

  private :
     virtual void set_cheb_base(Base_spectral&) const ;  
     virtual void set_legendre_base(Base_spectral&) const ;
     virtual void set_anti_cheb_base(Base_spectral&) const ;      
     virtual void set_anti_legendre_base(Base_spectral&) const ; 
     virtual void set_cheb_base_with_m(Base_spectral&, int m) const ;  
     virtual void set_legendre_base_with_m(Base_spectral&, int m) const ;
     virtual void set_anti_cheb_base_with_m(Base_spectral&, int m) const ;      
     virtual void set_anti_legendre_base_with_m(Base_spectral&, int m) const ;     

     virtual void do_coloc () ;
     virtual int give_place_var (char*) const ;
  public:     
     virtual Point get_center () const {return center ;} ;
     virtual bool is_in(const Point& xx, double prec=1e-13) const ;
     virtual const Point absol_to_num(const Point&) const;
     virtual void do_der_abs_from_der_var(Val_domain** der_var, Val_domain** der_abs) const ;
   
     virtual Base_spectral mult (const Base_spectral&, const Base_spectral&) const ;

  public:   
    
     virtual Val_domain mult_cos_theta (const Val_domain&) const ;
     virtual Val_domain mult_sin_theta (const Val_domain&) const ;
     virtual Val_domain div_sin_theta (const Val_domain&) const ;
     virtual Val_domain mult_r (const Val_domain&) const ;
     virtual Val_domain div_r (const Val_domain&) const ;
     virtual Val_domain laplacian (const Val_domain&, int) const ;
     virtual Val_domain laplacian2 (const Val_domain&, int) const ;
     virtual Val_domain der_r (const Val_domain&) const ;
     virtual Val_domain dt (const Val_domain&) const ;
     virtual Val_domain div_xp1 (const Val_domain&) const ; 
     virtual double integrale (const Val_domain&) const ;
     virtual double integ_volume (const Val_domain&) const ;

     virtual void find_other_dom (int, int, int&, int&) const ;     
     virtual Val_domain der_normal (const Val_domain&, int) const ;
     virtual double val_boundary (int, const Val_domain&, const Index&) const ;
     
     virtual int nbr_unknowns (const Tensor&, int) const ;
	/**
	* Computes the number of true unknowns of a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions to determine the precise number of degrees of freedom.
	* @param so : the field.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @returns the number of true unknowns.
	*/
     int nbr_unknowns_val_domain (const Val_domain& so, int mquant) const ;
     virtual Array<int> nbr_conditions (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
      /**
	* Computes number of discretized equations associated with a given tensorial equation in the bulk.
	* It takes into account the various Galerkin basis used.
	* @param so : the residual of the equation.
        * @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param order : order of the equation (i.e. 2 for a Laplacian for instance)
	* @returns the number of true unknowns.
	*/
     int nbr_conditions_val_domain (const Val_domain& so, int mquant, int order) const ;
     virtual Array<int> nbr_conditions_boundary (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given equation on a boundary.
	* It takes into account the various Galerkin basis used.
	* It is used for implementing boundary conditions and matching ones.
	* @param eq : the residual of the equation.
        * @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @returns the number of true conditions.
	*/
     int nbr_conditions_val_domain_boundary (const Val_domain& eq, int mquant) const ;
     virtual void export_tau (const Tensor&, int, int, Array<double>&, int&, const Array<int>&,int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports a residual equation in the bulk.
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param order : describes the order of the equation (2 for a Laplacian for instance).
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond :  the corresponding number of equations. It is used when the equation is null.
	*/
     void export_tau_val_domain (const Val_domain& eq, int mquant, int order, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void export_tau_boundary (const Tensor&, int, int, Array<double>&, int&, const Array<int>&, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports all the residual equations corresponding to a tensorial one on a given boundary
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param bound : the boundary at which the equation is enforced.
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond : the corresponding number of equations. It is used when the residual is null.
	*/
     void export_tau_val_domain_boundary (const Val_domain& eq, int mquant, int bound, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void affecte_tau (Tensor&, int, const Array<double>&, int&) const ;
	/**
	* Affects some coefficients to a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the field to be affected.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param cf : \c Array of the coefficients used.
	* @param pos_cf : current position in the array of coefficients.
	*/
     void affecte_tau_val_domain (Val_domain& so, int mquant, const Array<double>& cf, int& pos_cf) const ;
     virtual void affecte_tau_one_coef (Tensor&, int, int, int&) const ;
	/**
	* Sets at most one coefficient of a \c Val_domain to 1.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the \c Val_domain to be affected. It is set to zero if cc does not corresponds to another field.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param cc : location, in the overall system, of the coefficient to be set to 1.
	* @param pos_cf : current position.
	*/
    void affecte_tau_one_coef_val_domain (Val_domain& so, int mquant, int cc, int& pos_cf) const ;

public:
     virtual ostream& print (ostream& o) const ;
} ;

/**
* Class for a 2-dimensional spherical shell and a symmetry with respect to the plane \f$ z=0 \f$.
* It is intended to deal with systems with a symmetry with respect to \f$\varphi\f$.
* \li 2 dimensions.
* \li centered on the point \c center \f$(X_c, Z_c)\f$
* \li The numerical coordinates are :
*
* \f$ -1 \leq x \leq 1 \f$
*
* \f$ 0 \leq \theta^\star \leq \pi/2 \f$
*
* \li Standard spherical coordinates :
*
* \f$ \rho = \displaystyle\frac{1}{\alpha x -1}\f$
*
* \f$ \theta = \theta^\star \f$
*
* \ingroup domain
*/
class Domain_polar_compact : public Domain {

 private:
  double alpha ; ///< Relates the numerical to the physical radii.
  Point center ; ///< Absolute coordinates of the center.

 public:
  /**
  * Standard constructor :
  * @param num : index of the domain, used by \c Space.
  * @param ttype : Chebyshev or Legendre type of spectral expansion.
  * @param r_int [input] : inner radius of the shell.
  * @param cr [input] : center of the spherical coordinates.
  * @param nbr [nbr] : number of points in each dimension.
  */
  Domain_polar_compact (int num, int ttype, double r_int, const Point& cr, const Dim_array& nbr) ;
  Domain_polar_compact (const Domain_polar_compact& so) ; ///< Copy constructor.
 /**
  * Constructor from a file
  * @param num : number of the domain (used by the \c Space).
  * @param ff: file containd the domain, generated by the save function.
  */
  Domain_polar_compact (int num, FILE* ff) ;

  virtual ~Domain_polar_compact() ;
  virtual void save(FILE*) const ;

  private:
    virtual void do_radius () const ;
  
  private:
     virtual void set_cheb_base(Base_spectral&) const ;     
     virtual void set_legendre_base(Base_spectral&) const ;
     virtual void set_anti_cheb_base(Base_spectral&) const ;       
     virtual void set_anti_legendre_base(Base_spectral&) const ; 
     virtual void set_cheb_base_with_m(Base_spectral&, int m) const ;     
     virtual void set_legendre_base_with_m(Base_spectral&, int m) const ;
     virtual void set_anti_cheb_base_with_m(Base_spectral&, int m) const ;       
     virtual void set_anti_legendre_base_with_m(Base_spectral&, int m) const ;     
     virtual void do_coloc()  ;     

     virtual void do_absol ()  const ; 
     virtual void do_cart ()  const ;      
     virtual void do_cart_surr ()  const ; 
     virtual int give_place_var (char*) const ;
 public:      
    virtual Point get_center () const {return center ;} ;
    /**
     *  Returns the \f$\alpha\f$ of the mapping
     */
    double get_alpha() const {return alpha ;} ;

   
     virtual bool is_in(const Point& xx, double prec=1e-13) const ;    
     virtual const Point absol_to_num(const Point&) const;
     virtual void do_der_abs_from_der_var(Val_domain** der_var, Val_domain** der_abs) const ; 
    
     virtual Base_spectral mult (const Base_spectral&, const Base_spectral&) const ;

  public:      
    
     virtual Val_domain mult_cos_theta (const Val_domain&) const ;
     virtual Val_domain mult_sin_theta (const Val_domain&) const ;
     virtual Val_domain div_sin_theta (const Val_domain&) const ;
     virtual Val_domain div_xm1 (const Val_domain&) const ;      
    
     virtual Val_domain mult_xm1 (const Val_domain&) const ;          
     virtual Val_domain mult_r (const Val_domain&) const ;   
     virtual Val_domain div_r (const Val_domain&) const ;
     virtual Val_domain laplacian (const Val_domain&, int) const ;
     virtual Val_domain laplacian2 (const Val_domain&, int) const ;
     virtual Val_domain der_r (const Val_domain&) const ;   
     virtual Val_domain der_r_rtwo (const Val_domain&) const ;
     virtual Val_domain dt (const Val_domain&) const ; 
     virtual Val_domain div_xp1 (const Val_domain&) const ; 
     virtual double integrale (const Val_domain&) const ;
     virtual double integ_volume (const Val_domain&) const ;
     virtual void set_val_inf (Val_domain& so, double xx) const ;    

     virtual void find_other_dom (int, int, int&, int&) const ;     
     virtual Val_domain der_normal (const Val_domain&, int) const ;
     virtual double integ (const Val_domain&, int) const ;
     virtual double val_boundary (int, const Val_domain&, const Index&) const ;

     virtual int nbr_unknowns (const Tensor&, int) const ;
	/**
	* Computes the number of true unknowns of a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions to determine the precise number of degrees of freedom.
	* @param so : the field.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @returns the number of true unknowns.
	*/
     int nbr_unknowns_val_domain (const Val_domain& so, int mquant) const ;
     virtual Array<int> nbr_conditions (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given tensorial equation in the bulk.
	* It takes into account the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param order : order of the equation (i.e. 2 for a Laplacian for instance)
	* @returns the number of true unknowns.
	*/
     int nbr_conditions_val_domain (const Val_domain& eq, int mquant, int order) const ;
     virtual Array<int> nbr_conditions_boundary (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given equation on a boundary.
	* It takes into account the various Galerkin basis used.
	* It is used for implementing boundary conditions and matching ones.
	* @param eq : the residual of the equation.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @returns the number of true conditions.
	*/
     int nbr_conditions_val_domain_boundary (const Val_domain& eq, int mquant) const ;
     virtual void export_tau (const Tensor&, int, int, Array<double>&, int&, const Array<int>&,int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports a residual equation in the bulk.
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param order : describes the order of the equation (2 for a Laplacian for instance).
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond : the corresponding number of equations. It is used when the equation is null.
	*/
     void export_tau_val_domain (const Val_domain& eq, int mquant, int order, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void export_tau_boundary (const Tensor&, int, int, Array<double>&, int&, const Array<int>&, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports all the residual equations corresponding to a tensorial one on a given boundary
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param bound : the boundary at which the equation is enforced.
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond : the corresponding number of equations. It is used when the equation is null. 
	*/
     void export_tau_val_domain_boundary (const Val_domain& eq, int mquant, int bound, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void affecte_tau (Tensor&, int, const Array<double>&, int&) const ;
	/**
	* Affects some coefficients to a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the field to be affected.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param cf : \c Array of the coefficients used.
	* @param pos_cf : current position in the array of coefficients.
	*/
     void affecte_tau_val_domain (Val_domain& so, int mquant, const Array<double>& cf, int& pos_cf) const ;
     virtual void affecte_tau_one_coef (Tensor&, int, int, int&) const ;
	/**
	* Sets at most one coefficient of a \c Val_domain to 1.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the \c Val_domain to be affected. It is set to zero if cc does not corresponds to another field.
	* @param  mquant : the field is the \c mquant harmonic with respexcto to \f$\varphi\f$.
	* @param cc : location, in the overall system, of the coefficient to be set to 1.
	* @param pos_cf : current position.
	*/
     void affecte_tau_one_coef_val_domain (Val_domain& so, int mquant, int cc, int& pos_cf) const ;

public:
     virtual ostream& print (ostream& o) const ;
} ;


/**
 * The \c Space_polar class fills the space with one polar nucleus, several polar shells and a compactified polar domain, all centered on the same point.
 * \ingroup domain
 */
class Space_polar : public Space {
     public:
    	/**
     	* Standard constructor 
     	* @param ttype [input] : the type of basis.
	* @param cr [input] : absolute coordinates of the center.
	* @param nbr [input] : number of points in each domain.
	* @param bounds [input] : radii of the various shells (and also determines the total number of domains).
	*/
	Space_polar (int ttype, const Point& cr, const Dim_array& nbr, const Array<double>& bounds) ;
	Space_polar (FILE*) ; ///< Constructor from a file
	virtual ~Space_polar() ;      
	virtual void save(FILE*) const ;

	/**
	* Sets a boundary condition at the inner radius of the innermost polar shell.
	* Intended for systems where the nucleus is not used.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the boundary condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_inner_bc (System_of_eqs& syst, const char* eq, int nused=-1, Array<int>** pused=0x0)  const ;
	
	/**
	* Sets a boundary condition at the outer radius of the compactified polar domain.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the boundary condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_outer_bc (System_of_eqs& syst, const char* eq, int nused=-1, Array<int>** pused=0x0)  const ;

	/**
	* Adds a bulk equation in all the domains of a given system (equation is assumed to be second order)
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the equation.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_eq (System_of_eqs& syst, const char* eq, int nused=-1, Array<int>** pused=0x0)  const ;

	/**
	* Adds a bulk equation in all the domains of a given system (equation is assumed to be 0th order)
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the equations.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_eq_full (System_of_eqs& syst, const char* eq, int nused=-1, Array<int>** pused=0x0)  const ;

	/**
	* Adds a matching condition, at all the interface present in a given system.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the matching condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_matching (System_of_eqs& syst, const char* eq, int nused=-1, Array<int>** pused=0x0)  const ;
	
	/**
	* Adds a bulk equation and two matching conditions.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the bulk equation.
	* @param rac : the string describing the first matching condition.
	* @param rac_der : the string describing the second matching condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_eq (System_of_eqs& syst, const char* eq, const char* rac, const char* rac_der, int nused=-1, Array<int>** pused=0x0) const ; 
	
	/**
	* Adds an equation being a volume integral in the whole computational space.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the equation (should contain something like integvolume(f)=b)
	*/
	void add_eq_int_volume (System_of_eqs& syst, const char* eq) ;
	
	/**
	* Adds an equation being a surface integral at infinity.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the equation (should contain something like integ(f)=b)
	*/
	void add_eq_int_inf (System_of_eqs& syst, const char* eq) ;
	
	/**
	* Adds an equation being a surface integral in the innermost shell.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the equation (should contain something like integ(f)=b)
	*/
	void add_eq_int_inner (System_of_eqs& syst, const char* eq) ;
	
	/**
	* Adds an equation saying that one coefficient of a field is zero in a given domain
	* @param syst : the \c System_of_eqs.
	* @param f : the field
	* @param domtarget : the target domain.
	* @param itarget : the index \f$r\f$ of the mode that must vanish.
	* @param jtarget : the index \f$\theta\f$ of the mode that must vanish.
	*/
	void add_eq_mode (System_of_eqs& syst, const char* f, int domtarget, int itarget, int jtarget) ;


	/**
	* Adds an equation being the value of some field at the origin.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the quantity that must be zero at the origin
	*/
	void add_eq_ori (System_of_eqs& syst, const char* eq) ;

	/**
	* Adds an equation being the value of some field at a given point.
	* @param syst : the \c System_of_eqs.
	* @param pp : the point
	* @param eq : the string describing the quantity that must be zero at the \c pp
	*/
	void add_eq_point (System_of_eqs& syst, const Point& pp, const char* eq) ;
		
	/**
	* Computes the surface integral at infinity.
	* @param so : the field to be integrated.
	* @returns the surface integral.
	*/
	double int_inf (const Scalar& so) const ;
} ;
}
#endif
