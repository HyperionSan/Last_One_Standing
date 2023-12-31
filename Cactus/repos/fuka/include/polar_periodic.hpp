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

#ifndef __POLAR_PERIODIC_HPP_
#define __POLAR_PERIODIC_HPP_

#include "space.hpp"

#define TO_PI 0
namespace Kadath {

/**
* Class for a spherical nucleus with a symmetry in \f$\varphi\f$. The spacelike coordinates are thus 2d being \f$(r,\theta)\f$.
* There is a third dimension, being time. The fields are supposed to be periodic.
* The numerical time goes either to \f$\pi\f$.
*
* The period is an unknown, prescribed by additionnal constraints on the solution. 
*
* \li The numerical coordinates are :
*
* \f$ 0 \leq x \leq 1 \f$
*
* \f$ 0 \leq \theta^\star \leq \pi/2 \f$
*
* \f$ 0 \leq t^\star \leq \quad \pi\f$

* \li Standard spherical coordinates :
*
* \f$ r = \alpha x \f$
*
* \f$ \theta = \theta^\star \f$
*
* \f$ t = t^\star / \omega\f$
*
* \ingroup domain
*/
class Domain_polar_periodic_nucleus : public Domain {

 private:
  double alpha ; ///< Relates the numerical radius to the physical one.

  double ome ; ///< The pulsation.
  mutable Term_eq* ome_term_eq ; ///< Pointer on the \c Term_eq version of the pulsation.

  /**
   * Gives the type of time periodicity.
   * \li if 0 \f$ t^\star\f$ goes to \f$ \pi\f$.
   */
  int type_time ;
  double maxt ; ///< Upper bound of \f$ t^\star\f$ which is \f$\pi\f$ or \f$\pi/2\f$.

 public: 
   /**
  * Standard constructor :
  * @param num : number of the domain (used by the \c Space).
  * @param ttype : Chebyshev or Legendre type of spectral expansion.
  * @param radius : radius of the nucleus.
  * @param ome : pulsation.
  * @param nbr : number of points in each dimension.
  */
  Domain_polar_periodic_nucleus (int num, int ttype, double radius, double ome, const Dim_array& nbr) ;
  Domain_polar_periodic_nucleus (const Domain_polar_periodic_nucleus& so) ; ///< Copy constructor.
 /**
  * Constructor from a file
  * @param num : number of the domain (used by the \c Space).
  * @param ff: file containd the domain, generated by the save function.
  */
  Domain_polar_periodic_nucleus (int num, FILE* ff) ;

  virtual ~Domain_polar_periodic_nucleus() ;
  virtual void save (FILE*) const ;

  private:    
    virtual void do_absol ()  const ;

  private:
     virtual void set_cheb_base(Base_spectral&) const ;       
     virtual void set_legendre_base(Base_spectral&) const ;
     virtual void set_anti_cheb_base(Base_spectral&) const ;       
     virtual void set_anti_legendre_base(Base_spectral&) const ;
     virtual void set_cheb_base_with_m(Base_spectral&, int m) const ;       
     virtual void set_legendre_base_with_m(Base_spectral&, int m) const ;
     virtual void set_anti_cheb_base_with_m(Base_spectral&, int m) const ;       
     virtual void set_anti_legendre_base_with_m(Base_spectral&, int m) const ;
     virtual void set_cheb_base_r_spher(Base_spectral&) const ;
     virtual void set_cheb_base_t_spher(Base_spectral&) const ;
     virtual void set_cheb_base_p_spher(Base_spectral&) const ;
     virtual void set_legendre_base_r_spher(Base_spectral&) const ;
     virtual void set_legendre_base_t_spher(Base_spectral&) const ;
     virtual void set_legendre_base_p_spher(Base_spectral&) const ;
 

     virtual void do_coloc () ;
     virtual void do_radius () const ;

  public:
     virtual bool is_in(const Point&xx, double prec=1e-13) const ;
     virtual const Point absol_to_num(const Point&) const;
     virtual void do_der_abs_from_der_var(Val_domain** der_var, Val_domain** der_abs) const ;
     virtual Base_spectral mult (const Base_spectral&, const Base_spectral&) const ;
     /**
	* Returns omega
	*/ 
     double get_ome() const {return ome ;} ;

  public:
     virtual Val_domain mult_cos_theta (const Val_domain&) const ;
     virtual Val_domain mult_sin_theta (const Val_domain&) const ;
     virtual Val_domain div_sin_theta (const Val_domain&) const ;
     virtual Val_domain div_x (const Val_domain&) const ;     
     virtual Val_domain mult_r (const Val_domain&) const ;
     virtual Val_domain div_r (const Val_domain&) const ;
     virtual Val_domain der_r (const Val_domain&) const ; 
     virtual Val_domain dt (const Val_domain&) const ;
     virtual Val_domain srdr (const Val_domain&) const ;
     virtual Val_domain dtime (const Val_domain&) const ;
     virtual Val_domain mult_cos_time (const Val_domain&) const ;
     virtual Val_domain mult_sin_time (const Val_domain&) const ;
   
     virtual Term_eq partial_spher (const Term_eq&) const ;
     virtual Term_eq derive_flat_spher (int, char, const Term_eq&, const Metric*) const ;

     virtual double val_boundary (int, const Val_domain&, const Index&) const ;
     virtual void find_other_dom (int, int, int&, int&) const ;
     virtual Val_domain der_normal (const Val_domain&, int) const ;
 	

  //   virtual int nbr_unknowns_from_adapted() const ;
   //  virtual void vars_to_terms() const ;
   //  virtual void affecte_coef (int&, int, bool&) const ;
   //  virtual void xx_to_vars_from_adapted (double, const Array<double>&, int&) const ;
   //  virtual void xx_to_ders_from_adapted (const Array<double>&, int&) const ;
   //  virtual void update_term_eq (Term_eq*) const ;
   //  virtual void update_variable (double, const Scalar&, Scalar&) const ;   
   //  virtual void update_constante (double, const Scalar&, Scalar&) const ;   
   //  virtual void update_mapping(double) ;

     /**
     * Updates the quantities that depend on the frequency.
     */
    // void update() const ;
    
     virtual Term_eq dtime_term_eq (const Term_eq& so) const ;   
     virtual Term_eq ddtime_term_eq (const Term_eq& so) const ;   

         virtual int nbr_unknowns (const Tensor&, int) const ;
	/**
	* Computes the number of true unknowns of a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions to determine the precise number of degrees of freedom.
	* @param so : the field.
	* @param llim: limit for the regularity (quantum number wrt \f$\theta\f$).
	* @returns the number of true unknowns.
	*/
     int nbr_unknowns_val_domain (const Val_domain& so,  int llim) const ;

     virtual Array<int> nbr_conditions (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given tensorial equation in the bulk.
	* It takes into account the various Galerkin basis used.
	* @param so : the residual of the equation.
	* @param llim : limit for the regularity (quantum number wrt \f$\theta\f$).
	* @param order : order of the equation (i.e. 2 for a Laplacian for instance)
	* @returns the number of true unknowns.
	*/
     int nbr_conditions_val_domain (const Val_domain& so, int llim, int order) const ;
     virtual Array<int> nbr_conditions_boundary (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given equation on a boundary.
	* It takes into account the various Galerkin basis used.
	* It is used for implementing boundary conditions and matching ones.
	* @param eq : the residual of the equation.
	* @returns the number of true conditions.
	*/
     int nbr_conditions_val_domain_boundary (const Val_domain& eq) const ;
     virtual void export_tau (const Tensor&, int, int, Array<double>&, int&, const Array<int>&,int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports a residual equation in the bulk.
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param llim : limit for the regularity (quantum number wrt \f$\theta\f$).
	* @param order : describes the order of the equation (2 for a Laplacian for instance).
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond :  the corresponding number of equations. It is used when the equation is null.
	*/
     void export_tau_val_domain (const Val_domain& eq, int llim, int order, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void export_tau_boundary (const Tensor&, int, int, Array<double>&, int&, const Array<int>&, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports all the residual equations corresponding to a tensorial one on a given boundary
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param bound : the boundary at which the equation is enforced.
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond :  the corresponding number of equations. It is used when the residual is null.
	*/
     void export_tau_val_domain_boundary (const Val_domain& eq, int bound, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void affecte_tau (Tensor&, int, const Array<double>&, int&) const ;
	/**
	* Affects some coefficients to a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the field to be affected.
	* @param llim : limit for the regularity (quantum number wrt \f$\theta\f$).
	* @param cf : \c Array of the coefficients used.
	* @param pos_cf : current position in the array of coefficients.
	*/
     void affecte_tau_val_domain (Val_domain& so, int llim, const Array<double>& cf, int& pos_cf) const ;
     virtual void affecte_tau_one_coef (Tensor&, int, int, int&) const ;
	/**
	* Sets at most one coefficient of a \c Val_domain to 1.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the \c Val_domain to be affected. It is set to zero if cc does not corresponds to another field.
	* @param llim : limit for the regularity (quantum number wrt \f$\theta\f$).
	* @param cc : location, in the overall system, of the coefficient to be set to 1.
	* @param pos_cf : current position.
	*/
      void affecte_tau_one_coef_val_domain (Val_domain& so, int llim, int cc, int& pos_cf) const ;


public:
     virtual ostream& print (ostream& o) const ;
} ;

/**
* Class for a spherical shell with a symmetry in \f$\varphi\f$. The spacelike coordinates are thus 2d being \f$(r,\theta)\f$.
* There is a third dimension, being time. The fields are supposed to be periodic.
* The numerical time goes either to \f$\pi\f$.
*
* The period is an unknown, prescribed by additionnal constraints on the solution. 
*
* \li The numerical coordinates are :
*
* \f$ -1 \leq x \leq 1 \f$
*
* \f$ 0 \leq \theta^\star \leq \pi/2 \f$
*
* \f$ 0 \leq t^\star \leq \quad \pi\f$

* \li Standard spherical coordinates :
*
* \f$ r = \alpha x + \beta\f$
*
* \f$ \theta = \theta^\star \f$
*
* \f$ t = t^\star / \omega\f$
*
* \ingroup domain
*/
class Domain_polar_periodic_shell : public Domain {

 private:
  double alpha ; ///< Relates the numerical radius to the physical one.
  double beta ; ///< Relates the numerical radius to the physical one.


  double ome ; ///< The pulsation.
  mutable Term_eq* ome_term_eq ; ///< Pointer on the \c Term_eq version of the pulsation.

  /**
   * Gives the type of time periodicity.
   * \li if 0 \f$ t^\star\f$ goes to \f$ \pi\f$.
   */
  int type_time ;
  double maxt ; ///< Upper bound of \f$ t^\star\f$ which is \f$\pi\f$ or \f$\pi/2\f$.

 public: 
   /**
  * Standard constructor :
  * @param num : number of the domain (used by the \c Space).
  * @param ttype : Chebyshev or Legendre type of spectral expansion.
  * @param rin : inner radius of the shell.  
  * @param rout : outer radius of the shell.
  * @param ome : pulsation.
  * @param nbr : number of points in each dimension.
  */
  Domain_polar_periodic_shell (int num, int ttype, double rin, double rout, double ome, const Dim_array& nbr) ;
  Domain_polar_periodic_shell (const Domain_polar_periodic_shell& so) ; ///< Copy constructor.
 /**
  * Constructor from a file
  * @param num : number of the domain (used by the \c Space).
  * @param ff: file containd the domain, generated by the save function.
  */
  Domain_polar_periodic_shell (int num, FILE* ff) ;

  virtual ~Domain_polar_periodic_shell() ;
  virtual void save (FILE*) const ;

  private:    
    virtual void do_absol ()  const ;

  private:
     virtual void set_cheb_base(Base_spectral&) const ;       
     virtual void set_legendre_base(Base_spectral&) const ;
     virtual void set_anti_cheb_base(Base_spectral&) const ;       
     virtual void set_anti_legendre_base(Base_spectral&) const ;
     virtual void set_cheb_base_with_m(Base_spectral&, int m) const ;       
     virtual void set_legendre_base_with_m(Base_spectral&, int m) const ;
     virtual void set_anti_cheb_base_with_m(Base_spectral&, int m) const ;       
     virtual void set_anti_legendre_base_with_m(Base_spectral&, int m) const ;
     virtual void set_cheb_base_r_spher(Base_spectral&) const ;
     virtual void set_cheb_base_t_spher(Base_spectral&) const ;
     virtual void set_cheb_base_p_spher(Base_spectral&) const ;
     virtual void set_legendre_base_r_spher(Base_spectral&) const ;
     virtual void set_legendre_base_t_spher(Base_spectral&) const ;
     virtual void set_legendre_base_p_spher(Base_spectral&) const ;
 
     virtual void do_coloc () ;
     virtual void do_radius () const ;

  public:
     virtual bool is_in(const Point&xx, double prec=1e-13) const ;
     virtual const Point absol_to_num(const Point&) const;
     virtual void do_der_abs_from_der_var(Val_domain** der_var, Val_domain** der_abs) const ;
     virtual Base_spectral mult (const Base_spectral&, const Base_spectral&) const ;

  public:
     virtual Val_domain mult_cos_theta (const Val_domain&) const ;
     virtual Val_domain mult_sin_theta (const Val_domain&) const ;
     virtual Val_domain div_sin_theta (const Val_domain&) const ;
     virtual Val_domain mult_r (const Val_domain&) const ;
     virtual Val_domain div_r (const Val_domain&) const ;
     virtual Val_domain der_r (const Val_domain&) const ; 
     virtual Val_domain dt (const Val_domain&) const ; 
     virtual Val_domain div_1mrsL (const Val_domain&) const ;
     virtual Val_domain dtime (const Val_domain&) const ;
     virtual Val_domain mult_cos_time (const Val_domain&) const ;
     virtual Val_domain mult_sin_time (const Val_domain&) const ;
   
     virtual Term_eq partial_spher (const Term_eq&) const ;
     virtual Term_eq derive_flat_spher (int, char, const Term_eq&, const Metric*) const ;

     virtual double val_boundary (int, const Val_domain&, const Index&) const ;
     virtual void find_other_dom (int, int, int&, int&) const ;
     virtual Val_domain der_normal (const Val_domain&, int) const ;

  /*   virtual int nbr_unknowns_from_adapted() const ;
     virtual void vars_to_terms() const ;
     virtual void affecte_coef (int&, int, bool&) const ;
     virtual void xx_to_vars_from_adapted (double, const Array<double>&, int&) const ;
     virtual void xx_to_ders_from_adapted (const Array<double>&, int&) const ;
     virtual void update_term_eq (Term_eq*) const ;
     virtual void update_variable (double, const Scalar&, Scalar&) const ;   
     virtual void update_constante (double, const Scalar&, Scalar&) const ;   
     virtual void update_mapping(double) ;
*/
     /**
     * Updates the quantities that depend on the frequency.
     */
  //   void update() const ;
    
     virtual Term_eq dtime_term_eq (const Term_eq& so) const ;   
     virtual Term_eq ddtime_term_eq (const Term_eq& so) const ;   


         virtual int nbr_unknowns (const Tensor&, int) const ;
	/**
	* Computes the number of true unknowns of a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions to determine the precise number of degrees of freedom.
	* @param so : the field.
	* @returns the number of true unknowns.
	*/
     int nbr_unknowns_val_domain (const Val_domain& so) const ;

     virtual Array<int> nbr_conditions (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given tensorial equation in the bulk.
	* It takes into account the various Galerkin basis used.
	* @param so : the residual of the equation.
	* @param order : order of the equation (i.e. 2 for a Laplacian for instance)
	* @returns the number of true unknowns.
	*/
     int nbr_conditions_val_domain (const Val_domain& so, int order) const ;
     virtual Array<int> nbr_conditions_boundary (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given equation on a boundary.
	* It takes into account the various Galerkin basis used.
	* It is used for implementing boundary conditions and matching ones.
	* @param eq : the residual of the equation.
	* @returns the number of true conditions.
	*/
     int nbr_conditions_val_domain_boundary (const Val_domain& eq) const ;
     virtual void export_tau (const Tensor&, int, int, Array<double>&, int&, const Array<int>&,int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports a residual equation in the bulk.
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param order : describes the order of the equation (2 for a Laplacian for instance).
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond :  the corresponding number of equations. It is used when the equation is null.
	*/
     void export_tau_val_domain (const Val_domain& eq, int order, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void export_tau_boundary (const Tensor&, int, int, Array<double>&, int&, const Array<int>&, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports all the residual equations corresponding to a tensorial one on a given boundary
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param bound : the boundary at which the equation is enforced.
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond :  the corresponding number of equations. It is used when the residual is null.
	*/
     void export_tau_val_domain_boundary (const Val_domain& eq, int bound, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void affecte_tau (Tensor&, int, const Array<double>&, int&) const ;
	/**
	* Affects some coefficients to a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the field to be affected.
	* @param cf : \c Array of the coefficients used.
	* @param pos_cf : current position in the array of coefficients.
	*/
     void affecte_tau_val_domain (Val_domain& so, const Array<double>& cf, int& pos_cf) const ;
     virtual void affecte_tau_one_coef (Tensor&, int, int, int&) const ;
	/**
	* Sets at most one coefficient of a \c Val_domain to 1.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the \c Val_domain to be affected. It is set to zero if cc does not corresponds to another field.
	* @param cc : location, in the overall system, of the coefficient to be set to 1.
	* @param pos_cf : current position.
	*/
      void affecte_tau_one_coef_val_domain (Val_domain& so, int cc, int& pos_cf) const ;


public:
     virtual ostream& print (ostream& o) const ;
} ;


/**
 * The \c Space_polar_periodic class fills the space with one polar nucleus and several polar shells, with a periodicity in time.
 * \ingroup domain
 */
class Space_polar_periodic : public Space {
     public:
    	/**
     	* Standard constructor 
     	* @param ttype [input] : the type of basis.
	* @param omega [input] : the pulsation for the time dependance.
	* @param nbr [input] : number of points in each domain.
	* @param bounds [input] : radii of the various shells (and also determines the total number of domains).
	*/
	Space_polar_periodic (int ttype, double omega, const Dim_array& nbr, const Array<double>& bounds) ;
	Space_polar_periodic (FILE*) ; ///< Constructor from a file
	virtual ~Space_polar_periodic() ;      
	virtual void save(FILE*) const ;

	/*virtual int nbr_unknowns_from_variable_domains() const ;
	virtual void affecte_coef_to_variable_domains(int& , int, Array<int>&) const ;
	virtual void xx_to_ders_variable_domains(const Array<double>&, int&) const ;
	virtual void xx_to_vars_variable_domains(System_of_eqs*, const Array<double>&, int&) const ;
*/
	double get_omega() const ;
} ;

}
#endif
