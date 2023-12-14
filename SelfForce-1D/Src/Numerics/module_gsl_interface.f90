module gsl_interface
!! Module that defines the Fortran interface for a couple of C-routines from
!! the Gnu Scientific Libray (GSL).

  use kinds
  use iso_c_binding

  implicit none

  interface

    function legendre_sphPlm ( l, m, x ) bind(c, name='gsl_sf_legendre_sphPlm')
    !! Interface to the C-routine gsl_sf_legendre_sphPlm for calculating the
    !! associated Legendre polynomials used for spherical harmonics.
      use iso_c_binding
      integer(c_int), value :: l
      !! The l-value of the associated Legendre polynomial.
      integer(c_int), value :: m
      !! The m-value of the associated Legendre polynomial.
      real(c_double), value :: x
      !! The polar angle argument \(x=\cos(\theta)\).
      real(c_double) :: legendre_sphPlm
    end function legendre_sphPlm 

    subroutine legendre_sphPlm_d ( lmax, x, ylm, dylm ) bind(c, name='gsl_sf_legendre_sphPlm_deriv')
    !! Interface to the C-routine gsl_sf_legendre_sphPlm_deriv for calculating the
    !! associated Legendre polynomials and their first derivative.
      use iso_c_binding
      integer(c_int), value :: lmax
      !! The max l-value of the associated Legendre polynomial array.
      real(c_double), value :: x
      !! The polar angle argument \(x=\cos(\theta)\).
      real(c_double) :: ylm(*), dylm(*)
    end subroutine legendre_sphPlm_d

    subroutine legendre_sphPlm_d2 ( lmax, x, ylm, dylm, d2ylm ) bind(c, name='gsl_sf_legendre_sphPlm_deriv2')
    !! Interface to the C-routine gsl_sf_legendre_sphPlm_deriv2 for calculating the
    !! associated Legendre polynomials and their first two derivatives.
      use iso_c_binding
      integer(c_int), value :: lmax
      !! The max l-value of the associated Legendre polynomial array.
      real(c_double), value :: x
      !! The polar angle argument \(x=\cos(\theta)\).
      real(c_double) :: ylm(*), dylm(*), d2ylm(*)
    end subroutine legendre_sphPlm_d2

    function legendre_array_index ( l, m ) bind(c, name='gsl_sf_legendre_index')
    !! Interface to the C-routine gsl_sf_legendre_array_index for finding the
    !! index for an (l, m) mode of the calculated legendre derivatives.
      use iso_c_binding
      integer(c_int), value :: l
      !! The l-value of the associated Legendre polynomial.
      integer(c_int), value :: m
      !! The m-value of the associated Legendre polynomial.
      integer(c_int) :: legendre_array_index
      !! The array index for the given mode.
    end function legendre_array_index

    function legendre_array_n ( lmax ) bind(c, name='gsl_sf_legendre_n')
    !! Interface to the C-routine gsl_sf_legendre_array_index for finding the
    !! index for an (l, m) mode of the calculated legendre derivatives.
      use iso_c_binding
      integer(c_int), value :: lmax
      !! The maximum l-value needed for the legendre polynomials.
      integer(c_int) :: legendre_array_n
      !! The array length for the legendre objects.
    end function legendre_array_n

    subroutine multifit_linear ( n, p, xm, yv, cv, covm, chisq, ierr ) &
                    bind(c, name='c_multifit_linear')
    !! Interface to a C-routine that sets up the data structures needed for
    !! calling the GSL routine gsl_multifit_linear.
      use iso_c_binding
      integer(c_int), value :: n
      !! The number of data points.
      integer(c_int), value :: p
      !! The number of parameters to fit.
      real(c_double) :: xm(p*n)
      !! A real 1d array of size p*n containing all the fit models evaluated
      !! at all data points.
      real(c_double) :: yv(n)
      !! A real 1d array of size n that contains all the observed data values.
      real(c_double) :: cv(p)
      !! A real 1d array of size p that on output contains the best fit
      !! parameters.
      real(c_double) :: covm(p*p)
      !! A real 1d array of size p*p that on output contains the covariance
      !! matrix.
      real(c_double) :: chisq
      !! A real variable that on output contains the sum of squares of the
      !! residuals.
      integer(c_int) :: ierr
      !! An error code for the call to the C-routine.
    end subroutine multifit_linear

  end interface

end module gsl_interface
