 /*@@
   @file    InitSymBound.F
   @date    May 13 1999
   @author  Gabrielle Allen
   @desc
            Sets the symmetries for the Constraint grid functions
   @enddesc
   @version $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"

 /*@@
   @routine    InitSymBound
   @date       May 13 1999
   @author     Gabrielle Allen
   @desc 
               Sets the symmetries for the Constraint grid functions
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine ADMConstraint_InitSymBound(CCTK_ARGUMENTS)

      implicit none
      
      DECLARE_CCTK_ARGUMENTS 
      
      INTEGER one 
      PARAMETER (one=1)

      INTEGER sym(3)
      INTEGER ierr

      sym(1) = one
      sym(2) = one
      sym(3) = one
      call SetCartSymVN(ierr, cctkGH, sym, 'admconstraints::ham')
      if (ierr .lt. 0) then
        call CCTK_WARN(0,"Error setting symmetries for ADMConstraints::ham")
      end if
      sym(1) = one
      sym(2) = one
      sym(3) = one
      call SetCartSymVN(ierr, cctkGH, sym, 'admconstraints::hamnormalized')
      if (ierr .lt. 0) then
        call CCTK_WARN(0,"Error setting symmetries for ADMConstraints::hamnormalized")
      end if
      sym(1) = -one
      sym(2) = one
      sym(3) = one
      call SetCartSymVN(ierr, cctkGH, sym, 'admconstraints::momx')
      if (ierr .lt. 0) then
        call CCTK_WARN(0,"Error setting symmetries for ADMConstraints::momx")
      end if
      sym(1) = one
      sym(2) = -one
      sym(3) = one
      call SetCartSymVN(ierr, cctkGH, sym, 'admconstraints::momy')
      if (ierr .lt. 0) then
        call CCTK_WARN(0,"Error setting symmetries for ADMConstraints::momy")
      end if
      sym(1) = one
      sym(2) = one
      sym(3) = -one
      call SetCartSymVN(ierr, cctkGH, sym, 'admconstraints::momz')
      if (ierr .lt. 0) then
        call CCTK_WARN(0,"Error setting symmetries for ADMConstraints::momz")
      end if

      return
      end 

