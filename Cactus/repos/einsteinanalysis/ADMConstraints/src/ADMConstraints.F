/*@@
  @file      ADMConstraints.F77
  @date      August 98
  @desc 
             Calculate the ADM Constraints for output:
    
             Hamiltonian Constraint is:

             H = R - K^i_j K^j_i + trK^2 - 16 Pi rho

             Momentum Constraints are:

             M_i = Del_j K_i^j - Del_i trK - 8 Pi S_i

  @enddesc
  @version $Header$
@@*/

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "SpaceMask.h"

      subroutine ADMConstraints(CCTK_ARGUMENTS)
      
      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer i,j,k
      integer nx,ny,nz

c     Stencil width used for calculating constraints
c     (for outer boundary condition)

      integer sw(3)
      logical docalc

c     flags for excision

      integer ex_field, ex_type_excised

c     Various real variables.

      CCTK_REAL m_rho,m_sx,m_sy,m_sz
      CCTK_REAL pi,ialp,ialp2
      CCTK_REAL det,uxx,uyy,uzz,uxy,uxz,uyz

#include "EinsteinBase/ADMMacros/src/macro/ADM_Spacing_declare.h"

c     Temporaries for the Stress-Energy tensor.

      CCTK_REAL Ttt,Ttx,Tty,Ttz,Txx,Txy,Txz,Tyy,Tyz,Tzz

c     Matter declarations.

#include "CalcTmunu_temps.inc"

c     Macros from Standard Einstein.

#include "EinsteinBase/ADMMacros/src/macro/HAMADM_declare.h"
#include "EinsteinBase/ADMMacros/src/macro/MOMXADM_declare.h"
#include "EinsteinBase/ADMMacros/src/macro/MOMYADM_declare.h"
#include "EinsteinBase/ADMMacros/src/macro/MOMZADM_declare.h"
#include "EinsteinBase/ADMMacros/src/macro/DETG_declare.h"
#include "EinsteinBase/ADMMacros/src/macro/UPPERMET_declare.h"

c --------------------------------------------------------------

      if (check_excision_bitmask .ne. 0) then
          call SpaceMask_GetTypeBits(ex_field, excision_mask_name)
          call SpaceMask_GetStateBits(ex_type_excised, excision_mask_name, \
                                      excision_type_excised)
      end if

      sw(1) = 1
      sw(2) = 1
      sw(3) = 1

c     Grid parameters.

#include "EinsteinBase/ADMMacros/src/macro/ADM_Spacing.h"

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

c     Fill with zeros.

      do k=1,nz
         do j=1,ny
            do i=1,nx

              ham(i,j,k) = 0.0D0
              hamnormalized(i,j,k) = 0.0D0
        
              momx(i,j,k) = 0.0D0
              momy(i,j,k) = 0.0D0
              momz(i,j,k) = 0.0D0       

            end do
         end do
      end do

c     Calculate constraints.

      pi = acos(-1.0D0)

      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1

              if ( (i<3).or.(i>cctk_lsh(1)-2).or.
     .             (j<3).or.(j>cctk_lsh(2)-2).or.
     .             (k<3).or.(k>cctk_lsh(3)-2) ) then
                local_spatial_order = 2
              else
                local_spatial_order = spatial_order
              end if

               docalc = .TRUE.

               if (use_mask .eq. 1) then
                  if (abs(emask(i,j,k)-1) > 0.001) then
                     docalc = .FALSE.
                  elseif (check_excision_bitmask .ne. 0) then
                     if (SpaceMask_CheckStateBitsF90(space_mask, i, j, k,\
                          ex_field, ex_type_excised)) then
                        docalc = .FALSE.
                     end if
                  end if
               end if
                  
               if (docalc) then

                  ialp  = 1.0D0/alp(i,j,k)
                  ialp2 = ialp**2
                  
c                 Calculate the stress energy tensor at this point
c                 ------------------------------------------------

c                 This may be needed for CalcTmunu

#include "EinsteinBase/ADMMacros/src/macro/DETG_guts.h"

                  det = DETG_DETCG

#include "EinsteinBase/ADMMacros/src/macro/UPPERMET_guts.h"

                  uxx = UPPERMET_UXX; uxy = UPPERMET_UXY; uxz = UPPERMET_UXZ
                  uyy = UPPERMET_UYY; uyz = UPPERMET_UYZ; uzz = UPPERMET_UZZ

c                 Initialize stress-energy tensor components.

                  if (stress_energy_state .ne. 0) then

                     Ttt = eTtt(i,j,k)

                     Ttx = eTtx(i,j,k); Tty = eTty(i,j,k); Ttz = eTtz(i,j,k)

                     Txx = eTxx(i,j,k); Tyy = eTyy(i,j,k); Tzz = eTzz(i,j,k)

                     Txy = eTxy(i,j,k); Txz = eTxz(i,j,k); Tyz = eTyz(i,j,k)

                  else

                     Ttt = 0.0D0

                     Ttx = 0.0D0; Tty = 0.0D0; Ttz = 0.0D0

                     Txx = 0.0D0; Tyy = 0.0D0; Tzz = 0.0D0
                     Txy = 0.0D0; Txz = 0.0D0; Tyz = 0.0D0

c                    Include macro for stress energy tensor.

#include "CalcTmunu.inc"
                  
                  end if

c                 Calculate the hamiltonian constraint
c                 ------------------------------------

c                 Geometric piece.

#include "EinsteinBase/ADMMacros/src/macro/HAMADM_guts.h"

c                 Add matter terms: - 16*pi*rho
c
c                 with rho defined as:
c
c                 rho = n_a n_b T^{ab} = n^a n^b T_{ab}
c                     = (T_00 - 2 beta^i T_{i0} + beta^i beta^j T_{ij})/alpha^2

                  m_rho = ialp2*Ttt

                  if (shift_state .eq. 1) then

                     m_rho = m_rho + ialp2
     &                  *(betax(i,j,k)**2*Txx
     &                  + betay(i,j,k)**2*Tyy
     &                  + betaz(i,j,k)**2*Tzz
     &                  +(betax(i,j,k)*betay(i,j,k)*Txy
     &                  + betax(i,j,k)*betaz(i,j,k)*Txz
     &                  + betay(i,j,k)*betaz(i,j,k)*Tyz)*2.0D0
     &                  -(betax(i,j,k)*Ttx
     &                  + betay(i,j,k)*Tty
     &                  + betaz(i,j,k)*Ttz)*2.0D0)

                  end if

                  ham(i,j,k) = HAMADM_HAMADM - 16.0D0*pi*m_rho
                  if ((HAMADM_HAMADMABS + abs(16.0D0*pi*m_rho))==0) then
                     hamnormalized(i,j,k) = abs(ham(i,j,k))
                  else
                     hamnormalized(i,j,k) = abs(ham(i,j,k))/
     &              (HAMADM_HAMADMABS + abs(16.0D0*pi*m_rho))
                  end if


c                 Calculate the Momentum constraints
c                 ----------------------------------

c                 Geometric piece.

#include "EinsteinBase/ADMMacros/src/macro/MOMXADM_guts.h"
#include "EinsteinBase/ADMMacros/src/macro/MOMYADM_guts.h"
#include "EinsteinBase/ADMMacros/src/macro/MOMZADM_guts.h"

c                 Add matter terms: - 8*pi*S_i
c
c                 with S_i defined as:
c
c                 S_i = - g_{ia} n_b T^{ab} = - g_i^a n^b T_{ab} 
c                     = - (T_{i0} - beta^j T_{ij})/alpha

                  m_sx = - ialp*Ttx
                  m_sy = - ialp*Tty
                  m_sz = - ialp*Ttz

                  if (shift_state .eq. 1) then

                     m_sx = m_sx + ialp
     &                    *(betax(i,j,k)*Txx
     &                    + betay(i,j,k)*Txy
     &                    + betaz(i,j,k)*Txz)

                     m_sy = m_sy + ialp
     &                    *(betax(i,j,k)*Txy
     &                    + betay(i,j,k)*Tyy
     &                    + betaz(i,j,k)*Tyz)
                     
                     m_sz = m_sz + ialp
     &                    *(betax(i,j,k)*Txz
     &                    + betay(i,j,k)*Tyz
     &                    + betaz(i,j,k)*Tzz)
                     
                  end if
                  
                  momx(i,j,k) = MOMXADM_MOMXADM - 8.0D0*pi*m_sx
                  momy(i,j,k) = MOMYADM_MOMYADM - 8.0D0*pi*m_sy
                  momz(i,j,k) = MOMZADM_MOMZADM - 8.0D0*pi*m_sz

               else

                  ham(i,j,k) = excised_value
                  momx(i,j,k) = excised_value
                  momy(i,j,k) = excised_value
                  momz(i,j,k) = excised_value

               end if

            end do
         end do 
      end do

#include "EinsteinBase/ADMMacros/src/macro/DETG_undefine.h"
#include "EinsteinBase/ADMMacros/src/macro/UPPERMET_undefine.h"
#include "EinsteinBase/ADMMacros/src/macro/HAMADM_undefine.h"
#include "EinsteinBase/ADMMacros/src/macro/MOMXADM_undefine.h"
#include "EinsteinBase/ADMMacros/src/macro/MOMYADM_undefine.h"
#include "EinsteinBase/ADMMacros/src/macro/MOMZADM_undefine.h"

      return
      end

      

      subroutine ADMConstraints_Boundaries(CCTK_ARGUMENTS)
      
      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      CCTK_INT, parameter :: izero = 0
      integer,  parameter :: ik = kind (izero)

c     Return code from Cactus sync routine and boundary conditions.

      integer ierr

c     Apply flat boundary conditions at outer boundaries.

      if (CCTK_EQUALS(bound,"flat")) then
         ierr = Boundary_SelectGroupForBC(cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik,
     $        "admconstraints::hamiltonian", "Flat")
         ierr = Boundary_SelectGroupForBC(cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik,
     $        "admconstraints::normalized_hamiltonian", "Flat")
         ierr = Boundary_SelectGroupForBC(cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik,
     $        "admconstraints::momentum", "Flat")
      else
         ierr = Boundary_SelectGroupForBC(cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik,
     $        "admconstraints::hamiltonian", "None")
         ierr = Boundary_SelectGroupForBC(cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik,
     $        "admconstraints::normalized_hamiltonian", "None")
         ierr = Boundary_SelectGroupForBC(cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik,
     $        "admconstraints::momentum", "None")
      end if

      return
      end 
