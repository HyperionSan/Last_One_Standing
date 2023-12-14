!/*@@
!  @file      MoncriefQ.F90
!  @date      unknown
!  @author    unknown
!  @desc
!             Compute Regge Wheeler quantities and from them the Moncrief
!             Qeven, Qodd functions.
!  @enddesc
!  @@*/


#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

!/*@@
!  @routine    WavExtrCPM_MoncriefQ
!  @date       unknown
!  @author     unknown
!  @desc
!
!              Compute Regge Wheeler quantities and from them the Moncrief
!              Qeven, Qodd functions.
!  @enddesc
!@@*/
subroutine WavExtrCPM_MoncriefQ(CCTK_ARGUMENTS)

  use WavExtrCPMConstants

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, il, im
  CCTK_INT :: status, istat
  integer ::  ierr, sumhandle
! SHH CHANGE
  CCTK_REAL :: st, st2, inv_st, inv_st2,inv_st3, ist
  CCTK_INT,dimension(2) :: lsh
  integer :: factorial

! SHH CHANGE
!  CCTK_REAL :: fac_h1, fac_H2, fac_G, fac_K, &
!               fac_c1, fac_c2, fac_dG, fac_dK, fac_dc2

    CCTK_REAL :: fac_ht, fac_dr_ht, fac_dt_hr
    CCTK_REAL :: fac_hr, fac_h2, fac_dr_h2
    CCTK_REAL :: fac_hrr, fac_jr, fac_K, fac_G, fac_dr_K
    CCTK_REAL :: fac_htr, fac_jt, fac_dr_jt, fac_dt_jr, fac_dt_K, fac_dt_G

!  CCTK_REAL :: gththcomb, gthphicomb, gphiphicomb

  CCTK_REAL :: Lambda, la, mass

    CCTK_REAL,dimension(2) :: Ylm,Y1,Y2,Y3,Y4

! SHH CHANGE
! These are all quantities with lowered indices
    CCTK_REAL,dimension(2) :: Yth, Yphi, Ythth, Ythphi, Yphiphi
    CCTK_REAL,dimension(2) :: Xth, Xphi, Xthth, Xthphi, Xphiphi

! SHH CHANGE
!                            h1,H2,K,G,c1,c2,dG,dK,dc2, &
!                            ih1,iH2,iK,iG,ic1,ic2,idG,idK,idc2

    CCTK_REAL,dimension(2) :: ht, dr_ht, dt_hr
    CCTK_REAL,dimension(2) :: hr, h2, dr_h2
    CCTK_REAL,dimension(2) :: hrr, jr, K, G, dr_K
    CCTK_REAL,dimension(2) :: htr, jt, dr_jt, dt_jr, dt_K, dt_G



  CCTK_REAL :: dtheta, dphi, dTh_dPhi

  integer :: num_out_vals, num_in_fields, minus_one
!  CCTK_REAL,dimension(18) :: out_vals, local_reduced_vals(18)
! SHH CHANGE
    CCTK_REAL,dimension(34) :: out_vals, local_reduced_vals(34)

  CCTK_REAL :: lapse,cor_angle
  CCTK_INT :: marr

  CCTK_REAL :: h_fac

  character(len=80) :: infoline

! _________________________________________________________________

  if (verbose>4) &
    call CCTK_INFO("Compute Regge Wheeler quantities and from them Qeven, Qodd")

  if (do_nothing == 1) &
    return

  if (cctk_iteration .ne. 0) then
    if (mod(cctk_iteration,my_out_every_det(current_detector)).ne.0) then
      if (verbose>2) call CCTK_INFO("No time for this detector")
      return
    end if
  end if

  if (calc_when_necessary .eq. 1) then
    if (cctk_time .lt. current_detector_radius-50) then
      if (verbose>2) call CCTK_INFO("No time for this detector")
      return
    endif
    call CCTK_IsFunctionAliased(istat, "MergerHandler_WeHaveMerger")
    if (istat .eq. 1) then
      if (MergerHandler_WeHaveMerger() .eq. 1) then
        if (cctk_time .gt. MergerHandler_MergerTime()+current_detector_radius+ringdown_margin) then
          if (verbose>2) call CCTK_INFO("No time for this detector")
          return
        endif
      endif
    endif
  end if

  ! local shape of grid arrays on sphere
  call CCTK_GrouplshGN(status, cctkGH, 2, lsh, "WaveExtractCPM::surface_arrays")
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, "cannot get local size for surface arrays" )
  end if

  dtheta = ctheta(2,1) - ctheta(1,1)
  dphi = cphi(1,2) - cphi(1,1)

  if (cartoon .ne. 0) then
    dphi = two*pi
  end if

  dTh_dPhi= dtheta*dphi

  ! l-mode setup
  if (CCTK_EQUALS(mode_type,"specific mode")) then
    l_min = l_mode ; l_max = l_mode
  else if (CCTK_EQUALS(mode_type,"all modes")) then
    l_min = 2 ; l_max = l_mode
  end if


  if (cartoon .ne. 0) then
    l_step =2
    m_step =2
  else if (CCTK_EQUALS(domain,"full")) then
    l_step =1
    m_step =1
  else if (CCTK_EQUALS(domain,"octant")) then
    l_step =2
    m_step =2
  ! FIXME : look at symmetries
  else if (CCTK_EQUALS(domain,"quadrant")) then
    l_step =1
    m_step =1
  ! FIXME : REALLY ALL MODES?? Does the bitant symmetry kill some, ie odd modes ???
  else if (CCTK_EQUALS(domain,"bitant")) then
    l_step =1
    m_step =1
  end if

  if (verbose > 3) then
    write(infoline,'(A29,I2,I2,I2)') 'mode setup: l [min,max,step]:',l_min,l_max,l_step
    call CCTK_INFO(infoline)
  end if

  ! Loop over l-modes
  loop_l: do il = l_min,l_max,l_step

    ! m-mode setup (depends on l_mode, ie m_mode <= l_mode)
    if (CCTK_EQUALS(mode_type,"specific mode")) then
      m_min = -m_mode ; m_max = m_mode
    else if(CCTK_EQUALS(mode_type,"all modes")) then
      m_min = -m_mode ; m_max = m_mode
      if (m_max>il) then
        m_max=il
      end if
      if (m_min<-il) then
        m_min=-il
      end if
    end if

    ! Factors independent of angular coordinates and of m-mode, but not l_mode

! SHH CHANGE
    fac_ht  = one/dble(il*(il+1))
    fac_hr  = dri_drsch/dble(il*(il+1))
    fac_h2  = two/dble(il*(il+1)*(il-1)*(il+2))

    fac_dr_ht = fac_ht
    fac_dr_h2 = fac_h2
    fac_dt_hr = fac_hr

    fac_htr = dri_drsch
    fac_hrr = dri_drsch**2
    fac_jt  = one/dble(il*(il+1))
    fac_jr  = dri_drsch/dble(il*(il+1))
    fac_K   = half/rsch**2
    fac_G   = two/(rsch**2 * dble(il*(il+1)*(il-1)*(il+2)))

    fac_dr_jt   = fac_jt
    fac_dr_K    = fac_K
    fac_dt_jr   = fac_jr
    fac_dt_G    = fac_G
    fac_dt_K    = fac_K


    call CCTK_TimerStart(ierr,"MoncriefQ_CPM")

    if (verbose > 3) then
      write(infoline,'(A29,I2,I2,I2)') '            m [min,max,step]:',m_min,m_max,m_step
      call CCTK_INFO(infoline)
    end if

    ! Get sum reduction operator handle
    call CCTK_ReductionArrayHandle ( sum_handle, 'sum' )
    if ( sum_handle .lt. 0 ) then
      call CCTK_WARN(0,'Could not obtain a handle for sum reduction')
    end if

    ! Loop over m-modes
    loop_m: do im = m_min, m_max, m_step
      ! Real parts of the Regge-Wheeler variables
      loop_phi1: do j = 1, lsh(2)
        loop_theta1: do i = 1, lsh(1)

!   SHH CHANGE
            st      = sintheta(i,j)
            
           ist = one/st

!   In terms of Martel and Poisson's (2005) spherical harmonics:
!
!   Y1 = Y_\th = X_\phi / sin \th
!   Y2 = Y_\phi = -sin \th X_\th
!   Y3 = Y_{\th \th} = - Y_{\phi \phi} / sin^2 \th = X_{\th \phi} / sin \th
!   Y4 = Y_{\th \phi} = - sin \th X_{\th \th} = X_{\phi \phi} / sin \th

          call WavExtrCPM_spher_harm_combs(ctheta(i,j),cphi(i,j),il,im,Ylm,Y1,Y2,Y3,Y4)

            Yth(1)     = Y1(1)
            Yphi(1)    = Y2(1)

            Xth(1)     = - Yphi(1) / st
            Xphi(1)    = Yth(1) * st

            Ythth(1)   = Y3(1)
            Yphiphi(1) = - Ythth(1) * st**2
            Ythphi(1)  = Y4(1)

            Xthth(1)   = - Ythphi(1) / st
            Xphiphi(1) = - Xthth(1) * st**2
            Xthphi(1)  = Ythth(1) * st

!            if (verbose >4) then
            !        print*,'Real Ylm'
!            print*,'Ylm',Ylm(1)
!            print*,'Y1',Y1(1)
!            print*,'Y2',Y2(1)
!            print*,'Y3',Y3(1)
!            print*,'Y4',Y4(1)
!            print*,'sinth',st

            !        print*,'Img Ylm'
            !        print*,'Ylm',Ylm(2)
            !        print*,'Y1',Y1(2)
            !        print*,'Y2',Y2(2)
            !        print*,'Y3',Y3(2)
            !        print*,'Y4',Y4(2)

!            end if

!            if (verbose >4) then
!            print*,'Real Ylm'
!            print*,'Ylm',Ylm(1)
!            print*,'Yth',Yth(1)
!            print*,'Yphi',Yphi(1)
!            print*,'Xth',Xth(1)
!            print*,'Xphi',Xphi(1)
!            print*,'Ythth',Ythth(1)
!            print*,'Ythphi',Ythphi(1)
!            print*,'Yphiphi',Yphiphi(1)
!            print*,'Xthth',Xthth(1)
!            print*,'Xthphi',Xthphi(1)
!            print*,'Xphiphi',Xphiphi(1)
!            end if

            if (verbose >4 .and. i .eq. 25 .and. j .eq. 35) then
                print*,'Metric vals'
                print*,'grr',grr(i,j)
                print*,'grth',grth(i,j)
                print*,'grphi',grphi(i,j)
                print*,'gthth',gthth(i,j)
                print*,'gthphi',gthphi(i,j)
                print*,'gphiphi',gphiphi(i,j)
                print*,'dr_gthth',dr_gthth(i,j)
                print*,'dr_gthphi',dr_gthphi(i,j)
                print*,'dr_gphiphi',dr_gphiphi(i,j)
                print*,'dt_grth',dt_grth(i,j)
                print*,'dt_grphi',dt_grphi(i,j)
                print*,'dt_gthth',dt_gthth(i,j)
                print*,'dt_gpthphi',dt_gthphi(i,j)
                print*,'dt_gphiphi',dt_gphiphi(i,j)
            end if

            htRe_i(i,j) = st*(gtth(i,j) * Xth(1) + gtphi(i,j) * Xphi(1)/st**2)
            hrRe_i(i,j) = st*(grth(i,j) * Xth(1) + grphi(i,j) * Xphi(1)/st**2)
            h2Re_i(i,j) = st*(gthth(i,j) * Xthth(1) + 2.0 * gthphi(i,j) * Xthphi(1)/st**2 + gphiphi(i,j) * Xphiphi(1)/st**4)

            dr_htRe_i(i,j) = st*((dr_gtth(i,j) * Xth(1) + dr_gtphi(i,j) * Xphi(1)/st**2) * dri_drsch)
            dr_h2Re_i(i,j) = st*((  dr_gthth(i,j) * Xthth(1) &
                              + 2.0 * dr_gthphi(i,j) * Xthphi(1)/st**2 &
                              + dr_gphiphi(i,j) * Xphiphi(1)/st**4 ) * dri_drsch)

            dt_hrRe_i(i,j) = st*(dt_grth(i,j) * Xth(1) + dt_grphi(i,j) * Xphi(1)/st**2)

            htrRe_i(i,j) = st*(gtr(i,j) * Ylm(1))
            hrrRe_i(i,j) = st*(grr(i,j) * Ylm(1))
            jtRe_i(i,j)  = st*(gtth(i,j) * Yth(1) + gtphi(i,j) * Yphi(1)/st**2)
            jrRe_i(i,j)  = st*(grth(i,j) * Yth(1) + grphi(i,j) * Yphi(1)/st**2)

            G_Re_i(i,j)  = st*(gthth(i,j) * Ythth(1) + gphiphi(i,j) * Yphiphi(1) / st**4 + 2.0 * gthphi(i,j) * Ythphi(1)/st**2)

            K_Re_i(i,j)  = st*( ( gthth(i,j) + gphiphi(i,j) / st**2 ) * Ylm(1))

            dr_jtRe_i(i,j)  = st*(dr_gtth(i,j) * Yth(1) + dr_gtphi(i,j) * Yphi(1)/st**2) * dri_drsch
            dt_jrRe_i(i,j)  = st*(dt_grth(i,j) * Yth(1) + dt_grphi(i,j) * Yphi(1)/st**2)
            dt_G_Re_i(i,j)  = st*(dt_gthth(i,j) * Ythth(1) + 2.0 * dt_gthphi(i,j) * Ythphi(1)/st**2 + dt_gphiphi(i,j) * Yphiphi(1) / st**4)
            dt_K_Re_i(i,j)  = st*(( dt_gthth(i,j) + dt_gphiphi(i,j) / st**2 ) * Ylm(1))
            dr_K_Re_i(i,j)  = -two/rsch * K_Re_i(i,j) + st*(( dr_gthth(i,j) + dr_gphiphi(i,j) / st**2 ) * dri_drsch * Ylm(1))


          ! for m!=0 case we have imaginary part as well from the Ylm's
          if (im.ne.0) then
            ! switch signs - stupid convention  <==== Isn't is just that we integrate with the c.c. of Ylm?  SHH

            Ylm(2)  = -Ylm(2)
            Y1(2)   = -Y1(2)
            Y2(2)   = -Y2(2)
            Y3(2)   = -Y3(2)
            Y4(2)   = -Y4(2)

            Yth(2)     = Y1(2)
            Yphi(2)    = Y2(2)

            Xth(2)     = - Yphi(2) / st
            Xphi(2)    = Yth(2) * st

            Ythth(2)   = Y3(2)
            Yphiphi(2) = - Ythth(2) * st**2
            Ythphi(2)  = Y4(2)

            Xthth(2)   = - Ythphi(2) / st
            Xphiphi(2) = - Xthth(2) * st**2
            Xthphi(2)  = Ythth(2) * st

            htIm_i(i,j) = st*(gtth(i,j) * Xth(2) + gtphi(i,j) * Xphi(2)/st**2)
            hrIm_i(i,j) = st*(grth(i,j) * Xth(2) + grphi(i,j) * Xphi(2)/st**2)
            h2Im_i(i,j) = st*(gthth(i,j) * Xthth(2) + 2.0 * gthphi(i,j) * Xthphi(2)/st**2 + gphiphi(i,j) * Xphiphi(2)/st**4)

            dr_htIm_i(i,j) = st*((dr_gtth(i,j) * Xth(2) + dr_gtphi(i,j) * Xphi(2)/st**2) * dri_drsch)
            dr_h2Im_i(i,j) = st*((  dr_gthth(i,j) * Xthth(2) &
                            + 2.0 * dr_gthphi(i,j) * Xthphi(2)/st**2 &
                            + dr_gphiphi(i,j) * Xphiphi(2)/st**4 ) * dri_drsch)
!            dr_h2Im_i(i,j) = st*((  dr_gthth(i,j) * Xthth(2) &
!                                + 2.0 * dr_gthphi(i,j) * Xthphi(2)/st**2 &
!                                + dr_gphiphi(i,j) * Xphiphi(2)/st**4 ))

            dt_hrIm_i(i,j) = st*(dt_grth(i,j) * Xth(2) + dt_grphi(i,j) * Xphi(2)/st**2)

            htrIm_i(i,j) = st*(gtr(i,j) * Ylm(2))
            hrrIm_i(i,j) = st*(grr(i,j) * Ylm(2))
            jtIm_i(i,j)  = st*(gtth(i,j) * Yth(2) + gtphi(i,j) * Yphi(2)/st**2)
            jrIm_i(i,j)  = st*(grth(i,j) * Yth(2) + grphi(i,j) * Yphi(2)/st**2)

            G_Im_i(i,j)  = st*(gthth(i,j) * Ythth(2) + gphiphi(i,j) * Yphiphi(2) / st**4 + 2.0 * gthphi(i,j) * Ythphi(2)/st**2)

            K_Im_i(i,j)  = st*(( gthth(i,j) + gphiphi(i,j) / st**2 ) * Ylm(2))

            dr_jtIm_i(i,j)  = st*(dr_gtth(i,j) * Yth(2) + dr_gtphi(i,j) * Yphi(2)/st**2) * dri_drsch
            dt_jrIm_i(i,j)  = st*(dt_grth(i,j) * Yth(2) + dt_grphi(i,j) * Yphi(2)/st**2)
            dt_G_Im_i(i,j)  = st*(dt_gthth(i,j) * Ythth(2) + 2.0 * dt_gthphi(i,j) * Ythphi(2)/st**2 + dt_gphiphi(i,j) * Yphiphi(2) / st**4)
            dt_K_Im_i(i,j)  = st*(( dt_gthth(i,j) + dt_gphiphi(i,j) / st**2 ) * Ylm(2))
            dr_K_Im_i(i,j)  = st*(-two/rsch * K_Im_i(i,j) + ( dr_gthth(i,j) + dr_gphiphi(i,j) / st**2 ) * dri_drsch * Ylm(2))


          end if
        end do loop_theta1
      end do loop_phi1

      ! Integrations over the 2-sphere
      ! Note the abscence of sintheta which is already included
      ! in the above expressions

        int_tmp1 = sym_factor*weights*dTh_dPhi* htRe_i
        int_tmp2 = sym_factor*weights*dTh_dPhi* hrRe_i
        int_tmp3 = sym_factor*weights*dTh_dPhi* h2Re_i
        int_tmp4 = sym_factor*weights*dTh_dPhi* dr_htRe_i
        int_tmp5 = sym_factor*weights*dTh_dPhi* dr_h2Re_i
        int_tmp6 = sym_factor*weights*dTh_dPhi* dt_hrRe_i

        int_tmp7  = sym_factor*weights*dTh_dPhi* htrRe_i
        int_tmp8  = sym_factor*weights*dTh_dPhi* hrrRe_i
        int_tmp9  = sym_factor*weights*dTh_dPhi* jtRe_i
        int_tmp10 = sym_factor*weights*dTh_dPhi* jrRe_i
        int_tmp11 = sym_factor*weights*dTh_dPhi* G_Re_i
        int_tmp12 = sym_factor*weights*dTh_dPhi* K_Re_i
        int_tmp13 = sym_factor*weights*dTh_dPhi* dr_jtRe_i
        int_tmp14 = sym_factor*weights*dTh_dPhi* dt_jrRe_i
        int_tmp15 = sym_factor*weights*dTh_dPhi* dt_G_Re_i
        int_tmp16 = sym_factor*weights*dTh_dPhi* dt_K_Re_i
        int_tmp17 = sym_factor*weights*dTh_dPhi* dr_K_Re_i


        local_reduced_vals(1)  = sum(int_tmp1,weights.gt.1.0e-15)
        local_reduced_vals(2)  = sum(int_tmp2,weights.gt.1.0e-15)
        local_reduced_vals(3)  = sum(int_tmp3,weights.gt.1.0e-15)
        local_reduced_vals(4)  = sum(int_tmp4,weights.gt.1.0e-15)
        local_reduced_vals(5)  = sum(int_tmp5,weights.gt.1.0e-15)
        local_reduced_vals(6)  = sum(int_tmp6,weights.gt.1.0e-15)

        local_reduced_vals(7)  = sum(int_tmp7,weights.gt.1.0e-15)
        local_reduced_vals(8)  = sum(int_tmp8,weights.gt.1.0e-15)
        local_reduced_vals(9)  = sum(int_tmp9,weights.gt.1.0e-15)
        local_reduced_vals(10) = sum(int_tmp10,weights.gt.1.0e-15)
        local_reduced_vals(11) = sum(int_tmp11,weights.gt.1.0e-15)
        local_reduced_vals(12) = sum(int_tmp12,weights.gt.1.0e-15)
        local_reduced_vals(13) = sum(int_tmp13,weights.gt.1.0e-15)
        local_reduced_vals(14) = sum(int_tmp14,weights.gt.1.0e-15)
        local_reduced_vals(15) = sum(int_tmp15,weights.gt.1.0e-15)
        local_reduced_vals(16) = sum(int_tmp16,weights.gt.1.0e-15)
        local_reduced_vals(17) = sum(int_tmp17,weights.gt.1.0e-15)

      num_out_vals=1
      minus_one = -1
      sumhandle = sum_handle ! i.e., convert from CCTK_INT to integer

      if (im.eq.0) then
            num_in_fields=17

      else
            num_in_fields=34

        int_tmp18 = sym_factor*weights*dTh_dPhi* htIm_i
        int_tmp19 = sym_factor*weights*dTh_dPhi* hrIm_i
        int_tmp20 = sym_factor*weights*dTh_dPhi* h2Im_i
        int_tmp21 = sym_factor*weights*dTh_dPhi* dr_htIm_i
        int_tmp22 = sym_factor*weights*dTh_dPhi* dr_h2Im_i
        int_tmp23 = sym_factor*weights*dTh_dPhi* dt_hrIm_i

        int_tmp24 = sym_factor*weights*dTh_dPhi* htrIm_i
        int_tmp25 = sym_factor*weights*dTh_dPhi* hrrIm_i
        int_tmp26 = sym_factor*weights*dTh_dPhi* jtIm_i
        int_tmp27 = sym_factor*weights*dTh_dPhi* jrIm_i
        int_tmp28 = sym_factor*weights*dTh_dPhi* G_Im_i
        int_tmp29 = sym_factor*weights*dTh_dPhi* K_Im_i
        int_tmp30 = sym_factor*weights*dTh_dPhi* dr_jtIm_i
        int_tmp31 = sym_factor*weights*dTh_dPhi* dt_jrIm_i
        int_tmp32 = sym_factor*weights*dTh_dPhi* dt_G_Im_i
        int_tmp33 = sym_factor*weights*dTh_dPhi* dt_K_Im_i
        int_tmp34 = sym_factor*weights*dTh_dPhi* dr_K_Im_i


        local_reduced_vals(18) = sum(int_tmp18,weights.gt.1.0e-15)
        local_reduced_vals(19) = sum(int_tmp19,weights.gt.1.0e-15)
        local_reduced_vals(20) = sum(int_tmp20,weights.gt.1.0e-15)
        local_reduced_vals(21) = sum(int_tmp21,weights.gt.1.0e-15)
        local_reduced_vals(22) = sum(int_tmp22,weights.gt.1.0e-15)
        local_reduced_vals(23) = sum(int_tmp23,weights.gt.1.0e-15)

        local_reduced_vals(24) = sum(int_tmp24,weights.gt.1.0e-15)
        local_reduced_vals(25) = sum(int_tmp25,weights.gt.1.0e-15)
        local_reduced_vals(26) = sum(int_tmp26,weights.gt.1.0e-15)
        local_reduced_vals(27) = sum(int_tmp27,weights.gt.1.0e-15)
        local_reduced_vals(28) = sum(int_tmp28,weights.gt.1.0e-15)
        local_reduced_vals(29) = sum(int_tmp29,weights.gt.1.0e-15)
        local_reduced_vals(30) = sum(int_tmp30,weights.gt.1.0e-15)
        local_reduced_vals(31) = sum(int_tmp31,weights.gt.1.0e-15)
        local_reduced_vals(32) = sum(int_tmp32,weights.gt.1.0e-15)
        local_reduced_vals(33) = sum(int_tmp33,weights.gt.1.0e-15)
        local_reduced_vals(34) = sum(int_tmp34,weights.gt.1.0e-15)

      end if

      call CCTK_ReduceLocArrayToArray1D(ierr, cctkGH, minus_one,&
                        sumhandle, local_reduced_vals(1:num_in_fields),&
                        out_vals(1:num_in_fields), num_in_fields,&
                        CCTK_VARIABLE_REAL)

      if (ierr.ne.0) then
        call CCTK_WARN(1,"The reduction of the MoncriefQ integrands failed!")
      end if

        ht(1)    = out_vals(1)
        hr(1)    = out_vals(2)
        h2(1)    = out_vals(3)
        dr_ht(1) = out_vals(4)
        dr_h2(1) = out_vals(5)
        dt_hr(1) = out_vals(6)

        htr(1)   = out_vals(7)
        hrr(1)   = out_vals(8)
        jt(1)    = out_vals(9)
        jr(1)    = out_vals(10)
        G(1)     = out_vals(11)
        K(1)     = out_vals(12)
        dr_jt(1) = out_vals(13)
        dt_jr(1) = out_vals(14)
        dt_G(1)  = out_vals(15)
        dt_K(1)  = out_vals(16)
        dr_K(1)  = out_vals(17)

        ht(1)    = fac_ht * ht(1)
        hr(1)    = fac_hr * hr(1)
        h2(1)    = fac_h2 * h2(1)
        dr_ht(1) = fac_dr_ht * dr_ht(1)
        dr_h2(1) = fac_dr_h2 * dr_h2(1)
        dt_hr(1) = fac_dt_hr * dt_hr(1)

        htr(1)   = fac_htr * htr(1)
        hrr(1)   = fac_hrr * hrr(1)
        jt(1)    = fac_jt  * jt(1)
        jr(1)    = fac_jr  * jr(1)
        G(1)     = fac_G   * G(1)
        K(1)     = fac_K   * K(1)
        dr_jt(1) = fac_dr_jt * dr_jt(1)
        dt_jr(1) = fac_dt_jr * dt_jr(1)
        dt_G(1)  = fac_dt_G  * dt_G(1)
        dt_K(1)  = fac_dt_K  * dt_K(1)
        dr_K(1)  = fac_dr_K  * dr_K(1)



    if (im.ne.0) then

            ht(2)    = out_vals(18)
            hr(2)    = out_vals(19)
            h2(2)    = out_vals(20)
            dr_ht(2) = out_vals(21)
            dr_h2(2) = out_vals(22)
            dt_hr(2) = out_vals(23)

            htr(2)   = out_vals(24)
            hrr(2)   = out_vals(25)
            jt(2)    = out_vals(26)
            jr(2)    = out_vals(27)
            G(2)     = out_vals(28)
            K(2)     = out_vals(29)
            dr_jt(2) = out_vals(30)
            dt_jr(2) = out_vals(31)
            dt_G(2)  = out_vals(32)
            dt_K(2)  = out_vals(33)
            dr_K(2)  = out_vals(34)

            ht(2)    = fac_ht * ht(2)
            hr(2)    = fac_hr * hr(2)
            h2(2)    = fac_h2 * h2(2)
            dr_ht(2) = fac_dr_ht * dr_ht(2)
            dr_h2(2) = fac_dr_h2 * dr_h2(2)
            dt_hr(2) = fac_dt_hr * dt_hr(2)

            htr(2)   = fac_htr * htr(2)
            hrr(2)   = fac_hrr * hrr(2)
            jt(2)    = fac_jt  * jt(2)
            jr(2)    = fac_jr  * jr(2)
            G(2)     = fac_G   * G(2)
            K(2)     = fac_K   * K(2)
            dr_jt(2) = fac_dr_jt * dr_jt(2)
            dt_jr(2) = fac_dt_jr * dt_jr(2)
            dt_G(2)  = fac_dt_G  * dt_G(2)
            dt_K(2)  = fac_dt_K  * dt_K(2)
            dr_K(2)  = fac_dr_K  * dr_K(2)

      end if

        if (verbose >4) then

        print*,'Real Quantities'
        print*,'ht',ht(1)
        print*,'hr',hr(1)
        print*,'h2',h2(1)
        print*,'dr_ht', dr_ht(1)
        print*,'dt_hr', dt_hr(1)
        print*,'dr_h2', dr_h2(1)

        print*,'htr', htr(1)
        print*,'hrr', hrr(1)
        print*,'jt', jt(1)
        print*,'jr', jr(1)
        print*,'K', K(1)
        print*,'G', G(1)
        print*,'dr_jt', dr_jt(1)
        print*,'dt_jr', dt_jr(1)
        print*,'dt_K', dt_K(1)
        print*,'dt_G', dt_G(1)
        print*,'dr_K', dr_K(1)

        print*,'Imaginary Quantities'
        print*,'ht',ht(2)
        print*,'hr',hr(2)
        print*,'h2',h2(2)
        print*,'dr_ht', dr_ht(2)
        print*,'dt_hr', dt_hr(2)
        print*,'dr_h2', dr_h2(2)

        print*,'htr', htr(2)
        print*,'hrr', hrr(2)
        print*,'jt', jt(2)
        print*,'jr', jr(2)
        print*,'K', K(2)
        print*,'G', G(2)
        print*,'dr_jt', dr_jt(2)
        print*,'dt_jr', dt_jr(2)
        print*,'dt_K', dt_K(2)
        print*,'dt_G', dt_G(2)
        print*,'dr_K', dr_K(2)

        end if

      ! Moncrief Q
! SHH CHANGE - note that my Lambda is half this Lambda
!      Lambda = dble((il-1)*(il+2))+three*(one-S_factor)
        la = dble((il-1)*(il+2))/two
        mass = (one-S_factor) * rsch / two
        Lambda = la + three*(one-S_factor)/two

      ! m index into array : fortran: 1,2,3 index -> -m_max,-m_max+1,...
        marr=-m_min+im +1
        
      if (verbose>3) then
        print*,''
        print*,'il', il
        print*,'l_mode', l_mode
        print*,'marr', marr
        print*,'m_mode', m_mode
        print*,'2*m_mode+1', 2*m_mode+1
      end if
      if (il<1 .or. il>l_mode .or. marr<1 .or. marr>2*m_mode+1) call CCTK_WARN (0, "internal error")

! SHH CHANGE
! The purpose of all these changes is to compute Psi_odd, the Cunningham-Price-Moncrief 
! master function, which is twice the time integral of the original Regge-Wheeler function
! (equivalently Qodd).  We also print 2*Qodd, the time deriv of Psi_odd.
! The even-parity should be the same as before, except that we print out the
! time derivative as well.
      !
      
      if (verbose>2) then
         print*,'f', S_factor
         print*,'rsch', rsch
         print*,'Schw_M', Schwarzschild_Mass
         print*,'mass', mass
      end if

        Psi_odd_Re(il,marr)     = rsch/la * ( dr_ht(1) - dt_hr(1) - two/rsch * ht(1) )
        dt_Psi_odd_Re(il,marr)  = S_factor/rsch * ( two*hr(1) + two/rsch * h2(1) - dr_h2(1) )

        Psi_even_Re(il,marr)    =   rsch * G(1) - two*S_factor/Lambda * jr(1) &
                                  + rsch / (la+one) * (K(1) + S_factor/Lambda * (S_factor*hrr(1) - rsch*dr_K(1)))

        dt_Psi_even_Re(il,marr) =   rsch * dt_G(1) + one/Lambda * ( - S_factor * htr(1) &
                                                                    - two*mass/(rsch*rsch) * jt(1) &
                                                                    + S_factor * dr_jt(1) &
                                                                    + rsch * dt_K(1) &
                                                                    - S_factor * dt_jr(1) )

    ! For comparisons:

        Psi_odd_Re(il,marr)  = Psi_odd_Re(il,marr) * sqrt(two*dble((il+2)*(il+1)*il*(il-1))) / two
        dt_Psi_odd_Re(il,marr)  = dt_Psi_odd_Re(il,marr) * sqrt(two*dble((il+2)*(il+1)*il*(il-1))) / two

        Psi_even_Re(il,marr) = Psi_even_Re(il,marr) * (il+1)*il*sqrt(two*dble((il-1)*(il+2)) / dble(il*(il+1)))/two

        dt_Psi_even_Re(il,marr) = dt_Psi_even_Re(il,marr) * (il+1)*il*sqrt(two*dble((il-1)*(il+2)) / dble(il*(il+1)))/two

!      Qodd_Re(il,marr) = sqrt(two*dble((il+2)*(il+1)*il*(il-1)))* &
!                      S_factor/rsch*(c1(1)+half*(dc2(1)-two/rsch*c2(1)))


!      Qeven_Re(il,marr) = one/Lambda*sqrt(two*dble((il-1)*(il+2))/ &
!                          dble(il*(il+1)))*( dble(il*(il+1))*S_factor* &
!              (rsch**2*dG(1)-two*h1(1))+two*rsch*S_factor*(H2(1)-rsch*dK(1)) &
!                             +Lambda*rsch*K(1) )

      if (im.ne.0) then

        Psi_odd_Im(il,marr)     = rsch/la * ( dr_ht(2) - dt_hr(2) - two/rsch * ht(2) )
        dt_Psi_odd_Im(il,marr)  = S_factor/rsch * ( two*hr(2) + two/rsch * h2(2) - dr_h2(2) )

        Psi_even_Im(il,marr)    =   rsch * G(2) - two*S_factor/Lambda * jr(2) &
                                  + rsch / (la+one) * (K(2) + S_factor/Lambda * (S_factor*hrr(2) - rsch*dr_K(2)))

        dt_Psi_even_Im(il,marr) =   rsch * dt_G(2) + one/Lambda * ( - S_factor * htr(2) &
                                                                    - two*mass/(rsch*rsch) * jt(2) &
                                                                    + S_factor * dr_jt(2) &
                                                                    + rsch * dt_K(2) &
                                                                    - S_factor * dt_jr(2) )

! For Comparison purposes ... different normalization

        Psi_odd_Im(il,marr)  = Psi_odd_Im(il,marr) * sqrt(two*dble((il+2)*(il+1)*il*(il-1))) / two
        dt_Psi_odd_Im(il,marr)  = dt_Psi_odd_Im(il,marr) * sqrt(two*dble((il+2)*(il+1)*il*(il-1))) / two

        Psi_even_Im(il,marr) = Psi_even_Im(il,marr) * (il+1)*il*sqrt(two*dble((il-1)*(il+2)) / dble(il*(il+1)))/two

        dt_Psi_even_Im(il,marr) = dt_Psi_even_Im(il,marr) * (il+1)*il*sqrt(two*dble((il-1)*(il+2)) / dble(il*(il+1)))/two


!        Qodd_Im(il,marr) = sqrt(two*dble((il+2)*(il+1)*il*(il-1)))* &
!                        S_factor/rsch*(c1(2)+half*(dc2(2)-two/rsch*c2(2)))


 !       Qeven_Im(il,marr) = one/Lambda*sqrt(two*dble((il-1)*(il+2))/ &
 !                           dble(il*(il+1)))*( dble(il*(il+1))*S_factor* &
 !               (rsch**2*dG(2)-two*h1(2))+two*rsch*S_factor*(H2(2)-rsch*dK(2)) &
 !                              +Lambda*rsch*K(2) )

        ! Strain
        h_fac = -quarter * half * sqrt(dble(factorial(il+2))/dble(factorial(il-2)))
        
        h_Re(il,marr) = h_fac*(Psi_even_Re(il,marr) - Psi_odd_Im(il,marr))
        h_Im(il,marr) = h_fac*(Psi_even_Im(il,marr) + Psi_odd_Re(il,marr))

     end if

      if (verbose>3) then
        write(infoline,'(A9,I2,A1,I2,A1)') '  (l,m)=(',il,',',im,')'
        call CCTK_INFO(infoline)
        write(infoline,'(A10,G20.8,A7,G20.8)') '    Psi_even=', &
                         Psi_even_Re(il,marr),', Psi_odd=',Psi_odd_Re(il,marr)
        call CCTK_INFO(infoline)
        if (im.ne.0) then
          write(infoline,'(A16,G20.8,A7,G20.8)') '    imag: Psi_even=', &
                         Psi_even_Im(il,marr),', Psi_odd=',Psi_odd_Im(il,marr)
          call CCTK_INFO(infoline)
        end if
      end if

!if (verbose>3) then
!write(infoline,'(A9,I2,A1,I2,A1)') '  (l,m)=(',il,',',im,')'
!call CCTK_INFO(infoline)
!write(infoline,'(A10,G20.8,A7,G20.8)') '    Qeven=', &
!Qeven_Re(il,marr),', Qodd=',Qodd_Re(il,marr)
!call CCTK_INFO(infoline)
!if (im.ne.0) then
!write(infoline,'(A16,G20.8,A7,G20.8)') '    imag: Qeven=', &
!Qeven_Im(il,marr),', Qodd=',Qodd_Im(il,marr)
!call CCTK_INFO(infoline)
!end if
!end if

    if (current_detector<1 .or. current_detector>maximum_detector_number .or. il<1 .or. il>l_mode .or. marr<1 .or. marr>2*m_mode+1) call CCTK_WARN (0, "internal error")
        Psi_odd_Re_Array(current_detector,il,marr)      = Psi_odd_Re(il,marr)
        dt_Psi_odd_Re_Array(current_detector,il,marr)   = dt_Psi_odd_Re(il,marr)
        Psi_even_Re_Array(current_detector,il,marr)     = Psi_even_Re(il,marr)
        dt_Psi_even_Re_Array(current_detector,il,marr)  = dt_Psi_even_Re(il,marr)
       

!Qodd_Re_Array(current_detector,il,marr) = Qodd_Re(il,marr)
!Qeven_Re_Array(current_detector,il,marr) = Qeven_Re(il,marr)

      if (im.ne.0) then

        Psi_odd_Im_Array(current_detector,il,marr)      = Psi_odd_Im(il,marr)
        dt_Psi_odd_Im_Array(current_detector,il,marr)   = dt_Psi_odd_Im(il,marr)
        Psi_even_Im_Array(current_detector,il,marr)     = Psi_even_Im(il,marr)
        dt_Psi_even_Im_Array(current_detector,il,marr)  = dt_Psi_even_Im(il,marr)

!        Qodd_Im_Array(current_detector,il,marr) = Qodd_Im(il,marr)
!        Qeven_Im_Array(current_detector,il,marr) = Qeven_Im(il,marr)
      else

        Psi_odd_Im_Array(current_detector,il,marr)      = 0
        dt_Psi_odd_Im_Array(current_detector,il,marr)   = 0
        Psi_even_Im_Array(current_detector,il,marr)     = 0
        dt_Psi_even_Im_Array(current_detector,il,marr)  = 0
!        Qodd_Im_Array(current_detector,il,marr) = zero
!        Qeven_Im_Array(current_detector,il,marr) = zero
      end if
      h_Re_Array(current_detector,il,marr)  = h_Re(il,marr)
      h_Im_Array(current_detector,il,marr)  = h_Im(il,marr)
    end do loop_m
  end do loop_l

  call CCTK_TimerStop(ierr,"MoncriefQ_CPM")

  if (phicorotate .ne. 0) then
    call CCTK_INFO("reset corotation on phi at the end of MoncriefQ for next step")
    lapse=one
    cor_angle=cctk_time*rotation_omega/lapse

!   undo rotation transformation on phi itself.
    do j = 1, lsh(2)
      cphi(:,j) = cphi(:,j) + cor_angle
    end do

    sinphi = sin(cphi)
    cosphi = cos(cphi)
  end if

end subroutine WavExtrCPM_MoncriefQ



! ------------------------------------------------------------------
!
! Calculate the (l,m) spherical harmonic at given angular
! coordinates. This number is in general complex, and
!
! Ylm = Ylm(1) + i Ylm(2)
!
! where
!
!           a    ( 2 l + 1 (l-|m|)! )                      i m phi
! Ylm = (-1) SQRT( ------- -------  ) P_l|m| (cos(theta)) e
!                (   4 Pi  (l+|m|)! ) 
!
! and where
!
! a = m/2 (sign(m)+1)
!
! ------------------------------------------------------------------
subroutine WavExtrCPM_spherical_harmonic(l,m,theta,phi,Ylm)

  use WavExtrCPMConstants

  implicit none

! Input variables
  CCTK_INT :: l,m
  CCTK_REAL :: theta,phi

! Output variables
  CCTK_REAL :: Ylm(2)

! Local variables
  CCTK_INT :: i
  CCTK_REAL :: a,fac,WavExtrCPM_plgndr
! _________________________________________________________________

  fac = one
  do i = l-abs(m)+1,l+abs(m)
        fac = fac*dble(i)
  end do
  fac = one/fac

! a = (-one)**((m*ISIGN(m,1)/abs(m)+m)/2)*SQRT(dble(2*l+1)
! &    /four/Pi*fac)*WavExtrCPM_plgndr(l,abs(m),cos(theta))

  a = (-one)**max(m,0)*sqrt(dble(2*l+1)/ &
      four/Pi*fac)*WavExtrCPM_plgndr(l,abs(m),cos(theta))
  Ylm(1) = a*cos(dble(m)*phi)
  Ylm(2) = a*sin(dble(m)*phi)
end subroutine WavExtrCPM_spherical_harmonic


! __________________________________________________________________
!
! FIXME : CHECK THE LICENSE OF THIS ROUTINE
! From Numerical Recipes FIXME
!
!   Calculates the associated Legendre polynomial Plm(x).
!   Here m and l are integers satisfying 0 <= m <= l,
!   while x lies in the range -1 <= x <= 1
!
! __________________________________________________________________
function WavExtrCPM_plgndr(l,m,x)

  use WavExtrCPMConstants

  implicit none

! Input variables
  CCTK_INT,INTENT(IN) :: l,m
  CCTK_REAL,INTENT(IN) :: x

! Output variables
  CCTK_REAL :: WavExtrCPM_plgndr

! Local Variables
  CCTK_INT :: i,ll
  CCTK_REAL :: pmm,somx2,fact,pmmp1,pll

! __________________________________________________________________

  pmm = one

  if (m.gt.0) then
    somx2=sqrt((one-x)*(one+x))
    fact=one
    do i=1,m
      pmm  = -pmm*fact*somx2
      fact = fact+two
    end do
  end if

  if (l.eq.m) then
    WavExtrCPM_plgndr = pmm
  else
    pmmp1 = x*(two*m+one)*pmm
    if (l.eq.m+1) then
      WavExtrCPM_plgndr=pmmp1
    else
      do ll=m+2,l
        pll = ( x*dble(2*ll-1)*pmmp1-dble(ll+m-1)*pmm )/dble(ll-m)
        pmm = pmmp1
        pmmp1 = pll
      end do
      WavExtrCPM_plgndr = pll
    end if
  end if
end function WavExtrCPM_plgndr


! __________________________________________________________________
!
! Calculates the various combinations of spherical harmonics needed
! for the extraction (all are complex):     
!
!   Y  = Ylm
!   Y1 = Ylm,theta
!   Y2 = Ylm,phi
!   Y3 = Ylm,theta,theta + l(l+1)/2 * Y_lm      <== SHH CHANGE
!   Y4 = Ylm,theta,phi - cot theta Ylm,phi
!
!   In terms of Martel and Poisson's (2005) spherical harmonics:
!
!   Y1 = Y_\th = X_\phi / sin \th
!   Y2 = Y_\phi = -sin \th X_\th
!   Y3 = Y_{\th \th} = - Y_{\phi \phi} / sin^2 \th = X_{\th \phi} / sin \th
!   Y4 = Y_{\th \phi} = - sin \th X_{\th \th} = X_{\phi \phi} / sin \th
!
!   OLD COMMENTED OUT BELOW
!   Y3 = Ylm,theta,theta-cot theta Ylm,theta-Ylm,phi,phi/sin^2 theta
!   OLD VALUE
!
! The local variables Yplus is the spherical harmonic at (l+1,m)
!
! All the return values are 2 dim to account for real and complex part.
! __________________________________________________________________
subroutine WavExtrCPM_spher_harm_combs(theta,phi,l,m,Y,Y1,Y2,Y3,Y4)

  use WavExtrCPMConstants

  implicit none

! Input variables
  CCTK_INT :: l,m
  CCTK_REAL :: theta,phi

! Output variables
  CCTK_REAL,DIMENSION(2) :: Y,Y1,Y2,Y3,Y4

! Local variables
  CCTK_INT :: i
  CCTK_REAL ::  Yplus(2),rl,rm,cot_theta
! __________________________________________________________________

  rl = dble(l)
  rm = dble(m)

  cot_theta = cos(theta)/sin(theta)

  call WavExtrCPM_spherical_harmonic(l+1,m,theta,phi,Yplus)

  ! Find Y
  call WavExtrCPM_spherical_harmonic(l,m,theta,phi,Y)

  ! Find Y1
  do i = 1,2
    Y1(i) = -(rl+one)*cot_theta*Y(i)+Yplus(i)/sin(theta) *  sqrt(((rl+one)**2-rm**2)*(rl+half)/(rl+one+half))
  end do

  ! Find Y2
  Y2(1) = -rm*Y(2)
  Y2(2) =  rm*Y(1)

  ! Find Y3
    ! SHH CHANGE
  do i = 1,2
    Y3(i) = -cot_theta*Y1(i) + (rm*rm/(sin(theta)**2) - rl*(rl+one)/two)*Y(i)
!    Y3(i) = -two*cot_theta*Y1(i) + (two*rm*rm/(sin(theta)**2) - rl*(rl+one))*Y(i)
  end do

  ! Find Y4
  Y4(1) = rm*(cot_theta*Y(2)-Y1(2))
  Y4(2) = rm*(Y1(1)-cot_theta*Y(1))
end subroutine WavExtrCPM_spher_harm_combs

function factorial(n)
  implicit none
  integer factorial, n, p, i
  p = 1
  do i = 1, n
     p = p * i
  end do
  factorial = p
end function factorial
