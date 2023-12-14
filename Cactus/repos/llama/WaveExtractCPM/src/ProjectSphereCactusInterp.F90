! $Header: /numrelcvs/AEIDevelopment/WaveExtract/src/ProjectSphere.F90,v 1.12 2008/02/19 04:35:46 schnetter Exp $

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"



subroutine WavExtrCPM_ProjectSphereCactusInterp(CCTK_ARGUMENTS)

  use WavExtrCPMConstants

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS


  integer, save :: interpolation_timer = -1

  CCTK_INT :: it, ip, j
  CCTK_INT :: interp_handle, table_handle, coord_system_handle

  character(len=200) :: interp
  CCTK_INT :: interp_len
  character(len=20) :: interp_order

  CCTK_INT, dimension(2) :: lsh,lbnd

  CCTK_REAL :: dtheta, dphi, dthetainv, dphiinv
  CCTK_POINTER, dimension(3) :: interp_coords
  CCTK_POINTER, dimension(27) :: out_array
  CCTK_INT, dimension(27) :: in_array
  CCTK_INT, dimension(27) :: out_types

! temp shortcuts
  CCTK_REAL :: tgxx, dgxx, &
               tgxy, dgxy, &
               tgxz, dgxz, &
               tgyy, dgyy, &
               tgyz, dgyz, &
               tgzz, dgzz

  CCTK_INT :: di

  CCTK_INT,dimension(27) :: op_indices
  CCTK_INT,dimension(27) :: op_codes

! FIXME : make this the real lapse !
  CCTK_REAL :: lapse

  CCTK_INT :: ierr,num_in_arrays,num_out_arrays,status, istat
  CCTK_REAL :: rad,r2,ct,st,ct2,st2,cp,cp2,sp,sp2
  CCTK_REAL :: cor_angle,siO,coO,cosO,sinO
  character(len=200) :: infoline
  CCTK_INT :: metric_string_len, metric_group_index, metric_firstvar_index
  character(len=200) :: local_metric, warnline


  integer :: i,k
  integer :: interpolator_error_level
  integer :: sf_i, sf_j, sf_k

  CCTK_REAL :: g11,g12,g13,g21,g22,g23,g31,g32,g33
  CCTK_REAL :: g11_n,g12_n,g13_n,g21_n,g22_n,g23_n,g31_n,g32_n,g33_n
  CCTK_REAL :: d1_g11,d1_g12,d1_g13,d1_g21,d1_g22,d1_g23,d1_g31,d1_g32,d1_g33
  CCTK_REAL :: d2_g11,d2_g12,d2_g13,d2_g21,d2_g22,d2_g23,d2_g31,d2_g32,d2_g33
  CCTK_REAL :: d3_g11,d3_g12,d3_g13,d3_g21,d3_g22,d3_g23,d3_g31,d3_g32,d3_g33
  CCTK_REAL :: dx_g11,dx_g12,dx_g13,dx_g21,dx_g22,dx_g23,dx_g31,dx_g32,dx_g33
  CCTK_REAL :: dy_g11,dy_g12,dy_g13,dy_g21,dy_g22,dy_g23,dy_g31,dy_g32,dy_g33
  CCTK_REAL :: dz_g11,dz_g12,dz_g13,dz_g21,dz_g22,dz_g23,dz_g31,dz_g32,dz_g33

  CCTK_REAL :: jac(3,3), det_jac

  CCTK_REAL :: dxcdxs, dxcdys, dycdxs, dycdys

  CCTK_REAL :: xv,yv
! ________________________________________________________________________

  
  if (verbose>2) &
    call CCTK_INFO('Interpolating metric and derivatives onto sphere')

  if (do_nothing == 1) &
    return

  if (cctk_iteration .ne. 0) then
    if (mod(cctk_iteration,my_out_every_det(current_detector)).ne.0) then
      if (verbose>2) call CCTK_INFO("No time for this detector")
      return
    end if
  end if

  current_detector_radius = detector_radius(current_detector)

  di=current_detector

  ! Sanity Check
  if (current_detector_radius .lt. 1.d-10) then
    do_nothing =1
    call CCTK_WARN(1,"This should never happen: The detector radius is 0!")
    return
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
    endiF
  end if

  if (verbose>1) then
    write(infoline,'(a,i4,a,g16.8)') 'Analysing Detector No.: ',current_detector, &
                                        ' Radius ',current_detector_radius
    call CCTK_INFO(infoline)
  end if

  if (interpolation_timer < 0) then
     call CCTK_TimerCreate(interpolation_timer, "WaveExtractCPM_CactusInterp")
     if (interpolation_timer < 0) then
        call CCTK_WARN(1,"could not create WaveExtractCPM_CactusInterp timer")
     end if
  end if


  ! local shape of the 2D grid arrays
  call CCTK_GrouplshGN(ierr, cctkGH, 2, lsh, "WaveExtractCPM::surface_arrays")
  if ( ierr .lt. 0 ) then
    call CCTK_WARN(0, "cannot get local size for surface arrays")
  end if

  if (size(ctheta,1)<2) call CCTK_WARN (0, "internal error")
  if (size(cphi,2)<2) call CCTK_WARN (0, "internal error")
  dtheta = ctheta(2,1) - ctheta(1,1)
  dphi = cphi(1,2) - cphi(1,1)

  if (cartoon .ne. 0) then
    dphi = two*pi
  end if

  dthetainv = one / dtheta
  dphiinv = one / dphi

  call CCTK_FortranString(interp_len, interpolation_operator, interp)

  call CCTK_InterpHandle(interp_handle,interp(1:interp_len))
  if ( interp_handle .lt. 0 ) then
    call CCTK_WARN( 0, "Cannot get handle for interpolation. Forgot to activate an implementation providing interpolation operators?" )
  end if

  write(interp_order,'(a6,i1)') 'order=',interpolation_order
  call Util_TableCreateFromString(table_handle,interp_order)
  if ( table_handle .lt. 0 ) then
    call CCTK_WARN( 0, "Cannot create parameter table for interpolator" )
  end if

  if (make_interpolator_warnings_fatal .ne. 0) then
    interpolator_error_level = 0
  else
    interpolator_error_level = 1
  end if

  ! For observers we use it's own coord system
  if (observers .ne. 0) then
    ! FIXME : add observers coordinate system handle !
    call CCTK_WARN(1,"IMPLEMENT COORDINATE SYSTEM HANDLE FOR OBSERVERS")
    call CCTK_CoordSystemHandle ( coord_system_handle, "cart3d" )
    if ( coord_system_handle .lt. 0) then
      call CCTK_WARN( 0, "Cannot get handle for observers coordinate system. Forgot to activate an implementation providing coordinates?" )
    end if
  else 
    call CCTK_CoordSystemHandle ( coord_system_handle, "cart3d" )
    if ( coord_system_handle .lt. 0) then
      call CCTK_WARN( 0, "Cannot get handle for cart3d coordinate system. Forgot to activate an implementation providing coordinates?" )
    end if
  end if

  call CCTK_TimerStart(ierr,"WaveExtractCPM_CactusInterp")

  interp_coords(1) = CCTK_PointerTo(interp_x)
  interp_coords(2) = CCTK_PointerTo(interp_y)
  interp_coords(3) = CCTK_PointerTo(interp_z)


  out_array(1)  = CCTK_PointerTo(gxxi)
  out_array(2)  = CCTK_PointerTo(gxyi)
  out_array(3)  = CCTK_PointerTo(gxzi)
  out_array(4)  = CCTK_PointerTo(gyyi)
  out_array(5)  = CCTK_PointerTo(gyzi)
  out_array(6)  = CCTK_PointerTo(gzzi)
  out_array(7)  = CCTK_PointerTo(dr_gxxi)
  out_array(8)  = CCTK_PointerTo(dr_gxyi)
  out_array(9)  = CCTK_PointerTo(dr_gxzi)
  out_array(10)  = CCTK_PointerTo(dr_gyyi)
  out_array(11)  = CCTK_PointerTo(dr_gyzi)
  out_array(12)  = CCTK_PointerTo(dr_gzzi)

  ! SHH CHANGE
  out_array(13)  = CCTK_PointerTo(gtxi)
  out_array(14)  = CCTK_PointerTo(gtyi)
  out_array(15)  = CCTK_PointerTo(gtzi)
  out_array(16)  = CCTK_PointerTo(dt_gtxi)
  out_array(17)  = CCTK_PointerTo(dt_gtyi)
  out_array(18)  = CCTK_PointerTo(dt_gtzi)

  out_array(19)  = CCTK_PointerTo(dt_gxxi)
  out_array(20)  = CCTK_PointerTo(dt_gxyi)
  out_array(21)  = CCTK_PointerTo(dt_gxzi)
  out_array(22)  = CCTK_PointerTo(dt_gyyi)
  out_array(23)  = CCTK_PointerTo(dt_gyzi)
  out_array(24)  = CCTK_PointerTo(dt_gzzi)

  out_array(25)  = CCTK_PointerTo(dr_gtxi)
  out_array(26)  = CCTK_PointerTo(dr_gtyi)
  out_array(27)  = CCTK_PointerTo(dr_gtzi)

  num_out_arrays=27

  ! set all out_types(1:27) to type real
  out_types=CCTK_VARIABLE_REAL

  ! Find out the lower bounds of the distributed integration grid arrays.
  call CCTK_GrouplbndGN(ierr, cctkGH,2,lbnd,"waveextractcpm::surface_integrands")
  if ( ierr .lt. 0 ) then
    call CCTK_WARN(0, "cannot get lower bounds for surface integrands")
  end if

  ! find the cartesian coordinates for the interpolation points
  if (use_spherical_surface .ne. 0) then
    if(verbose>2) then
      call CCTK_INFO("go for spherical surface")
    end if


    !if (int_ntheta(di).ne.sf_ntheta(surface_index(current_detector)+1)) then
    !  call CCTK_WARN(1,"different theta surface for spherical surface")
    !end if

    !if (int_nphi(di).ne.sf_nphi(surface_index(current_detector)+1)) then
    !  call CCTK_WARN(1,"different phi surface for spherical surface")
    !end if
    
    
    do ip = 1, lsh(2)
      do it = 1, lsh(1)
        ! FIXME : check the index computation for parallel case
        
        sf_k = surface_index(current_detector)+1
        
        sf_i = it+lbnd(1) + nghoststheta(sf_k)
        sf_j = ip+lbnd(2) + nghostsphi(sf_k)
        
        ! get proper sf_radius-index by taking into account symmetries
        if (symmetric_z(sf_k)/=0) then
           if (sf_i > sf_ntheta(sf_k)-nghoststheta(sf_k)) then
              sf_i = sf_ntheta(sf_k)-nghoststheta(sf_k) - (sf_i - (sf_ntheta(sf_k)-nghoststheta(sf_k)))
           end if
        end if
        if (symmetric_x(sf_k)/=0) then
           if (symmetric_y(sf_k)/=0) then
              if (sf_j > sf_nphi(sf_k)-nghostsphi(sf_k) .and. sf_j < 2*(sf_nphi(sf_k)-nghostsphi(sf_k))) then
                 sf_j = sf_nphi(sf_k)-nghostsphi(sf_k) - (sf_j - (sf_nphi(sf_k)-nghostsphi(sf_k)))
              else if (sf_j > 2*(sf_nphi(sf_k)-nghostsphi(sf_k))) then
                 sf_j = sf_j - 2*(sf_nphi(sf_k)-nghostsphi(sf_k))
              end if
           else
              if (sf_j > sf_nphi(sf_k)-nghostsphi(sf_k)) then
                 sf_j = sf_nphi(sf_k)-nghostsphi(sf_k) - (sf_j - (sf_nphi(sf_k)-nghostsphi(sf_k)))
              end if
           end if
        end if
        if (symmetric_y(sf_k)/=0) then
           if (sf_j > sf_nphi(sf_k)-nghostsphi(sf_k)) then
              sf_j = sf_nphi(sf_k)-nghostsphi(sf_k) - (sf_j - (sf_nphi(sf_k)-nghostsphi(sf_k)))
           end if
        end if
        
        rad=sf_radius(sf_i,sf_j,sf_k)
        if ((it+lbnd(1) <= int_ntheta(current_detector) .and. ip+lbnd(2) <= int_nphi(current_detector)) .and. rad .eq. 0) then 
           !print*,'FIXME rad is ',rad,it,ip,lbnd(1),lbnd(2),sf_i,sf_j,int_ntheta(current_detector),int_nphi(current_detector),current_detector,surface_index(current_detector)+1
           call CCTK_WARN(1, "Detector radius as given by SphericalSurface is 0!")
        endif
        
        if (it+lbnd(1)>int_ntheta(current_detector) .or. &
            ip+lbnd(2)>int_nphi(current_detector) ) then
          interp_x(it,ip)=zero
          interp_y(it,ip)=zero
          interp_z(it,ip)=zero
        else
          interp_x(it,ip) = origin_x + rad * sintheta(it,ip) * cosphi(it,ip)
          interp_y(it,ip) = origin_y + rad * sintheta(it,ip) * sinphi(it,ip)
          interp_z(it,ip) = origin_z + rad * costheta(it,ip)
        end if
      end do
    end do
    rad=current_detector_radius  ! set to mean radius because we need that later on  FIXME: May want to set this to areal radius!
 else
    rad=current_detector_radius
    do ip = 1, lsh(2)
      do it = 1, lsh(1)
        interp_x(it,ip) = origin_x + rad * sintheta(it,ip) * cosphi(it,ip)
        interp_y(it,ip) = origin_y + rad * sintheta(it,ip) * sinphi(it,ip)
        interp_z(it,ip) = origin_z + rad * costheta(it,ip)
      end do
    end do
 end if

  !print*,'interp_x',interp_x
  !print*,'interp_y',interp_y
  !print*,'interp_z',interp_z


  call CCTK_FortranString(metric_string_len, &
       detector_metric(current_detector), &
       local_metric)
  call CCTK_GroupIndex(metric_group_index, local_metric(1:metric_string_len))
  call CCTK_ActiveTimeLevelsGI(ierr, cctkGH, metric_group_index)
  if (ierr < 1) then
    write(warnline,*) "The metric ", &
         local_metric(1:metric_string_len), " for detector number", &
         current_detector, " does not have storage on"
    call CCTK_WARN(0, warnline)
  end if
  call CCTK_FirstVarIndexI(metric_firstvar_index, metric_group_index)

  if (verbose > 1) then
    write(infoline,*) "Extracting with respect to metric ",&
         local_metric(1:metric_string_len)
    call CCTK_INFO(infoline)
  end if

!  in_array(1) = metric_firstvar_index + 0

  call CCTK_VarIndex(in_array(1), "admbase::gxx")
  call CCTK_VarIndex(in_array(2), "admbase::gxy")
  call CCTK_VarIndex(in_array(3), "admbase::gxz")
  call CCTK_VarIndex(in_array(4), "admbase::gyy")
  call CCTK_VarIndex(in_array(5), "admbase::gyz")
  call CCTK_VarIndex(in_array(6), "admbase::gzz")
  call CCTK_VarIndex(in_array(7), "admderivatives::gxx_dr")
  call CCTK_VarIndex(in_array(8), "admderivatives::gxy_dr")
  call CCTK_VarIndex(in_array(9), "admderivatives::gxz_dr")
  call CCTK_VarIndex(in_array(10), "admderivatives::gyy_dr")
  call CCTK_VarIndex(in_array(11), "admderivatives::gyz_dr")
  call CCTK_VarIndex(in_array(12), "admderivatives::gzz_dr")
  ! SHH CHANGE
  call CCTK_VarIndex(in_array(13), "admbase::betax")
  call CCTK_VarIndex(in_array(14), "admbase::betay")
  call CCTK_VarIndex(in_array(15), "admbase::betaz")
  call CCTK_VarIndex(in_array(16), "admbase::dtbetax")
  call CCTK_VarIndex(in_array(17), "admbase::dtbetay")
  call CCTK_VarIndex(in_array(18), "admbase::dtbetaz")

  call CCTK_VarIndex(in_array(19), "admderivatives::gxx_dt")
  call CCTK_VarIndex(in_array(20), "admderivatives::gxy_dt")
  call CCTK_VarIndex(in_array(21), "admderivatives::gxz_dt")
  call CCTK_VarIndex(in_array(22), "admderivatives::gyy_dt")
  call CCTK_VarIndex(in_array(23), "admderivatives::gyz_dt")
  call CCTK_VarIndex(in_array(24), "admderivatives::gzz_dt")

  call CCTK_VarIndex(in_array(25), "admderivatives::betax_dr")
  call CCTK_VarIndex(in_array(26), "admderivatives::betay_dr")
  call CCTK_VarIndex(in_array(27), "admderivatives::betaz_dr")
  
  num_in_arrays=27

  ! table entries for interpolator

  do i = 1, num_out_arrays
     op_indices(i) = i-1
     op_codes(i) = 0
  end do
  
  call Util_TableSetIntArray ( status, table_handle, num_out_arrays, &
       op_indices(1:num_out_arrays), "operand_indices" )
  if(status.lt. 0) then
     call CCTK_WARN(0, "Cannot set operand indices array in parameter table")
  end if
  call Util_TableSetIntArray ( status, table_handle, num_out_arrays, &
       op_codes(1:num_out_arrays), "operation_codes" )
  if (status.lt. 0) then
     call CCTK_WARN(0, "Cannot set operation codes array in parameter table")
  end if

  call CCTK_InterpGridArrays(ierr, cctkGH, 3, interp_handle, &
       table_handle, coord_system_handle, &
       lsh(1) * lsh(2), CCTK_VARIABLE_REAL, &
       interp_coords, &
       num_in_arrays, &
       in_array(1:num_in_arrays), &
       num_out_arrays, &
       out_types(1:num_out_arrays), &
       out_array(1:num_out_arrays))
  if (ierr.ne.0) then
     call CCTK_WARN(interpolator_error_level,"Interpolation of metric failed")
  end if

  call CCTK_TimerStop(ierr,"WaveExtractCPM_CactusInterp")
  call Util_TableDestroy ( status, table_handle )
          
end subroutine WavExtrCPM_ProjectSphereCactusInterp

