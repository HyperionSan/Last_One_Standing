module flux_observer
!! Module that defines the interface for a metric flux observer class for a
!! point mass in orbit around a Schwarzschild black hole, i.e. it provides
!! extraction of the flux from the variables in [[rwz_schw]]).
!!
!! The implementation is found in
!! [[submodule_flux_observer_implementation.f90]]

  use rwz_pert_schw
  use observers

  implicit none

  type, extends(observer) :: flx_observer
  !! A class that defines an observer of the metric flux for a point mass
  !! in orbit around a Schwarzschild black hole.
    type(rwz_schw), pointer :: p
    !! A pointer to the [[rwz_schw]] equation that provides the data.
    complex(wp) :: Eflux_hori
    !! The summed energy flux at the horizon.
    complex(wp) :: Jflux_hori
    !! The summed angular momentum flux at the horizon.
    complex(wp) :: Eflux_scri
    !! The summed energy flux at Scri.
    complex(wp) :: Jflux_scri
    !! The summed angular momentum flux at Scri.
    complex(wp), dimension(:), allocatable :: Eflux_hori_l
    !! The \(\ell\)-modes of the energy flux at the horizon.
    complex(wp), dimension(:), allocatable :: Jflux_hori_l
    !! The \(\ell\)-modes of the angular momentum flux at the horizon.
    complex(wp), dimension(:), allocatable :: Eflux_scri_l
    !! The \(\ell\)-modes of the energy flux at Scri.
    complex(wp), dimension(:), allocatable :: Jflux_scri_l
    !! The \(\ell\)-modes of the angular momentum flux at Scri.
    complex(wp), dimension(:), allocatable :: Eflux_hori_lm
    !! The \(\ell\)m-modes of the energy flux at the horizon.
    complex(wp), dimension(:), allocatable :: Jflux_hori_lm
    !! The \(\ell\)m-modes of the angular momentum flux at Scri.
    complex(wp), dimension(:), allocatable :: Eflux_scri_lm
    !! The \(\ell\)m-modes of the energy flux at Scri.
    complex(wp), dimension(:), allocatable :: Jflux_scri_lm
    !! The \(\ell\)m-modes of the angular momentum flux at Scri.
  contains
    procedure :: init => flux_init
    !! The [[observer:init]] routine is provided by [[flux_init]]
    procedure :: extract => flux_extract
    !! The [[observer:extract]] routine is provided by [[flux_extract]]
    procedure :: output => flux_output 
    !! The [[observer:output]] routine is provided by [[flux_output]]
    final :: close_flux_observer
    !! The finalizer.
  end type flx_observer

  interface
    module subroutine flux_init ( this, rad, coord, object )
    !! The interface for the [[flux_observer]] version of [[observer:init]].
    !! This interface is consistent with [[obs_init_interface]].
      class(flx_observer), intent(inout) :: this
      !! The flux observer that is being initialized.
      real(wp), dimension(:), intent(in) :: rad
      !! A 1d-array containing the radii where the observations have to be
      !! performed. Obviously it only makes sense to pass in the particle
      !! location.
      type(rgf), intent(in) :: coord
      !! A real grid function containing the coordinates of the grid.
      class(*), target, intent(in) :: object
      !! The object on which observations have to be done. This has to be of
      !! type [[rwz_schw]].
    end subroutine flux_init

    module subroutine flux_extract ( this, effs )
    !! The interface for the [[flux_observer]] version of [[observer:extract]].
    !! This interface is consistent with [[obs_extract_interface]].
      class(flx_observer), intent(inout) :: this
      !! The routine is called on this flux observer.
      class(eff_source), optional, pointer, intent(in) :: effs
      !! Optionally pass in an effective source object when regularization has
      !! to be done.
    end subroutine flux_extract

    module subroutine flux_output ( this )
    !! The interface for the [[flux_output]] version of [[observer:output]].
    !! This interface is consistent with [[obs_output_interface]].
      class(flx_observer), intent(inout) :: this
      !! The routine is called on this flux observer.
    end subroutine flux_output

    module subroutine close_flux_observer ( this )
    !! The interface for the finalizer.
      type(flx_observer) :: this
      !! The flux observer to be finalized.
    end subroutine close_flux_observer
  end interface
end module flux_observer
