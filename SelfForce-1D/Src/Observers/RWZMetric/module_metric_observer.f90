module metric_observer
!! Module that defines the interface for a metric perturbation observer class for a
!! point mass in orbit around a Schwarzschild black hole, i.e. it provides
!! extraction of the perturbation from the variables in [[rwz_schw]]).
!!
!! The implementation is found in
!! [[submodule_metric_observer_implementation.f90]]

  use rwz_pert_schw
  use observers

  implicit none

  type, extends(observer) :: met_observer
  !! A class that defines an observer of the metric perturbation for a point mass
  !! in orbit around a Schwarzschild black hole.
    type(rwz_schw), pointer :: p
    !! A pointer to the [[rwz_schw]] equation that provides the data.
    real(wp), dimension(:), allocatable :: htt
    !! The \(\ell\)-modes of the htt component of the metric perturbation.
    real(wp), dimension(:), allocatable :: htr
    !! The \(\ell\)-modes of the htr component of the metric perturbation.
    real(wp), dimension(:), allocatable :: hrr
    !! The \(\ell\)-modes of the hrr component of the metric perturbation.
    real(wp), dimension(:), allocatable :: ht
    !! The \(\ell\)-modes of the ht component of the metric perturbation.
    real(wp), dimension(:), allocatable :: hr
    !! The \(\ell\)-modes of the hr component of the metric perturbation.
  contains
    procedure :: init => metric_init
    !! The [[observer:init]] routine is provided by [[metric_init]]
    procedure :: extract => metric_extract
    !! The [[observer:extract]] routine is provided by [[metric_extract]]
    procedure :: output => metric_output 
    !! The [[observer:output]] routine is provided by [[metric_output]]
    final :: close_metric_observer
    !! The finalizer.
  end type met_observer

  interface
    module subroutine metric_init ( this, rad, coord, object )
    !! The interface for the [[metric_observer]] version of [[observer:init]].
    !! This interface is consistent with [[obs_init_interface]].
      class(met_observer), intent(inout) :: this
      !! The metric observer that is being initialized.
      real(wp), dimension(:), intent(in) :: rad
      !! A 1d-array containing the radii where the observations have to be
      !! performed. Obviously it only makes sense to pass in the particle
      !! location.
      type(rgf), intent(in) :: coord
      !! A real grid function containing the coordinates of the grid.
      class(*), target, intent(in) :: object
      !! The object on which observations have to be done. This has to be of
      !! type [[rwz_schw]].
    end subroutine metric_init

    module subroutine metric_extract ( this, effs )
    !! The interface for the [[metric_observer]] version of [[observer:extract]].
    !! This interface is consistent with [[obs_extract_interface]].
      class(met_observer), intent(inout) :: this
      !! The routine is called on this metric observer.
      class(eff_source), optional, pointer, intent(in) :: effs
      !! Optionally pass in an effective source object when regularization has
      !! to be done.
    end subroutine metric_extract

    module subroutine metric_output ( this )
    !! The interface for the [[metric_output]] version of [[observer:output]].
    !! This interface is consistent with [[obs_output_interface]].
      class(met_observer), intent(inout) :: this
      !! The routine is called on this metric observer.
    end subroutine metric_output

    module subroutine close_metric_observer ( this )
    !! The interface for the finalizer.
      type(met_observer) :: this
      !! The metric observer to be finalized.
    end subroutine close_metric_observer
  end interface
end module metric_observer
