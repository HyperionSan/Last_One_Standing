module strain_observer
!! Module that defines the interface for a metric strain observer class for a
!! point mass in orbit around a Schwarzschild black hole, i.e. it provides
!! extraction of the strain from the variables in [[rwz_schw]]).
!!
!! The implementation is found in
!! [[submodule_strain_observer_implementation.f90]]

  use rwz_pert_schw
  use observers

  implicit none

  type, extends(observer) :: str_observer
  !! A class that defines an observer of the metric strain for a point mass
  !! in orbit around a Schwarzschild black hole.
    type(rwz_schw), pointer :: p
    !! A pointer to the [[rwz_schw]] equation that provides the data.
    real(wp), dimension(:), allocatable :: hplus
    !! The \(\ell\)-modes of the plus component of the strain.
    real(wp), dimension(:), allocatable :: hcross
    !! The \(\ell\)-modes of the cross component of the strain.
  contains
    procedure :: init => strain_init
    !! The [[observer:init]] routine is provided by [[strain_init]]
    procedure :: extract => strain_extract
    !! The [[observer:extract]] routine is provided by [[strain_extract]]
    procedure :: output => strain_output 
    !! The [[observer:output]] routine is provided by [[strain_output]]
    final :: close_strain_observer
    !! The finalizer.
  end type str_observer

  interface
    module subroutine strain_init ( this, rad, coord, object )
    !! The interface for the [[strain_observer]] version of [[observer:init]].
    !! This interface is consistent with [[obs_init_interface]].
      class(str_observer), intent(inout) :: this
      !! The strain observer that is being initialized.
      real(wp), dimension(:), intent(in) :: rad
      !! A 1d-array containing the radii where the observations have to be
      !! performed. Obviously it only makes sense to pass in the particle
      !! location.
      type(rgf), intent(in) :: coord
      !! A real grid function containing the coordinates of the grid.
      class(*), target, intent(in) :: object
      !! The object on which observations have to be done. This has to be of
      !! type [[rwz_schw]].
    end subroutine strain_init

    module subroutine strain_extract ( this, effs )
    !! The interface for the [[strain_observer]] version of [[observer:extract]].
    !! This interface is consistent with [[obs_extract_interface]].
      class(str_observer), intent(inout) :: this
      !! The routine is called on this strain observer.
      class(eff_source), optional, pointer, intent(in) :: effs
      !! Optionally pass in an effective source object when regularization has
      !! to be done.
    end subroutine strain_extract

    module subroutine strain_output ( this )
    !! The interface for the [[strain_output]] version of [[observer:output]].
    !! This interface is consistent with [[obs_output_interface]].
      class(str_observer), intent(inout) :: this
      !! The routine is called on this strain observer.
    end subroutine strain_output

    module subroutine close_strain_observer ( this )
    !! The interface for the finalizer.
      type(str_observer) :: this
      !! The strain observer to be finalized.
    end subroutine close_strain_observer
  end interface
end module strain_observer
