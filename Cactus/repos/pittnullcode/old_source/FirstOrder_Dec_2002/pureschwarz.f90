module pureschwarzsboundary
  use null_params
  use null_grid
  use null_vars
  implicit none
contains
  subroutine schwarzschild_driver
    implicit none
    double precision :: SMass  ! here mass will be inferred from rwt
    double precision :: RR
    integer :: i
    SMass = .5d0*rwt 
    !  "boundary" will consist of the points i=1, and i=2
    maskn = 2
    ! hence first evolution point will be at i=3

    do i=1, 2
      RR = rb(i)
      jnn(:,:,i) = (0.0d0, 0.0d0)
   
      bnn(:,:,i) = .5d0 * log(-4.0d0*SMass / time)

      cbnn(:,:,i) = (0.0d0, 0.0d0)
         
      cknn(:,:,i) = (0.0d0, 0.0d0)
 
      nunn(:,:,i) = (0.0d0, 0.0d0)
 
      ! remeber that U lies on the half grid
      unn(:,:,i) = 0.0d0

      wnn(:,:,i) = (-4.0d0 * SMass / time * ( 1.0d0 - 2.0d0*SMass / RR) - 1.0d0 ) /RR  
    end do

    
  end subroutine schwarzschild_driver

end module pureschwarzsboundary
