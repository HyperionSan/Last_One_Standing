module null_params

  implicit none

  integer, save :: nn = 31, &  ! number of grid points in p/q
                   nx = 65, &  ! number of radial grid points
                   nt = 32, &  ! total number of time steps
             it_start =  0, &  ! where to start time loop iteration count
                 iot0 =  1, &  ! delta it for 0D output
                 iot1 =  1, &  ! delta it for 1D output
                 iot2 = 20, &  ! delta it for 2D output
                 iot3 = 50, &  ! delta it for 3D output
               it_cfl =  4, &  ! delta it for resetting dt by CFL-cond.
              id_case =  1, &  ! type of initial data
              ONE, ZERO, HALF, EDGE, HFED, PNTT,&
              NIntPts = 3, it_checkpoint=10000, &
              cdumpj  = 1, &
              cdumpit = 10,&
              stdang = 13, &
              stdrad = 31

  double precision, save :: rwt = 1.0d0,  &
                            cfl = 0.5d0,  &
                          qsize = 1.0d0,  &
                          disip = 0.05d0, &
                           time = 0.0d0,  &
                        R_left  = 4.,     &
	                R_right = 8.,     &
	                ID_AMP  = 0.01,   &
                        dt_fix = 1.0d-4,  &
                        circ_rad = .5d0,  &
                        time_real_start,  &
                        real_start_dt

  logical, save ::  horizon_output  = .true.,   &
       &            scri_output     = .true.,   &
       &            boundary_output = .true.,   &
       &            cfl_compute     = .true.,   &
       &            chkrcvr         = .false.,  &
       &            do_checkpoint   = .false.,  &
       &            pure_schwarzs   = .false.,  &
       &            calc_news       = .true.,   &
       &            horizon_on_grid = .false.,  &
       &            docdump         = .false.

end module null_params
