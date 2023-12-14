module null_params

  implicit none

  integer,          save :: it_skip, nn, nx, nt
  logical,          save :: null_output
  double precision, save :: mass, rwt, cfl, time, disip

end module null_params
