submodule(metric_observer) metric_observer_implementation
!! The implementation of the interface in [[metric_observer]].

  implicit none

contains

  module procedure metric_init

    implicit none

    select type(object)
    type is (rwz_schw)
      this%p => object
    class default
      print*,'Metric Perturbation observer initialized with invalid type'
      stop
    end select

    allocate ( this%vname, source = this%p%ename )

    if (size(rad)/=1) then
      print*,'Metric Perturbation observer should be initialized with 1 extraction radius'
      stop
    end if
    this%nradii = size(rad)
    this%radii = rad
    this%ioo_id = -1
    
    call find_indices ( rad, coord, this%elem_index, this%node_index )

    associate ( nmodes => this%p%nmodes)
      allocate ( this%htt(nmodes), this%htr(nmodes), this%hrr(nmodes), &
                 this%ht(nmodes), this%hr(nmodes) )
    end associate

  end procedure metric_init


  module procedure metric_extract

    use orbit_base
    use parameters, only : mass, lmax, order
    use grid, only : drdrho

    implicit none

    real(wp) :: r, phi, ur, En, Lz
    real(wp) ::twommr, m2, m3, mr, m2r, mr2, r2, r3, longexp
    complex(wp) :: dpsiedr, d2psiedr2, &
                   d2psiedtdr, dpsiodr
    real(wp), dimension(:), allocatable :: drho, dpi
    integer(ip) :: i, j, k, nm
    integer(ip) :: lfac1, lfac2, lfac22, lfac3

    i = this%elem_index(1)
    j = this%node_index(1)
    this%htt = 0.0_wp
    this%htr = 0.0_wp
    this%hrr = 0.0_wp
    this%ht = 0.0_wp
    this%hr = 0.0_wp

    call orbit_info%get_orbit ( r, phi, ur, En, Lz )    

    m2 = mass**2
    m3 = mass*m2
    mr = mass*r
    m2r = m2*r
    r2 = r**2
    mr2 = mass*r2
    r3 = r*r2
    twommr = 2.0_wp*mass - r

    associate ( nmodes => this%p%nmodes, &
                mm => this%p%mm, &
                ll => this%p%ll, &
                htt => this%htt, &
                htr => this%htr, &
                hrr => this%hrr, &
                ht => this%ht, &
                hr => this%hr )

      allocate ( drho(order+1), dpi(order+1) )

      do k = 1, nmodes
        drho(:) = drdrho%elems(i)%var * &
                    matmul ( this%p%refelem%dr, &
                             this%p%eq_data(2,k)%elems(i)%var )
        dpi(:) = drdrho%elems(i)%var * &
                    matmul ( this%p%refelem%dr, &
                             this%p%eq_data(3,k)%elems(i)%var )

        associate ( psi => this%p%eq_data(1,k)%elems(i)%var(j), &
                    dpsidt => this%p%eq_data(2,k)%elems(i)%var(j), &
                    dpsidr => this%p%eq_data(3,k)%elems(i)%var(j), &
                    d2psidtdr => drho(j), &
                    d2psidr2 => dpi(j) )
          lfac1 = ll(k)*(ll(k)+1)
          lfac2 = (ll(k)-1)*(ll(k)+2)
          lfac22 = lfac2**2
          longexp = 72.0_wp*m3 + 36.0_wp*lfac2*m2r + 6.0_wp*lfac22*mr2 + lfac3*r3

          if (mod(ll(k)+mm(k),2) == 0) then
            d2psiedr2 = ( r**2*d2psidr2 - 2.0_wp*mass*dpsidr ) / twommr**2
            d2psiedtdr = -r * d2psidtdr / twommr
            dpsiedr = -r * dpsidr / twommr

            htt(k) = twommr*longexp*psi / ( 2.0_wp*r3 * ( 6.0_wp*mass + lfac2*r )**2 ) &
                 -twommr * ( 6.0_wp*m2 + lfac2*(r2 - mr) ) * dpsiedr &
                 /( r2 * (6.0_wp*mass + lfac2*r ) )  &
                 +twommr**2 * d2psiedr2 / r

            htr(k) = ( 6.0_wp*m2 + lfac2 * ( 3.0_wp*mr - r2 ) ) * dpsidt &
                 /( twommr * ( 6.0_wp * mass + lfac2*r ) ) &
                 +r*d2psiedtdr

            hrr(k) = longexp*psi / ( 2.0_wp*twommr*r * ( 6.0_wp*mass +lfac2*r)**2 ) &
                 +( -6.0_wp*m2 + lfac2 * ( mr - r2 ) ) * dpsiedr  &
                 /( twommr * ( 6.0_wp*mass +lfac2*r ) ) &
                 +r*d2psiedr2

            ht(k) = 0.0_wp
            hr(k) = 0.0_wp
          else
            htt(k) = 0.0_wp
            htr(k) = 0.0_wp
            hrr(k) = 0.0_wp
            dpsiodr = -r * dpsidr / twommr
            ht(k) = ( 0.5_wp - mass/r ) * psi - 0.5_wp * ( twommr * dpsiodr )
            hr(k) = - r2*dpsidt / ( 2.0_wp*twommr )
          end if
        end associate
      end do 

    end associate
  end procedure metric_extract


  module procedure metric_output

    use output_base
    use time_info

    implicit none

    integer(ip) :: ioo_id, tmp_id
    character(len=len_trim(this%vname)+10) :: filename

    ioo_id = this%ioo_id
    if (ioo_id<0) then
      filename = trim(this%vname) // '.htt.asc'
      ioo_id = next_available_io_id ()
      this%ioo_id = ioo_id
      print*,'Opening ', filename, ' with id ', ioo_id
      open(ioo_id, file=filename, status='replace', action='write')
      filename = trim(this%vname) // '.htr.asc'
      tmp_id = next_available_io_id ()
      print*,'Opening ', filename, ' with id ', tmp_id
      open(tmp_id, file=filename, status='replace', action='write')
      filename = trim(this%vname) // '.hrr.asc'
      tmp_id = next_available_io_id ()
      print*,'Opening ', filename, ' with id ', tmp_id
      open(tmp_id, file=filename, status='replace', action='write')
      filename = trim(this%vname) // '.ht.asc'
      tmp_id = next_available_io_id ()
      print*,'Opening ', filename, ' with id ', tmp_id
      open(tmp_id, file=filename, status='replace', action='write')
      filename = trim(this%vname) // '.hr.asc'
      tmp_id = next_available_io_id ()
      print*,'Opening ', filename, ' with id ', tmp_id
      open(tmp_id, file=filename, status='replace', action='write')
    end if
    write(ioo_id,'(*(es23.15e3,1x))') get_current_time ( ), this%htt(:) 
    write(ioo_id+1,'(*(es23.15e3,1x))') get_current_time ( ), this%htr(:) 
    write(ioo_id+2,'(*(es23.15e3,1x))') get_current_time ( ), this%hrr(:) 
    write(ioo_id+3,'(*(es23.15e3,1x))') get_current_time ( ), this%ht(:) 
    write(ioo_id+4,'(*(es23.15e3,1x))') get_current_time ( ), this%hr(:) 
    
  end procedure metric_output


  module procedure close_metric_observer
      deallocate ( this%htt, this%htr, this%hrr, this%ht, this%hr )
  end procedure close_metric_observer

end submodule metric_observer_implementation 
