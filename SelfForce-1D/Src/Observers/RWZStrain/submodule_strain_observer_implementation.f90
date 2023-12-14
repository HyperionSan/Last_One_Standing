submodule(strain_observer) strain_observer_implementation
!! The implementation of the interface in [[strain_observer]].

  implicit none

contains

  module procedure strain_init

    implicit none

    select type(object)
    type is (rwz_schw)
      this%p => object
    class default
      print*,'Strain observer initialized with invalid type'
      stop
    end select

    allocate ( this%vname, source = this%p%ename )

    if (size(rad)/=1) then
      print*,'Strain observer should be initialized with 1 extraction radius'
      stop
    end if
    this%nradii = size(rad)
    this%radii = rad
    this%ioo_id = -1
    
    call find_indices ( rad, coord, this%elem_index, this%node_index )

    associate ( nmodes => this%p%nmodes)
      allocate ( this%hplus(nmodes), this%hcross(nmodes) )
    end associate

  end procedure strain_init


  module procedure strain_extract

    use orbit_base
    use gsl_interface, only : legendre_sphPlm_d2, legendre_array_n, &
                              legendre_array_index
    use parameters, only : mass, lmax, q_charge

    implicit none

    real(wp) :: fac1, fac2
    real(wp), dimension(:), allocatable :: ylm, dylm, d2ylm2
    integer(ip) :: i, j, k, nm, lmindex

    i = this%elem_index(1)
    j = this%node_index(1)
    this%hplus = 0.0_wp
    this%hcross = 0.0_wp

    associate ( nmodes => this%p%nmodes, &
                mm => this%p%mm, &
                ll => this%p%ll, &
                hplus => this%hplus, &
                hcross => this%hcross )

      nm = legendre_array_n ( lmax )
      allocate ( ylm(nm), dylm(nm), d2ylm2(nm) )
      call legendre_sphPlm_d2 ( lmax, 0.0_wp, ylm, dylm, d2ylm2 )

      do k = 1, nmodes
        associate ( psi => this%p%eq_data(1,k)%elems(i)%var(j) )
          lmindex = legendre_array_index ( ll(k), mm(k) )
          fac1 = 0.5_wp*ll(k)*(ll(k)+1)*ylm(lmindex) + d2ylm2(lmindex)
          fac2 = zi*mm(k)*dylm(lmindex)

          if (mod(ll(k)+mm(k),2) == 0) then
            hplus(k) = psi*fac1
            hcross(k) = psi*fac2
          else
            hplus(k) = -psi*fac2
            hcross(k) = psi*fac1
          end if
        end associate
      end do 

    end associate
  end procedure strain_extract


  module procedure strain_output

    use output_base
    use time_info

    implicit none

    integer(ip) :: ioo_id, tmp_id
    character(len=len_trim(this%vname)+10) :: filename

    ioo_id = this%ioo_id
    if (ioo_id<0) then
      filename = trim(this%vname) // '.hpl.asc'
      ioo_id = next_available_io_id ()
      this%ioo_id = ioo_id
      print*,'Opening ', filename, ' with id ', ioo_id
      open(ioo_id, file=filename, status='replace', action='write')
      filename = trim(this%vname) // '.hcr.asc'
      tmp_id = next_available_io_id ()
      print*,'Opening ', filename, ' with id ', tmp_id
      open(tmp_id, file=filename, status='replace', action='write')
    end if
    write(ioo_id,'(*(es23.15e3,1x))') get_current_time ( ), this%hplus(:) 
    write(ioo_id+1,'(*(es23.15e3,1x))') get_current_time ( ), this%hcross(:) 
    
  end procedure strain_output

  module procedure close_strain_observer
      deallocate ( this%hplus, this%hcross )
  end procedure close_strain_observer

end submodule strain_observer_implementation 
