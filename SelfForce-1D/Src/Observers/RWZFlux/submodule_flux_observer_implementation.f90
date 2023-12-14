submodule(flux_observer) flux_observer_implementation
!! The implementation of the interface in [[flux_observer]].

  implicit none

contains

  module procedure flux_init

    use parameters, only : lmin, lmax
    implicit none

    select type(object)
    type is (rwz_schw)
      this%p => object
    class default
      print*,'Flux observer initialized with invalid type'
      stop
    end select

    allocate ( this%vname, source = this%p%ename )

    if (size(rad)/=2) then
      print*,'Flux observer should be initialized with 2 extraction radii'
      stop
    end if
    this%nradii = size(rad)
    this%radii = rad
    this%ioo_id = -1

    call find_indices ( rad, coord, this%elem_index, this%node_index )

    associate ( nmodes => this%p%nmodes)
      allocate ( this%Eflux_scri_lm(nmodes), this%Jflux_scri_lm(nmodes), &
                 this%Eflux_hori_lm(nmodes), this%Jflux_hori_lm(nmodes) )
      
      allocate ( this%Eflux_scri_l(lmin:lmax), this%Jflux_scri_l(lmin:lmax), &
                 this%Eflux_hori_l(lmin:lmax), this%Jflux_hori_l(lmin:lmax) )
    end associate

  end procedure flux_init


  module procedure flux_extract

    use orbit_base

    implicit none

    real(wp) :: lr, lfactor
    real(wp), parameter :: pi = acos(-1.0_wp)
    real(wp), parameter :: pifactor = 1.0_wp/(64.0_wp*pi)
    integer(ip) :: i, j, k

    this%Eflux_scri_l = 0; this%Jflux_scri_l = 0; this%Eflux_hori_l = 0; this%Jflux_hori_l = 0

! if loops are paralellized, consider combining the two loops over modes
    i = this%elem_index(1)
    j = this%node_index(1)

    do k = 1, this%p%nmodes
      associate ( l => this%p%ll(k), &
                  m => this%p%mm(k), &
                  psi => this%p%eq_data(1,k)%elems(i)%var(j), &
                  dpsidt => this%p%eq_data(2,k)%elems(i)%var(j) )
        lr = real(l,wp)
        lfactor = pifactor*(lr-1.0_wp)*lr*(lr+1.0_wp)*(lr+2.0_wp)
        this%Eflux_hori_lm(k) = lfactor * conjg(dpsidt)*dpsidt
        this%Jflux_hori_lm(k) = lfactor * zi*m*conjg(psi)*dpsidt
        this%Eflux_hori_l(l) = this%Eflux_hori_l(l) + this%Eflux_hori_lm(k)
        this%Jflux_hori_l(l) = this%Jflux_hori_l(l) + this%Jflux_hori_lm(k)
      end associate
    end do 

    i = this%elem_index(2)
    j = this%node_index(2)

    do k = 1, this%p%nmodes
      associate ( l => this%p%ll(k), &
                  m => this%p%mm(k), &
                  psi => this%p%eq_data(1,k)%elems(i)%var(j), &
                  dpsidt => this%p%eq_data(2,k)%elems(i)%var(j) )
        lr = real(l,wp)
        lfactor = pifactor*(lr-1.0_wp)*lr*(lr+1.0_wp)*(lr+2.0_wp)
        this%Eflux_scri_lm(k) = lfactor * conjg(dpsidt)*dpsidt
        this%Jflux_scri_lm(k) = lfactor * zi*m*conjg(psi)*dpsidt
        this%Eflux_scri_l(l) = this%Eflux_scri_l(l) + this%Eflux_scri_lm(k)
        this%Jflux_scri_l(l) = this%Jflux_scri_l(l) + this%Jflux_scri_lm(k)
      end associate
    end do 

    this%Eflux_hori = sum ( this%Eflux_hori_l )
    this%Jflux_hori = sum ( this%Jflux_hori_l )
    this%Eflux_scri = sum ( this%Eflux_scri_l )
    this%Jflux_scri = sum ( this%Jflux_scri_l )

  end procedure flux_extract


  module procedure flux_output

    use output_base
    use time_info
    use parameters, only : lmin, lmax, out_by_l, out_by_lm

    implicit none

    integer(ip) :: ioo_id, tmp_id, cnt, k
    integer(ip) :: momentum_output
    character(7) :: EName, JName, Hor, Inf, l, m
    character(3) :: l, m
    character(100) :: filename

    momentum_output = 1

    ioo_id = this%ioo_id
    if (ioo_id<0) then
      EName = 'EFlux_'
      JName = 'JFlux_'
      Hor = 'Hor.asc'
      Inf = 'Inf.asc'

      filename = trim(EName) // trim(Hor)
      ioo_id = next_available_io_id ()
      this%ioo_id = ioo_id
      print*,'Opening ', filename, ' with id ', ioo_id
      open(ioo_id, file=filename, status='replace', action='write')

      filename = trim(EName) // trim(Inf)
      tmp_id = next_available_io_id ()
      print*,'Opening ', filename, ' with id ', tmp_id
      open(tmp_id, file=filename, status='replace', action='write')

      if (momentum_output==1) then
        filename = trim(JName) // trim(Hor)
        tmp_id = next_available_io_id ()
        print*,'Opening ', filename, ' with id ', tmp_id
        open(tmp_id, file=filename, status='replace', action='write')
  
        filename = trim(JName) // trim(Inf)
        tmp_id = next_available_io_id ()
        print*,'Opening ', filename, ' with id ', tmp_id
        open(tmp_id, file=filename, status='replace', action='write')
      endif

      if (out_by_l) then
        do k = lmin, lmax
          write(l,"(I0)") k
  
          filename = trim(EName) // trim(l) // "_" // trim(Hor)
          tmp_id = next_available_io_id ()
          print*,'Opening ', filename, ' with id ', tmp_id
          open(tmp_id, file=filename, status='replace', action='write')
  
          filename = trim(EName) // trim(l) // "_" // trim(Inf)
          tmp_id = next_available_io_id ()
          print*,'Opening ', filename, ' with id ', tmp_id
          open(tmp_id, file=filename, status='replace', action='write')
  
          if (momentum_output==1) then
            filename = trim(JName) // trim(l) // "_" // trim(Hor)
            tmp_id = next_available_io_id ()
            print*,'Opening ', filename, ' with id ', tmp_id
            open(tmp_id, file=filename, status='replace', action='write')
  
            filename = trim(JName) // trim(l) // "_" // trim(Inf)
            tmp_id = next_available_io_id ()
            print*,'Opening ', filename, ' with id ', tmp_id
            open(tmp_id, file=filename, status='replace', action='write')
          endif
        end do
      endif

      if (out_by_lm) then
        do k = 1, this%p%nmodes
          write(l,"(I0)") this%p%ll(k)
          write(m,"(I0)") this%p%mm(k)

          filename = trim(EName) // trim(l) // "_" // trim(m) // "_" // trim(Hor)
          tmp_id = next_available_io_id ()
          print*,'Opening ', filename, ' with id ', tmp_id
          open(tmp_id, file=filename, status='replace', action='write')

          filename = trim(EName) // trim(l) // "_" // trim(m) // "_" // trim(Inf)
          tmp_id = next_available_io_id ()
          print*,'Opening ', filename, ' with id ', tmp_id
          open(tmp_id, file=filename, status='replace', action='write')

          if (momentum_output==1) then
            filename = trim(JName) // trim(l) // "_" // trim(m) // "_" // trim(Hor)
            tmp_id = next_available_io_id ()
            print*,'Opening ', filename, ' with id ', tmp_id
            open(tmp_id, file=filename, status='replace', action='write')

            filename = trim(JName) // trim(l) // "_" // trim(m) // "_" // trim(Inf)
            tmp_id = next_available_io_id ()
            print*,'Opening ', filename, ' with id ', tmp_id
            open(tmp_id, file=filename, status='replace', action='write')
          end if
        end do
      end if
    end if

    write(ioo_id,'(*(es23.15e3,1x))') get_current_time ( ), this%Eflux_hori
    write(ioo_id+1,'(*(es23.15e3,1x))') get_current_time ( ), this%Eflux_scri
    if (momentum_output==1) then
      write(ioo_id+2,'(*(es23.15e3,1x))') get_current_time ( ), this%Jflux_hori
      write(ioo_id+3,'(*(es23.15e3,1x))') get_current_time ( ), this%Jflux_scri
    endif

    cnt = 2*(1+momentum_output)

    if(out_by_l) then
      do k = lmin, lmax
        write(ioo_id+cnt,'(*(es24.16e3,1x))') get_current_time ( ), this%Eflux_hori_l(k)
        write(ioo_id+cnt+1,'(*(es23.15e3,1x))') get_current_time ( ), this%Eflux_scri_l(k)
        if (momentum_output==1) then
          write(ioo_id+cnt+2,'(*(es23.15e3,1x))') get_current_time ( ), this%Jflux_hori_l(k)
          write(ioo_id+cnt+3,'(*(es23.15e3,1x))') get_current_time ( ), this%Jflux_scri_l(k)
        endif
        cnt = cnt + 2*(1+momentum_output)
      end do
    endif

    if(out_by_lm) then
      do k = 1, this%p%nmodes
        write(ioo_id+cnt,'(*(es23.15e3,1x))') get_current_time ( ), this%Eflux_hori_lm(k)
        write(ioo_id+cnt+1,'(*(es23.15e3,1x))') get_current_time ( ), this%Eflux_scri_lm(k)
        if (momentum_output==1) then
          write(ioo_id+cnt+2,'(*(es23.15e3,1x))') get_current_time ( ), this%Jflux_hori_lm(k)
          write(ioo_id+cnt+3,'(*(es23.15e3,1x))') get_current_time ( ), this%Jflux_scri_lm(k)
        endif
        cnt = cnt + 2*(1+momentum_output)
      end do
    endif

  end procedure flux_output

  module procedure close_flux_observer
    deallocate ( this%Eflux_scri_lm, this%Jflux_scri_lm, this%Eflux_hori_lm, this%Jflux_hori_lm, &
                 this%Eflux_scri_l, this%Jflux_scri_l, this%Eflux_hori_l, this%Jflux_hori_l )
  end procedure close_flux_observer

end submodule flux_observer_implementation 
