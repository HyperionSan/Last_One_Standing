module cauchydump
  use null_params
  use null_vars
  implicit none
  integer, save :: uunit = 601
  integer, save :: wunit = 602
  integer, save :: junit = 603
  integer, save :: bunit = 604
contains
  subroutine cauchy_dump
    use null_grid
    use checkpoint_defs
    implicit none
    logical, save :: neverbeencalled  = .true.
    integer       :: uunit = 601
    integer       :: wunit = 603
    integer       :: junit = 605
    integer       :: bunit = 607
    integer       :: tunit = 609
    integer       :: atunit = 610
    double precision, dimension(stdang,stdang,stdrad) :: realarray
    double complex, dimension(stdang,stdang,stdrad)  :: complxarray
    double complex, dimension(nn,nn,nx) :: uint  
    double complex, dimension(nn,nn,nx) :: cmplxcheck
    double precision, dimension(nn,nn,nx) :: realcheck
    integer :: kk

    realarray = 0.0d0 
    complxarray = (0.0d0, 0.0d0) 
    cmplxcheck=0.0d0
    realcheck=0.0d0
    if (alreadydumped) then   ! if we are recovering from a terminated run
      neverbeencalled = .false.   ! then *.cauchy will contain run data so
    endif                         ! don't delete
    if ( neverbeencalled ) then
      neverbeencalled = .false.
      open(unit=uunit, file='u.cauchy', &
                       status='replace', FORM='unformatted')   
      open(unit=wunit, file='w.cauchy', &
                       status='replace', FORM='unformatted')   
      open(unit=junit, file='j.cauchy', &
                       status='replace', FORM='unformatted')   
      open(unit=bunit, file='b.cauchy', &
                       status='replace', FORM='unformatted')

      open(unit=tunit, file='time.cauchy', &
                       status='replace', FORM='unformatted')
      open(unit=atunit, file='atime.cauchy', &
                       status='replace')
    else
      open(unit=uunit, file='u.cauchy', &
                       status='old', FORM='unformatted', position='APPEND')   
      open(unit=wunit, file='w.cauchy', &
                       status='old', FORM='unformatted', position='APPEND')   
      open(unit=junit, file='j.cauchy', &
                       status='old', FORM='unformatted', position='APPEND')   
      open(unit=bunit, file='b.cauchy', &
                       status='old', FORM='unformatted', position='APPEND')

      open(unit=tunit, file='time.cauchy', &
                       status='old', FORM='unformatted', position='APPEND')
      open(unit=atunit, file='atime.cauchy', &
                       status='old', position='APPEND')
     end if 
    
    call dump_interp_u (unn, uint)

    call trimcmplxarray(uint, complxarray)
    write(uunit) complxarray

       
    call trimrealarray(wnn, realarray)
    write(wunit) realarray 


    call trimcmplxarray(jnn, complxarray)
    write(junit) complxarray 

    call trimrealarray(bnn, realarray)
    write(bunit) realarray 

    write(tunit) time
    write(atunit, *) time

   
    close(uunit)
    close(wunit)
    close(junit)
    close(bunit)
    close(tunit)
    close(atunit)

    alreadydumped = .true.  ! tell the checkpointer that we have saved some cauchy stuff
  end subroutine cauchy_dump
 
  subroutine   dump_interp_u (unn, uint)
    use null_params
    use null_interp
    use null_grid
    implicit none
    double complex, dimension(nn,nn,nx), intent(in) :: unn
    double complex, dimension(nn,nn,nx), intent(out) :: uint
    double precision :: xdesire, xh1, xh2, xh3, xh4, c1, c2, c3, c4
    integer :: kk
    uint(:,:,1) = 0.0d0  ! dont extrapolate

    kk = 2
    xdesire = x(2)
    xh1 = xh(1)
    xh2 = xh(2)
    xh3 = xh(3)
    xh4 = xh(4)


    c1 = (xdesire - xh2) * (xdesire - xh3) * (xdesire - xh4) /  (  &
         (xh1 - xh2) * (xh1 - xh3) * (xh1 - xh4)    )

    c2 = (xdesire - xh1) * (xdesire - xh3) * (xdesire - xh4) /  (  &
         (xh2 - xh1) * (xh2 - xh3) * (xh2 - xh4)    )

    c3 = (xdesire - xh1) * (xdesire - xh2) * (xdesire - xh4) /  (  &
         (xh3 - xh1) * (xh3 - xh2) * (xh3 - xh4)    )

    c4 = (xdesire - xh1) * (xdesire - xh2) * (xdesire - xh3) /  (  &
         (xh4 - xh1) * (xh4 - xh2) * (xh4 - xh3)    )

    uint(:,:,2) = c1*unn(:,:, 1) + c2 * unn(:,:,2) + c3 * unn(:,:,3) + &
                  c4 * unn(:,:,4)
    do kk = 3, nx - 1

      xdesire = x(kk)
      xh1 = xh(kk - 2) 
      xh2 = xh(kk - 1) 
      xh3 = xh(kk    ) 
      xh4 = xh(kk + 1) 
 
      c1 = (xdesire - xh2) * (xdesire - xh3) * (xdesire - xh4) /  (  &
         (xh1 - xh2) * (xh1 - xh3) * (xh1 - xh4)    )

      c2 = (xdesire - xh1) * (xdesire - xh3) * (xdesire - xh4) /  (  &
         (xh2 - xh1) * (xh2 - xh3) * (xh2 - xh4)    )

      c3 = (xdesire - xh1) * (xdesire - xh2) * (xdesire - xh4) /  (  &
         (xh3 - xh1) * (xh3 - xh2) * (xh3 - xh4)    )

      c4 = (xdesire - xh1) * (xdesire - xh2) * (xdesire - xh3) /  (  &
         (xh4 - xh1) * (xh4 - xh2) * (xh4 - xh3)    )

      uint(:,:, kk) = unn(:,:, kk -2) * c1 + unn(:,:, kk -1) * c2 + &
               unn(:,:, kk) * c3 + unn(:,:, kk +1) * c4
    end do
    xh1 = xh(nx - 3)
    xh2 = xh(nx - 2)
    xh3 = xh(nx - 1)
    xh4 = xh(nx    )

    xdesire = x(nx)

    c1 = (xdesire - xh2) * (xdesire - xh3) * (xdesire - xh4) /  (  &
         (xh1 - xh2) * (xh1 - xh3) * (xh1 - xh4)    )

    c2 = (xdesire - xh1) * (xdesire - xh3) * (xdesire - xh4) /  (  &
         (xh2 - xh1) * (xh2 - xh3) * (xh2 - xh4)    )

    c3 = (xdesire - xh1) * (xdesire - xh2) * (xdesire - xh4) /  (  &
         (xh3 - xh1) * (xh3 - xh2) * (xh3 - xh4)    )

    c4 = (xdesire - xh1) * (xdesire - xh2) * (xdesire - xh3) /  (  &
         (xh4 - xh1) * (xh4 - xh2) * (xh4 - xh3)    )

    uint(:,:,nx) = unn(:,:, nx - 3) * c1 + unn(:,:, nx - 2) * c2 + &
               unn(:,:, nx - 1) * c3 + unn(:,:, nx) * c4
  end subroutine   dump_interp_u
  subroutine trimrealarray(bin, bout)
    use null_params
    implicit none
    double precision, dimension(nn,nn,nx), intent(in) :: bin
    double precision, dimension(stdang,stdang,stdrad), intent(out) :: bout
    integer :: ll, mm, rr, ii, jj, kk

    do rr = 1, stdrad
      kk  = (rr-1)* cdumpj +1+NIntPts   ! here we dump the first radial point
      do ll = 1, stdang
        do mm = 1, stdang
           ii = (ll-1) * cdumpj + 3
           jj = (mm-1) * cdumpj + 3
           bout(ll, mm, rr) = bin(ii,jj,kk)
        end do
      end do
    end do
  end subroutine trimrealarray
  subroutine trimcmplxarray(jin, jout)
    use null_params
    implicit none
    double complex, dimension(nn,nn,nx), intent(in) :: jin
    double complex, dimension(stdang,stdang,stdrad), intent(out) :: jout
    integer :: ll, mm, rr, ii, jj, kk

    do rr = 1, stdrad
      kk  = (rr-1)* cdumpj +1+NIntPts   ! here we dump the first radial point
      do ll = 1, stdang
        do mm = 1, stdang
           ii = (ll-1) * cdumpj + 3
           jj = (mm-1) * cdumpj + 3
           jout(ll, mm, rr) = jin(ii,jj,kk)
        end do
      end do
    end do
  end subroutine trimcmplxarray
   
end module cauchydump

