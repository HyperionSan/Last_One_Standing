program null

  use null_params, only: time, nt, iot0, iot2, it_start, pure_schwarzs
  use null_grid
  use null_analytic
  use null_code
  use null_io
  !use flexio
  use null_vars
  use ascii_io
  use null_cfl_test

  ! to couple the null code to the boundary & expansion routines
  use model, only: mass, null_evolution_switch, p_switch, eps
  use hdrive
  use boundary
  use pureschwarzsboundary

  ! for the news
  use news2
  use StereoIntegration
  use particle
  use checkpoint

  implicit none

  ! news output
  integer poutx, pouty, it_begin
  double precision :: cfl_dt(3)
  integer:: LN
  double precision, dimension(:), allocatable :: LX, LY
  logical, dimension(:,:,:), allocatable :: LMask
  print *,'<<< PITTcode $Id: null.f90,v 1.xx 2002/05/21 17:05:57 yosef Exp $ >>>'
  print *

!call errset(12, 257, 257, 0, 1, 1)
! call errset(176, 257, 257, 0, 1, 1)   !
 ! call errset(1398, 257, 257, 0, 1, 1)
  ! initialize parameters and allocate variables
  call null_read_params
  call null_allocate
  call boundary_allocate
  call null_setup_grid
  call anSetup
  if (calc_news) then
    call news_allocate
  end if

  jnn = (0., 0.)
  nunn = (0., 0.)
  bnn = 0.
  cbnn = (0., 0.)
  cknn = (0., 0.)
  unn = (0., 0.)
  wnn = 0.

  if (do_checkpoint) then
    checkpointrecover=chkrcvr
    call check_allocate   ! this must preceed any recovery step
  else
    checkpointrecover=.false.
  endif

  time_real_start = time   ! this is the time that null.in says was the first time step
  real_start_dt = dt_fix   ! this will have to changed in the cvs code
  if (checkpointrecover) then
    ! horizon module will unset checkpointrecover
    call recover   ! this is only a partial recovery, rest is done in hdriver
    it_begin = it  ! it and dt were read in from the recovery files
    dt_fix = dt     ! if we supress cfl then this is usefull
    cfl_dt(1) = dt / cfl
    cfl_dt(2) = dt / cfl
  else
    it_begin = it_start
    dt = dt_fix
  endif

  ! for News output
  ! poutx = nn-2      ! alternative: poutx = int( 1.5/dd ) + 3
  ! pouty = (nn+1)/2  ! alternative: pouty = poutx 

  ! Compute the indices of the point for which the News will be written.
  ! Take the closest point to, and inside the equator, along the q=p line.

  open (unit = 15, file = 'pout.dat',  status = 'unknown')
  write (unit = 15, fmt = '(a)') '# I  SQRT[ qs(I,I)**2 + ps(I,I)**2 ]'
  do poutx = 1, nn
     pouty = poutx
     write (unit = 15, fmt = '(i3,1x,e16.5)') &
        poutx, sqrt(qs(poutx,pouty)**2+ps(poutx,pouty)**2)
     if (sqrt(qs(poutx,pouty)**2+ps(poutx,pouty)**2) .le. 1.) exit ! close enough to 1.0 ...
  end do
  close (unit = 15)

  print *, "News poutx: ", poutx
  print *, "News pouty: ", pouty
  print *, sqrt(qs(poutx,pouty)**2+ps(poutx,pouty)**2)

  open (unit = 91, file = 'newseq.dat', status = 'unknown')
  open (unit = 92, file = 'mloss.dat',  status = 'unknown')

  EDGE  = nn - 2 
  ONE  = nint((1.0d0 + qsize) / dd)  + 3    ! the real equator, provided that
  ZERO = 3 + (nn - 5) / 2             ! (1+qsize)/(2 qsize) * (nn -5 ) is an int
  HALF = nint((.5d0 + qsize) / dd) + 3      ! the real halfway point provided (.5+qsize) ...
  HFED =  3 + 3 * ( ( nn - 5 ) / 4 )
  PNTT = nint((.2d0 + qsize) / dd) + 3
  write(*,*)  qs(ONE,ZERO), qs(ZERO, ZERO), qs(HALF,ZERO)
  write(*,*)  ps(ZERO, ONE), ps(ZERO, ZERO), ps(ZERO, HALF)
  write(*,*)  qs(HFED, ZERO), qs(EDGE, ZERO)
  write(*,*)  dx, dd, dt

  ! here is some pseudo code explaining the main loop
  !
  ! the variable time is the time coresponding to the new level
  ! 
  ! On the first time slice the 'new' level is filled in (top filled)
  !
  !
  ! DO
  !   time = time of level to be filled
  !   ConformalModel(time) get horizon quantities at time 'time'  
  !   HorizonExpansion (fill in beta, w, U, J) on 'new' level
  !      InitialDATA   Fill in J   ( first time only) on 'new' level
  !   HypersurfaceAndEvolution (no evolution on first iteration)
  !   News  (only 2nd iteration and on)
  !   cfl_calc (get dt to next level)
  !   Output
  !   NullCopy
  !   time = time + dt
  ! LOOP
  
  print*, 'start time = ', time
  print*, 'delta t    = ', dt

  print*
  print*, '<<< starting time loop >>>'
  print*
  if (pure_schwarzs) then
    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print*
    print*, '<<<WARNING HORIZON IS SCHWARZCHILD >>>'
    print*, 'ALL OTHER MODEL PARAMETERS ARE IGNORED'
    print*
    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  end if
  do it = it_begin, nt
     if ( do_checkpoint ) then 
       CheckTimeLvls(it) = time     ! lets hope no one wants it_start < 0
       if ( (it > it_begin .AND. (.NOT. checkpointrecover)) .AND. &
                           mod(it, it_checkpoint) == 0) then
         call check_point
       endif 
     end if

     ! 'time' is the time of the new level (the one not yet filled). 
     print*, '*** step it =', it, '@ t =', time, ': dt=', dt, '***'

     if (pure_schwarzs) then 
       call schwarzschild_driver
       checkpointrecover = .false.   ! we have to uset this since hdriver 
     else                            ! is not called
       call hdriver (time, dt, nn, it, nt)
     end if

     if (null_evolution_switch == 1) then
        if(.NOT.pure_schwarzs) then
          call boundary_expansion ( mass )
        end if

        if (it == it_start) then
           call null_data
           call null_evolve     
           if(p_switch) call particle_initialize
           print*, '<<< initializions done >>>'
        else
           call null_evolve
           if (calc_news) then
             call getnews


             ! mass loss between two consecutive time slices
             mloss(1) = StereoInt (nn, dMdOmega(:,:,1) ) / (4.0d0 * pi) ! north 
             mloss(2) = StereoInt (nn, dMdOmega(:,:,2) ) / (4.0d0 * pi) ! south 

             mloss(0) = mloss(1) + mloss(2)               ! both patches
             mloss(3) = mloss(3) + mloss(0) * dt          ! cumulative dM
          end if
         ! Particle motion, if required
          if (p_switch) call particle_evolve(it - it_start)
        end if

        ! compute CFL on 1st level too
        if (cfl_compute) then
               if (mod(it - it_start, it_cfl) == 0) then
                  if(it.eq.it_start) then
                    cfl_dt(1) = dt / cfl
	          end if
                  call null_cfl0(jnn,  bnn, unn, wnn, cfl_dt(1), maskn)

                  cfl_dt(3) = cfl * cfl_dt(1)
                  print *, "CFL: ", it, time, dt, cfl_dt(3)
                  dt = min(dt, cfl_dt(3))
               end if
         else
            dt = dt_fix
        endif

        ! Bondi news at one point  
        if (calc_news) then
          if (mod(it, iot0) == 0 .and. it > it_start + 1) then
             write(unit = 91, fmt = '(4e23.13)')               &
                  uBondi(poutx, pouty, 1),                     &
                  dble (NewsB(poutx, pouty, 1)),               &
                  dimag(NewsB(poutx, pouty, 1)), time_of_news
             call flush (91)
             write(unit = 92, fmt = '(3e23.13)')               &
                &                   time_of_news, mloss(3), mloss(0)
             call flush (92)
!             call LinfNews(NewsB, time_of_news)
          end if
        end if        
       ! output and copy timelevels new -> old
        call null_io_ascii
        call point_io
       !call write_flexio

       if (mod(it - it_start, iot2) ==0 ) then
          call mid_dump
       end if
        call null_copy
     end if
       
     time = time + dt
     print*
  end do

  print*, '<<< Null code finished - graceful exit >>>'
contains
  
  subroutine mid_dump
    use ascii_io
    implicit none
    integer :: imid, ret

    imid = (nx + 1 + NIntPts) /2

    ret = gft_write('midJJb', time, dble(jnn(:,:,imid)*conjg(jnn(:,:,imid))))
    ret = gft_write('midUUb', time, dble(unn(:,:,imid)*conjg(unn(:,:,imid))))
    ret = gft_write('midiU', time, dimag(unn(:,:,imid)))
    ret = gft_write('midB', time, dble(bnn(:,:,imid)))
    ret = gft_write('midW', time, dble(wnn(:,:,imid)))
    
  end subroutine mid_dump
  subroutine point_io
    use point_dump
    use cauchydump
    implicit none

    double complex, dimension(nn,nn) :: jxscri, jxmid, uscri, umid
    integer :: imid
    
!    imid = (nx + 1 + NIntPts) /2


!    jxmid(:,:) = (.5d0 / dx ) * (jnn(:,:,imid+1) - jnn(:,:,imid-1))
!    jxscri(:,:)=(1.0d0/dx)*(1.5d0*jnn(:,:,nx) -2.0d0*jnn(:,:,nx-1) +0.5d0*jnn(:,:,nx-2) )

!    call interp_u(imid, uscri, umid)
    
!    call pointdump(time, dble(jnn(:,:,nx)), 'rJScri')
!    call pointdump(time, dimag(jnn(:,:,nx)), 'iJScri')
!    call pointdump(time, dble(jnn(:,:,imid)), 'rJMid')
!    call pointdump(time, dimag(jnn(:,:,imid)), 'iJMid')
!    call pointdump(time, dble(jxscri(:,:)), 'rJxScri')
!    call pointdump(time, dimag(jxscri(:,:)), 'iJxScri')
!    call pointdump(time, dble(jxmid(:,:)), 'rJxMid')
!    call pointdump(time, dimag(jxmid(:,:)), 'iJxMid')
!    call pointdump(time, dble(bnn(:,:,nx)), 'BScri')
!    call pointdump(time, dble(bnn(:,:,imid)), 'BMid')
!    call pointdump(time, dble(cbnn(:,:,nx)), 'rcBScri')
!    call pointdump(time, dimag(cbnn(:,:,nx)), 'icBScri')
!    call pointdump(time, dble(cbnn(:,:,imid)), 'rcBMid')
!    call pointdump(time, dimag(cbnn(:,:,imid)), 'icBMid')
!    call pointdump(time, dble(cknn(:,:,nx)), 'rcKScri')
!    call pointdump(time, dimag(cknn(:,:,nx)), 'icKScri')
!    call pointdump(time, dble(cknn(:,:,imid)), 'rcKMid')
!    call pointdump(time, dimag(cknn(:,:,imid)), 'icKMid')
!    call pointdump(time, dble(nunn(:,:,nx)), 'rnuScri')
!    call pointdump(time, dimag(nunn(:,:,nx)), 'inuScri')
!    call pointdump(time, dble(nunn(:,:,imid)), 'rnuMid')
!    call pointdump(time, dimag(nunn(:,:,imid)), 'inuMid')
!    call pointdump(time, dble(uscri(:,:)), 'rUScri')
!    call pointdump(time, dimag(uscri(:,:)), 'iUScri')
!    call pointdump(time, dble(umid(:,:)), 'rUMid')
!    call pointdump(time, dimag(umid(:,:)), 'iUMid')
!    call pointdump(time, dble(wnn(:,:,nx)), 'WScri')
!    call pointdump(time, dble(wnn(:,:,imid)), 'WMid')
!    call pointdump(time, dble(qs(:,:)), 'qq')
!    call pointdump(time, dble(ps(:,:)), 'pp')
    if (it > it_start .AND. calc_news ) then
      call pointdump(time_of_news, dble(NewsB(:,:,1)), 'reNEWS')
      call pointdump(time_of_news, dimag(NewsB(:,:,1)), 'imNEWS')
    endif
!    if (mod(it - it_start, cdumpj*cdumpit) == 0 .and. docdump ) then
!        call cauchy_dump
!    end if
  end subroutine point_io
  
  subroutine interp_u(imid, uscri, umid)
    implicit none
    integer, intent(in) :: imid
    double complex, dimension(nn, nn),  intent(out) :: uscri, umid 
    double precision :: xdesire, xhmid, xh1, xh2, xh3, xh4
    double precision :: c1, c2, c3, c4


    xh1 = xh(imid - 2)
    xh2 = xh(imid - 1)
    xh3 = xh(imid    )
    xh4 = xh(imid + 1)

    xdesire = x(imid)

    c1 = (xdesire - xh2) * (xdesire - xh3) * (xdesire - xh4) /  (  &
         (xh1 - xh2) * (xh1 - xh3) * (xh1 - xh4)    )

    c2 = (xdesire - xh1) * (xdesire - xh3) * (xdesire - xh4) /  (  &
         (xh2 - xh1) * (xh2 - xh3) * (xh2 - xh4)    )

    c3 = (xdesire - xh1) * (xdesire - xh2) * (xdesire - xh4) /  (  &
         (xh3 - xh1) * (xh3 - xh2) * (xh3 - xh4)    )

    c4 = (xdesire - xh1) * (xdesire - xh2) * (xdesire - xh3) /  (  &
         (xh4 - xh1) * (xh4 - xh2) * (xh4 - xh3)    )

    umid(:,:) = unn(:,:, imid -2) * c1 + unn(:,:, imid -1) * c2 + &
               unn(:,:, imid) * c3 + unn(:,:, imid +1) * c4

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

    uscri(:,:) = unn(:,:, nx - 3) * c1 + unn(:,:, nx - 2) * c2 + &
               unn(:,:, nx - 1) * c3 + unn(:,:, nx) * c4
  end subroutine interp_u
  subroutine LinfNews( NB, x )
    implicit none
    double complex, dimension(nn,nn,2), intent(in) :: NB
    double precision, intent(in) :: x
    
    double complex, dimension(nn, nn,2 ) :: anNews
    integer :: i, j
    integer, save :: count = 0
    double precision :: amp, linf
    double precision :: x1, x2, x3, x4, y1, y2, y3, y4
    do
      count = count + 1
      if (count > LN -1 ) stop 'Out of News'
      if (LX(count) > x) EXIT
    end do
    if ( count > 2 ) then   
      ! interp anNews onto time
       x1 = LX(count-2)
       x2 = LX(count-1)
       x3 = LX(count)
       x4 = LX(count+1)

       y1 = LY(count-2)
       y2 = LY(count-1)
       y3 = LY(count)
       y4 = LY(count+1)

       amp = y1 * ( x - x2 ) * ( x - x3 ) * ( x - x4 ) / ( &
                  (x1 - x2 ) * (x1 - x3 ) * (x1 - x4 ) ) + &
             y2 * ( x - x1 ) * ( x - x3 ) * ( x - x4 ) / ( &
                  (x2 - x1 ) * (x2 - x3 ) * (x2 - x4 ) ) + &
             y3 * ( x - x1 ) * ( x - x2 ) * ( x - x4 ) / ( &
                  (x3 - x1 ) * (x3 - x2 ) * (x3 - x4 ) ) + &
             y4 * ( x - x1 ) * ( x - x2 ) * ( x - x3 ) / ( &
                  (x4 - x1 ) * (x4 - x2 ) * (x4 - x3 ) )
     else 
       amp = 0.0d0
     endif
     amp = amp * eps !/ 16.0d0
     anNews(:,:,1) = amp * 4.0*(z*zb)/(1+z*zb)**2 
     anNews(:,:,2) = amp * 4.0*(z*zb)/(1+z*zb)**2 
     linf =  maxval(abs (NB - anNews) , LMask)
     write(101, *) x, linf, amp
     write(102,*)  x, anNews(poutx, pouty,1), NB(poutx, pouty, 1)
     write(103, *) x, maxloc(abs (NB - anNews), LMask )
     call flush(101)
     call flush(102)
     call flush(103)
  end subroutine LinfNews

  subroutine anSetup
    implicit none
    integer :: i, j
    double precision :: x, y
    open (unit=110, status='old', file = 'annews', IOSTAT=i)
    if (i == -1) stop 'No News'
    j = 0
    do
      read(110, *, IOSTAT=i) x, y
      if ( i == -1) EXIT
      j = j + 1
    end do
    close(110)
    LN = j
    allocate(LX(LN), LY(LN))
    allocate( LMask(nn,nn,2) )
    open (unit=110, status='old', file = 'annews', IOSTAT=i)  
    do j=1, LN
      read(110, *, IOSTAT=i) x, y
      LX(j)  = x
      LY(j) = y
    end do
    LMask = .false.
    do i=1, nn
      do j=1, nn
        if ( sqrt((qs(i,j)**2 + ps(i,j)**2))  < 1.00000000001d0) then 
          LMask(i,j,1) = .true.
          LMask(i,j,2) = .true.
        end if
      end do
    end do
    close (110)
  end subroutine anSetup


end program null
