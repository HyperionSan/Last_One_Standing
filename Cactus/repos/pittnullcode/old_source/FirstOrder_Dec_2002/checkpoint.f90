module checkpoint
  use checkpoint_defs 
  use horizon
  use null_vars
  use null_params
  use cauchydump
  use news2
  implicit none
contains
  subroutine check_allocate
    use null_params
    implicit none
    allocate(checkTimeLvls(0:nt+1))   ! we will assume it_start = 0 or 1
    checkTimeLvls = 0.0d0
  end subroutine check_allocate
  subroutine check_point
!    use news2
    implicit none 
    
    ! save  global params ( see null.f90
    !   time, dt, it, nt ONE, ZERO, HALF, EDGE, HFED .....
    ! save horizon
    ! boundaryu (sep from hor, esp mask)
    ! save null code
    ! save news
    checkpointcount  = mod( checkpointcount+1, 2)
    if (checkpointcount .eq. 0) then
      folderu = folder0
    else
      folderu = folder1
    endif

    open(unit=979, file=folderu // 'success', status='replace')
    write(979,*) 'bad'
    close(979) 
   
    call main_check
    call horizon_check
    if (calc_news) then
      call news_check   !located in news module
    end if
    if (alreadydumped) then
      call cauchy_check
    endif
    open(unit=979, file=folderu // 'success', status='replace')
    write(979,*) 'good'
    close(979) 

  end subroutine check_point
  subroutine cauchy_check
    implicit none
    integer :: test, ll
    double precision  :: tm 
! here we wil simply count the number of time steps saved in the cauchy data
! and record it. The recover routine will then extract this many timesteps
! from the old *.cauchy files to make the new ones. It is up to the user
! to place the old cauchy files into the 'recover' directory. We dont do it
! here to save disk space.
    ll = 0 
    open(unit=770, file='time.cauchy', FORM='unformatted', status='old')

    do 
      read(770, IOSTAT=test) tm
      if (test == -1) EXIT
      ll = ll + 1
   end do

    close(770)

   open(unit=770, file=folderu // 'CauchyLevels', &
                       status='replace', FORM='unformatted')   
   write(770) ll
   close(770)
  end subroutine cauchy_check
  subroutine ccauchy_copy(iunit, ounit, complx, ll)
    implicit none
    integer, intent(in) :: iunit, ounit
    double complex, dimension(stdang, stdang, stdrad), intent(inout) :: complx
    integer, intent(in) :: ll
    integer :: test, i
    do i=1, ll
      read(iunit, IOSTAT=test) complx
      if (test.eq.-1) then
        stop 'FILE ERROR'
      end if
      write(ounit) complx
    end do 
  end subroutine ccauchy_copy
  subroutine rcauchy_copy(iunit, ounit, real, ll)
    implicit none
    integer, intent(in) :: iunit, ounit
    double precision, dimension(stdang, stdang, stdrad), intent(inout) :: real
    integer, intent(in) :: ll
    integer :: test, i
    do i=1, ll
      read(iunit, IOSTAT=test) real
      if (test.eq.-1) then
        stop 'FILE ERROR'
      end if
      write(ounit) real
    end do 
  end subroutine rcauchy_copy
subroutine main_check
  use null_vars
  use null_grid
  implicit none
 
  open(unit=MJnFile, file=folderu // 'MainJn', &
                       status='replace', FORM='unformatted')   
  open(unit=MnunFile, file=folderu // 'Mainnun', &
                       status='replace', FORM='unformatted')   
  open(unit=McknFile, file=folderu // 'Mainckn', &
                       status='replace', FORM='unformatted')   
  open(unit=MbnFile, file=folderu // 'Mainbn', &
                       status='replace', FORM='unformatted')   
  open(unit=McbnFile, file=folderu // 'Maincbn', &
                       status='replace', FORM='unformatted')   
  open(unit=MunFile, file=folderu // 'Mainun', &
                       status='replace', FORM='unformatted')   
  open(unit=MwnFile, file=folderu // 'Mainwn', &
                       status='replace', FORM='unformatted')   
  open(unit=MMlossFile, file=folderu // 'MainMloss', &
                       status='replace', FORM='unformatted')   
  open(unit=MVarsFile, file=folderu // 'MainVars', &
                       status='replace', FORM='unformatted')   
  open(unit=MTimeLevels, file=folderu // 'TimeLevels', &
                       status='replace', FORM='unformatted')   


  open (unit=MVarsAsciiFile, file=folderu // 'MainVarsAscii',  status='replace')
  write(MVarsAsciiFile, *) time, it, dt, alreadydumped


   write(MJnFile) jon

   write(MnunFile) nuon 

   write(McknFile) ckon

   write(MbnFile) bon

   write(McbnFile) cbon

   write(MunFile) uon

   write(MwnFile) won

   write(MMlossFile) mloss
   write(MVarsFile) time, it, dt, alreadydumped
   write(MTimeLevels) checkTimeLvls

   close( MJnFile  )
   close( MnunFile )
   close( McknFile )
   close( MbnFile  )
   close( McbnFile )
   close( MunFile  )
   close( MwnFile  )
   close( MMlossFile )
   close (MVarsFile)
   close (MTimeLevels)

  close (MVarsAsciiFile)
end subroutine main_check

subroutine horizon_check
  use horizon
  implicit none
  open(unit=HJFile, file=folderu // 'HJ', &
                       status='replace', FORM='unformatted')
  open(unit=HrhoFile, file=folderu // 'Hrho', &
                       status='replace', FORM='unformatted')
  open(unit=HrhodotFile, file=folderu // 'Hrhodot', &
                       status='replace', FORM='unformatted')
  open(unit=HomegaFile, file=folderu // 'Homega', &
                       status='replace', FORM='unformatted')
  open(unit=HrholFile, file=folderu // 'Hrhol', &
                       status='replace', FORM='unformatted')
  open(unit=HJlFile, file=folderu // 'HJl', &
                       status='replace', FORM='unformatted')


  write(HJFile) j
  write(HrhoFile) rho
  write(HrhodotFile) rhodot
  write(HomegaFile) omega
  write(HrholFile) rhol
  write(HJlFile) jl

  close(HJFile)
  close(HrhoFile)
  close(HrhodotFile)
  close(HomegaFile)
  close(HrholFile)
  close(HJlFile)

end subroutine horizon_check

  subroutine recover
    use news2
    implicit none
    ! recover global params
    ! re-allocate horizon variables and fill them in
    ! re-allocate null variables and fill them in
    ! re-allocate news and fill in
    call main_recover
    if (alreadydumped) then
      call cauchy_recover
    endif
    if (calc_news) then
      call recover_news  ! locates in news module (because stuff is private there)
    end if
    ! horizon_recover will be called  by horizon module
  end subroutine


subroutine main_recover
  use null_vars
  use null_grid
  use null_code
 
  implicit none


  open(unit=MJnFile, file=folderr // 'MainJn', &
                       status='old', FORM='unformatted')   
  open(unit=MnunFile, file=folderr // 'Mainnun', &
                       status='old', FORM='unformatted')   
  open(unit=McknFile, file=folderr // 'Mainckn', &
                       status='old', FORM='unformatted')   
  open(unit=MbnFile, file=folderr // 'Mainbn', &
                       status='old', FORM='unformatted')   
  open(unit=McbnFile, file=folderr // 'Maincbn', &
                       status='old', FORM='unformatted')   
  open(unit=MunFile, file=folderr // 'Mainun', &
                       status='old', FORM='unformatted')   
  open(unit=MwnFile, file=folderr // 'Mainwn', &
                       status='old', FORM='unformatted')   
  open(unit=MMlossFile, file=folderr // 'MainMloss', &
                       status='old', FORM='unformatted')   
  open (unit=MVarsFile, file=folderr // 'MainVars',  &  
                       status='old', FORM='unformatted')
  open(unit=MTimeLevels, file=folderr // 'TimeLevels', &
                       status='old', FORM='unformatted')   


   read(MJnFile) jon

   read(MnunFile) nuon 

   read(McknFile) ckon

   read(MbnFile) bon

   read(McbnFile) cbon

   read(MunFile) uon

   read(MwnFile) won

   read(MMlossFile) mloss
   read(MVarsFile) time, it, dt, alreadydumped
   read(MTimeLevels) checkTimeLvls

   close( MJnFile  )
   close( MnunFile )
   close( McknFile )
   close( MbnFile  )
   close( McbnFile )
   close( MunFile  )
   close( MwnFile  )
   close( MMlossFile )
   close (MVarsFile)
   close (MTimeLevels)
 
end subroutine main_recover

 subroutine recover_horizon
  use horizon
  implicit none
  open(unit=HJFile, file=folderr // 'HJ', &
                       status='old', FORM='unformatted')
  open(unit=HrhoFile, file=folderr // 'Hrho', &
                       status='old', FORM='unformatted')
  open(unit=HrhodotFile, file=folderr // 'Hrhodot', &
                       status='old', FORM='unformatted')
  open(unit=HomegaFile, file=folderr // 'Homega', &
                       status='old', FORM='unformatted')
  open(unit=HrholFile, file=folderr // 'Hrhol', &
                       status='old', FORM='unformatted')
  open(unit=HJlFile, file=folderr // 'HJl', &
                       status='old', FORM='unformatted')


  read(HJFile) j
  read(HrhoFile) rho
  read(HrhodotFile) rhodot
  read(HomegaFile) omega
  read(HrholFile) rhol
  read(HJlFile) jl


  close(HJFile)
  close(HrhoFile)
  close(HomegaFile)
  close(HrholFile)
  close(HJlFile)

end subroutine recover_horizon
  subroutine cauchy_recover
    implicit none
    double precision, dimension(stdang, stdang, stdrad) :: real
    double complex, dimension(stdang, stdang, stdrad) :: complx
    integer :: test, ll, i
    double precision :: tm
  
    open(unit=770, file=folderr // 'CauchyLevels', FORM='unformatted', status='old')
    read(770) ll
    close(770)
 
    open(unit=770, file='j.cauchy', FORM='unformatted', status='replace')
    open(unit=880, file=folderr // 'j.cauchy', FORM='unformatted', status='old')
    open(unit=771, file='u.cauchy', FORM='unformatted', status='replace')
    open(unit=881, file=folderr // 'u.cauchy', FORM='unformatted', status='old')
    open(unit=772, file='b.cauchy', FORM='unformatted', status='replace')
    open(unit=882, file=folderr // 'b.cauchy', FORM='unformatted', status='old')
    open(unit=773, file='w.cauchy', FORM='unformatted', status='replace')
    open(unit=883, file=folderr // 'w.cauchy', FORM='unformatted', status='old')
    open(unit=774, file='time.cauchy', FORM='unformatted', status='replace')
    open(unit=884, file=folderr // 'time.cauchy', FORM='unformatted',status='old')
    open(unit=775, file='atime.cauchy', status='replace')
    open(unit=885, file=folderr // 'atime.cauchy', status='old')

    do i = 1, ll
      read(884, IOSTAT=test) tm
      if (test == -1) then
        stop 'FILE ERROR'
      endif
      write(774) tm

      read(885, *, IOSTAT=test) tm
      if (test == -1) then
        stop 'FILE ERROR'
      endif
      write(775, *) tm
    end do

 
    call ccauchy_copy(880, 770, complx, ll)
    call ccauchy_copy(881, 771, complx, ll)
    call rcauchy_copy(882, 772, real, ll)
    call rcauchy_copy(883, 773, real, ll)

    close(770)
    close(880)
    close(771)
    close(881)
    close(772)
    close(882)
    close(773)
    close(883)
    close(774)
    close(784)
    close(775)
    close(785)

  end subroutine cauchy_recover

end module checkpoint
