module hio

   implicit none
   
   integer :: index = 0
   double precision, dimension (:,:,:), allocatable :: rtemp
   double complex,   dimension (:,:,:), allocatable :: ctemp

contains

subroutine io (nn, it, time)

   use horizon
   use model
   use null_params, only : iot0, iot2, it_start
   implicit none
   integer, intent (in) :: nn, it
   double precision, intent (in) :: time
   logical :: reduced_output = .false., &
              race_output = .false., &
              what_the_hell_is_this_output = .false.
   integer, save :: times_called

   double precision :: r_m, r_ml

   times_called = times_called + 1

! compute the expansion (do it here in case we don't do a 3-d evolution)
   r_m = 2. * mass
   r_ml = - time / (4. * mass)
   exp_out = (r_M * rholnew + r_Ml * rhonew) / (r_M * rhonew)

   if (times_called == 1) then

      allocate (rtemp(nn,nn,2), ctemp(nn,nn,2))

      open (unit=20, status="unknown", file="maxH.dat")
      open (unit=25, status="unknown", file="wminH.dat")
      open (unit=26, status="unknown", file="betaminH.dat")
      open (unit=27, status="unknown", file="expminH.dat")
      open (unit=28, status="unknown", file="rminH.dat")

      ! write a header line for each of these files
      write(20, *) &
           "# time, max(|rho - 1|), max(|omega|), max(|rhol|), max(|jl|)" 

      write(25, *) &
           "# time min(|omeganew(:,:,1)|), max(|omeganew(:,:,1)|)"

      write(26, *) "# time, max(|betanew(:,:,1)|)"
  
      write(27, *) &
           "# time, min(|exp_out(:,:,1)|), max(|exp_out(:,:,1)|)"
 
      write(28, *) "# time, min(|rho(:,:,1)|), max(|rho(:,:,1)|)" 

   end if

! put this so we dump every time step at the end - watch for the hardcoded no.
!  if (race_output .and. mod(it - 1,iot2) == 0) then
!  if ((race_output .and. mod(it - 1,iot2) == 0) .or. (it > 1100)) then
!  if ((race_output .and. mod(it - 1,iot2) == 0) .or. (it > 2200)) then

!! if ((race_output .and. mod(it - 1,iot2) == 0)) then
!!    call convio (nn, time)
!! end if

   if (mod(it - (it_start + 1),iot2) == 0) then
      write(*,*) 'write grid functions in hio @ time step', it, 'plot #', index

      if (reduced_output) then
        call outc(nn, index, 'jnH',  jnew(:,:,1))
        call outr(nn, index, 'rnH',  rhonew(:,:,1) - 1.)
        call outc(nn, index, 'omeganH',  omeganew(:,:,1))
        call outr(nn, index, 'rlnH', rholnew(:,:,1))
        call outc(nn, index, 'jlnH', jlnew(:,:,1))
        call outr(nn, index, 'bnH',  betanew(:,:,1))
        call outr(nn, index, 'brnH', betarnew(:,:,1))
        call outc(nn, index, 'unH',  unew(:,:,1))
        call outr(nn, index, 'thnH', t_hat(:,:,1))
        call outc(nn, index, 'jrnH', jrnew(:,:,1))
        call outc(nn, index, 'urnH', urnew(:,:,1))
        call outr(nn, index, 'vnH',  vnew(:,:,1))
        call outr(nn, index, 'vrnH', vrnew(:,:,1))
        call outr(nn, index, 'exp_outnH', exp_out(:,:,1))
      else if (what_the_hell_is_this_output) then
         call outc(nn, index, 'jnH',  jnew(:,:,1))
         call outc(nn, index, 'jsH',  jnew(:,:,2))
         
         call outr(nn, index, 'rnH',  rhonew(:,:,1) - 1.)
         call outr(nn, index, 'rsH',  rhonew(:,:,2) - 1.)
         
         call outc(nn, index, 'omeganH', omeganew(:,:,1))
         call outc(nn, index, 'omegasH', omeganew(:,:,2))
      
         call outr(nn, index, 'rlnH', rholnew(:,:,1))
         call outr(nn, index, 'rlsH', rholnew(:,:,2))
      
         call outc(nn, index, 'jlnH', jlnew(:,:,1))
         call outc(nn, index, 'jlsH', jlnew(:,:,2))

         call outr(nn, index, 'bnH', betanew(:,:,1))
         call outr(nn, index, 'bsH', betanew(:,:,2))

         call outr(nn, index, 'brnH', betarnew(:,:,1))
         call outr(nn, index, 'brsH', betarnew(:,:,2))

         call outc(nn, index, 'unH', unew(:,:,1))
         call outc(nn, index, 'usH', unew(:,:,2))

         call outr(nn, index, 'thnH', t_hat(:,:,1))
         call outr(nn, index, 'thsH', t_hat(:,:,2))

         call outc(nn, index, 'jrnH', jrnew(:,:,1))
         call outc(nn, index, 'jrsH', jrnew(:,:,2))

         call outc(nn, index, 'urnH', urnew(:,:,1))
         call outc(nn, index, 'ursH', urnew(:,:,2))
         
         call outr(nn, index, 'vnH', vnew(:,:,1))
         call outr(nn, index, 'vsH', vnew(:,:,2))

         call outr(nn, index, 'vrnH', vrnew(:,:,1))
         call outr(nn, index, 'vrsH', vrnew(:,:,2))

         call outr(nn, index, 'exp_outnH', exp_out(:,:,1))
         call outr(nn, index, 'exp_outsH', exp_out(:,:,2))

         call outr(nn, index, 'nmodjsH', dble(jnew(:,:,1)*conjg(jnew(:,:,1)) ))
         call outr(nn, index, 'smodjsH', dble(jnew(:,:,2)*conjg(jnew(:,:,2)) ))

         call outr(nn, index, 'e2bnH', e2beta(:,:,1))
         call outr(nn, index, 'e2bsH', e2beta(:,:,2))
         
         call outc(nn, index, 'qnH', qnew(:,:,1))
         call outc(nn, index, 'qsH', qnew(:,:,2))
      end if
      index = index + 1
   end if

   if (mod(it - (it_start+1), iot0) == 0) then
      write (unit = 20, fmt = '(5e13.4)') time, &
            maxval(abs(rho - 1.0)), maxval(abs(omega)), &
            maxval(abs(rhol)),      maxval(abs(jl))
      call flush (20)

     write (unit = 25, fmt = '(3E17.8)') time, &
           minval(abs(omeganew(:,:,1))),  maxval(abs(omeganew(:,:,1))) 
     call flush (25)

     write (unit = 26, fmt = '(2E17.8)') time, maxval(betanew(:,:,1))
     call Flush (26)
     
     write (unit = 27, fmt = '(3E17.8)') time, minval(exp_out(:,:,1)), &
                                               maxval(exp_out(:,:,1))
     call Flush (27)
     
     write (unit = 28, fmt = '(3E17.8)') time, minval(rho(:,:,1)), &
                                               maxval(rho(:,:,1)) 
     call Flush (28) 

   end if

end subroutine io

subroutine outr_singlefiles(nn, index, filnam, north)

   use null_grid, only : qs, ps
   implicit none

   integer, intent (in) :: nn, index
   character*(*) filnam
   double precision, dimension (nn,nn), intent (in) :: north

   integer :: i, j, iunit = 11
   character number*4

   write (number,'(i4)') 1000 + index

   open (unit = iunit, file = filnam // '.gnu', status = 'unknown', &
         position = 'append')
   write (unit = iunit, fmt = '(a)')  &
         'splot "' // filnam // '_' // number(2:4) // '.dat"'
   close (unit = iunit)

   open (unit = iunit, file = filnam // '_' // number(2:4) // '.dat', &
         status = 'unknown')

   do j = 1, nn
      do i = 1, nn
         write (unit = iunit, fmt = '(f10.6,1x,f10.6,1x,e17.8)') &
               qs(i,j), ps(i,j), north(i,j)
      end do
      write (unit = iunit, fmt = '("")')
   end do

   close (unit = iunit)
   
end subroutine outr_singlefiles

subroutine outr(nn, index, filnam, north)

   use null_grid, only : qs, ps
   implicit none

   integer, intent (in) :: nn, index
   character*(*) filnam
   double precision, dimension (nn,nn), intent (in) :: north

   integer :: i, j, iunit = 11

   open (unit = iunit, file = filnam // '.gnu', status = 'unknown', &
         position = 'append')
   write (unit = iunit, fmt = '(a,i4)')  &
         'splot "' // filnam // '.dat" index ', index
   close (unit = iunit)

   open (unit = iunit, file = filnam // '.dat', status = 'unknown', &
         position = 'append')

!  write (unit = iunit, fmt = '("# ",e21.13e3)') time

   do j = 1, nn
      do i = 1, nn
         write (unit = iunit, fmt = '(f10.6,1x,f10.6,1x,e17.8)') &
               qs(i,j), ps(i,j), north(i,j)
      end do
      write (unit = iunit, fmt = '("")')
   end do
   write (unit = iunit, fmt = '("")')

   close (unit = iunit)
   
end subroutine outr

subroutine outc_reim(nn, index, filnam, north)

   use null_grid, only : qs, ps
   implicit none

   integer, intent (in) :: nn, index
   character*(*) filnam
   double complex, dimension (nn,nn), intent (in) :: north

   integer :: i, j, iunit = 11

   open (unit = iunit + 1, file = filnam // 'r.gnu', status = 'unknown', &
         position = 'append')
   open (unit = iunit + 2, file = filnam // 'i.gnu', status = 'unknown', &
         position = 'append')
   write (unit = iunit + 1, fmt = '(a,i4)')  &
         'splot "' // filnam // 'r.dat" index ', index
   write (unit = iunit + 2, fmt = '(a,i4)')  &
         'splot "' // filnam // 'i.dat" index ', index
   close (unit = iunit + 1)
   close (unit = iunit + 2)

   open (unit = iunit + 1, file = filnam // 'r.dat', status = 'unknown', &
         position = 'append')
   open (unit = iunit + 2, file = filnam // 'i.dat', status = 'unknown', &
         position = 'append')

   do j = 1, nn
      do i = 1, nn
         write (unit = iunit + 1, fmt = '(f10.6,1x,f10.6,1x,e17.8)') &
               qs(i,j), ps(i,j), dble(north(i,j))

         write (unit = iunit + 2, fmt = '(f10.6,1x,f10.6,1x,e17.8)') &
               qs(i,j), ps(i,j), dimag(north(i,j))


      end do
      do i = 1, 2
         write (unit = iunit + i, fmt = '("")')
      end do
   end do

   do i = 1, 2
      write (unit = iunit + i, fmt = '("")')
      close (unit = iunit + i)
   end do
   
end subroutine outc_reim

subroutine outc(nn, index, filnam, north)
! dumps real part and modulus instead of re and im
   use null_grid, only : qs, ps
   implicit none

   integer, intent (in) :: nn, index
   character*(*) filnam
   double complex, dimension (nn,nn), intent (in) :: north

   integer :: i, j, iunit = 11

   open (unit = iunit + 1, file = filnam // 'r.gnu', status = 'unknown', &
         position = 'append')
   open (unit = iunit + 2, file = filnam // 'mod.gnu', status = 'unknown', &
         position = 'append')
   write (unit = iunit + 1, fmt = '(a,i4)')  &
         'splot "' // filnam // 'r.dat" index ', index
   write (unit = iunit + 2, fmt = '(a,i4)')  &
         'splot "' // filnam // 'mod.dat" index ', index
   close (unit = iunit + 1)
   close (unit = iunit + 2)

   open (unit = iunit + 1, file = filnam // 'r.dat', status = 'unknown', &
         position = 'append')
   open (unit = iunit + 2, file = filnam // 'mod.dat', status = 'unknown', &
         position = 'append')

   do j = 1, nn
      do i = 1, nn
         write (unit = iunit + 1, fmt = '(f10.6,1x,f10.6,1x,3e13.4)') &
               qs(i,j), ps(i,j), dble(north(i,j))

         write (unit = iunit + 2, fmt = '(f10.6,1x,f10.6,1x,3e13.4)') &
               qs(i,j), ps(i,j), sqrt(dble(north(i,j)*conjg(north(i,j))))

      end do
      do i = 1, 2
         write (unit = iunit + i, fmt = '("")')
      end do
   end do

   do i = 1, 2
      write (unit = iunit + i, fmt = '("")')
      close (unit = iunit + i)
   end do
   
end subroutine outc

subroutine gio (nn, it, time)

   use null_params, only: iot2, it_start
   use horizon
   use model
   use geomview

   implicit none
   integer, intent (in) :: nn, it
   double precision, intent (in) :: time

   double precision scale, offset

   if ( mod(it - (it_start+1), iot2 ) == 0) then

      scale = 3.0/maxval(t_hat(:,:,2))
      offset = 0.0d0
      call zmeshr(nn, it, scale, offset, 'th', t_hat(:,:,2))

      scale = 3.0/maxval(dble(jnew(:,:,2)))
      offset = 0.0d0
      call zmeshr(nn, it, scale, offset, 'modjs', &
                         &  dble(jnew(:,:,2)*conjg(jnew(:,:,2))))

      scale = 3.0/maxval(dble(jnew(:,:,2)))
      offset = 0.0d0
      call zmeshc(nn, it, scale, offset, 'js', jnew(:,:,2))

      scale = 3.0/maxval(dble(jrnew(:,:,2)))
      offset = 0.0d0
      call zmeshc(nn, it, scale, offset, 'jrs', jrnew(:,:,2))

      scale = 6.0/maxval(dble(omeganew(:,:,2)))
      offset = 0.0d0
      call zmeshc(nn, it, scale, offset, 'oms', omeganew(:,:,2))

      scale = 3.0/maxval(betanew(:,:,2))
      offset = -1.6d0
      call zmeshr(nn, it, scale, offset, 'bs', betanew(:,:,2))

      scale = 3.0/maxval(betarnew(:,:,2))
      offset = 0.0d0
      call zmeshr(nn, it, scale, offset, 'brs', betarnew(:,:,2))

      scale = 3.0/maxval(dble(unew(:,:,2)))
      offset = 0.0d0
      call zmeshc(nn, it, scale, offset, 'us', unew(:,:,2))

      scale =  3.0/maxval(dble(qnew(:,:,2)))
      offset = 0.0d0
      call zmeshc(nn, it, scale, offset, 'qs', qnew(:,:,2))

      scale = 3.0/maxval(dble(urnew(:,:,2)))
      offset = 0.0d0
      call zmeshc(nn, it, scale, offset, 'urs', urnew(:,:,2))

      scale =  3.0/maxval(vnew(:,:,2))
      offset = 0.0d0
      call zmeshr(nn, it, scale, offset, 'vs', vnew(:,:,2))

      scale =  3.0/maxval(vrnew(:,:,2))
      offset = 0.0d0
      call zmeshr(nn, it, scale, offset, 'vr', vrnew(:,:,2))

   end if

end subroutine gio

subroutine raceio (nn, time)

   use null_grid, only : qs, ps
   use horizon
   implicit none

   integer nn
   double precision time

   integer ki, kj 
   double precision normxi2, costheta, theta, phi, pi

   pi = 4. * atan(1.0d0)

   write (*,*) 'twist and shear at: ', time

   open (unit=71, status='unknown', position='append', file='exp.dat')
   open (unit=72, status='unknown', position='append', file='rho.dat')
   open (unit=73, status='unknown', position='append', file='ricci.dat')
   open (unit=74, status='unknown', position='append', file='twist.dat')
   open (unit=75, status='unknown', position='append', file='shearp.dat')
   open (unit=76, status='unknown', position='append', file='omega.dat')
   open (unit=77, status='unknown', position='append', file='rhol.dat')
   open (unit=78, status='unknown', position='append', file='jl.dat')

   kj = (nn+1)/2
   do ki = (nn+1)/2, nn
      normxi2 = qs(ki,kj)**2 + ps(ki,kj)**2
      costheta = (1. - normxi2) / (1. + normxi2)
      theta = acos(costheta)
      phi = 0.5 * pi - theta	! the GNU polar angle

      write (unit=71, fmt='(3(1x,e23.13))') theta, time, exp_out(ki,kj,1)
      write (unit=72, fmt='(3(1x,e23.13))') theta, time, rho(ki,kj,1) - 1.
      write (unit=73, fmt='(3(1x,e23.13))') theta, time, ricci(ki,kj,1)
      write (unit=74, fmt='(3(1x,e23.13))') theta, time, twist_theta(ki,kj,1)
      write (unit=75, fmt='(3(1x,e23.13))') theta, time, shearplus(ki,kj,1)
      write (unit=76, fmt='(6(1x,e23.13))') theta, time, &
            dble(omeganew(ki,kj,1)), dimag(omeganew(ki,kj,1)), &
            dble(debug_omega(ki,kj,1)), dimag(debug_omega(ki,kj,1))
      write (unit=77, fmt='(3(1x,e23.13))') theta, time, rhol(ki,kj,1)
      write (unit=78, fmt='(3(1x,e23.13))') theta, time, abs(jl(ki,kj,1))

   end do

   do ki = 71, 78
!     write (unit = ki, fmt = '("&",e17.8)') time
      write (unit = ki, fmt = '("")')
      close (unit = ki)
   end do

end subroutine raceio

subroutine convio (nn, time)

   use null_grid, only : qs, ps
   use horizon
   use model, only : io_conv_skip
   implicit none

   integer nn
   double precision time

   integer ki, kj 
   double precision normxi2, costheta, theta, phi, pi

   pi = 4. * atan(1.0d0)

   write (*,*) 'twist and shear at: ', time

   open (unit=71, status='unknown', position='append', file='exp.dat')
   open (unit=72, status='unknown', position='append', file='rho.dat')
   open (unit=73, status='unknown', position='append', file='ricci.dat')
   open (unit=74, status='unknown', position='append', file='twist.dat')
   open (unit=75, status='unknown', position='append', file='shearp.dat')
   open (unit=76, status='unknown', position='append', file='omega.dat')
   open (unit=77, status='unknown', position='append', file='rhol.dat')
   open (unit=78, status='unknown', position='append', file='jl.dat')

   kj = (nn+1)/2
   do ki = (nn+1)/2, nn - 2, io_conv_skip
      normxi2 = qs(ki,kj)**2 + ps(ki,kj)**2
      costheta = (1. - normxi2) / (1. + normxi2)
      theta = acos(costheta)
      phi = 0.5 * pi - theta	! the GNU polar angle

      write (unit=71, fmt='(2(1x,e23.13))') theta, exp_out(ki,kj,1)
      write (unit=72, fmt='(2(1x,e23.13))') theta, rho(ki,kj,1) - 1.
      write (unit=73, fmt='(2(1x,e23.13))') theta, ricci(ki,kj,1)
      write (unit=74, fmt='(2(1x,e23.13))') theta, twist_theta(ki,kj,1)
      write (unit=75, fmt='(2(1x,e23.13))') theta, shearplus(ki,kj,1)
      write (unit=76, fmt='(5(1x,e23.13))') theta, &
            dble(omeganew(ki,kj,1)), dimag(omeganew(ki,kj,1)), &
            dble(debug_omega(ki,kj,1)), dimag(debug_omega(ki,kj,1))
      write (unit=77, fmt='(2(1x,e23.13))') theta, rhol(ki,kj,1)
      write (unit=78, fmt='(2(1x,e23.13))') theta, abs(jl(ki,kj,1))

   end do

   do ki = 71, 78
      write (unit = ki, fmt = '("&")')
      close (unit = ki)
   end do

end subroutine convio

end module hio
