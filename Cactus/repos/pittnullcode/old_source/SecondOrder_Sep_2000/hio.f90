module hio

   implicit none

   integer :: index = 0
   double precision, dimension (:,:,:), allocatable :: rtemp
   double complex,   dimension (:,:,:), allocatable :: ctemp

contains

subroutine io (nn, it, time)

   use horizon
   use model
   implicit none
   integer, intent (in) :: nn, it
   double precision, intent (in) :: time

   logical :: reduced_output = .true.

   if (it == 1) then

      allocate (rtemp(nn,nn,2), ctemp(nn,nn,2))

      Open ( unit = 20, file = "max.dat",      status = "unknown" )
 !     Open ( unit = 21, file = "omega.dat",    status = "unknown" )
 !     Open ( unit = 22, file = "rhol.dat",     status = "unknown" )
 !     Open ( unit = 23, file = "jl.dat",       status = "unknown" )
 !     Open ( unit = 24, file = "jl.data",      status = "unknown" )
      Open ( unit = 25, file = "wminmax.dat",  status = "unknown" )
      Open ( unit = 26, file = "betamax.dat",  status = "unknown" )
      Open ( unit = 27, file = "expmin.dat",   status = "unknown" )
      Open ( unit = 28, file = "rmin.dat",     status = "unknown" )
   end if

   if (mod(it - 1,iot2) == 0) then
      write(*,*) 'write grid functions in hio @ time step', it, 'plot #', index

      if (reduced_output) then
         call outc(nn, index, 'jn',  jnew(:,:,1))
         call outr(nn, index, 'rn',  rhonew(:,:,1) - 1.)
         call outc(nn, index, 'wn',  omeganew(:,:,1))
         call outc(nn, index, 'ws',  omeganew(:,:,2))
         call outr(nn, index, 'rln', rholnew(:,:,1))
         call outr(nn, index, 'rls', rholnew(:,:,2))
         call outc(nn, index, 'jln', jlnew(:,:,1))
!         call outr(nn, index, 'bn',  betanew(:,:,1))
!         call outr(nn, index, 'brn', betarnew(:,:,1))
!         call outc(nn, index, 'un',  unew(:,:,1))
         call outr(nn, index, 'thn', t_hat(:,:,1))
!         call outc(nn, index, 'jrn', jrnew(:,:,1))
!         call outc(nn, index, 'urn', urnew(:,:,1))
!         call outr(nn, index, 'vn',  vnew(:,:,1))
!         call outr(nn, index, 'vrn', vrnew(:,:,1))
         call outr(nn, index, 'exp_outn', exp_out(:,:,1))
      else
         call outc(nn, index, 'jn',  jnew(:,:,1))
         call outc(nn, index, 'js',  jnew(:,:,2))
         
         call outr(nn, index, 'rn',  rhonew(:,:,1) - 1.)
         call outr(nn, index, 'rs',  rhonew(:,:,2) - 1.)
         
         call outc(nn, index, 'wn',  omeganew(:,:,1))
         call outc(nn, index, 'ws',  omeganew(:,:,2))
      
         call outr(nn, index, 'rln', rholnew(:,:,1))
         call outr(nn, index, 'rls', rholnew(:,:,2))
      
         call outc(nn, index, 'jln', jlnew(:,:,1))
         call outc(nn, index, 'jls', jlnew(:,:,2))

         call outr(nn, index, 'bn', betanew(:,:,1))
         call outr(nn, index, 'bs', betanew(:,:,2))

         call outr(nn, index, 'brn', betarnew(:,:,1))
         call outr(nn, index, 'brs', betarnew(:,:,2))

         call outc(nn, index, 'un', unew(:,:,1))
         call outc(nn, index, 'us', unew(:,:,2))

         call outr(nn, index, 'thn', t_hat(:,:,1))
         call outr(nn, index, 'ths', t_hat(:,:,2))

         call outc(nn, index, 'jrn', jrnew(:,:,1))
         call outc(nn, index, 'jrs', jrnew(:,:,2))

         call outc(nn, index, 'urn', urnew(:,:,1))
         call outc(nn, index, 'urs', urnew(:,:,2))
         
         call outr(nn, index, 'vn', vnew(:,:,1))
         call outr(nn, index, 'vs', vnew(:,:,2))

         call outr(nn, index, 'vrn', vrnew(:,:,1))
         call outr(nn, index, 'vrs', vrnew(:,:,2))

         call outr(nn, index, 'exp_outn', exp_out(:,:,1))
         call outr(nn, index, 'exp_outs', exp_out(:,:,2))

         !     call model_twist (nn, time, j, ctemp)
         !     call outc(nn, index, 'awn', ctemp(:,:,1))
         !     call outc(nn, index, 'aws', ctemp(:,:,2))

         !     call model_rho_lambda (nn, time, rtemp)
         !     call outr(nn, index, 'arln', rtemp(:,:,1))
         !     call outr(nn, index, 'arls', rtemp(:,:,2))
         !     call outr(nn, index, 'erln', rhol(:,:,1) - rtemp(:,:,1))
         !     call outr(nn, index, 'erls', rhol(:,:,2) - rtemp(:,:,2))
         
         !     call model_j_lambda (nn, time, ctemp)
         !     call outc(nn, index, 'ajln', ctemp(:,:,1))
         !     call outc(nn, index, 'ajls', ctemp(:,:,2))
         
         !      call outc(nn, index, 'ejln', jl(:,:,1) - ctemp(:,:,1))
         !      call outc(nn, index, 'ejls', jl(:,:,2) - ctemp(:,:,2))

         call outr(nn, index, 'nmodjs', & 
              & dble(jnew(:,:,1)*conjg(jnew(:,:,1)) ))
         call outr(nn, index, 'smodjs', &
              & dble(jnew(:,:,2)*conjg(jnew(:,:,2)) ))

         !     call model_j_r (nn, time, ctemp)
         !     call outc(nn, index, 'ejrn', jrnew(:,:,1) - ctemp(:,:,1))
         !     call outc(nn, index, 'ejrs', jrnew(:,:,2) - ctemp(:,:,2))
         
         call outr(nn, index, 'e2bn', e2beta(:,:,1))
         call outr(nn, index, 'e2bs', e2beta(:,:,2))
         
         call outc(nn, index, 'qn', qnew(:,:,1))
         call outc(nn, index, 'qs', qnew(:,:,2))
      end if
      index = index + 1
   end if

   if (mod(it - 1,iot0) == 0) then
      write (unit = 20, fmt = '(5e13.4)') time, &
            maxval(abs(rho - 1.0)), maxval(abs(omega)), &
            maxval(abs(rhol)),      maxval(abs(jl))
      call flush (20)

!     call model_twist (nn, time, j, ctemp)
!     write (unit = 21, fmt = '(4e13.4)') time, &
!           maxval(abs(omega)), maxval(abs(ctemp)), maxval(abs(omega - ctemp))
!     call flush (21)

!     call model_rho_lambda (nn, time, rtemp)
!     write (unit = 22, fmt = '(4e13.4)') time, &
!           maxval(abs(rhol)), maxval(abs(rtemp)), maxval(abs(rhol - rtemp))
!     call flush (22)

!     call model_j_lambda (nn, time, ctemp)
!     write (unit = 23, fmt = '(4e13.4)') time, &
!           maxval(abs(jl)), maxval(abs(ctemp)), maxval(abs(jl - ctemp))
!     call flush (23)

!     write (unit = 24, fmt = '(3d23.13e3)') time,    &
!           dble (jl(nn-2,nn-2,2)), dimag(jl(nn-2,nn-2,2))
!     call flush (24)

     write (unit = 25, fmt = '(3E17.8)') time, &
           minval(abs(omeganew(:,:,1))),  maxval(abs(omeganew(:,:,1))) 
     call flush (25)

     Write (unit = 26, fmt = '(2E17.8)') time, maxval(betanew(:,:,1))
     Call Flush (26)
     
     Write (unit = 27, fmt = '(3E17.8)') time, minval(exp_out(:,:,1)), &
                                               maxval(exp_out(:,:,1))
     Call Flush (27)
     
     Write (unit = 28, fmt = '(3E17.8)') time, minval(rho(:,:,1)), &
                                               maxval(rho(:,:,1)) 
     Call Flush (28) 

   end if

end subroutine io

subroutine outr(nn, index, filnam, north)

   implicit none

   integer, intent (in) :: nn, index
   character*(*) filnam
   double precision, dimension (nn,nn), intent (in) :: north

   integer :: i, j, iunit = 11
   double precision :: dd, qs, ps

   open (unit = iunit, file = filnam // '.gnu', status = 'unknown', &
         position = 'append')
   write (unit = iunit, fmt = '(a,i4)')  &
         'splot "' // filnam // '.dat" index ', index
   close (unit = iunit)

   open (unit = iunit, file = filnam // '.dat', status = 'unknown', &
         position = 'append')

!  write (unit = iunit, fmt = '("# ",e21.13e3)') time

   dd = 2.0 / dble(nn - 5)

   do j = 1, nn
      ps = -1.0 + (j-3) * dd
      do i = 1, nn
         qs = -1.0 + (i-3) * dd
         write (unit = iunit, fmt = '(f10.6,1x,f10.6,1x,e17.8)') &
               qs, ps, north(i,j)
      end do
      write (unit = iunit, fmt = '("")')
   end do
   write (unit = iunit, fmt = '("")')

   close (unit = iunit)
   
end subroutine outr

subroutine outc_reim(nn, index, filnam, north)

   implicit none

   integer, intent (in) :: nn, index
   character*(*) filnam
   double complex, dimension (nn,nn), intent (in) :: north

   integer :: i, j, iunit = 11
   double precision :: dd, qs, ps

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

   dd = 2.0 / dble(nn - 5)

   do j = 1, nn
      ps = -1.0 + (j-3) * dd
      do i = 1, nn
         qs = -1.0 + (i-3) * dd
         write (unit = iunit + 1, fmt = '(f10.6,1x,f10.6,1x,e17.8)') &
               qs, ps, dble(north(i,j))

         write (unit = iunit + 2, fmt = '(f10.6,1x,f10.6,1x,e17.8)') &
               qs, ps, dimag(north(i,j))


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
   implicit none

   integer, intent (in) :: nn, index
   character*(*) filnam
   double complex, dimension (nn,nn), intent (in) :: north

   integer :: i, j, iunit = 11
   double precision :: dd, qs, ps

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

   dd = 2.0 / dble(nn - 5)

   do j = 1, nn
      ps = -1.0 + (j-3) * dd
      do i = 1, nn
         qs = -1.0 + (i-3) * dd
         write (unit = iunit + 1, fmt = '(f10.6,1x,f10.6,1x,3e13.4)') &
               qs, ps, dble(north(i,j))

         write (unit = iunit + 2, fmt = '(f10.6,1x,f10.6,1x,3e13.4)') &
               qs, ps, sqrt(dble(north(i,j)*conjg(north(i,j))))

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

   use horizon
   use model
   use geomview

   implicit none
   integer, intent (in) :: nn, it
   double precision, intent (in) :: time

   double precision scale, offset

   if (mod(it,iot2) == 0) then

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

end module hio
