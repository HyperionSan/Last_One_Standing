module numservicef90

 interface

  subroutine vintegrate(invect, outvect, a, b, debug)
   ! purpose: integrate a vector invect of funtion values
   !          at equally spaced sbscissae between and b

   implicit none

   double precision, intent(in),  dimension(:)            :: invect
   double precision, intent(out), dimension(size(invect)) :: outvect
   double precision, intent(in)                :: a, b
   integer,          intent(in)                :: debug
  end subroutine vintegrate

  subroutine vintegratex(invect, outvect, xvect, debug)
    ! purpose: integrate a vector invect of funtion values
    !          at non-equally space abscissae between and b
    implicit none
    double precision, intent(in),  dimension(:)            :: invect, xvect
    double precision, intent(out), dimension(size(invect)) :: outvect
    integer,          intent(in)                           :: debug
  end subroutine vintegratex

  subroutine vdiff(invect, outvect, a, b, debug)
   ! purpose: integrate a vector invect of funtion values
   !          at equally spaced sbscissae between and b

   implicit none

   double precision, intent(in),  dimension(:)            :: invect
   double precision, intent(out), dimension(size(invect)) :: outvect
   double precision, intent(in)                :: a, b
   integer,          intent(in)                :: debug
  end subroutine vdiff

  subroutine vdif2(invect, outvect, a, b, debug)
   ! purpose: integrate a vector invect of funtion values
   !          at equally spaced sbscissae between and b

   implicit none

   double precision, intent(in),  dimension(:)            :: invect
   double precision, intent(out), dimension(size(invect)) :: outvect
   double precision, intent(in)                :: a, b
   integer,          intent(in)                :: debug
  end subroutine vdif2     

  subroutine pvdiff(par, invect, outvect, a, b, debug)
   ! purpose: integrate a vector invect of funtion values
   !          at equally spaced sbscissae between and b

   implicit none

   double precision, intent(in),  dimension(:)            :: invect
   double precision, intent(out), dimension(size(invect)) :: outvect
   double precision, intent(in)                :: a, b
   integer,          intent(in)                :: par, debug
  end subroutine pvdiff

  subroutine pvdif2(par, invect, outvect, a, b, debug)
   ! purpose: integrate a vector invect of funtion values
   !          at equally spaced sbscissae between and b

   implicit none

   double precision, intent(in),  dimension(:)            :: invect
   double precision, intent(out), dimension(size(invect)) :: outvect
   double precision, intent(in)                :: a, b
   integer,          intent(in)                :: par, debug
  end subroutine pvdif2     

  subroutine vdiffx(invect, outvect, xvect, debug)
    ! purpose: differentiate a vector invect of funtion values,
    !          with non-equally spaced abscissae xvect
    implicit none
    double precision, intent(in),  dimension(:)            :: invect, xvect
    double precision, intent(out), dimension(size(invect)) :: outvect
    integer,          intent(in)                           :: debug
  end subroutine vdiffx

  subroutine xyz2cmesh(x, y, z, r, g, b, o, nx, nz, name, unit)
   ! purpose:         
   use prec
   implicit none   
   real(kind=wp),    intent(in), dimension (:,:)   :: x, y, z, r, g, b, o
   integer,          intent(in)                    :: nx, nz, unit
   character(len=*), intent(in)                    :: name
  end subroutine xyz2cmesh

  subroutine a_xyz2cmesh(x, y, z, r, g, b, o, nx, nz, name, unit)
   ! purpose:         
   use prec
   implicit none   
   real(kind=wp),    intent(in), dimension (:,:)   :: x, y, z, r, g, b, o
   integer,          intent(in)                    :: nx, nz, unit
   character(len=*), intent(in)                    :: name
  end subroutine a_xyz2cmesh
 
  subroutine add2cmesh_2Dto3D(x, y, r, g, b, o, nx, nphi, unit)
   ! purpose:  write one picture in a movie in Geomview's CMESH format        
   use prec
   implicit none
   real(kind=wp),    intent(in), dimension (:)     :: x, y, r, g, b, o
   integer,          intent(in)                    :: nx, nphi, unit
  end subroutine add2cmesh_2Dto3D


  subroutine addslice2cmesh(x, y, r, g, b, o, nx, index, nt, fname)
   !============================================================
   ! purpose:  write one slice of an object Geomview's CMESH format
   !           this is particularly intended for stacking a pair-of-pants   
   use prec, only: wp
   implicit none
   real(kind=wp),    intent(in), dimension (:)     :: x, y ! x, y coords
   real(kind=wp),    intent(in), dimension (:)     :: r, g, b, o
   integer,          intent(in)                    :: nx, index, nt
   character(len=*), intent(in)                    :: fname
  end subroutine addslice2cmesh

end interface
end module numservicef90

subroutine vintegrate(invect, outvect, a, b, debug)
! purpose: integrate a vector invect of funtion values
!          at equally space abscissae between and b
! method: simpson's rule, see e.g. Numerical Recipes/F77 p.126         

implicit none

double precision, intent(in),  dimension(:)            :: invect
double precision, intent(out), dimension(size(invect)) :: outvect
double precision, intent(in)                           :: a, b
integer,          intent(in)                           :: debug

integer :: steps, intervals
double precision :: h, integral3points
integer :: i

!! executable statemets

steps = size(invect)
intervals = steps-1
h = (b-a)/dble(intervals) 

outvect(1) = 0.0d0
outvect(2) = 0.5d0*h*(invect(1) + invect(2))

do i = 3, steps
   integral3points = h*(invect(i-2) + 4.0d0*invect(i-1) + invect(i) )/3.0d0
   outvect (i) = outvect(i-2) + integral3points
end do


if(debug==1) then
  write(*,*) 'here comes array outvect from subroutine vintegrate'
  do i = 1, steps
      write(*,*) outvect(i)
  end do
  write(*,*)   'this was array outvect from subroutine vintegrate'
end if

return

end subroutine vintegrate

subroutine vintegratex(invect, outvect, xvect, debug)
! purpose: integrate a vector invect of funtion values
!          at non-equally space abscissae between and b
! method:  

implicit none

double precision, intent(in),  dimension(:)            :: invect, xvect
double precision, intent(out), dimension(size(invect)) :: outvect
integer,          intent(in)                           :: debug

integer :: steps, intervals
integer :: i

!! executable statemets

steps = size(invect) ! intervals = steps-1

outvect(1) = 0.0d0

do i = 2, steps
   outvect (i) = outvect(i-1) + 0.5*( invect(i) + invect(i-1) ) &
                              * ( xvect(i) - xvect(i-1) )
end do 

if(debug==1) then
  write(*,*) 'here comes array outvect from subroutine vintegrate'
  do i = 1, steps
      write(*,*) outvect(i)
  end do
  write(*,*)   'this was array outvect from subroutine vintegrate'
end if

return

end subroutine vintegratex


subroutine vdiff(invect, outvect, a, b, debug)
! purpose: differentiate a vector invect of funtion values
!          at equally spaced abscissae between and b
! method:  centered differences, upwind/downwind at boundary 
implicit none

double precision, intent(in),  dimension(:)            :: invect
double precision, intent(out), dimension(size(invect)) :: outvect
double precision, intent(in)                           :: a, b
integer,          intent(in)                           :: debug

integer :: steps, intervals
double precision :: h, delta2
integer :: i
!! executable statements

steps = size(invect)
intervals = steps-1
h = (b-a)/dble(intervals) 
delta2 = 0.5d0/h

outvect(1)     = (-3.0d0*invect(1) + 4.0d0*invect(2) - invect(3))*delta2 
outvect(steps) =-(-3.0d0*invect(steps) + 4.0d0*invect(steps-1) &
                       - invect(steps-2))*delta2

do i = 2, steps - 1
   outvect(i) = ( invect(i+1) - invect(i-1) ) * delta2
end do
end subroutine vdiff

subroutine vdif2(invect, outvect, a, b, debug)
! purpose: differentiate a vector invect of funtion values
!          _twice_ at equally spaced abscissae between and b
! method:  centered differences, upwind/downwind at boundary
implicit none

double precision, intent(in),  dimension(:)            :: invect
double precision, intent(out), dimension(size(invect)) :: outvect
double precision, intent(in)                           :: a, b
integer,          intent(in)                           :: debug

integer :: steps, intervals
double precision :: h, delta2
integer :: i
!! executable statemets

steps = size(invect)
intervals = steps-1
h = (b-a)/dble(intervals)
delta2 = 1.0d0/h**2

do i = 2, steps - 1
   outvect(i) = ( invect(i-1) - 2.0d0*invect(i) + invect(i+1) ) * delta2
end do

outvect(1)     = (-invect(4)       + 4.0d0*invect(3)         &
                      -5.0d0*invect(2)       + 2.0d0*invect(1))     * delta2
outvect(steps) = (-invect(steps-3) + 4.0d0*invect(steps-2)   &
      &               -5.0d0*invect(steps-1) + 2.0d0*invect(steps)) * delta2

end subroutine vdif2      

subroutine pvdiff(par, invect, outvect, a, b, debug)
! purpose: differentiate a vector invect of funtion values
!          at equally spaced abscissae between and b
! method:  centered differences, upwind/downwind at boundary 
implicit none

double precision, intent(in),  dimension(:)            :: invect
double precision, intent(out), dimension(size(invect)) :: outvect
double precision, intent(in)                           :: a, b
integer,          intent(in)                           :: par, debug

integer :: steps, intervals
double precision :: h, delta2
integer :: i
!! executable statements

steps = size(invect)
intervals = steps-1
h = (b-a)/dble(intervals) 
delta2 = 0.5d0/h

outvect(1)     = (invect(2) - invect(steps-1) ) * delta2 
outvect(steps) = - outvect(1) * par 

do i = 2, steps - 1
   outvect(i) = ( invect(i+1) - invect(i-1) ) * delta2
end do
end subroutine pvdiff

subroutine pvdif2(par, invect, outvect, a, b, debug)
! purpose: differentiate a vector invect of funtion values
!          _twice_ at equally spaced abscissae between and b
! method:  centered differences, upwind/downwind at boundary
implicit none

double precision, intent(in),  dimension(:)            :: invect
double precision, intent(out), dimension(size(invect)) :: outvect
double precision, intent(in)                           :: a, b
integer,          intent(in)                           :: par, debug

integer :: steps, intervals
double precision :: h, delta2
integer :: i
!! executable statemets

steps = size(invect)
intervals = steps-1
h = (b-a)/dble(intervals)
delta2 = 1.0d0/h**2

do i = 2, steps - 1
   outvect(i) = ( invect(i-1) - 2.0d0*invect(i) + invect(i+1) ) * delta2
end do

outvect(1)     =  ( invect(steps-1) - 2.0d0*invect(1) + invect(2) )  * delta2
outvect(steps) =  outvect(1) * par

end subroutine pvdif2      

subroutine vdiffx(invect, outvect, xvect, debug)
! purpose: differentiate a vector invect of funtion values
!          at abscissae values xvect
! method:  centered differences, upwind/downwind at boundary 
! details: we compute the derivative from three unequally
!          spaced data points by second order finite differencing
! 
!          f'(x) = c1 f(x + a1) + c2 f(x) + c3 f(x - a2)
! 
!          where: f(x + a1) = f1, f(x) = f2, f(x - a2) = f3


implicit none

double precision, intent(in),  dimension(:)            :: invect, xvect
double precision, intent(out), dimension(size(invect)) :: outvect
integer,          intent(in)                           :: debug

integer          :: steps
double precision :: a1, a2, c1, c2, c3, f1, f2, f3
integer          :: i

!! executable statemets !!

steps = size(invect)

outvect(1)     = (invect(2) - invect(1))/(xvect(2) - xvect(1))
outvect(steps) = (invect(steps) - invect(steps-1))/ &
                 (xvect(steps)  - xvect(steps-1))
!FIXEME: make the above 2nd order!

do i = 2, steps - 1
   a1 = xvect(i+1) - xvect(i) 
   a2 = xvect(i)   - xvect(i-1)
   c1 =  a2/(a1 * (a1 + a2))
   c2 =  (a1 - a2)/(a1 * a2)
   c3 = -a1/(a2 * (a1 + a2))
   f1 = invect(i+1) 
   f2 = invect(i)
   f3 = invect(i-1)

   outvect(i) = c1 * f1 + c2 * f2 + c3 * f3
end do


if(debug == 1) then
  write(*,*) 'here comes array outvect from subroutine vdiff'
  do i = 1, steps
      write(*,*) outvect(i)
  end do
  write(*,*)   'this was array outvect from subroutine vdiff'
end if

return

end subroutine vdiffx


subroutine tdiffint(n)

use numservicef90
implicit none

integer, intent(in) :: n
double precision, dimension(:), allocatable  :: invect, outvect, temp, delta, &
                                                xvect
double precision :: h, x
integer :: i

allocate(xvect(n))
allocate(invect(n))
allocate(outvect(n))
allocate(temp(n))
allocate(delta(n))

h = 5.0d0/dble(n-1)

do i = 1, n
   x = dble(i)*h
   invect(i) = sin(x)
   xvect(i) = x
end do

temp = invect
call vintegrate(invect,outvect,0.0d0,5.0d0,0)

do i = 1, n
   write(11,*) xvect(i), outvect(i)
end do

invect = outvect
call vdiff(invect,outvect,0.0d0,5.0d0,0)

do i = 1, n
   write(12,*) xvect(i), outvect(i)
end do

delta = outvect - temp

do i = 1, n
   write(13,*) xvect(i), delta(i)
end do

deallocate(xvect)
deallocate(invect)
deallocate(outvect)
deallocate(temp)
deallocate(delta)

end subroutine tdiffint

integer function fmod(i,j)

integer, intent(in)  :: i, j

fmod = mod(i,j)

return
end function fmod


subroutine xyz2cmesh(x, y, z, r, g, b, o, nx, nz, name, unit)
! purpose:         
use prec
implicit none

real(kind=wp),    intent(in), dimension (:,:)   :: x, y, z, r, g, b, o
character(len=*), intent(in)                    :: name
integer,          intent(in)                    :: nx, nz, unit
integer                                         :: ix, iz

open(unit, file = name)

write(unit, *) '{CMESH'
write(unit, 998) nx, nz

do iz = 1, nz
   do ix = 1, nx
     write(unit,999) x(ix, iz), y(ix, iz), z(ix, iz),                &
               &     r(ix, iz), g(ix, iz), b(ix, iz), o(ix, iz)
   enddo
enddo

write(unit,*) '}'
close(unit)
return
998 format(2I3)
999 format(3E13.5,4F6.2) 
end subroutine xyz2cmesh


subroutine a_xyz2cmesh(x, y, z, r, g, b, o, nx, nz, name, unit)
! purpose:         
use prec
implicit none

real(kind=wp),    intent(in), dimension (:,:)   :: x, y, z, r, g, b, o
character(len=*), intent(in)                    :: name
integer,          intent(in)                    :: nx, nz, unit
integer                                         :: ix, iz

open(unit, file = name)

write(unit, *) '{CMESH'
write(unit, 998) 2*nx-1, nz

do iz = 1, nz
   do ix = 1, nx
     write(unit,999)  x(ix, iz), y(ix, iz), z(ix, iz),               &
          &           r(ix, iz), g(ix, iz), b(ix, iz), o(ix, iz)
   enddo
enddo

do iz = nz, 1, -1
   do ix = nx-1, 1, -1
     write(unit,999) -x(ix, iz), y(ix, iz), z(ix, iz),               &
          &           r(ix, iz), g(ix, iz), b(ix, iz), o(ix, iz)
   enddo
enddo

write(unit,*) '}'
close(unit)

return
998 format(2I3)
999 format(3E13.5,4F6.2) 
end subroutine a_xyz2cmesh


subroutine add2cmesh_2Dto3D(x, y, r, g, b, o, nx, nphi, unit)
!============================================================
! purpose:  write one picture in a movie in Geomview's CMESH format        
use prec
implicit none

real(kind=wp),    intent(in), dimension (:)     :: x, y, r, g, b, o
integer,          intent(in)                    :: nx, nphi, unit

integer                                         :: ix, iphi
real(kind=wp)                                   :: dphi, phi

!! === executable statemets === !!

dphi = 6.2831853072_wp/dble(nphi-1)

write(unit, *) '{CMESH'
write(unit, 998) nphi, modulo(nx,20)

do ix = 1, nx, 20  ! just an ad-hoc reduction of output, need a better scheme!
   do iphi = 1, nphi
      phi = dphi*(iphi-1.0_wp)
      write(unit,999) x(ix)*cos(phi), &
           &          x(ix)*sin(phi), &
           &          y(ix),          &
           &          r(ix), g(ix), b(ix), o(ix)
   enddo
enddo

write(unit,*) '}'

return
998 format(2I5)
999 format(3E13.5,4F6.2) 
end subroutine  add2cmesh_2Dto3D

subroutine addslice2cmesh(x, y, r, g, b, o, nx, index, nt, fname)
!============================================================
! purpose:  write one slice of an object Geomview's CMESH format
!           this is particularly intended for stacking a pair-of-pants pic   
use prec, only: wp
implicit none

real(kind=wp),    intent(in), dimension (:)     :: x, y ! x, y coords
real(kind=wp),    intent(in), dimension (:)     :: r, g, b, o ! rgb & opacity
integer,          intent(in)                    :: nx, index, nt
character(len=*), intent(in)                    :: fname
integer                                         :: unit = 96
integer                                         :: ix

!! === executable statemets === !!

open (unit, file = fname, status = 'unknown', position = 'append')

if (index == 1) then
   write(unit, *) '{CMESH'
   write(unit, 998) nx, nt
endif

do ix = 1, nx
      write(unit, 999) x(ix), y(ix), dble(index)/dble(nt),     &
           &           r(ix), g(ix), b(ix), o(ix)
enddo

if (index == nt) then
   write(unit,*) '}'
endif

close(unit)

return
998 format(2I5)
999 format(3E13.5,4F6.2) 
end subroutine  addslice2cmesh
