module point_dump
  use null_params
  use null_grid
  implicit none

contains
 
  subroutine  pointdump(time, j, name)
    implicit none
    double precision, intent(in) :: time
    double precision, dimension(nn, nn), intent(in) :: j
    character*(*), intent(in) :: name

    logical, save :: alreadycalled=.false.
    integer, save :: i
    double precision :: xi, xim1, xip1, xip2, intde
    double precision, save :: ci, cip1, cip2, cim1
    integer :: iunit1 = 991
    integer :: iunit2 = 992
    integer :: iunit3 = 993
    integer :: iunit4 = 994
    integer :: iunit5 = 995
    integer :: iunit6 = 996
    integer :: iunit7 = 997
    integer :: iunit8 = 998
    integer :: iunit9 = 999
    integer :: iunit10 =  1000
    integer :: iunit11 =  1001
    integer :: iunit12 =  1002

     if (.NOT.alreadycalled) then
      alreadycalled=.true.
      i  = (.5d0 * sqrt(2.0d0) + qsize) / dd + 3
      xi   = sqrt(qs(i  ,i  )**2 + ps(i  ,i  )**2)
      xim1 = sqrt(qs(i-1,i-1)**2 + ps(i-1,i-1)**2)
      xip1 = sqrt(qs(i+1,i+1)**2 + ps(i+1,i+1)**2)
      xip2 = sqrt(qs(i+2,i+2)**2 + ps(i+2,i+2)**2)

      ci = (1.0d0 - xim1)*(1.0d0 - xip1)*(1.0d0 - xip2) / ( &
           (xi    - xim1)*(xi    - xip1)*(xi    - xip2) )
     
      cim1 = (1.0d0 - xi)*(1.0d0 - xip1)*(1.0d0 - xip2) / ( &
           (xim1    - xi)*(xim1    - xip1)*(xim1    - xip2) )
     
      cip1 = (1.0d0 - xim1)*(1.0d0 - xi)*(1.0d0 - xip2) / ( &
           (xip1    - xim1)*(xip1    - xi)*(xip1    - xip2) )
     
      cip2 = (1.0d0 - xim1)*(1.0d0 - xi)*(1.0d0 - xip1) / ( &
           (xip2    - xim1)*(xip2    - xi)*(xip2    - xip1) )
    end if 

   ! WARNING! 
   ! notation : note this is a new naming convention
   ! H = 1/2 (not 1/2 qsize)
   ! E = Edge (qsize)
   ! 1 = 1 (not edge)
   ! DE = diagonal equator (interpolated value)
   ! G = qsize /2
   ! P = .2
    open(unit = iunit1,  file = name // '.00.dat', position = 'append')
    open(unit = iunit2,  file = name // '.H0.dat', position = 'append')
    open(unit = iunit3,  file = name // '.G0.dat', position = 'append')
    open(unit = iunit4,  file = name // '.10.dat', position = 'append')
    open(unit = iunit5,  file = name // '.E0.dat', position = 'append')
    open(unit = iunit6,  file = name // '.HH.dat', position = 'append')
    open(unit = iunit7,  file = name // '.DE.dat', position = 'append')
    open(unit = iunit8,  file = name // '.GG.dat', position = 'append')
    open(unit = iunit9,  file = name // '.11.dat', position = 'append')
    open(unit = iunit10, file = name // '.EE.dat', position = 'append')
    open(unit = iunit11, file = name // '.P0.dat', position = 'append')
    open(unit = iunit12, file = name // '.PP.dat', position = 'append')
 
    intde = ci * j(i,i)  + cim1 * j(i-1, i-1) + cip1*j(i+1, i+1) + cip2*j(i+2, i+2)
    write( unit = iunit1, fmt = * )  time, j(ZERO , ZERO)
    write( unit = iunit2, fmt = * )  time, j(HALF , ZERO)
    write( unit = iunit3, fmt = * )  time, j(HFED , ZERO)
    write( unit = iunit4, fmt = * )  time, j(ONE  , ZERO)
    write( unit = iunit5, fmt = * )  time, j(EDGE , ZERO)
    write( unit = iunit6, fmt = * )  time, j(HALF , HALF)
    write( unit = iunit7, fmt = * )  time, intde 
    write( unit = iunit8, fmt = * )  time, j(HFED , HFED)
    write( unit = iunit9, fmt = * )  time, j(ONE  , ONE)
    write( unit = iunit10, fmt = * ) time, j(EDGE , EDGE) 
    write( unit = iunit11, fmt = * ) time, j(PNTT , ZERO) 
    write( unit = iunit12, fmt = * ) time, j(PNTT , PNTT) 

 
    close(unit = iunit1)
    close(unit = iunit2)
    close(unit = iunit3)
    close(unit = iunit4)
    close(unit = iunit5)
    close(unit = iunit6)
    close(unit = iunit7)
    close(unit = iunit8)
    close(unit = iunit9)
    close(unit = iunit10)
    close(unit = iunit11)
    close(unit = iunit12)
  end subroutine  pointdump
   
 
end module point_dump
