module StereoIntegration

integer, parameter, private :: rp = kind(1.0d0)
integer, parameter, private :: ip = kind(1)

contains

  ! StereoInt computes an integral of a real function over a 
  ! stereographic coordinate patch. The routine treats the
  ! integral as a weighted sum, computing the value of the
  ! integration weight at each point only at the initial
  ! function call. The value of the weight is stored in memory
  ! for use in subesquent calls.

  function StereoInt( n, field ) result( integral )
    integer(IP), intent(IN)              :: n
    real(RP), dimension(:,:), intent(IN) :: field
    real(RP)                             :: integral

    logical, save                     :: initial = .true.
    real(RP), dimension(:,:), &
         &    allocatable, save       :: weight
    integer(IP)                       :: i,j
    real(RP)                          :: dq, da, area, delq, delp
    real(RP), dimension(2,2)          :: q, p, rho2, wt
    logical,  dimension(2,2)          :: mask
    integer(IP), dimension(2)         :: k


    ! Compute the integration weight at the initial function call only
    if( initial ) then

       allocate( weight(n,n) )

       ! This routine assumes the patch extends beyond the hemisphere
       !   by two zones in each direction.
       dq = 2._RP/(n-5)
       da = dq*dq

       Weight = 0._RP

       ! Loop over each cell, q & p contain stereographic coordinate
       !   values for the corners of the cell
       do j = 3, n-3
          p(1,1) = -1.0_RP + (j-3)*dq
          p(2,1) = p(1,1)
          p(1,2) = p(1,1)+dq
          p(2,2) = p(1,2)

          do i = 1, n-1
             q(1,1) = -1.0_RP + (i-3)*dq
             q(2,1) = q(1,1) + dq
             q(1,2) = q(1,1)
             q(2,2) = q(2,1)

             rho2 = q*q + p*p

             ! Elements of mask are true if the corresponding corner
             !   of a cell lies inside the hemisphere, rho^2 <= 1
             mask = (rho2 <= 1.0_RP)

             ! Operate based on how many corners lie within the 
             !   hemisphere
             select case( count(mask) )
             case(0)
                area = 0._RP
             case(1)
                k = minloc(rho2)
                delq = sqrt( 1._RP - p(k(1),k(2))*p(k(1),k(2)))&
                     & - abs(q(k(1),k(2)))
                delp = sqrt( 1._RP - q(k(1),k(2))*q(k(1),k(2)))&
                     & - abs(p(k(1),k(2)))
                area = delq*delp/2._RP
             case(2)
                if(    (mask(1,1).and.mask(1,2)) .or.&
                     & (mask(2,1).and.mask(2,2))     ) then
                   delp = dq
                   if(mask(1,1)) then
                      delq = ( sqrt(1._RP-p(1,1)*p(1,1))-abs(q(1,1))&
                           & + sqrt(1._RP-p(1,2)*p(1,2))-abs(q(1,2)) )/2._RP
                   else
                      delq = ( sqrt(1._RP-p(2,1)*p(2,1))-abs(q(2,1))&
                           & + sqrt(1._RP-p(2,2)*p(2,2))-abs(q(2,2)) )/2._RP
                   end if
                else
                   delq = dq
                   if(mask(1,1)) then
                      delp = ( sqrt(1._RP-q(1,1)*q(1,1))-abs(p(1,1))&
                           & + sqrt(1._RP-q(2,1)*q(2,1))-abs(p(2,1)) )/2._RP
                   else             
                      delp = ( sqrt(1._RP-q(1,2)*q(1,2))-abs(p(1,2))&
                           & + sqrt(1._RP-q(2,2)*q(2,2))-abs(p(2,2)) )/2._RP
                   end if
                end if
                area = delq*delp
             case(3)
                k = maxloc(rho2)
                delq = abs(q(k(1),k(2)))-&
                     & sqrt( 1._RP - p(k(1),k(2))*p(k(1),k(2)))
                delp = abs(p(k(1),k(2)))-&
                     & sqrt( 1._RP - q(k(1),k(2))*q(k(1),k(2)))
                area = da - delq*delp/2._RP
             case(4)
                area = da
             end select

             wt = area/((1._RP + rho2)*(1._RP + rho2))

             Weight(i,j)     = Weight(i,j)     + wt(1,1)
             Weight(i+1,j)   = Weight(i+1,j)   + wt(2,1)
             Weight(i,j+1)   = Weight(i,j+1)   + wt(1,2)
             Weight(i+1,j+1) = Weight(i+1,j+1) + wt(2,2)
          end do
       end do

    initial = .false.

    end if

    ! Compute integral as weighted sum of field*weight
    integral = sum ( field * weight )

  end function StereoInt


  subroutine intster(field, q, p, dd, nn, integral)

    implicit none
    integer nn
    real*8 field(nn,nn), integral 
    real*8 q(nn), p(nn), dd

    integer i, l, k, m, orig
    real*8 area
    complex*16 sum,v(5)
    real*8 zz(5), cv(5), ee(5), aa(2)
    real*8 xx(5), yy(5), jac(5)

    orig = int( (nn-5)/2.0d0 ) + 3

!!$ Partial upper result
!!$ Look at all points except the ones at the upper and right
!!$ corners. Those are taken care of.

    integral = 0.0d0

    do l = 3, nn-3
       do k = 3, nn-3

!!$ Store the coordinates in array form

          if(k.ge.orig.and.l.ge.orig) then
             xx(1) = q(k)
             xx(2) = q(k+1)
             xx(3) = q(k+1)
             xx(4) = q(k)
             xx(5) = q(k)
             yy(1) = p(l)
             yy(2) = p(l)
             yy(3) = p(l+1)
             yy(4) = p(l+1)
             yy(5) = p(l)
          else if (k.lt.orig.and.l.ge.orig) then
             xx(1) = q(k+1)
             xx(2) = q(k+1)
             xx(3) = q(k)
             xx(4) = q(k)
             xx(5) = q(k+1)
             yy(1) = p(l)
             yy(2) = p(l+1)
             yy(3) = p(l+1)
             yy(4) = p(l)
             yy(5) = p(l)
          else if (k.lt.orig.and.l.lt.orig) then
             xx(1) = q(k+1)
             xx(2) = q(k)
             xx(3) = q(k)
             xx(4) = q(k+1)
             xx(5) = q(k+1)
             yy(1) = p(l+1)
             yy(2) = p(l+1)
             yy(3) = p(l)
             yy(4) = p(l)
             yy(5) = p(l+1)
          else if (k.ge.orig.and.l.lt.orig) then
             xx(1) = q(k)
             xx(2) = q(k)
             xx(3) = q(k+1)
             xx(4) = q(k+1)
             xx(5) = q(k)
             yy(1) = p(l+1)
             yy(2) = p(l)
             yy(3) = p(l)
             yy(4) = p(l+1)
             yy(5) = p(l+1)
          end if

!!$ Store the field values in a convenient array

          v(1)  = field(k,l)
          v(2)  = field(k+1,l)
          v(3)  = field(k+1,l+1)
          v(4)  = field(k,l+1)

!!$ Calculate the radial distance for each cell edge
!!$ Calculate the signed radial distance from the boundary
!!$ Create an array that signals if a vertex crossed or not

          do i = 1, 5
             zz(i) = xx(i)*xx(i) + yy(i)*yy(i)
             ee(i) = 1. - zz(i)
          end do
          do i = 1, 4
             cv(i) = ee(i) * ee(i+1)
          end do

!!$ Interior cells

          if(ee(3).ge.0) then
             area = dd * dd

!!$ Exterior cells

          else if(ee(1).le.0) then
             area = 0.

!!$ Mixed cells

          else 

!!$ Go around the cell points  1 ->  2 ->   3 ->    4 ->   1
!!$                           k,l   k+1,l  k+1,l+1 k,l+1  k,l
!!$ And figure out the crossing points and the associated lengths

             m = 0
             do i = 1, 4
                if(cv(i).lt.0) then 
                   m = m + 1
                   if((k.ge.orig.and.l.ge.orig).or.  &
                        (k.lt.orig.and.l.lt.orig)) then             
                      if(i.eq.1) then
                         aa(m) = sqrt(1.- yy(1)**2) - abs(xx(1))
                      else if(i.eq.2) then 
                         aa(m) = sqrt(1.- xx(2)**2) - abs(yy(2))
                      else if(i.eq.3) then
                         aa(m) = sqrt(1.- yy(4)**2) - abs(xx(4))
                      else if(i.eq.4) then
                         aa(m) = sqrt(1.- xx(1)**2) - abs(yy(1))
                      end if
                   else if((k.ge.orig.and.l.lt.orig).or.   &
                        (k.lt.orig.and.l.ge.orig)) then 
                      if(i.eq.1) then
                         aa(m) = sqrt(1.- xx(1)**2) - abs(yy(1))
                      else if(i.eq.2) then 
                         aa(m) = sqrt(1.- yy(2)**2) - abs(xx(2))
                      else if(i.eq.3) then
                         aa(m) = sqrt(1.- xx(4)**2) - abs(yy(4))
                      else if(i.eq.4) then
                         aa(m) = sqrt(1.- yy(1)**2) - abs(xx(1))
                      end if
                   end if
                else if(ee(i).eq.0) then
                   m = m + 1
                   aa(m) = dd
                end if
             end do

             if((cv(1).le.0).and.(cv(4).le.0)) then
                area = 0.5 * aa(1) * aa(2)
             else if((cv(1).le.0).and.(cv(3).le.0)) then
                area = 0.5 * (aa(1) + aa(2)) * dd
             else if((cv(2).le.0).and.(cv(3).le.0)) then
                area = dd**2 - 0.5 * (dd - aa(1))*(dd -aa(2))
             else if((cv(2).le.0).and.(cv(4).le.0)) then
                area = 0.5 * (aa(1) + aa(2)) * dd
             end if

          end if

          sum = (0.,0.)
          do i=1, 4
             jac(i) = 4.D0/((1.D0+zz(i))*(1.D0+zz(i)))
             sum = sum + jac(i) * v(i)
          end do
          integral = integral +  0.25 * sum * area
       end do
    end do

    return

  end subroutine intster

end module StereoIntegration
