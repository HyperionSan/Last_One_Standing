	module coordtrans
        implicit none

      contains 

! compute the value of the news in the intertial frame at scri...(north patch)

	subroutine null_coord(newsn, newss, un, us, bmidn, bmids)

      use null_grid
      use null_eth
      use null_interp
      use null_params
      use null_coortranvars
      use null_newsvarsouth
      use null_newsvarnorth

      implicit none

	   double complex, dimension(nn,nn) :: newss, newsn, un, us	
	   double precision, dimension(nn,nn) :: bmidn, bmids	

! position at scri
        integer px1, px2, py1, py2,  ll, kk

!temporary arrays
!*********************************************
!coord transf
	double precision, allocatable, dimension(:,:) :: y1arr, y2arr
	double precision, allocatable, dimension(:,:) :: xint, yint
	double precision, allocatable, dimension(:,:) :: u1n, u1s, u2n, u2s

!bondi time transf 
	double precision, allocatable, dimension(:,:) :: asfactnorth, asfactsouth

	double complex, allocatable, dimension(:,:) :: zeta, zetab

!*********************************************

	double precision :: xs, ys,xd,yd

	integer i,j, ite, cran


	allocate(u1n(nn,nn), u1s(nn,nn), u2n(nn,nn), u2s(nn,nn), &
     &             zeta(nn,nn), zetab(nn,nn),xint(nn,nn), yint(nn,nn), &
     &             asfactnorth(nn,nn), asfactsouth(nn,nn))


	zeta = qs + ii * ps
	zetab = qs - ii * ps
!	pp = (1. + qs**2 + ps**2)


!************** get the values of u1, u2 from u and \bar u


	u1n = 0.5 * pp * real(un)
	u1s = 0.5 * pp * real(us)

	u2n = 0.5 * pp * aimag(un)
	u2s = 0.5 * pp * aimag(us)


!***************************************************
! define the mask as follows: mask=1 --> north patch inside circle
! mask = 0. outside circle, dead points!, mask=1 points originally
! inside the circle that evolved to the other patch
!****************************************************


!*****************************************
! at u=u_o, set y^a = x^a
!***************************************
	if (it.eq.1) then
	   xno = qs
	   yno = ps

!**************************
! masking points. xno, yno temporary array that have the values of the
! points at both interveining patches 1 one array by virtue of keeping
! the mask....
!**************************
     
         where((xno**2+yno**2)>=(1.- dd)**2) 
	       xno = xno /(xno**2+yno**2)
	       yno = -yno /(xno**2+yno**2)
	       maskscri = 0 
	    elsewhere
	       maskscri = 1
	    end where

	end if
!*********************************
! since after evolving the points shift one inerpolates
! to get the values of u where the shifted point is, xint and yint
! are temporary variables that will be used to do the interpolation
!***********************************

	xint = xno
	yint = yno

! crank nicholson integration 
	do cran = 1, 40

!bracket the position on the array... and obtain de interpolated u (intu)
	  do i = 1, nn
	   do j = 1, nn

	    px2 = int( (xint(i,j)+1)/dd) + 3
	    py2 = int( (yint(i,j)+1)/dd) + 3

	    kk = px2
	    ll = py2

            xs =  xint(i,j)
            ys =  yint(i,j)
	
	   if(maskscri(i,j).eq.-1) then

	      call totalrinte(nn, dd,  u1s, intu1n(i,j), xs, ys)
	      call totalrinte(nn, dd,  u2s, intu2n(i,j), xs, ys)

	    else if(maskscri(i,j).eq.1) then

	      call totalrinte(nn, dd,  u1n, intu1n(i,j), xs, ys)
	      call totalrinte(nn, dd,  u2n, intu2n(i,j), xs, ys)

	    else

	   intu1n(i,j) = 0.
	   intu2n(i,j) = 0.

	end if

	 end do
	end do


! evolve
	  xnn = xno + dt * intu1n 
	  ynn = yno + dt * intu2n 

	  if (maxval(abs(0.5 * (xno + xnn) - xint)).lt.dt*dd**2.and.  &
     &		maxval(abs(0.5 * (yno + ynn) - yint)).lt.dt*dd**2) exit

	  xint = 0.5 * (xno + xnn)
	  yint = 0.5 * (yno + ynn)

	    where((xint**2+yint**2)>=(1.-dd)**2) 
	       maskscri = -maskscri 
	    end where


	  end do

	  
!remask
	  xint = 0.5 * (xno + xnn)
	  yint = 0.5 * (yno + ynn)

	    where((xint**2+yint**2)>=(1.-dd)**2) 
	       maskscri = -maskscri 
	    end where

! update
          xno = xnn
	  yno = ynn


!get the values for the asym conformal factor to be used in the ubondi....

	asfactnorth = 0.5 * (nwn+nwo) * exp(2.*bmidn)
	asfactsouth = 0.5 * (swn+swo) * exp(2.*bmids)


!interpolate and get the news on the inertial north patch

	  do i = 1, nn
	   do j = 1, nn

            xs =  xint(i,j)
            ys =  yint(i,j)
	
	   if(maskscri(i,j).eq.-1) then

	      call totalcinte(nn, dd,  newss, newsout(i,j), xs, ys)
	      call totalrinte(nn, dd,  asfactsouth, asfac(i,j), xs, ys)


	    else if(maskscri(i,j).eq.1) then

	      call totalcinte(nn, dd,  newsn, newsout(i,j), xs, ys)
	      call totalrinte(nn, dd,  asfactnorth, asfac(i,j), xs, ys)

	    else

	   newsout(i,j) = 0.
           asfac(i,j) = 1.

	end if

	 end do
	end do


!*****************************
! get bondi time
!***********************

! bondi time integration 
 
	    utimebon =  dt * asfac

!********* remask

	    where((xnn**2+ynn**2)>=(1.-dd)**2) 
	       maskscri = -maskscri 
	    end where


	deallocate(u1n, u1s, u2n, u2s, zeta, zetab, xint, yint,  &
     &              asfactnorth, asfactsouth)


	end subroutine null_coord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!_
!!!!!!!!!!!!!!!  south patch !!!!!!!!!!!!!!!!!!!!!!!!!!!!

! compute the value of the news in the intertial frame at scri...(south patch)

	subroutine null_coord_s(newsn, newss, un, us, bmidn, bmids)

      use null_grid
      use null_eth
      use null_interp
      use null_params
      use null_coortranvars
      use null_newsvarsouth
      use null_newsvarnorth

      implicit none

	   double complex, dimension(nn,nn) :: newss, newsn, un, us	
	   double precision, dimension(nn,nn) :: bmidn, bmids	

! position at scri
        integer px1, px2, py1, py2,  ll, kk

!temporary arrays
!*********************************************
!coord transf
	double precision, allocatable, dimension(:,:) :: y1arr, y2arr
	double precision, allocatable, dimension(:,:) :: xint, yint
	double precision, allocatable, dimension(:,:) :: u1n, u1s, u2n, u2s

!bondi time transf 
	double precision, allocatable, dimension(:,:) :: asfactnorth, asfactsouth

	double complex, allocatable, dimension(:,:) :: zeta, zetab

!*********************************************

	double precision :: xs, ys,xd,yd

	integer i,j, ite, cran


	allocate(u1n(nn,nn), u1s(nn,nn), u2n(nn,nn), u2s(nn,nn), &
     &             zeta(nn,nn), zetab(nn,nn),xint(nn,nn), yint(nn,nn), &
     &             asfactnorth(nn,nn), asfactsouth(nn,nn))


	zeta = qs + ii * ps
	zetab = qs - ii * ps
!	pp = (1. + qs**2 + ps**2)


!************** get the values of u1, u2 from u and \bar u


	u1n = 0.5 * pp * real(un)
	u1s = 0.5 * pp * real(us)

	u2n = 0.5 * pp * aimag(un)
	u2s = 0.5 * pp * aimag(us)


!***************************************************
! define the mask as follows: mask=1 --> south patch inside circle
! mask = 0. outside circle, dead points!, mask=-1 points originally
! inside the circle that evolved to the other patch
!****************************************************


!*****************************************
! at u=u_o, set y^a = x^a
!***************************************
	if (it.eq.1) then
	   xso = qs
	   yso = ps

!**************************
! masking points. xso, yso temporary array that have the values of the
! points at both interveining patches 1 one array by virtue of keeping
! the mask....
!**************************
     
         where((xso**2+yso**2)>=(1.+ 0.2 * dd)**2) 
	       xso = xso /(xso**2+yso**2)
	       yso = -yso /(xso**2+yso**2)
	       maskscri_s = 0 
	    elsewhere
	       maskscri_s = 1
	    end where

	end if
!*********************************
! since after evolving the points shift one inerpolates
! to get the values of u where the shifted point is, xint and yint
! are temporary variables that will be used to do the interpolation
!***********************************

	xint = xso
	yint = yso

! crank nicholson integration 
	do cran = 1, 40

!bracket the position on the array... and obtain de interpolated u (intu)
	  do i = 1, nn
	   do j = 1, nn

	    px2 = int( (xint(i,j)+1)/dd) + 3
	    py2 = int( (yint(i,j)+1)/dd) + 3

	    kk = px2
	    ll = py2

            xs =  xint(i,j)
            ys =  yint(i,j)
	
	   if(maskscri_s(i,j).eq.-1) then

	      call totalrinte(nn, dd,  u1n, intu1s(i,j), xs, ys)
	      call totalrinte(nn, dd,  u2n, intu2s(i,j), xs, ys)

	    else if(maskscri(i,j).eq.1) then

	      call totalrinte(nn, dd,  u1s, intu1s(i,j), xs, ys)
	      call totalrinte(nn, dd,  u2s, intu2s(i,j), xs, ys)

	    else

	   intu1s(i,j) = 0.
	   intu2s(i,j) = 0.

	end if

	 end do
	end do


! evolve
	  xsn = xso + dt * intu1s 
	  ysn = yso + dt * intu2s 

	  if (maxval(abs(0.5 * (xso + xsn) - xint)).lt.dt*dd**2.and.  &
     &		maxval(abs(0.5 * (yso + ysn) - yint)).lt.dt*dd**2) exit

	  xint = 0.5 * (xso + xsn)
	  yint = 0.5 * (yso + ysn)

	    where((xint**2+yint**2)>=(1.+0.2*dd)**2) 
	       maskscri_s = -maskscri_s 
	    end where


	  end do

	  
!remask
	  xint = 0.5 * (xso + xsn)
	  yint = 0.5 * (yso + ysn)

	    where((xint**2+yint**2)>=(1.+0.2*dd)**2) 
	       maskscri_s = -maskscri_s 
	    end where

! update
          xso = xsn
	  yso = ysn


!get the values for the asym conformal factor to be used in the ubondi....

	asfactnorth = 0.5 * (nwn+nwo) * exp(2.*bmidn)
	asfactsouth = 0.5 * (swn+swo) * exp(2.*bmids)


!interpolate and get the news on the inertial south patch

	  do i = 1, nn
	   do j = 1, nn

            xs =  xint(i,j)
            ys =  yint(i,j)
	
	   if(maskscri_s(i,j).eq.1) then

	      call totalcinte(nn, dd,  newss, newsout_s(i,j), xs, ys)
	      call totalrinte(nn, dd,  asfactsouth, asfac_s(i,j), xs, ys)


	    else if(maskscri_s(i,j).eq.-1) then

	      call totalcinte(nn, dd,  newsn, newsout_s(i,j), xs, ys)
	      call totalrinte(nn, dd,  asfactnorth, asfac_s(i,j), xs, ys)

	    else

	   newsout_s(i,j) = 0.
           asfac_s(i,j) = 1.

	end if

	 end do
	end do


!*****************************
! get bondi time
!***********************

! bondi time integration 
 
	    utimebon_s =  dt * asfac_s

!********* remask

	    where((xsn**2+ysn**2)>=(1.+0.2*dd)**2) 
	       maskscri_s = -maskscri_s 
	    end where


	deallocate(u1n, u1s, u2n, u2s, zeta, zetab, xint, yint,  &
     &              asfactnorth, asfactsouth)


	end subroutine null_coord_s


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!************************************************************
! interpolating real routine
!************************************************************
	 subroutine totalrinte(nn, dd,  in, out, xs, ys)

   	implicit none

   	integer :: nn
   	double precision, dimension (nn,nn) :: in
        double precision out
	double precision xs, ys, dd

	double precision :: x1, x2, x3, x4, y1, y2,y3,y4, xk, yk
	double precision :: xd, yd, factor
	integer kk, ll


	   kk = int( (xs+1)/dd) + 3
	   ll = int( (ys+1)/dd) + 3

	   xk = -1. + (kk-3) * dd
	   yk = -1. + (ll-3) * dd


         factor = 1. / (36. * (dd * dd) * (dd * dd) * (dd * dd))

            xd = xs - xk
            x1 = - xd * (xd - 2 * dd) * (xd - dd)
            x2 = 3 * (xd - 2 * dd) * (xd - dd) * (xd + dd)
            x3 = - 3 * xd * (xd -2 * dd) * (xd + dd)
            x4 = xd * (xd - dd) * (xd + dd)
      
            yd = ys - yk
            y1 = - yd * (yd - 2 * dd) * (yd - dd)
            y2 = 3 * (yd - 2 * dd) * (yd - dd) * (yd + dd)
            y3 = - 3 * yd * (yd - 2 * dd) * (yd + dd)
            y4 = yd * (yd - dd) * (yd + dd)
      

            out = factor *                                &
     &        (y1 * ( x1 * in(kk-1,ll-1) + x2 * in(kk,ll-1)         &   
     &              + x3 * in(kk+1,ll-1) + x4 * in(kk+2,ll-1))      &  
     &       + y2 * ( x1 * in(kk-1,ll)   + x2 * in(kk,ll)            &  
     &              + x3 * in(kk+1,ll)   + x4 * in(kk+2,ll))         & 
     &       + y3 * ( x1 * in(kk-1,ll+1) + x2 * in(kk,ll+1)           &
     &              + x3 * in(kk+1,ll+1) + x4 * in(kk+2,ll+1))        &
     &       + y4 * ( x1 * in(kk-1,ll+2) + x2 * in(kk,ll+2)           &
     &              + x3 * in(kk+1,ll+2) + x4 * in(kk+2,ll+2)))


	   return

            end subroutine totalrinte


!*****************************************************
!************************************************************
! interpolating complex routine
!************************************************************
	 subroutine totalcinte(nn, dd,  in, out, xs, ys)

   	implicit none

   	integer :: nn
   	double complex, dimension (nn,nn) :: in
        double complex out
	double precision xs, ys, dd

	double precision :: x1, x2, x3, x4, y1, y2,y3,y4, xk, yk
	double precision :: xd, yd, factor
	integer kk, ll


	   kk = int( (xs+1)/dd) + 3
	   ll = int( (ys+1)/dd) + 3

	   xk = -1. + (kk-3) * dd
	   yk = -1. + (ll-3) * dd


         factor = 1. / (36. * (dd * dd) * (dd * dd) * (dd * dd))

            xd = xs - xk
            x1 = - xd * (xd - 2 * dd) * (xd - dd)
            x2 = 3 * (xd - 2 * dd) * (xd - dd) * (xd + dd)
            x3 = - 3 * xd * (xd -2 * dd) * (xd + dd)
            x4 = xd * (xd - dd) * (xd + dd)
      
            yd = ys - yk
            y1 = - yd * (yd - 2 * dd) * (yd - dd)
            y2 = 3 * (yd - 2 * dd) * (yd - dd) * (yd + dd)
            y3 = - 3 * yd * (yd - 2 * dd) * (yd + dd)
            y4 = yd * (yd - dd) * (yd + dd)
      

            out = factor *                                &
     &        (y1 * ( x1 * in(kk-1,ll-1) + x2 * in(kk,ll-1)       &     
     &              + x3 * in(kk+1,ll-1) + x4 * in(kk+2,ll-1))    &    
     &       + y2 * ( x1 * in(kk-1,ll)   + x2 * in(kk,ll)          &    
     &              + x3 * in(kk+1,ll)   + x4 * in(kk+2,ll))       &   
     &       + y3 * ( x1 * in(kk-1,ll+1) + x2 * in(kk,ll+1)        &   
     &              + x3 * in(kk+1,ll+1) + x4 * in(kk+2,ll+1))     &   
     &       + y4 * ( x1 * in(kk-1,ll+2) + x2 * in(kk,ll+2)        &   
     &              + x3 * in(kk+1,ll+2) + x4 * in(kk+2,ll+2)))


	   return

            end subroutine totalcinte


	end module coordtrans
