      module calc_facsouth
	implicit none
      contains

	subroutine calcconfactsouth(un, uo, betan, betao)

      use null_grid
      use null_eth
      use null_newsvarsouth
      use null_vars

! news <-- ::wn, wo , wom1, alphan, alphao, alphaom1
! jscri <-- :: jn, jo, jom1, j_ln, j_lo, j_lom1

! actually alpha stands for delta = alpha - gamma

      implicit none

	double complex, dimension(nn,nn) :: j, un, uo, u
	double precision, dimension(nn,nn) :: k,k_u, k_l, k_l_u,alpha, betan, betao, beta

!temporary arrays

	double complex, dimension(nn,nn) :: eth_j, ethb_j, eth_k, eth_u, eth_ub,&
     			 eth_k_l, eth_ethb_w,eth_beta,&
     			  eth_w, ethb_w, ethb_u,&
                         eth_j_l, ethb_j_l, eth2_w,eth2_beta,eth_ethb_beta,&
                         eth_alpha

	double complex, dimension(nn,nn) :: t1, t2
	double complex, dimension(nn,nn) :: j_u, zeta

	double complex, dimension(nn,nn):: j_l_u,  j_l

	double complex, dimension(nn,nn,nx):: temp

	double complex, dimension(nn,nn) :: v11, v12, v21, v22, d1f1, d1f2,& 
     &                                d2f1, d2f2,& 
     &                                f1, f2


	double precision, dimension(nn,nn) :: w,u1, u2, d1u1, d1u2, d2u1, d2u2, t3
	integer :: il,jl, cran


!************** eval de derivatives of j!


	zeta = qs + ii * ps

	j_los  =  - .5 *(3. * jos(:,:,nx) - 4. * jos(:,:,nx-1) &
     &                           + jos(:,:,nx-2) ) / dx

	j_lns  =  - .5 *(3. * jns(:,:,nx) - 4. * jns(:,:,nx-1) &
     &                           + jns(:,:,nx-2) ) / dx

!	call wt_csmooth(nn,j_los)
!	!call wt_csmooth(nn,j_lns)

	if(it.eq.1) then
	   swo = 1. 
	   alphaos = 0.
	end if

      w = swo
      u = 0.5 * (uo+un)
      alpha = alphaos

!************* crank nikolson
	do cran = 1, 60 


      call null_d1 (eth_w, dcmplx(w), 0, 1)
      call null_d1 (eth_ub, conjg(u), -1, 1)

	swn = swo - dt * real(eth_w * conjg(u)+ .5 * w * eth_ub )

	if(cran.gt.3.and.maxval(abs(.5*(swo+swn)-w)).lt..5*dt*dd**2) exit
	
	w = .5 * (swn+swo)

	end do
!*************end  crank nikolson

!        call wt_rsmooth(nn,swo)

      j_u = (jns(:,:,nx) - jos(:,:,nx)) / dt
      j_l = 0.5 * (j_los + j_lns)
      j_l_u = (j_lns - j_los)/ dt
      j = 0.5 * (jos(:,:,nx) + jns(:,:,nx) )
      k = sqrt(1.+j*conjg(j) )

      call null_d1 (eth_j, j, 2, 1)
      call null_d1 (ethb_j, j, 2, -1)
      call null_d1 (eth_u, u, 1, 1)
      call null_d1 (eth_ub, conjg(u), -1, 1)

	 ethb_u = conjg(eth_ub)

	f1 = 0.5 * pp * ( sqrt(.5 * (k+1.)) - j / sqrt(2. * (k+1.)) )
	f2 = ii * 0.5 * pp * ( sqrt(.5 * (k+1.)) + j / sqrt(2. * (k+1.)) )

	u1 = 0.5 * pp * real(u)
	u2 = 0.5 * pp * aimag(u)

	call srdx(u1, nn, dd, d1u1)
	call srdx(u2, nn, dd, d1u2)
	call scdx(f1, nn, dd, d1f1)
	call scdx(f2, nn, dd, d1f2)

	call srdy(u1, nn, dd, d2u1)
	call srdy(u2, nn, dd, d2u2)
	call scdy(f1, nn, dd, d2f1)
	call scdy(f2, nn, dd, d2f2)



	v11 = (1./pp) * sqrt(2*(k+1) ) * ( u1 * d1f1 + u2 * d2f1 &
     &                                   - f1 * d1u1 - f2 * d2u1 )

	v12 = - ii * (1./pp) * sqrt(2*(k+1) ) * ( u1 * d1f2 + u2 * d2f2 &
     &                                   - f1 * d1u2 - f2 * d2u2 )


	v21 = (1./pp) * conjg(j) * sqrt(2/(k+1) ) *         &
     &           ( u1 * d1f1 + u2 * d2f1  - f1 * d1u1 - f2 * d2u1 )


	v22 = ii *  (1./pp) * conjg(j) * sqrt(2/(k+1) ) *   &
     &           ( u1 * d1f2 + u2 * d2f2  - f1 * d1u2 - f2 * d2u2 )


	t3 = real( - 0.5 * ii * ( v11 + v12 + v21 + v22 + d1u1 + d2u2) )

!	call wt_rsmooth(nn,t3)

	t2 = 0.5 * ( conjg(j_u) * j ) / ( k + 1.)


!************* crank nikolson

	do cran = 1, 60

           call null_d1 (eth_alpha, dcmplx(alpha), 0, 1)

	t1 = 2.* conjg(u) * eth_alpha
	alphans = alphaos +  &
     &              .5 * dt * ( -real(t1) + aimag(t2) + t3 )

if((cran.gt.3).and.maxval(abs(.5*(alphaos+alphans)-alpha)).lt..5*dt*dd**2) exit

	alpha = .5 * (alphaos + alphans)

	end do

!************* end crank nikolson

!        call wt_rsmooth(nn,alphaos)

	j_lom1s  = j_los
	jom1s = jos(:,:,nx)
	jjs = jos(:,:,:)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!_

	return
	end subroutine calcconfactsouth



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine scdx(in, nn, dd, out)

	implicit none
      integer, intent (in) :: nn
      double complex, intent (inout) :: out(nn,nn)
      double complex, intent (in) :: in(nn,nn)
	double precision, intent(in) :: dd

      out(2:nn-1,:) = (in(3:nn,:) - in(1:nn-2,:)) / (2. * dd)
      out(1,:) = 0.5*(-3.*in(1,:)+ 4.*in(2,:) - in(3,:) ) /dd
      out(nn,:) =0.5*(3.*in(nn,:)-4.*in(nn-1,:)+in(nn-2,:) )/dd


	return
	end  subroutine scdx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine scdy(in, nn, dd, out)

	implicit none
      integer, intent (in) :: nn
      double complex, intent (inout) :: out(nn,nn)
      double complex, intent (in) :: in(nn,nn)
	double precision, intent(in) :: dd

      out(:,2:nn-1) = (in(:,3:nn) - in(:,1:nn-2)) / (2. * dd)
      out(:,1)=0.5*(-3.*in(:,1)+4.*in(:,2)-in(:,3))/dd
      out(:,nn)=0.5*(3.*in(:,nn)-4.*in(:,nn-1)+in(:,nn-2))/dd

	return
	end subroutine scdy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine srdx(in, nn, dd, out)

	implicit none
      integer, intent (in) :: nn
      double precision, intent (inout) :: out(nn,nn)
      double precision, intent (in) :: in(nn,nn)
	double precision, intent(in) :: dd

      out(2:nn-1,:) = (in(3:nn,:) - in(1:nn-2,:)) / (2. * dd)
      out(1,:) = 0.5*(-3.*in(1,:)+ 4.*in(2,:) - in(3,:) ) /dd
      out(nn,:) =0.5*(3.*in(nn,:)-4.*in(nn-1,:)+in(nn-2,:) )/dd


	return
	end subroutine srdx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine srdy(in, nn, dd, out)

	implicit none
      integer, intent (in) :: nn
      double precision, intent (inout) :: out(nn,nn)
      double precision, intent (in) :: in(nn,nn)
	double precision, intent(in) :: dd

      out(:,2:nn-1) = (in(:,3:nn) - in(:,1:nn-2)) / (2. * dd)
      out(:,1)=0.5*(-3.*in(:,1)+4.*in(:,2)-in(:,3))/dd
      out(:,nn)=0.5*(3.*in(:,nn)-4.*in(:,nn-1)+in(:,nn-2))/dd

	return
	end subroutine srdy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     end module calc_facsouth
