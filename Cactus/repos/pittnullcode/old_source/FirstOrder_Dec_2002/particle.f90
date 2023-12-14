module particle

  double precision, dimension(3), save :: p_zn, p_zo, &
                                          p_dvdo, p_dvdoo
  double precision, dimension(4), save :: p_vdn, p_vdo, &
                                          p_vuo, p_vuoo
  double precision, dimension (:,:,:), allocatable, save :: density, density_o
  double precision, save :: p_mass,p_mass_0,p_mass_r,p_mass_q,p_mass_p, pp_zn, U_2, &
                            bh_mass
  logical, save :: p_patchN, p_patchS, p_motion
  integer, save :: ir,iq,ip,irc,iqc,ipc,io_particle
  double complex, save :: p_van, p_vano

!Notes: the letters n o and oo at the end of an identifier refer to the new time-level,
!the previous time level or the one before that. 

contains

subroutine particle_initialize

   use null_grid
   implicit none
   double precision p1,p2,p3,v1,v2,v3

   namelist /particle_input/ p_mass, p1,p2,p3,v1,v2,v3,p_patchN,io_particle,U_2, &
                             bh_mass, p_motion

   open (unit = 10, file = "particle.in", status = "old" )
   read (unit = 10, nml = particle_input)
   close(unit = 10)

     print*,'particle evolution switched ON'
     allocate (density(nn,nn,nx), density_o(nn,nn,nx))  
     p_zn(1)=p1
     p_zn(2)=p2
     p_zn(3)=p3
     p_vdn(1)=v1
     p_vdn(2)=v2
     p_vdn(3)=v3
     p_patchS = .not.(p_patchN)
!!
!! N.B. Present code assumes that particle is always in patch N
!!
     call where_am_i
     p_zo  = p_zn
     p_vdo = p_vdn

     open (unit = 93, status = 'unknown', file = 'particle_p_zn.out')
     open (unit = 94, status = 'unknown', file = 'particle_p_vdo.out')
     open (unit = 95, status = 'unknown', file = 'particle_p_dvdo.out')
     open (unit = 96, status = 'unknown', file = 'particle_p_vuo.out')
     write(unit=93, fmt ='(4(1x,e23.13))')time,p_zn(1),p_zn(2),p_zn(3)
     write(unit=94, fmt ='(5(1x,e23.13))')time,p_vdo(1),p_vdo(2),p_vdo(3),p_vdo(4)
     write(unit=95, fmt ='(4(1x,e23.13))')time,p_dvdo(1),p_dvdo(2),p_dvdo(3)
     write(unit=96, fmt ='(5(1x,e23.13))')time,p_vuo(1),p_vuo(2),p_vuo(3),p_vuo(4)
     close (unit = 93)
     close (unit = 94)
     close (unit = 95)
     close (unit = 96)


end subroutine particle_initialize

subroutine where_am_i

  use null_grid
  implicit none
  integer :: i

!T
  integer j,k
  double precision ar,a0
!T

  double precision :: factor
  double precision dr,dq,dp,delta_r,delta_q,delta_p,det

   pp_zn=1+p_zn(2)**2+p_zn(3)**2

   do i=1,nx-1
     if((rb(i).le.p_zn(1)).and.(p_zn(1).le.rb(i+1))) then
        ir=i
     endif
   enddo
   if((p_zn(1)-rb(ir)).lt.(rb(ir+1)-p_zn(1))) then
       irc=ir+1
   else
       irc=ir
       ir=irc+1
   endif 

   do i=1,nn-1
     if((qs(i,i).le.p_zn(2)).and.(p_zn(2).le.qs(i+1,i+1))) then
        iq=i
     endif
     if((ps(i,i).le.p_zn(3)).and.(p_zn(3).le.ps(i+1,i+1))) then
        ip=i
     endif
   enddo
   if(p_zn(2)-qs(iq,iq).lt.qs(iq+1,iq+1)-p_zn(2)) then
       iqc=iq+1
   else
       iqc=iq
       iq=iqc+1
   endif 
   if(p_zn(3)-ps(ip,ip).lt.ps(ip+1,ip+1)-p_zn(3)) then
       ipc=ip+1
   else
       ipc=ip
       ip=ipc+1
   endif 

   dr = dabs(p_zn(1)-rb(ir))
   dq = dabs(p_zn(2)-qs(iq,iq))
   dp = dabs(p_zn(3)-ps(ip,ip))
   delta_r = dabs(rb(irc)-rb(ir))
   delta_q = dabs(qs(iqc,iqc)-qs(iq,iq))
   delta_p = dabs(ps(ipc,ipc)-ps(ip,ip))
   det = dp*delta_q*delta_r+dq*delta_p*delta_r+dr*delta_p*delta_q &
             +delta_p*delta_q*delta_r
   p_mass_r = (dr*delta_p*delta_q)/det
   p_mass_q = (dq*delta_p*delta_r)/det
   p_mass_p = (dp*delta_q*delta_r)/det 
   p_mass_0=1.0d0-(p_mass_r+p_mass_q+p_mass_p)
   density=0.0d0
   factor=-p_mass*pp_zn**2*rwt/ &
        (rwt+p_zn(1))**2/(4*p_vdn(1)*p_zn(1)**2*dd*dd*dx) 
   density(iq,ip,ir)=p_mass_0*factor
   density(iq,ip,irc)=p_mass_r*factor
   density(iqc,ip,ir)=p_mass_q*factor
   density(iq,ipc,ir)=p_mass_p*factor

!T For setting the density to a specified value at one point
!T   density=0.0d0
!T   density(:,:,ir)=0.2d-6
!T   print *, 'ir = ',ir, 'rb(ir) = ',rb(ir)
!T

!T For distributing the mass over a cube of grid-points
!T   a0=density(iq,ip,ir)/27
!T   do j=iq-1,iq+1
!T     do k=ip-1,ip+1
!T       do i=ir-1,ir+1
!T         density(j,k,i)=a0
!T       end do
!T     end do
!T   end do

!T For setting the density to a Gaussian
!T   a0=4.5d0
!T   do i=1,nx-1
!T     do j=1,nn
!T       do k=1,nn
!T         ar=dsqrt( (rb(i)-p_zn(1))**2+(p_zn(1)*datan(dsqrt(qs(j,j)**2+ps(k,k)**2)))**2)
!T         if (ar<a0)then
!T           density(j,k,i)=p_mass*dexp(-ar*4.5d0/a0)
!T           print*,i,j,k, density(j,k,i)
!T         end if
!T       end do
!T     end do
!T   end do
!T

!print*,'p_zn',p_zn
!print*,'ir,irc,iq,iqc,ip,ipc',ir,irc,iq,iqc,ip,ipc
!print*,'qs(iq,iq),qs(iqc,iqc)',qs(iq,iq),qs(iqc,iqc)
!print*,'ps(ip,ip),ps(ipc,ipc)',ps(ip,ip),ps(ipc,ipc)
!print*,'rb(ir),rb(irc)',rb(ir),rb(irc)
!print*,'p_mass,p_mass_0,p_mass_r,p_mass_q,p_mass_p',&
!        p_mass,p_mass_0,p_mass_r,p_mass_q,p_mass_p

end subroutine where_am_i

subroutine particle_evolve (iter)

   use null_grid
   implicit none
   integer,        intent (in) :: iter
   integer i
   double precision, dimension(4) ::  p_vuo_0, p_vuo_r, p_vuo_q, p_vuo_p
   double precision, dimension(3) ::  p_dvdo_0, p_dvdo_r, p_dvdo_q, p_dvdo_p
   double precision               ::  p_vdo4_0, p_vdo4_r, p_vdo4_q, p_vdo4_p

   p_vano = p_van
   density_o   = density
   call find_vu_dvd(iq,ip,ir,p_vuo_0,p_dvdo_0,p_vdo4_0)
   call find_vu_dvd(iq,ip,irc,p_vuo_r,p_dvdo_r,p_vdo4_r)
   call find_vu_dvd(iq,ipc,ir,p_vuo_p,p_dvdo_p,p_vdo4_p)
   call find_vu_dvd(iqc,ip,ir,p_vuo_q,p_dvdo_q,p_vdo4_q)
   p_vuo=(p_mass_0*p_vuo_0+p_mass_r*p_vuo_r+p_mass_q*p_vuo_q &
          +p_mass_p*p_vuo_p)
   p_vdo(4)=(p_mass_0*p_vdo4_0+p_mass_r*p_vdo4_r+p_mass_q*p_vdo4_q &
          +p_mass_p*p_vdo4_p)
   p_dvdo=(p_mass_0*p_dvdo_0+p_mass_r*p_dvdo_r+p_mass_q*p_dvdo_q &
          +p_mass_p*p_dvdo_p)

!T
!print*,'p_vuo',p_vuo
!print*,'p_vuo_0',p_vuo_0
!print*,'p_vuo_r',p_vuo_r
!print*,'p_vuo_q',p_vuo_q
!print*,'p_vuo_p',p_vuo_p
!print*,'p_dvdo',p_dvdo
!print*,'p_dvdo_0',p_dvdo_0
!print*,'p_dvdo_r',p_dvdo_r
!print*,'p_dvdo_q',p_dvdo_q
!print*,'p_dvdo_p',p_dvdo_p
!T

   if (iter.gt.1) then

!T For test with J set <>0, then put m to 0 to see how initial
!T radiation disperses.
!T     p_mass=0.0d0
!T
     do i=1,3
       if(p_motion) then ! p_motion switches particle motion on/off
         p_zn(i)  = p_zo(i)+dt*(1.5d0*p_vuo(i)-0.5d0*p_vuoo(i))
       end if
       p_vdn(i) = p_vdo(i)+dt*(1.5d0*p_dvdo(i)-0.5d0*p_dvdoo(i))
     enddo
   else
     do i=1,3
       if(p_motion) then
         p_zn(i)  = p_zo(i)+dt*p_vuo(i)
       end if
       p_vdn(i) = p_vdo(i)+dt*p_dvdo(i)
     enddo
   endif
   call where_am_i
   p_van = dcmplx(pp_zn*p_vdo(2)/2., pp_zn*p_vdo(3)/2.)
   if (iter==1) then
     p_vano = p_van
     density_o   = density
   endif

   if (mod(iter,io_particle)==0) then
     print*,'p_ev: p_zn',p_zn
     print*,'p_ev: p_vdn',p_vdn
     print*,'p_ev: p_dvdo',p_dvdo
     print*,'p_ev: p_vuo',p_vuo
     print*,'p_ev: p_vdo',p_vdo
     open (unit = 93, status = 'unknown', file = 'particle_p_zn.out', &
            position = 'append')
     open (unit = 94, status = 'unknown', file = 'particle_p_vdo.out', &
            position = 'append')
     open (unit = 95, status = 'unknown', file = 'particle_p_dvdo.out', &
            position = 'append')
     open (unit = 96, status = 'unknown', file = 'particle_p_vuo.out', &
            position = 'append')
     write(unit=93, fmt ='(4(1x,e23.13))')time,p_zn(1),p_zn(2),p_zn(3)
     write(unit=94, fmt ='(5(1x,e23.13))')time,p_vdo(1),p_vdo(2),p_vdo(3),p_vdo(4)
     write(unit=95, fmt ='(4(1x,e23.13))')time,p_dvdo(1),p_dvdo(2),p_dvdo(3)
     write(unit=96, fmt ='(5(1x,e23.13))')time,p_vuo(1),p_vuo(2),p_vuo(3),p_vuo(4)
     close (unit = 93)
     close (unit = 94)
     close (unit = 95)
     close (unit = 96)
   endif
 
   call particle_copy

end subroutine particle_evolve
  
subroutine particle_copy

   implicit none

   p_zo   = p_zn
   p_dvdoo= p_dvdo
   p_vdo  = p_vdn
   p_vuoo = p_vuo

end subroutine particle_copy

subroutine find_vu_dvd(jq,jp,jr,t_vuo,t_dvdo,t_vdo4)

   use null_grid
   use null_vars
   use null_eth
   implicit none
   integer,                        intent (in)  :: jq,jp,jr
   double precision, dimension(4), intent (out) ::  t_vuo
   double precision, dimension(3), intent (out) ::  t_dvdo
   double precision                             ::  t_vdo4

   double complex vdco, vdbo, Ub, U, J, Jb, p_vuoc, K, B, Va, r, &
                  t0, s1,s2,s3,s4, Jb_r,J_r,B_r,U_r,Ub_r,Va_r,q,p,&
                  ethK,ethJb,ethJ,ethVa,ethB,ethU,ethUb,dfac
   double complex, dimension (nn,nn) :: temp_eth, bnn_t, wnn_t

      r=dcmplx(rb(jr),0)
      q=dcmplx(qs(jq,jp),0)
      p=dcmplx(ps(jq,jp),0)
      Va=dcmplx(1.,0.)+dcmplx(wnn(jq,jp,jr),0.)*r
      B=dcmplx(bnn(jq,jp,jr),0.0d0)
      J=jnn(jq,jp,jr)
      U=unn(jq,jp,jr)
      Jb=conjg(J)
      K = dcmplx( dsqrt(1.0d0+ dble(J*Jb)),0.0d0)
      Ub=conjg(U)

      dfac = dcmplx(rwt/2/dx/(rb(jr)+rwt)**2,0.0d0)
      J_r = (jnn(jq,jp,jr+1) &
             -jnn(jq,jp,jr-1)) *dfac
      Jb_r = conjg(J_r)
      B_r =  dcmplx(bnn(jq,jp,jr+1) &
             -bnn(jq,jp,jr-1)) *dfac
      U_r = (unn(jq,jp,jr+1) &
             -unn(jq,jp,jr-1)) *dfac
      Ub_r = conjg(U_r)
      Va_r = r*(wnn(jq,jp,jr+1) &
             -wnn(jq,jp,jr-1)) *dfac &
                  +dcmplx(wnn(jq,jp,jr),0.0d0)

!T
!       Va=dcmplx(1.0-2.0/p_zo(1),0.0)
!       Va_r=dcmplx(2.0/p_zo(1)**2,0.0)
!      print*,r,wnn(jq,jp,jr),Va_r
!T

      call null_d1(temp_eth,jnn(:,:,jr),2,1)
      ethJ = temp_eth(jq,jp)
      call null_d1(temp_eth,jnn(:,:,jr),2,-1)
      ethJb = conjg(temp_eth(jq,jp))
      ethK = (ethJ*Jb+ethJb*J)/(2*K)
      wnn_t = dcmplx(wnn(:,:,jr),0.0d0)
      call null_d1(temp_eth,wnn_t,0,1)
      ethVa = r*temp_eth(jq,jp)
      bnn_t = dcmplx(bnn(:,:,jr),0.0d0)
      call null_d1(temp_eth,bnn_t,0,1)
      ethB = temp_eth(jq,jp)
      call null_d1(temp_eth,unn(:,:,jr),1,1)
      ethU = temp_eth(jq,jp)
      call null_d1(temp_eth,unn(:,:,jr),1,-1)
      ethUb = conjg(temp_eth(jq,jp))

      vdco = dcmplx(pp_zn*p_vdo(2)/2.,pp_zn*p_vdo(3)/2.)
      vdbo = conjg(vdco)
      p_vdo(4) = dble( (2*vdbo*vdco*K-J*vdbo**2-Jb*vdco**2 &
           -2*exp(-2*B)*p_vdo(1)*r**2*Ub*vdco-2*exp(-2*B) &
            *p_vdo(1)*r**2*U*vdbo+2*exp(-2*B)*Va &
            *p_vdo(1)**2*r**2+2*r**2)*exp(2*B)/p_vdo(1)/r**2/4 )
      t_vdo4=p_vdo(4)
!T
!     print*,'find ..1/0',K,J,Jb,B,U,Ub,Jb_r,J_r,B_r,U_r,Ub_r,&
!                      ethK,ethJb,ethJ,ethVa,ethB,ethU,ethUb   
!     print*,'find ..not0',Va,Va_r,vdco,vdbo,p_vdo(1),r,q,p,pp_zn
!      print*,' ..rb...',rb
!T
      t_vuo(4) = dble(-exp(-2*B)*p_vdo(1))

      t_vuo(1) = (exp(-2*B)*Va*p_vdo(1)-exp(-2*B)*U*vdbo/2 &
                -exp(-2*B)*Ub*vdco/2-exp(-2*B)*p_vdo(4))/t_vuo(4)

      p_vuoc= (-exp(-2*B)*p_vdo(1)*r**2*U+K*vdco-J*vdbo)/r**2
      t_vuo(2) = pp_zn* dble(p_vuoc)/2/t_vuo(4)
      t_vuo(3) = pp_zn* dimag(p_vuoc)/2/t_vuo(4)

      t0 = -(r*Jb_r*J*vdco*vdbo+r*J_r*Jb*vdco*vdbo &
      -4*vdco*vdbo-4*vdco*vdbo*Jb*J+2*J*vdbo**2*K &
      +2*Jb*vdco**2*K-4*exp(-2*B)*p_vdo(1)**2*B_r &
      *Va*r**3*K+8*B_r*exp(-2*B) &
      *p_vdo(1)*r**3*p_vdo(4)*K-2*exp(-2*B)*p_vdo(1) &
      *r**3*U_r*vdbo*K+4*B_r*exp(-2*B)*p_vdo(1)*r**3*U &
      *vdbo*K+4*B_r*exp(-2*B)*p_vdo(1)*r**3*Ub &
      *vdco*K-2*exp(-2*B)*p_vdo(1)*r**3*Ub_r*vdco*K &
      +2*exp(-2*B)*p_vdo(1)**2*Va_r*r**3*K &
      -r*J_r*vdbo**2*K-r*Jb_r*vdco**2*K)/K/r**3/4
      t_dvdo(1) = dble(t0)/t_vuo(4)

      s1 = 1.D0/4.D0
      s4 = 4*q*Jb*vdco**2+2*ethK*vdco*vdbo+ethJb*J**2* &
      vdbo**2+ethJ*Jb**2*vdco**2+ethJb*vdco**2+ethJ*vdbo**2 &
     -2*exp(-2*B)*p_vdo(1)**2*ethVa*r**2+ethJb*vdco**2*Jb &
      *J-4*q*vdco*vdbo*K-8*exp(-2*B)*p_vdo(1) &
      *ethB*r**2*p_vdo(4)+dcmplx(0.D0,4.D0)*p*Jb*vdco**2 &
      +ethJ*vdbo**2*Jb*J
      s3 = s4-2*ethK*J*K*vdbo**2-2*ethK*Jb*K &
      *vdco**2+dcmplx(0.D0,4.D0)*r**2*exp(-2*B)*p_vdo(1) &
      *p*Ub*vdco+2*r**2*exp(-2*B)*p_vdo(1)*ethU*vdbo+2 &
      *r**2*exp(-2*B)*p_vdo(1)*ethUb*vdco+4*exp(-2*B)* &
      p_vdo(1)**2*Va*ethB*r**2+4*r**2*exp(-2*B) &
      *p_vdo(1)*q*Ub*vdco+4*ethK*J*Jb*vdco*vdbo-2*ethJb &
      *K*J*vdco*vdbo-2*ethJ*K*Jb* &
      vdco*vdbo+dcmplx(0.D0,-4.D0)*p*vdco*vdbo*K-4*exp(-2*B) &
      *p_vdo(1)*ethB*r**2*Ub*vdco-4*exp(-2*B)*p_vdo(1) &
      *ethB*r**2*U*vdbo
      s4 = 1/r**2
      s2 = s3*s4
      t0 = s1*s2
      t_dvdo(2) = dble(t0)*2/pp_zn/t_vuo(4)
      t_dvdo(3) = dimag(t0)*2/pp_zn/t_vuo(4)

end subroutine find_vu_dvd


end module particle
