        o1 = 2.00000000000000d0*phi(i,j,k)
        o2 = cos(o1)
        o3 = -eta0
        o4 = eta(i,j,k) + o3
        o5 = o4**2
        o6 = sigma**2
        o7 = 1/o6
        o8 = -(o5*o7)
        o9 = exp(o8)
        o10 = eta(i,j,k) + eta0
        o11 = o10**2
        o12 = -(o11*o7)
        o13 = exp(o12)
        o14 = o13 + o9
        o15 = sin(q(i,j,k))
        o16 = o15**n
        o17 = 2.00000000000000d0*amp*o14*o16
        o18 = exp(o17)
        o19 = -1.00000000000000d0 + o18
        o20 = eta(i,j,k)**2
        o21 = 2.00000000000000d0*o20
        o22 = eta0**2
        o23 = 2.00000000000000d0*o22
        o24 = eta(i,j,k)*o6
        o25 = o21 + o23 + o24
        o26 = -(o25*o7)
        o27 = exp(o26)
        o28 = 1/o15
        o29 = o20 + o22
        o30 = 2.00000000000000d0*o29*o7
        o31 = exp(o30)
        o32 = sin(phi(i,j,k))
        o33 = sin(o1)
        o34 = o11*o7
        o35 = exp(o34)
        o36 = cos(phi(i,j,k))
        o37 = o32**2
        o38 = cos(q(i,j,k))
        o39 = o38**2
        o40 = o5*o7
        o41 = exp(o40)
        o42 = o15**2
        o43 = o36**2
        o44 = -2.00000000000000d0*o19*o31*o6
        o45 = -2.00000000000000d0*eta(i,j,k)
        o46 = -2.00000000000000d0*eta0
        o47 = n*o6
        o48 = 4.0000000000000d0*eta(i,j,k)*eta0*o7
        o49 = exp(o48)
        o50 = -2.00000000000000d0*eta(i,j,k)*o49
        o51 = 2.00000000000000d0*eta0*o49
        o52 = n*o49*o6
        o53 = 2.00000000000000d0*q(i,j,k)
        o54 = cos(o53)
        o55 = 2.00000000000000d0*eta(i,j,k)
        o56 = 2.00000000000000d0*eta0
        o57 = 2.00000000000000d0*eta(i,j,k)*o49
        o58 = -2.00000000000000d0*eta0*o49
        o59 = o47 + o52 + o55 + o56 + o57 + o58
        o60 = o54*o59
        o61 = o45 + o46 + o47 + o50 + o51 + o52 + o60
        o62 = o17 + o40
        o63 = exp(o62)
        o64 = amp*o16*o61*o63
        o65 = o44 + o64
        o66 = -1.00000000000000d0 + o49
        o67 = -2.00000000000000d0*eta0*o66
        o68 = 1.00000000000000d0 + o49
        o69 = 2.00000000000000d0*eta(i,j,k)*o68
        o70 = n*o6*o68
        o71 = o67 + o69 + o70
        o72 = o56 + o6
        o73 = eta(i,j,k)*o72
        o74 = o20 + o22 + o73
        o75 = -(o7*o74)
        o76 = o17 + o75
        o77 = exp(o76)
        o78 = -eta(i,j,k)
        o79 = exp(o78)
        o80 = -2.00000000000000d0*o29*o7
        o81 = o17 + o80
        o82 = exp(o81)
        o83 = -1.00000000000000d0 + n
        o84 = o15**o83
        o85 = o35 + o41
        o86 = -(n*o39*o6*o85)
        o87 = eta(i,j,k)*o49
        o88 = -(eta0*o49)
        o89 = eta(i,j,k) + eta0 + o87 + o88
        o90 = 2.00000000000000d0*o41*o42*o89
        o91 = o86 + o90
        o92 = o17 + o26
        o93 = exp(o92)
        o94 = n*o39*o6*o85
        o95 = -2.00000000000000d0*o41*o42*o89
        o96 = o94 + o95
        gxx(i,j,k) = 5.0000000000000d-1*(1.00000000000000d0 + o18 + o19*
     &  o2)
c        dxgxx(i,j,k) = 5.0000000000000d-1*o27*o28*o7*(o36*(-2.0000000000
c     &  0000d0*o31*o37*o6 - amp*o16*o18*(2.00000000000000d0*eta(i,j,k)*o35*o42 
c     &  - 3.00000000000000d0*eta0*o35*o42 + 2.00000000000000d0*eta(i,j,k)*o2*o3
c     &  5*o42 + 2.00000000000000d0*eta(i,j,k)*o41*o42 + 2.00000000000000d0*eta0
c     &  *o41*o42 + 2.00000000000000d0*eta(i,j,k)*o2*o41*o42 + 2.00000000000000d
c     &  0*eta0*o2*o41*o42 - n*o35*o39*o6 - n*o2*o35*o39*o6 - n*o39*o41*o
c     &  6 - n*o2*o39*o41*o6)) + o18*(o31*o32*o33*o6 + amp*eta0*o15**(2.0
c     &  0000000000000d0 + n)*o35*cos(3.00000000000000d0*phi(i,j,k))))
c        dygxx(i,j,k) = 5.0000000000000d-1*o27*o28*o32*o43*o65*o7
c        dzgxx(i,j,k) = -(amp*o16*o38*o43*o7*o71*o77)
        gxy(i,j,k) = 5.0000000000000d-1*o19*o33
c        dxgxy(i,j,k) = 5.0000000000000d-1*o7*o79*(-(o19*o2*o28*o32*o6) -
c     &   amp*o33*o36*o82*o84*o91)
c        dygxy(i,j,k) = 5.0000000000000d-1*o7*o79*(1.00000000000000d0*o19
c     &  *o2*o28*o36*o6 - amp*o32*o33*o82*o84*o91)
c        dzgxy(i,j,k) = -5.0000000000000d-1*amp*o16*o33*o38*o7*o71*o77
        gxz(i,j,k) = 0
c        dxgxz(i,j,k) = 0
c        dygxz(i,j,k) = 0
c        dzgxz(i,j,k) = 0
        gyy(i,j,k) = 5.0000000000000d-1*(1.00000000000000d0 + o18 - o19*
     &  o2)
c        dxgyy(i,j,k) = 5.0000000000000d-1*o27*o28*o36*o37*o65*o7
c        dygyy(i,j,k) = 5.0000000000000d-1*o27*o28*o7*(-2.00000000000000d
c     &  0*o31*o32*o43*o6 - 2.00000000000000d0*amp*o16*o18*o32*o37*o91 + 
c     &  o33*o36*o6*exp(2.00000000000000d0*(amp*o14*o16 + o29*o7)))
c        dzgyy(i,j,k) = -(amp*o16*o37*o38*o7*o71*o77)
        gyz(i,j,k) = 0
c        dxgyz(i,j,k) = 0
c        dygyz(i,j,k) = 0
c        dzgyz(i,j,k) = 0
        gzz(i,j,k) = o18
c        dxgzz(i,j,k) = amp*o36*o7*o84*o93*o96
c        dygzz(i,j,k) = amp*o32*o7*o84*o93*o96
c        dzgzz(i,j,k) = -(amp*o16*o38*o7*o71*o77)