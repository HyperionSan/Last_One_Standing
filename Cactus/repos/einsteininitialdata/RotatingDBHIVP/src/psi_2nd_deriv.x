        o1 = 5.0000000000000d-1*eta(i,j,k)
        o2 = exp(o1)
        o3 = psi3d(i,j,k)
        o4 = 1/o3
        o5 = cos(phi(i,j,k))
        o6 = o5**2
        o7 = cos(q(i,j,k))
        o8 = o7**2
        o9 = detapsi3d(i,j,k)
        o10 = -2.50000000000000d0*eta(i,j,k)
        o11 = exp(o10)
        o12 = dqqpsi3d(i,j,k)
        o13 = detaphipsi3d(i,j,k)
        o14 = sin(phi(i,j,k))
        o15 = dphipsi3d(i,j,k)
        o16 = o14**2
        o17 = sin(q(i,j,k))
        o18 = o17**2
        o19 = 1/o18
        o20 = dphiphipsi3d(i,j,k)
        o21 = detaqpsi3d(i,j,k)
        o22 = dqpsi3d(i,j,k)
        o23 = detaetapsi3d(i,j,k)
        o24 = tan(q(i,j,k))
        o25 = o24**2
        o26 = 1/o25
        o27 = dqphipsi3d(i,j,k)
        o28 = 1/o24
        psixx(i,j,k) = o2*o4*(1.00000000000000d0*o11*o16*o19*o20 + 1.000
     &  00000000000d0*o11*o16*o22*o28 - 5.0000000000000d-1*o11*o16*o3 - 
     &  2.00000000000000d0*o11*o13*o14*o5 + 2.00000000000000d0*o11*o14*o
     &  15*o5 + 1.00000000000000d0*o11*o14*o15*o19*o5 + 1.00000000000000
     &  d0*o11*o14*o15*o26*o5 - 2.00000000000000d0*o11*o14*o27*o28*o5 + 
     &  o11*o18*o23*o6 + 7.5000000000000d-1*o11*o18*o3*o6 + 2.0000000000
     &  0000d0*o11*o17*o21*o6*o7 - 3.00000000000000d0*o11*o17*o22*o6*o7 
     &  + o11*o12*o6*o8 - 5.0000000000000d-1*o11*o3*o6*o8 + o11*o16*o9 -
     &   2.00000000000000d0*o11*o18*o6*o9 + o11*o6*o8*o9)
        psixy(i,j,k) = o2*o4*(-(o11*o13*o16) + 1.50000000000000d0*o11*o1
     &  5*o16 + 1.00000000000000d0*o11*o15*o16*o26 - o11*o16*o27*o28 - o
     &  11*o14*o19*o20*o5 + o11*o14*o18*o23*o5 - o11*o14*o22*o28*o5 + 5.
     &  0000000000000d-1*o11*o14*o3*o5 + 7.5000000000000d-1*o11*o14*o18*
     &  o3*o5 + o11*o13*o6 - 5.0000000000000d-1*o11*o15*o6 - o11*o15*o19
     &  *o6 + 1.00000000000000d0*o11*o27*o28*o6 + 2.00000000000000d0*o11
     &  *o14*o17*o21*o5*o7 - 3.00000000000000d0*o11*o14*o17*o22*o5*o7 + 
     &  o11*o12*o14*o5*o8 - 5.0000000000000d-1*o11*o14*o3*o5*o8 - o11*o1
     &  4*o5*o9 - 2.00000000000000d0*o11*o14*o18*o5*o9 + o11*o14*o5*o8*o
     &  9)
        psixz(i,j,k) = o2*o4*(o11*o14*o27 - o11*o13*o14*o28 + 5.00000000
     &  00000d-1*o11*o14*o15*o28 - o11*o18*o21*o5 + 1.50000000000000d0*o
     &  11*o18*o22*o5 - o11*o12*o17*o5*o7 + o11*o17*o23*o5*o7 + 1.250000
     &  00000000d0*o11*o17*o3*o5*o7 + o11*o21*o5*o8 - 1.50000000000000d0
     &  *o11*o22*o5*o8 - 3.00000000000000d0*o11*o17*o5*o7*o9)
        psiyy(i,j,k) = o2*o4*(o11*o16*o18*o23 + 7.5000000000000d-1*o11*o
     &  16*o18*o3 + 2.00000000000000d0*o11*o13*o14*o5 - 2.00000000000000
     &  d0*o11*o14*o15*o5 - o11*o14*o15*o19*o5 - o11*o14*o15*o26*o5 + 2.
     &  00000000000000d0*o11*o14*o27*o28*o5 + 1.00000000000000d0*o11*o19
     &  *o20*o6 + 1.00000000000000d0*o11*o22*o28*o6 - 5.0000000000000d-1
     &  *o11*o3*o6 + 2.00000000000000d0*o11*o16*o17*o21*o7 - 3.000000000
     &  00000d0*o11*o16*o17*o22*o7 + o11*o12*o16*o8 - 5.0000000000000d-1
     &  *o11*o16*o3*o8 - 2.00000000000000d0*o11*o16*o18*o9 + o11*o6*o9 +
     &   o11*o16*o8*o9)
        psiyz(i,j,k) = o2*o4*(-(o11*o14*o18*o21) + 1.50000000000000d0*o1
     &  1*o14*o18*o22 - o11*o27*o5 + 1.00000000000000d0*o11*o13*o28*o5 -
     &   5.0000000000000d-1*o11*o15*o28*o5 - o11*o12*o14*o17*o7 + o11*o1
     &  4*o17*o23*o7 + 1.25000000000000d0*o11*o14*o17*o3*o7 + o11*o14*o2
     &  1*o8 - 1.50000000000000d0*o11*o14*o22*o8 - 3.00000000000000d0*o1
     &  1*o14*o17*o7*o9)
        psizz(i,j,k) = o2*o4*(o11*o12*o18 - 5.0000000000000d-1*o11*o18*o
     &  3 - 2.00000000000000d0*o11*o17*o21*o7 + 3.00000000000000d0*o11*o
     &  17*o22*o7 + o11*o23*o8 + 7.5000000000000d-1*o11*o3*o8 + o11*o18*
     &  o9 - 2.00000000000000d0*o11*o8*o9)
