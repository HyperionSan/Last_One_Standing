#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth1(u) ((-KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1o2dx)
#else
#  define PDstandardNth1(u) (PDstandardNth1_impl(u,p1o2dx,cdj,cdk))
static CCTK_REAL PDstandardNth1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1o2dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth2(u) ((-KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1o2dy)
#else
#  define PDstandardNth2(u) (PDstandardNth2_impl(u,p1o2dy,cdj,cdk))
static CCTK_REAL PDstandardNth2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1o2dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth3(u) ((-KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1))*p1o2dz)
#else
#  define PDstandardNth3(u) (PDstandardNth3_impl(u,p1o2dz,cdj,cdk))
static CCTK_REAL PDstandardNth3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNth2_impl(u, p1o2dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth11(u) ((-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1odx2)
#else
#  define PDstandardNth11(u) (PDstandardNth11_impl(u,p1odx2,cdj,cdk))
static CCTK_REAL PDstandardNth11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1odx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth22(u) ((-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1ody2)
#else
#  define PDstandardNth22(u) (PDstandardNth22_impl(u,p1ody2,cdj,cdk))
static CCTK_REAL PDstandardNth22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1ody2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth33(u) ((-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1))*p1odz2)
#else
#  define PDstandardNth33(u) (PDstandardNth33_impl(u,p1odz2,cdj,cdk))
static CCTK_REAL PDstandardNth33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNth22_impl(u, p1odz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth12(u) ((KRANC_GFOFFSET3D(u,-1,-1,0) - KRANC_GFOFFSET3D(u,-1,1,0) - KRANC_GFOFFSET3D(u,1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1o4dxdy)
#else
#  define PDstandardNth12(u) (PDstandardNth12_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL PDstandardNth12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,-1,-1,0) - KRANC_GFOFFSET3D(u,-1,1,0) - KRANC_GFOFFSET3D(u,1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1o4dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth13(u) ((KRANC_GFOFFSET3D(u,-1,0,-1) - KRANC_GFOFFSET3D(u,-1,0,1) - KRANC_GFOFFSET3D(u,1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1))*p1o4dxdz)
#else
#  define PDstandardNth13(u) (PDstandardNth13_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL PDstandardNth13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNth12_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth21(u) ((KRANC_GFOFFSET3D(u,-1,-1,0) - KRANC_GFOFFSET3D(u,-1,1,0) - KRANC_GFOFFSET3D(u,1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1o4dxdy)
#else
#  define PDstandardNth21(u) (PDstandardNth21_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL PDstandardNth21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNth12_impl(u, p1o4dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth23(u) ((KRANC_GFOFFSET3D(u,0,-1,-1) - KRANC_GFOFFSET3D(u,0,-1,1) - KRANC_GFOFFSET3D(u,0,1,-1) + KRANC_GFOFFSET3D(u,0,1,1))*p1o4dydz)
#else
#  define PDstandardNth23(u) (PDstandardNth23_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL PDstandardNth23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,0,-1,-1) - KRANC_GFOFFSET3D(u,0,-1,1) - KRANC_GFOFFSET3D(u,0,1,-1) + KRANC_GFOFFSET3D(u,0,1,1))*p1o4dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth31(u) ((KRANC_GFOFFSET3D(u,-1,0,-1) - KRANC_GFOFFSET3D(u,-1,0,1) - KRANC_GFOFFSET3D(u,1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1))*p1o4dxdz)
#else
#  define PDstandardNth31(u) (PDstandardNth31_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL PDstandardNth31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNth12_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth32(u) ((KRANC_GFOFFSET3D(u,0,-1,-1) - KRANC_GFOFFSET3D(u,0,-1,1) - KRANC_GFOFFSET3D(u,0,1,-1) + KRANC_GFOFFSET3D(u,0,1,1))*p1o4dydz)
#else
#  define PDstandardNth32(u) (PDstandardNth32_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL PDstandardNth32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNth23_impl(u, p1o4dydz, cdj, cdk);
}
#endif

