/*
 * Copyright 2016 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include <math.h>
#include <cassert>
#include <iostream>
#include "EffectiveSource.h"
#include "WignerDMatrix.h"

extern "C"
{
#include <gsl/gsl_sf_ellint.h>
}

using std::complex;

/* signum function - returns -1, 0, 1 depending on the sign of val*/
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

RWZ_EffectiveSource::RWZ_EffectiveSource(const int l, const int m, const double M, const parity_t parity)
  : M(M), l(l), m(m), parity(parity), T(1.0), dT_dt(0.0), d2T_dt2(0.0),
    coeffs{}, dcoeffs_dt{}, d2coeffs_dt2{},
    abscoeffs{}, dabscoeffs_dt{}, d2abscoeffs_dt2{}
{
  if (parity == EvenParity)
    assert((l+m) % 2 == 0);
  else if (parity == OddParity)
    assert((l+m) % 2 == 1);

  /* Pre-compute Wigner-D used to rotate the coordinates from
   * a system where the particle is on the pole to the evolved system. */
  WignerDMatrix WignerD;

  /* The Wigner-D matrices */
  w0 = WignerD(l, 0, m);
  w1p = l>=1 ? WignerD(l, 1, m) : 0.0;
  w1m = l>=1 ? WignerD(l, -1, m) : 0.0;
  w2p = l>=2 ? WignerD(l, 2, m) : 0.0;
  w2m = l>=2 ? WignerD(l, -2, m) : 0.0;
  w3p = l>=3 ? WignerD(l, 3, m) : 0.0;
  w3m = l>=3 ? WignerD(l, -3, m) : 0.0;
}

void RWZ_EffectiveSource::set_rwz_time_window(const double T, const double dT_dt, const double d2T_dt2)
{
  this->T = T;
  this->dT_dt = dT_dt;
  this->d2T_dt2 = d2T_dt2;
}

void RWZ_EffectiveSource::set_rwz_particle(const double r, const double phi,
  const double ur, const double E0, const double L, const double ar, const double aphi,
  const double dardt, const double daphidt, const double d2ardt2, const double d2aphidt2)
{
  /* Coordinate position of the particle */
  r0 = r;
  phi0 = phi;

  /* Elliptic integrals appearing in the coefficients */
  const double gamma = sqrt(L*L/(L*L + r0*r0));
  const double ellE = gsl_sf_ellint_Ecomp(gamma, GSL_PREC_DOUBLE);
  const double ellK = gsl_sf_ellint_Kcomp(gamma, GSL_PREC_DOUBLE);

  if (parity == EvenParity) {
    /* Coefficients of powers of r-r0: m'=0 */
    coeffs[0][0] = (8*((ellK*(-12*(72275 - 13800*l - 13016*pow(l,2) + 1568*pow(l,3) + 784*pow(l,4))*pow(M,4) + 12*(90475 - 9080*l - 13536*pow(l,2) - 8416*pow(l,3) - 2968*pow(l,4) + 1488*pow(l,5) + 496*pow(l,6))*pow(M,3)*r0 + (-479325 + 650*l + 83322*pow(l,2) + 146704*pow(l,3) + 27904*pow(l,4) - 51312*pow(l,5) - 11728*pow(l,6) + 4608*pow(l,7) + 1152*pow(l,8))*pow(M,2)*pow(r0,2) + (88200 + 11985*l - 23791*pow(l,2) - 61152*pow(l,3) - 5344*pow(l,4) + 28128*pow(l,5) + 5792*pow(l,6) - 3072*pow(l,7) - 768*pow(l,8))*M*pow(r0,3) + (-5775 - 2165*l + 2683*pow(l,2) + 8080*pow(l,3) + 128*pow(l,4) - 4336*pow(l,5) - 848*pow(l,6) + 512*pow(l,7) + 128*pow(l,8))*pow(r0,4)))/((3*M - r0)*sqrt((1 + 2*l)*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) + (ellE*(48*(8225 - 2820*l - 2552*pow(l,2) + 536*pow(l,3) + 268*pow(l,4))*pow(M,5) - 12*(57225 - 4075*l - 15435*pow(l,2) - 21432*pow(l,3) - 7496*pow(l,4) + 3864*pow(l,5) + 1288*pow(l,6))*pow(M,4)*r0 - 4*(-154350 - 11885*l + 69691*pow(l,2) + 144188*pow(l,3) + 25836*pow(l,4) - 52284*pow(l,5) - 12052*pow(l,6) + 4608*pow(l,7) + 1152*pow(l,8))*pow(M,3)*pow(r0,2) + 3*(-107625 - 7705*l + 74719*pow(l,2) + 140128*pow(l,3) + 10056*pow(l,4) - 66992*pow(l,5) - 13968*pow(l,6) + 7168*pow(l,7) + 1792*pow(l,8))*pow(M,2)*pow(r0,3) - 2*(-43050 - 925*l + 37579*pow(l,2) + 63800*pow(l,3) - 96*pow(l,4) - 35528*pow(l,5) - 7064*pow(l,6) + 4096*pow(l,7) + 1024*pow(l,8))*M*pow(r0,4) + (-8925 + 85*l + 8637*pow(l,2) + 13920*pow(l,3) - 744*pow(l,4) - 8528*pow(l,5) - 1648*pow(l,6) + 1024*pow(l,7) + 256*pow(l,8))*pow(r0,5)))/(pow(-3*M + r0,2.5)*sqrt((1 + 2*l)*(-2*M + r0)))))/(1.*l*(1 + l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*sqrt(M_PI)*r0*(6*M + (-2 + l + pow(l,2))*r0));
    coeffs[0][1] = (4*(-((ellK*(192*(-46200 + 3230*l + 3259*pow(l,2) + 58*pow(l,3) + 29*pow(l,4))*pow(M,5) + 4*(8446725 - 2292035*l - 2346939*pow(l,2) - 87568*pow(l,3) + 11816*pow(l,4) + 66720*pow(l,5) + 22240*pow(l,6))*pow(M,4)*r0 + 2*(-16305975 + 6702590*l + 6617979*pow(l,2) - 351038*pow(l,3) - 615203*pow(l,4) - 486024*pow(l,5) - 92680*pow(l,6) + 59424*pow(l,7) + 14856*pow(l,8))*pow(M,3)*pow(r0,2) + (13025250 - 7057975*l - 6846016*pow(l,2) + 785382*pow(l,3) + 1240751*pow(l,4) + 864296*pow(l,5) + 39384*pow(l,6) - 199360*pow(l,7) - 32560*pow(l,8) + 11520*pow(l,9) + 2304*pow(l,10))*pow(M,2)*pow(r0,3) - 2*(1138725 - 760705*l - 737381*pow(l,2) + 112280*pow(l,3) + 206724*pow(l,4) + 143680*pow(l,5) - 11504*pow(l,6) - 46304*pow(l,7) - 5816*pow(l,8) + 3840*pow(l,9) + 768*pow(l,10))*M*pow(r0,4) + (139650 - 108415*l - 107344*pow(l,2) + 18070*pow(l,3) + 44935*pow(l,4) + 32360*pow(l,5) - 6312*pow(l,6) - 13120*pow(l,7) - 1360*pow(l,8) + 1280*pow(l,9) + 256*pow(l,10))*pow(r0,5)))/(3*M - r0)) + (ellE*(-384*(144375 - 39580*l - 36722*pow(l,2) + 5716*pow(l,3) + 2858*pow(l,4))*pow(M,6) - 8*(-16690275 + 5650450*l + 5185578*pow(l,2) - 924688*pow(l,3) - 449704*pow(l,4) + 15168*pow(l,5) + 5056*pow(l,6))*pow(M,5)*r0 + 8*(-16039275 + 6458980*l + 6055062*pow(l,2) - 907540*pow(l,3) - 693070*pow(l,4) - 259272*pow(l,5) - 39944*pow(l,6) + 39840*pow(l,7) + 9960*pow(l,8))*pow(M,4)*pow(r0,2) + 2*(31991925 - 14925345*l - 14804282*pow(l,2) + 1114174*pow(l,3) + 2597471*pow(l,4) + 2063536*pow(l,5) + 65312*pow(l,6) - 495584*pow(l,7) - 76376*pow(l,8) + 31680*pow(l,9) + 6336*pow(l,10))*pow(M,3)*pow(r0,3) + (-17508750 + 9350455*l + 10155292*pow(l,2) + 328794*pow(l,3) - 2789339*pow(l,4) - 2862992*pow(l,5) + 138864*pow(l,6) + 851776*pow(l,7) + 106384*pow(l,8) - 71040*pow(l,9) - 14208*pow(l,10))*pow(M,2)*pow(r0,4) + 2*(1238475 - 743790*l - 913525*pow(l,2) - 144374*pow(l,3) + 373529*pow(l,4) + 419816*pow(l,5) - 43928*pow(l,6) - 141856*pow(l,7) - 15784*pow(l,8) + 13120*pow(l,9) + 2624*pow(l,10))*M*pow(r0,5) + (-139650 + 89515*l + 130294*pow(l,2) + 39254*pow(l,3) - 76357*pow(l,4) - 88448*pow(l,5) + 13152*pow(l,6) + 32704*pow(l,7) + 3376*pow(l,8) - 3200*pow(l,9) - 640*pow(l,10))*pow(r0,6)))/pow(-3*M + r0,2)))/(1.*l*(1 + l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*sqrt(M_PI)*pow(r0,2)*pow(6*M + (-2 + l + pow(l,2))*r0,2)*sqrt((1 + 2*l)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))));
    coeffs[0][2] = ((ellK*(3456*(1125075 - 543560*l - 497776*pow(l,2) + 91568*pow(l,3) + 45784*pow(l,4))*pow(M,7) + 864*(-14417025 + 6694690*l + 5881550*pow(l,2) - 1600648*pow(l,3) - 736244*pow(l,4) + 76896*pow(l,5) + 25632*pow(l,6))*pow(M,6)*r0 - 24*(-615967275 + 336923305*l + 286269061*pow(l,2) - 101208816*pow(l,3) - 50104604*pow(l,4) + 1301512*pow(l,5) + 1603416*pow(l,6) + 1002496*pow(l,7) + 250624*pow(l,8))*pow(M,5)*pow(r0,2) - 36*(249272275 - 170582485*l - 135936869*pow(l,2) + 68895880*pow(l,3) + 32959376*pow(l,4) - 3152360*pow(l,5) - 3224296*pow(l,6) - 1656416*pow(l,7) - 155864*pow(l,8) + 172160*pow(l,9) + 34432*pow(l,10))*pow(M,4)*pow(r0,3) - 6*(-516518100 + 450609750*l + 320486413*pow(l,2) - 251262961*pow(l,3) - 100881506*pow(l,4) + 35688859*pow(l,5) + 20622497*pow(l,6) + 5031520*pow(l,7) - 1661708*pow(l,8) - 1727800*pow(l,9) - 126968*pow(l,10) + 119232*pow(l,11) + 19872*pow(l,12))*pow(M,3)*pow(r0,4) - 9*(69255200 - 77135980*l - 45979493*pow(l,2) + 57707903*pow(l,3) + 16942582*pow(l,4) - 15192301*pow(l,5) - 5934487*pow(l,6) + 429072*pow(l,7) + 1400316*pow(l,8) + 596624*pow(l,9) - 131440*pow(l,10) - 118144*pow(l,11) - 4160*pow(l,12) + 7168*pow(l,13) + 1024*pow(l,14))*pow(M,2)*pow(r0,5) + 6*(11592000 - 16194670*l - 7570411*pow(l,2) + 15341464*pow(l,3) + 2850105*pow(l,4) - 5824002*pow(l,5) - 1698406*pow(l,6) + 766096*pow(l,7) + 759856*pow(l,8) + 186344*pow(l,9) - 140632*pow(l,10) - 78400*pow(l,11) + 2464*pow(l,12) + 7168*pow(l,13) + 1024*pow(l,14))*M*pow(r0,6) - pow(-2 + l + pow(l,2),2)*(861525 - 582495*l - 451875*pow(l,2) + 285928*pow(l,3) + 191948*pow(l,4) + 24144*pow(l,5) - 46608*pow(l,6) - 40704*pow(l,7) - 2496*pow(l,8) + 5120*pow(l,9) + 1024*pow(l,10))*pow(r0,7)))/(6*pow(M,2) - 5*M*r0 + pow(r0,2)) + (ellE*(13824*(907725 - 407480*l - 373648*pow(l,2) + 67664*pow(l,3) + 33832*pow(l,4))*pow(M,8) + 432*(-83978475 + 41028220*l + 35689244*pow(l,2) - 10482496*pow(l,3) - 4752608*pow(l,4) + 586368*pow(l,5) + 195456*pow(l,6))*pow(M,7)*r0 - 96*(-472721025 + 261746530*l + 221038258*pow(l,2) - 80784252*pow(l,3) - 38664596*pow(l,4) + 2484076*pow(l,5) + 1513092*pow(l,6) + 587200*pow(l,7) + 146800*pow(l,8))*pow(M,6)*pow(r0,2) - 12*(2689505175 - 1742388520*l - 1425403993*pow(l,2) + 642840258*pow(l,3) + 337024955*pow(l,4) + 725612*pow(l,5) - 28544220*pow(l,6) - 22245376*pow(l,7) - 2525824*pow(l,8) + 2023680*pow(l,9) + 404736*pow(l,10))*pow(M,5)*pow(r0,3) - 12*(-1191831375 + 927701910*l + 717398398*pow(l,2) - 429596461*pow(l,3) - 227720171*pow(l,4) + 9998977*pow(l,5) + 42174851*pow(l,6) + 26119264*pow(l,7) - 2076908*pow(l,8) - 5177080*pow(l,9) - 474680*pow(l,10) + 305856*pow(l,11) + 50976*pow(l,12))*pow(M,4)*pow(r0,4) - 9*(450014600 - 427873920*l - 299624859*pow(l,2) + 256855302*pow(l,3) + 121126825*pow(l,4) - 29868476*pow(l,5) - 39797012*pow(l,6) - 15420448*pow(l,7) + 7691888*pow(l,8) + 5954192*pow(l,9) - 472432*pow(l,10) - 804736*pow(l,11) - 48704*pow(l,12) + 39424*pow(l,13) + 5632*pow(l,14))*pow(M,3)*pow(r0,5) + 3*(239689800 - 279697960*l - 167011459*pow(l,2) + 216924065*pow(l,3) + 79392886*pow(l,4) - 54263579*pow(l,5) - 42168881*pow(l,6) - 5411952*pow(l,7) + 15019884*pow(l,8) + 6922864*pow(l,9) - 2336912*pow(l,10) - 1685120*pow(l,11) + 6464*pow(l,12) + 132608*pow(l,13) + 18944*pow(l,14))*pow(M,2)*pow(r0,6) - (74006100 - 105453360*l - 49440897*pow(l,2) + 102837304*pow(l,3) + 24632389*pow(l,4) - 40688738*pow(l,5) - 21748486*pow(l,6) + 3947264*pow(l,7) + 12173096*pow(l,8) + 3707312*pow(l,9) - 2708560*pow(l,10) - 1499776*pow(l,11) + 68416*pow(l,12) + 146944*pow(l,13) + 20992*pow(l,14))*M*pow(r0,7) + pow(-2 + l + pow(l,2),2)*(861525 - 620295*l - 462675*pow(l,2) + 384376*pow(l,3) + 333476*pow(l,4) + 83760*pow(l,5) - 107376*pow(l,6) - 100608*pow(l,7) - 5952*pow(l,8) + 12800*pow(l,9) + 2560*pow(l,10))*pow(r0,8)))/(pow(-3*M + r0,2)*(-2*M + r0)))/(3.*l*(1 + l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*sqrt(M_PI)*pow(r0,3)*pow(6*M + (-2 + l + pow(l,2))*r0,3)*sqrt((1 + 2*l)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))));
    coeffs[0][3] = ((ellK*(331776*(258825 - 69860*l - 64846*pow(l,2) + 10028*pow(l,3) + 5014*pow(l,4))*pow(M,9) + 27648*(-11719050 + 6323755*l + 5166852*pow(l,2) - 2245306*pow(l,3) - 951403*pow(l,4) + 205500*pow(l,5) + 68500*pow(l,6))*pow(M,8)*r0 + 1728*(310279725 - 282606345*l - 199702121*pow(l,2) + 156549088*pow(l,3) + 55470464*pow(l,4) - 26400800*pow(l,5) - 7193440*pow(l,6) + 1377280*pow(l,7) + 344320*pow(l,8))*pow(M,7)*pow(r0,2) + 48*(-11001545100 + 14453642755*l + 8421247523*pow(l,2) - 11055283624*pow(l,3) - 3062962552*pow(l,4) + 2792476616*pow(l,5) + 656517080*pow(l,6) - 233276800*pow(l,7) - 56013280*pow(l,8) + 1537280*pow(l,9) + 307456*pow(l,10))*pow(M,6)*pow(r0,3) - 24*(-14210171850 + 24072570095*l + 10981108330*pow(l,2) - 23237610956*pow(l,3) - 4503086267*pow(l,4) + 7845433690*pow(l,5) + 1459325758*pow(l,6) - 988827008*pow(l,7) - 240180464*pow(l,8) + 11968480*pow(l,9) + 9677984*pow(l,10) + 3973248*pow(l,11) + 662208*pow(l,12))*pow(M,5)*pow(r0,4) - 12*(12287956800 - 24865323910*l - 8321955827*pow(l,2) + 28411692349*pow(l,3) + 3036966388*pow(l,4) - 11950017707*pow(l,5) - 1564975985*pow(l,6) + 2057285776*pow(l,7) + 472290028*pow(l,8) - 62079584*pow(l,9) - 44681680*pow(l,10) - 15317664*pow(l,11) - 651408*pow(l,12) + 877632*pow(l,13) + 125376*pow(l,14))*pow(M,4)*pow(r0,5) + 2*(20950524000 - 48143274000*l - 10699254638*pow(l,2) + 62239585437*pow(l,3) + 1224577048*pow(l,4) - 31045995029*pow(l,5) - 2331950777*pow(l,6) + 6756747268*pow(l,7) + 1427632555*pow(l,8) - 350309756*pow(l,9) - 231850004*pow(l,10) - 70196208*pow(l,11) + 3403720*pow(l,12) + 7035168*pow(l,13) + 1051104*pow(l,14) + 18432*pow(l,15) + 2304*pow(l,16))*pow(M,3)*pow(r0,6) + (-7420140000 + 18734826480*l + 2275593976*pow(l,2) - 26538516246*pow(l,3) + 1677647935*pow(l,4) + 15127194934*pow(l,5) + 372724852*pow(l,6) - 3933488660*pow(l,7) - 798878591*pow(l,8) + 269101588*pow(l,9) + 187998388*pow(l,10) + 57714336*pow(l,11) - 3714992*pow(l,12) - 6627264*pow(l,13) - 1284672*pow(l,14) - 135168*pow(l,15) - 16896*pow(l,16))*pow(M,2)*pow(r0,7) + 2*pow(-2 + l + pow(l,2),2)*(91372050 - 156753405*l - 96945094*pow(l,2) + 116243514*pow(l,3) + 48494227*pow(l,4) - 14824348*pow(l,5) - 10166724*pow(l,6) - 4009696*pow(l,7) - 381640*pow(l,8) + 467360*pow(l,9) + 146976*pow(l,10) + 29184*pow(l,11) + 4864*pow(l,12))*M*pow(r0,8) - pow(-2 + l + pow(l,2),3)*(-3735900 + 5161125*l + 4594281*pow(l,2) - 1334632*pow(l,3) - 1159756*pow(l,4) - 561616*pow(l,5) - 133744*pow(l,6) + 55040*pow(l,7) + 25280*pow(l,8) + 7680*pow(l,9) + 1536*pow(l,10))*pow(r0,9)))/((3*M - r0)*pow(-2*M + r0,2)*sqrt((1 + 2*l)*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) - (ellE*(663552*(311325 - 87760*l - 81356*pow(l,2) + 12808*pow(l,3) + 6404*pow(l,4))*pow(M,10) + 6912*(-107606625 + 51072910*l + 42719112*pow(l,2) - 16250020*pow(l,3) - 6981070*pow(l,4) + 1372728*pow(l,5) + 457576*pow(l,6))*pow(M,9)*r0 + 864*(1487034150 - 1153650505*l - 824561485*pow(l,2) + 617692104*pow(l,3) + 209422508*pow(l,4) - 114292624*pow(l,5) - 29738160*pow(l,6) + 7165184*pow(l,7) + 1791296*pow(l,8))*pow(M,8)*pow(r0,2) + 48*(-28960614900 + 33055637650*l + 19151145275*pow(l,2) - 25147911346*pow(l,3) - 6144064579*pow(l,4) + 7098407540*pow(l,5) + 1356103868*pow(l,6) - 827599744*pow(l,7) - 159222496*pow(l,8) + 31784960*pow(l,9) + 6356992*pow(l,10))*pow(M,7)*pow(r0,3) + 24*(42594013350 - 64733366165*l - 28746050395*pow(l,2) + 62646380909*pow(l,3) + 9151668737*pow(l,4) - 23446196821*pow(l,5) - 2722384015*pow(l,6) + 4021131536*pow(l,7) + 578850620*pow(l,8) - 278135920*pow(l,9) - 49474928*pow(l,10) + 3355776*pow(l,11) + 559296*pow(l,12))*pow(M,6)*pow(r0,4) - 12*(43609299300 - 81795124040*l - 25671446997*pow(l,2) + 93936762476*pow(l,3) + 4173392634*pow(l,4) - 43225408134*pow(l,5) - 1515323472*pow(l,6) + 9755940972*pow(l,7) + 850929813*pow(l,8) - 1009320924*pow(l,9) - 150681748*pow(l,10) + 30208272*pow(l,11) + 6943528*pow(l,12) + 880992*pow(l,13) + 125856*pow(l,14))*pow(M,5)*pow(r0,5) - 2*(-93079224000 + 204201868020*l + 38877348284*pow(l,2) - 265817345019*pow(l,3) + 16034682743*pow(l,4) + 143252624327*pow(l,5) - 7525605883*pow(l,6) - 39817305376*pow(l,7) - 838092328*pow(l,8) + 5550283016*pow(l,9) + 577456544*pow(l,10) - 297054384*pow(l,11) - 54484216*pow(l,12) - 1780128*pow(l,13) + 114336*pow(l,14) + 147456*pow(l,15) + 18432*pow(l,16))*pow(M,4)*pow(r0,6) + (-45014911200 + 111319308480*l + 8268975764*pow(l,2) - 159068146824*pow(l,3) + 26638609805*pow(l,4) + 96906648080*pow(l,5) - 14286304462*pow(l,6) - 31706412292*pow(l,7) + 1766429537*pow(l,8) + 5558149100*pow(l,9) + 233227196*pow(l,10) - 447085056*pow(l,11) - 48370240*pow(l,12) + 11464320*pow(l,13) + 1207680*pow(l,14) - 172032*pow(l,15) - 21504*pow(l,16))*pow(M,3)*pow(r0,7) + (7032085200 - 19120302600*l + 659576532*pow(l,2) + 29299166132*pow(l,3) - 8098923737*pow(l,4) - 19623712236*pow(l,5) + 4926813434*pow(l,6) + 7303982020*pow(l,7) - 1067841825*pow(l,8) - 1534279596*pow(l,9) + 61990308*pow(l,10) + 168912576*pow(l,11) + 2666272*pow(l,12) - 10243072*pow(l,13) - 377856*pow(l,14) + 434176*pow(l,15) + 54272*pow(l,16))*pow(M,2)*pow(r0,8) - pow(-2 + l + pow(l,2),2)*(159043500 - 308912580*l - 125519963*pow(l,2) + 319039982*pow(l,3) + 46595229*pow(l,4) - 118023980*pow(l,5) - 11799972*pow(l,6) + 20477824*pow(l,7) + 1374688*pow(l,8) - 2237440*pow(l,9) - 188416*pow(l,10) + 141312*pow(l,11) + 23552*pow(l,12))*M*pow(r0,9) + pow(-2 + l + pow(l,2),3)*(-3131100 + 5182725*l + 3089313*pow(l,2) - 3777688*pow(l,3) - 919444*pow(l,4) + 1016720*pow(l,5) + 103856*pow(l,6) - 183040*pow(l,7) - 22720*pow(l,8) + 15360*pow(l,9) + 3072*pow(l,10))*pow(r0,10)))/(sqrt(1 + 2*l)*pow(6*pow(M,2) - 5*M*r0 + pow(r0,2),2.5)))/(6.*l*(1 + l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*sqrt(M_PI)*pow(r0,4)*pow(6*M + (-2 + l + pow(l,2))*r0,4));

    abscoeffs[0][0] = (8*sqrt(M_PI)*(2*M - r0)*sqrt((1 + 2*l)*r0*(-3*M + r0)))/(1.*l*(1 + l)*(3*M - r0)*(6*M + (-2 + l + pow(l,2))*r0));
    abscoeffs[0][1] = (-4*sqrt(M_PI + 2*l*M_PI)*(24*pow(M,2) + 6*(-2 + l + pow(l,2))*M*r0 + l*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(r0,2)))/(1.*l*(1 + l)*sqrt(r0*(-3*M + r0))*pow(6*M + (-2 + l + pow(l,2))*r0,2));
    abscoeffs[0][2] = (2*sqrt(M_PI + 2*l*M_PI)*(12*(-16 + 5*l + 5*pow(l,2))*pow(M,3) + 12*(16 - 10*l - 9*pow(l,2) + 2*pow(l,3) + pow(l,4))*pow(M,2)*r0 + 3*(-4 + l + pow(l,2))*pow(-2 + l + pow(l,2),2)*M*pow(r0,2) - 2*l*(1 + l)*pow(-2 + l + pow(l,2),2)*pow(r0,3)))/(1.*l*(1 + l)*(2*M - r0)*sqrt(r0*(-3*M + r0))*pow(6*M + (-2 + l + pow(l,2))*r0,3));
    abscoeffs[0][3] = ((1 + 2*l)*sqrt(M_PI)*(-1382400*pow(M,7) - 432*(-9408 + 1277*l + 1277*pow(l,2))*pow(M,6)*r0 - 1152*(4328 - 1277*l - 1239*pow(l,2) + 76*pow(l,3) + 38*pow(l,4))*pow(M,5)*pow(r0,2) + 8*(413056 - 195864*l - 181650*pow(l,2) + 29293*pow(l,3) + 16809*pow(l,4) + 2595*pow(l,5) + 865*pow(l,6))*pow(M,4)*pow(r0,3) + 8*(-159296 + 106232*l + 92858*pow(l,2) - 28007*pow(l,3) - 17041*pow(l,4) - 3337*pow(l,5) - 599*pow(l,6) + 440*pow(l,7) + 110*pow(l,8))*pow(M,3)*pow(r0,4) + 3*(94720 - 81904*l - 66720*pow(l,2) + 32160*pow(l,3) + 20172*pow(l,4) + 3831*pow(l,5) - 501*pow(l,6) - 1482*pow(l,7) - 318*pow(l,8) + 35*pow(l,9) + 7*pow(l,10))*pow(M,2)*pow(r0,5) - 2*pow(-2 + l + pow(l,2),2)*(4176 - 248*l - 466*pow(l,2) - 427*pow(l,3) - 191*pow(l,4) + 27*pow(l,5) + 9*pow(l,6))*M*pow(r0,6) + 4*pow(-2 + l + pow(l,2),3)*(-48 - 12*l - 11*pow(l,2) + 2*pow(l,3) + pow(l,4))*pow(r0,7)))/(6.*l*(1 + l)*pow(2*M - r0,3)*pow(r0,2.5)*sqrt(-3*(1 + 2*l)*M + r0 + 2*l*r0)*pow(6*M + (-2 + l + pow(l,2))*r0,4));

    dcoeffs_dt[0][0] = 0.0;
    dcoeffs_dt[0][1] = 0.0;
    dcoeffs_dt[0][2] = 0.0;
    dcoeffs_dt[0][3] = 0.0;

    dabscoeffs_dt[0][0] = 0.0;
    dabscoeffs_dt[0][1] = 0.0;
    dabscoeffs_dt[0][2] = 0.0;
    dabscoeffs_dt[0][3] = 0.0;

    d2coeffs_dt2[0][0] = 0.0;
    d2coeffs_dt2[0][1] = 0.0;
    d2coeffs_dt2[0][2] = 0.0;
    d2coeffs_dt2[0][3] = 0.0;

    d2abscoeffs_dt2[0][0] = 0.0;
    d2abscoeffs_dt2[0][1] = 0.0;
    d2abscoeffs_dt2[0][2] = 0.0;
    d2abscoeffs_dt2[0][3] = 0.0;

    /* Coefficients of powers of r-r0: m'=1 */
    coeffs[1][0] = 0.0;
    coeffs[1][1] = 0.0;
    coeffs[1][2] = 0.0;
    coeffs[1][3] = 0.0;

    abscoeffs[1][0] = 0.0;
    abscoeffs[1][1] = 0.0;
    abscoeffs[1][2] = 0.0;
    abscoeffs[1][3] = 0.0;

    dcoeffs_dt[1][0] = 0.0;
    dcoeffs_dt[1][1] = 0.0;
    dcoeffs_dt[1][2] = 0.0;
    dcoeffs_dt[1][3] = 0.0;

    dabscoeffs_dt[1][0] = 0.0;
    dabscoeffs_dt[1][1] = 0.0;
    dabscoeffs_dt[1][2] = 0.0;
    dabscoeffs_dt[1][3] = 0.0;

    d2coeffs_dt2[1][0] = 0.0;
    d2coeffs_dt2[1][1] = 0.0;
    d2coeffs_dt2[1][2] = 0.0;
    d2coeffs_dt2[1][3] = 0.0;

    d2abscoeffs_dt2[1][0] = 0.0;
    d2abscoeffs_dt2[1][1] = 0.0;
    d2abscoeffs_dt2[1][2] = 0.0;
    d2abscoeffs_dt2[1][3] = 0.0;

    /* Coefficients of powers of r-r0: m'=2 */
    coeffs[2][0] = (-32*(-(ellK*(3*M - r0)*(144*(-2797830*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1593820*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1405373*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 370378*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 168899*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 19548*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 6516*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,5) + 6*(40821921*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 44621664*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 37787164*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 13331064*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 5824916*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 996912*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 312592*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 16896*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 4224*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,4)*r0 + 3*(68739300*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 5016459*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 329295*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 8447468*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1990952*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2445796*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 434092*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 310592*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 57488*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 13440*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2688*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,3)*pow(r0,2) - 3*(73865358*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 28516407*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 27000431*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3587264*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3106360*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1366288*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 114800*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 276608*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 49952*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 12800*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2560*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,2)*pow(r0,3) + 4*(16817220*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 9119880*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 8287287*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1771349*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1133470*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 248649*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3531*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 64368*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 11532*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3040*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 608*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))))*M*pow(r0,4) + 16*(-425250*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 278559*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 248664*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 62182*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 36611*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 5352*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 288*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1680*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 300*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 80*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 16*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(r0,5))) + ellE*(2*M - r0)*(18*(-26799507*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 15279760*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 13467404*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3561976*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1624148*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 188208*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 62736*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,5) - 3*(-85334634*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 99637635*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 84466291*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 29614356*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 12994460*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2152548*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 679660*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 32448*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 8112*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,4)*r0 + 6*(42924357*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 138399*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2065847*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3909446*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 754636*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1308562*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 221558*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 174752*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 32168*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 7680*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1536*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,3)*pow(r0,2) - 3*(82396818*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 32915049*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 30968633*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 4487588*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3647852*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1455756*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 111844*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 303168*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 54672*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 14080*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2816*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,2)*pow(r0,3) + 4*(17667720*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 9676998*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 8784615*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1895713*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1206692*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 259353*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2955*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 67728*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 12132*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3200*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 640*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))))*M*pow(r0,4) - 16*(425250*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 278559*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 248664*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 62182*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 36611*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 5352*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 288*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1680*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 300*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 80*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 16*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(r0,5))))/(3.*(-1 + l)*l*(1 + l)*(2 + l)*(-7 + 2*l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(9 + 2*l)*M*sqrt(M_PI)*r0*pow(-3*M + r0,2)*(6*M + (-2 + l + pow(l,2))*r0)*sqrt((1 + 2*l)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))));
    coeffs[2][1] = (16*(-(ellK*(3*M - r0)*(36*(-291820095*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 163955408*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 142958588*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 41241752*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 18741156*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2255664*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 751888*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,6) + 12*(1781255196*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1471793483*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1197792439*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 529790680*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 219775236*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 53000560*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 15760912*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1633664*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 408416*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,5)*r0 - 3*(5665212504*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 5934866068*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 4449888447*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2839561594*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1097154169*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 377726736*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 110127496*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 13672368*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3599892*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 121200*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 24240*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,4)*pow(r0,2) + 6*(1072557360*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1344698690*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 912907841*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 813985892*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 286005593*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 136956170*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 32463690*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 10271512*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1320322*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 764120*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 85240*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 36864*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 6144*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,3)*pow(r0,3) + 3*(-358481844*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 542962012*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 317207515*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 414755450*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 119840741*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 93240088*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 12495056*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 13790096*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 857756*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1594160*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 186480*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 72192*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 12032*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,2)*pow(r0,4) + 4*(7708932*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 22843944*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 6942711*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 27028945*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2537223*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 10585947*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 510761*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2927916*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 84753*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 399100*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 47436*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 17664*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2944*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))))*M*pow(r0,5) - 32*(-223020*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 216768*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 197793*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 74417*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 116208*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 61983*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 30311*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 36372*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 228*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 5470*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 654*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 240*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 40*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(r0,6))) + ellE*(576*(-44590329*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 25083740*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 21866978*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 6318284*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2871042*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 345720*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 115240*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,7) + 12*(5235742197*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 4091148512*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3356138380*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1422860344*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 594543876*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 137425744*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 41078896*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 4054016*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1013504*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,6)*r0 - 12*(5262684966*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 5115724863*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3929715007*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2273472116*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 892876636*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 285591076*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 83422988*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 10172672*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2643968*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 67200*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 13440*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,5)*pow(r0,2) + 3*(10970708112*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 12581418188*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 8964586711*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 6858720698*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2510039369*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1053519472*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 270861448*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 64002608*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 10164740*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3565360*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 387824*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 177408*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 29568*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,4)*pow(r0,3) + 12*(-759743838*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1010345693*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 654215204*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 664747551*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 217752285*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 126204281*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 24092105*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 13587628*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1194997*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1353540*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 156308*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 62400*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 10400*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,3)*pow(r0,4) + (1154303892*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1859099052*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1022646363*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1515219130*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 385890285*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 385078488*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 32668576*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 70325232*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3422580*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 8722480*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1027824*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 390912*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 65152*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,2)*pow(r0,5) - 12*(1082844*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 6169528*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 995617*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 8513535*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 71021*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3941869*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 372327*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1218452*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 29771*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 169500*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 20172*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 7488*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1248*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))))*M*pow(r0,6) + 32*(-223020*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 216768*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 197793*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 74417*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 116208*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 61983*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 30311*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 36372*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 228*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 5470*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 654*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 240*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 40*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(r0,7))))/(3.*(-1 + l)*l*(1 + l)*(2 + l)*(-7 + 2*l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(9 + 2*l)*M*sqrt(M_PI)*pow(r0,2)*pow(-3*M + r0,2)*pow(6*M + (-2 + l + pow(l,2))*r0,2)*sqrt((1 + 2*l)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))));
    coeffs[2][2] = (4*(-(ellK*(3*M - r0)*(8640*(-293354145*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 224303056*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 193748812*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 60002680*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 27236820*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3317424*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1105808*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,8) + 1080*(13478853969*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 9647115560*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 7742848076*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3660565160*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1464312876*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 428090160*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 124240912*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 15819264*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3954816*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,7)*r0 + 60*(-424545122358*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 344787236273*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 252695234745*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 173245115212*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 59883082716*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 30398523564*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 7356614084*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2302589760*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 479355600*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 64194560*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 12838912*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,6)*pow(r0,2) + 6*(3672770182200*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3470369959400*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2301880442729*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2151267436830*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 626530023955*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 497040637280*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 97685690232*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 54686639920*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 9244453780*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2849165520*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 467527824*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 55802880*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 9300480*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,5)*pow(r0,3) - 3*(3660533679600*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 4007672124200*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2384823630050*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2929073408749*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 704969253469*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 823415884025*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 132215675395*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 112942931032*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 17104036732*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 7249017820*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1267737900*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 111975744*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 29218624*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 4872000*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 696000*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,4)*pow(r0,4) + 6*(553763512260*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 694907145040*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 367758395135*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 580500122381*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 114052228994*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 188578863109*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 25873918947*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 29465257668*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 4542735972*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1921014462*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 415048406*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 6918784*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 7101904*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3810016*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 544288*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,3)*pow(r0,5) + (-606652328520*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 860911889580*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 402695146470*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 802442401363*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 129286866031*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 290187575255*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 37706174125*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 48080098616*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 8636753852*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2619126700*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 857570620*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 144189376*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 7273408*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 14190400*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1842880*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 73728*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 9216*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,2)*pow(r0,6) - 4*(-15489283320*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 24451295220*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 10126691970*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 24888994561*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3435998755*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 9709298588*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1395060934*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1581820250*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 378787124*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 42622999*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 37843507*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 13228768*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 58072*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1001392*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 112336*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 12288*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1536*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))))*M*pow(r0,7) + 16*(-171460800*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 294726600*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 110160468*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 320542422*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 42464561*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 130636991*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 24364757*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 18765954*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 7011869*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 669047*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 683737*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 370264*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 316*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 26544*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2512*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 512*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 64*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(r0,8))) + ellE*(138240*(-44390871*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 33954580*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 29332522*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 9076876*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 4120338*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 501720*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 167240*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,9) + 2160*(17511920811*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 12596617976*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 10165550456*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 4678408544*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1884678432*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 532017888*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 154984096*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 19161600*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 4790400*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,8)*r0 + 24*(-3191782852755*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2531709912155*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1886663616507*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1217123643100*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 429992967540*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 203570245740*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 50235184676*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 14639395200*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3078845760*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 387335360*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 77467072*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,7)*pow(r0,2) + 12*(6618431063790*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 5965049821325*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 4079985555640*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3489822209582*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1065514409177*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 756232998580*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 156033057100*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 77637477536*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 13628743736*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3733248160*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 626147360*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 65728512*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 10954752*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,6)*pow(r0,3) - 6*(8175303387720*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 8391453687960*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 5235182049979*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 5740583854217*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1494270034178*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1502763140023*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 260634228021*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 191729367416*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 30310595984*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 11480876244*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2015484676*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 170586432*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 42999808*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 6724032*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 960576*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,5)*pow(r0,4) + 3*(6397325232480*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 7437023283160*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 4203107921310*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 5789147801681*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1270166943353*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1750070363965*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 261013905175*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 256866750616*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 39826595200*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 16393428100*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3354993940*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 30757120*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 65083296*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 27285440*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3621440*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 110592*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 13824*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,4)*pow(r0,5) + 2*(-2413899370260*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3151707180060*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1603185586755*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2734982963282*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 503168613977*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 924405537994*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 124837438322*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 147567981244*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 24526589212*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 8882521202*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2347843466*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 226070624*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 32970080*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 31897376*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 4049888*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 202752*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 25344*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,3)*pow(r0,6) - (-760962966120*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1104397481340*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 503632565166*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1050178689287*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 163619449943*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 386765512147*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 51692107529*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 63781905304*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 12459239200*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3021239876*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1247386964*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 282654464*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 8922656*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 25043648*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2993984*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 233472*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 29184*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,2)*pow(r0,7) + 4*(-17203891320*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 27398561220*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 11228296650*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 28094418781*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3860644365*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 11015668498*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1638708504*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1769479790*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 448905814*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 35932529*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 44680877*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 16931408*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 54912*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1266832*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 137456*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 17408*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2176*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))))*M*pow(r0,8) - 16*(-171460800*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 294726600*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 110160468*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 320542422*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 42464561*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 130636991*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 24364757*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 18765954*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 7011869*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 669047*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 683737*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 370264*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 316*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 26544*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2512*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 512*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 64*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(r0,9))))/(15.*(-1 + l)*l*(1 + l)*(2 + l)*(-7 + 2*l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(9 + 2*l)*M*sqrt(M_PI)*(2*M - r0)*pow(r0,3)*pow(-3*M + r0,2)*pow(6*M + (-2 + l + pow(l,2))*r0,3)*sqrt((1 + 2*l)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))));
    coeffs[2][3] = (2*((ellK*(-207360*(-439350975*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 261840976*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 227847388*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 66766936*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 30332868*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3660720*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1220240*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,10) + 34560*(-4555914237*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3126184420*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2836020554*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 590760452*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 320369706*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 26928880*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3878800*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 4369280*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1092320*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,9)*r0 + 2160*(-67786030557*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 81810660748*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 36604932240*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 81059455160*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 17936181180*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 24932688576*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 4745832256*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2920673280*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 561298560*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 112579840*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 22515968*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,8)*pow(r0,2) + 96*(6443145835560*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 7745635098575*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 4072432034369*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 6531111489410*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1314560374660*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2101195307630*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 313152484082*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 306659660080*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 45622610020*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 19855652720*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3131913264*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 457754880*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 76292480*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,7)*pow(r0,3) - 24*(30373979929650*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 41050649236100*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 19419978095637*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 37698008164154*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 5668142518959*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 13826139895260*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1427556543756*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2456155599312*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 283613598972*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 208206412320*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 29560166032*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 6594784224*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1103367664*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1955520*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 279360*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,6)*pow(r0,4) - 24*(-19528779609360*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 29641659339600*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 11980627175490*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 30085758532015*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2758730150209*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 12542554495277*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 695689389391*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2634024420620*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 227206143782*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 271202153446*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 37291401238*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 9996833840*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2241012424*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 209630848*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 9835136*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 15912960*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1989120*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,5)*pow(r0,5) + (-186130073298960*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 314434118161920*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 104123157733760*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 349902728404840*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 10584694939163*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 163464540639568*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1273411935294*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 39519036315948*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2408395414095*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 4737961005204*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 675899084832*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 183372097064*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 57752442956*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 10246903072*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 201967904*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 686065152*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 98096064*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 4354560*sqrt(pow(l,35)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 483840*sqrt(pow(l,37)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,4)*pow(r0,6) - 12*(-3917175423840*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 7290831668920*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1883394956980*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 8790363694570*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 271726842073*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 4530751125943*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 168486294704*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1232103619016*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 51853150142*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 166063323697*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 26604545910*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 5928700158*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2802884245*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 704198208*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2923008*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 41311104*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 6564144*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 494208*sqrt(pow(l,35)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 54912*sqrt(pow(l,37)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,3)*pow(r0,7) + (-7376403419280*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 14967614544480*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2828943152648*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 19315753870440*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1706957870823*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 10818880887904*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 818237502038*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3240499802212*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 111921283677*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 473189835756*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 93298978608*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 11919568376*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 10792426500*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3494814624*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 78366304*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 181580288*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 30961984*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2916864*sqrt(pow(l,35)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 324096*sqrt(pow(l,37)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,2)*pow(r0,8) - 4*(-164376535680*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 359406047760*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 45687574944*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 490540900008*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 67268314448*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 294159851985*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 29744739644*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 94976935564*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 4435146909*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 14331575161*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3868186890*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 19500678*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 451105521*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 180502112*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 6327264*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 8478080*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1493872*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 153216*sqrt(pow(l,35)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 17024*sqrt(pow(l,37)*(-2 - l + 2*pow(l,2) + pow(l,3))))*M*pow(r0,9) + 80*(-318608640*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 740170944*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 58072656*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1055582368*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 178352888*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 667212132*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 68020787*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 226123990*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 21417669*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 32625956*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 13885839*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1712058*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1542587*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 743372*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 31764*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 32384*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 5680*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 576*sqrt(pow(l,35)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 64*sqrt(pow(l,37)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(r0,10)))/((3*M - r0)*pow(-2*M + r0,2)*sqrt((1 + 2*l)*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) + (ellE*(3317760*(-66398409*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 39553228*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 34415278*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 10091428*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 4584534*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 553416*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 184472*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,11) + 4320*(109751701563*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 73374182684*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 66069648040*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 14747172328*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 7697856564*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 330360320*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 12179840*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 83948800*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 20987200*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,10)*r0 - 432*(-453101035275*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 738335721480*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 214219387004*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 934264099920*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 191798498100*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 303825911000*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 57799778632*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 35634615040*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 6871002880*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1358433920*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 271686784*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,9)*pow(r0,2) - 48*(34454275273215*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 41405109728650*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 21491631346387*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 35441830089214*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 7217086589759*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 11341165401040*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1739519978176*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1620782889632*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 247137755672*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 101314730080*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 16205698272*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2213044224*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 368840704*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,8)*pow(r0,3) + 24*(100156741597710*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 131409680972200*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 63852747836591*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 118271262622300*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 19152149040343*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 42165033026894*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 4748100129650*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 7216751750920*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 874580730184*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 586847231612*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 84476236012*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 17944682880*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2993223648*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1127616*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 161088*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,7)*pow(r0,4) + 6*(-313607938545720*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 455498270629520*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 195606069144988*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 446186911520036*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 49911346541537*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 178319080533556*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 12574751064022*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 35626729326884*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3360430321133*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3487641102508*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 481852524512*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 124714976184*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 26188861908*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1964176704*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 97639488*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 151294464*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 18911808*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,6)*pow(r0,5) - (-928232469082800*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1487023029456960*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 544884702847280*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1584507473147440*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 92767008223529*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 703544838937720*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 19981216752378*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 160418583862164*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 11380831308021*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 18151785758412*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2549122248072*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 692244682472*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 194466395516*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 29667791104*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 700439552*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2027367936*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 285852096*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 11446272*sqrt(pow(l,35)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1271808*sqrt(pow(l,37)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,5)*pow(r0,6) + (-302601087396720*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 531070801787040*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 160203241041536*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 610950278704000*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 2637313732037*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 297918230492356*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3682631092422*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 76050563506692*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3947853066549*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 9657756683736*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1461950599128*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 359501487776*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 140357414336*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 30926865760*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 103570016*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1895912448*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 291510912*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 19243008*sqrt(pow(l,35)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2138112*sqrt(pow(l,37)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,4)*pow(r0,7) - (-65343127716720*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 124683051547200*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 29643236025512*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 153466841872464*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 7486977596205*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 81233382964516*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 4048662347930*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 22827881847964*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 899944871571*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3166836996912*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 550410025560*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 100846383440*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 60384333288*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 17096863392*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 225022048*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 946239488*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 155332864*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 13077504*sqrt(pow(l,35)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1453056*sqrt(pow(l,37)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,3)*pow(r0,8) + (-9015750712080*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 18551477703840*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 3284949650888*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 24206700204360*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 2377284634463*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 13751429378164*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1114777747798*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 4187266434932*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 156024995817*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 616094249536*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 131832323928*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 12137351216*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 15294880920*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 5292887264*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 139937824*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 266487808*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 45818624*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 4414464*sqrt(pow(l,35)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 490496*sqrt(pow(l,37)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,2)*pow(r0,9) - 4*(-180306967680*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 396414594960*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 48591207744*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 543320018408*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 76185958848*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 327520458585*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 33145778994*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 106283135064*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 5506030359*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 15962872961*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 4562478840*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 66102222*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 528234871*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 217670712*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 7915464*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 10097280*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1777872*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 182016*sqrt(pow(l,35)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 20224*sqrt(pow(l,37)*(-2 - l + 2*pow(l,2) + pow(l,3))))*M*pow(r0,10) + 80*(-318608640*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 740170944*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 58072656*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1055582368*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 178352888*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 667212132*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 68020787*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 226123990*sqrt(pow(l,15)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 21417669*sqrt(pow(l,17)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 32625956*sqrt(pow(l,19)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 13885839*sqrt(pow(l,21)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 1712058*sqrt(pow(l,23)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 1542587*sqrt(pow(l,25)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 743372*sqrt(pow(l,27)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 31764*sqrt(pow(l,29)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 32384*sqrt(pow(l,31)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 5680*sqrt(pow(l,33)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 576*sqrt(pow(l,35)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 64*sqrt(pow(l,37)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(r0,11)))/(sqrt(1 + 2*l)*pow(6*pow(M,2) - 5*M*r0 + pow(r0,2),2.5))))/(15.*(-1 + l)*l*(1 + l)*(2 + l)*(-7 + 2*l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(9 + 2*l)*M*sqrt(M_PI)*pow(r0,4)*pow(6*M + (-2 + l + pow(l,2))*r0,4));

    abscoeffs[2][0] = 0.0;
    abscoeffs[2][1] = (2*M*sqrt(M_PI)*sqrt(-((l*(-2 - 5*l + 5*pow(l,3) + 2*pow(l,4))*r0)/(3*M - r0))))/(1.*(-1 + l)*l*(1 + l)*(2 + l)*(2*M - r0)*r0);
    abscoeffs[2][2] = (M*sqrt((l*(-2 - l + 2*pow(l,2) + pow(l,3))*(M_PI + 2*l*M_PI))/(-3*M + r0))*(12*pow(M,2) + (-2 + l + pow(l,2))*pow(r0,2)))/(1.*(-1 + l)*l*(1 + l)*(2 + l)*pow(r0,1.5)*pow(-2*M + r0,2)*(6*M + (-2 + l + pow(l,2))*r0));
    abscoeffs[2][3] = (M*sqrt(-((M_PI + 2*l*M_PI)/(3*M - r0)))*(288*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3)))*pow(M,4) + 36*(-40*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 13*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 13*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,3)*r0 + 8*(152*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 98*sqrt(pow(l,3)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 87*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 22*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 11*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(M,2)*pow(r0,2) + (-360*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 340*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) + 246*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 181*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 73*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 21*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 7*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))))*M*pow(r0,3) + (32*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 44*l*sqrt(l*(-2 - l + 2*pow(l,2) + pow(l,3))) - 24*sqrt(pow(l,5)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 37*sqrt(pow(l,7)*(-2 - l + 2*pow(l,2) + pow(l,3))) + 11*sqrt(pow(l,9)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 9*sqrt(pow(l,11)*(-2 - l + 2*pow(l,2) + pow(l,3))) - 3*sqrt(pow(l,13)*(-2 - l + 2*pow(l,2) + pow(l,3))))*pow(r0,4)))/(6.*(-1 + l)*l*(1 + l)*(2 + l)*pow(r0,2.5)*pow(-2*M + r0,3)*pow(6*M + (-2 + l + pow(l,2))*r0,2));

    dcoeffs_dt[2][0] = 0.0;
    dcoeffs_dt[2][1] = 0.0;
    dcoeffs_dt[2][2] = 0.0;
    dcoeffs_dt[2][3] = 0.0;

    dabscoeffs_dt[2][0] = 0.0;
    dabscoeffs_dt[2][1] = 0.0;
    dabscoeffs_dt[2][2] = 0.0;
    dabscoeffs_dt[2][3] = 0.0;

    d2coeffs_dt2[2][0] = 0.0;
    d2coeffs_dt2[2][1] = 0.0;
    d2coeffs_dt2[2][2] = 0.0;
    d2coeffs_dt2[2][3] = 0.0;

    d2abscoeffs_dt2[2][0] = 0.0;
    d2abscoeffs_dt2[2][1] = 0.0;
    d2abscoeffs_dt2[2][2] = 0.0;
    d2abscoeffs_dt2[2][3] = 0.0;

    /* Coefficients of powers of r-r0: m'=3 */
    coeffs[3][0] = 0.0;
    coeffs[3][1] = 0.0;
    coeffs[3][2] = 0.0;
    coeffs[3][3] = 0.0;

    abscoeffs[3][0] = 0.0;
    abscoeffs[3][1] = 0.0;
    abscoeffs[3][2] = 0.0;
    abscoeffs[3][3] = 0.;

    dcoeffs_dt[3][0] = 0.0;
    dcoeffs_dt[3][1] = 0.0;
    dcoeffs_dt[3][2] = 0.0;
    dcoeffs_dt[3][3] = 0.0;

    dabscoeffs_dt[3][0] = 0.0;
    dabscoeffs_dt[3][1] = 0.0;
    dabscoeffs_dt[3][2] = 0.0;
    dabscoeffs_dt[3][3] = 0.0;

    d2coeffs_dt2[3][0] = 0.0;
    d2coeffs_dt2[3][1] = 0.0;
    d2coeffs_dt2[3][2] = 0.0;
    d2coeffs_dt2[3][3] = 0.0;

    d2abscoeffs_dt2[3][0] = 0.0;
    d2abscoeffs_dt2[3][1] = 0.0;
    d2abscoeffs_dt2[3][2] = 0.0;
    d2abscoeffs_dt2[3][3] = 0.0;
  } else if (parity == OddParity) {
    /* Coefficients of powers of r-r0: m'=0 */
    coeffs[0][0] = 0.0;
    coeffs[0][1] = 0.0;
    coeffs[0][2] = 0.0;
    coeffs[0][3] = 0.0;

    abscoeffs[0][0] = 0.0;
    abscoeffs[0][1] = 0.0;
    abscoeffs[0][2] = 0.0;
    abscoeffs[0][3] = 0.0;

    dcoeffs_dt[0][0] = 0.0;
    dcoeffs_dt[0][1] = 0.0;
    dcoeffs_dt[0][2] = 0.0;
    dcoeffs_dt[0][3] = 0.0;

    dabscoeffs_dt[0][0] = 0.0;
    dabscoeffs_dt[0][1] = 0.0;
    dabscoeffs_dt[0][2] = 0.0;
    dabscoeffs_dt[0][3] = 0.0;

    d2coeffs_dt2[0][0] = 0.0;
    d2coeffs_dt2[0][1] = 0.0;
    d2coeffs_dt2[0][2] = 0.0;
    d2coeffs_dt2[0][3] = 0.0;

    d2abscoeffs_dt2[0][0] = 0.0;
    d2abscoeffs_dt2[0][1] = 0.0;
    d2abscoeffs_dt2[0][2] = 0.0;
    d2abscoeffs_dt2[0][3] = 0.0;

    /* Coefficients of powers of r-r0: m'=1 */
    coeffs[1][0] = (16*(-(ellE*((2*pow(-2*M + r0,2)*sqrt((l*(1 + l))/(1.*(1 + 2*l)*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(17*pow(M,2) - 13*M*r0 + 2*pow(r0,2)))/(1.*(-1 + 2*l)*(3 + 2*l)*pow(r0,1.5)*(-3*M + r0)) - (210*sqrt((l*(1 + l))/(M + 2*l*M))*pow((-2*M + r0)/(r0*(-3*M + r0)),1.5)*(46*pow(M,2) - 28*M*r0 + 6*pow(r0,2)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)) + (6*pow(-2*M + r0,2)*sqrt((l*(1 + l))/(1.*(1 + 2*l)*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(29*pow(M,2) - 45*M*r0 + 14*pow(r0,2)))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*pow(r0,1.5)*(-3*M + r0)) + 4*pow(-2*M + r0,2)*sqrt((l*(1 + l))/(1.*(1 + 2*l)*M*(6*pow(M,2)*r0 - 5*M*pow(r0,2) + pow(r0,3)))) + ((2*M - r0)*(-6*(-7595 - 270*l - 206*pow(l,2) + 128*pow(l,3) + 64*pow(l,4))*pow(M,4) + (-132825 + 114734*l + 77614*pow(l,2) - 71360*pow(l,3) - 28480*pow(l,4) + 8640*pow(l,5) + 2880*pow(l,6))*pow(M,3)*r0 - 2*(-48195 + 50378*l + 34186*pow(l,2) - 31136*pow(l,3) - 12448*pow(l,4) + 3744*pow(l,5) + 1248*pow(l,6))*pow(M,2)*pow(r0,2) + (-25515 + 28538*l + 19386*pow(l,2) - 17600*pow(l,3) - 7040*pow(l,4) + 2112*pow(l,5) + 704*pow(l,6))*M*pow(r0,3) - 4*(-525 + 646*l + 438*pow(l,2) - 400*pow(l,3) - 160*pow(l,4) + 48*pow(l,5) + 16*pow(l,6))*pow(r0,4)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*M*pow(-3*M + r0,2)*sqrt(((1 + 2*l)*pow(r0,3)*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))/((l + pow(l,2))*M))))) - ellK*((210*(13*M - 6*r0)*(2*M - r0)*sqrt((l*(1 + l)*(-2*M + r0))/(1.*(1 + 2*l)*M*(-3*M + r0))))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*pow(r0,1.5)) + (4*(2*M - r0)*(3*pow(M,2) - 6*M*r0 + pow(r0,2))*sqrt((l*(1 + l))/(1.*(1 + 2*l)*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))))/((-3 + 4*l + 4*pow(l,2))*pow(r0,1.5)) + (6*(2*M - r0)*sqrt((l*(1 + l))/(1.*(1 + 2*l)*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(20*pow(M,2) - 35*M*r0 + 14*pow(r0,2)))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*pow(r0,1.5)) - 4*(-3*M + r0)*(-2*M + r0)*sqrt((l*(1 + l))/(1.*(1 + 2*l)*M*(6*pow(M,2)*r0 - 5*M*pow(r0,2) + pow(r0,3)))) + (-6*(7945 - 542*l - 534*pow(l,2) + 16*pow(l,3) + 8*pow(l,4))*pow(M,4) + (118965 - 93686*l - 63886*pow(l,2) + 57296*pow(l,3) + 22888*pow(l,4) - 6912*pow(l,5) - 2304*pow(l,6))*pow(M,3)*r0 + 2*(-43365 + 42670*l + 28982*pow(l,2) - 26320*pow(l,3) - 10520*pow(l,4) + 3168*pow(l,5) + 1056*pow(l,6))*pow(M,2)*pow(r0,2) + (23835 - 25914*l - 17602*pow(l,2) + 15984*pow(l,3) + 6392*pow(l,4) - 1920*pow(l,5) - 640*pow(l,6))*M*pow(r0,3) + 4*(-525 + 646*l + 438*pow(l,2) - 400*pow(l,3) - 160*pow(l,4) + 48*pow(l,5) + 16*pow(l,6))*pow(r0,4))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*M*(3*M - r0)*sqrt(((1 + 2*l)*pow(r0,3)*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))/((l + pow(l,2))*M))))))/(l*(-2 - l + 2*pow(l,2) + pow(l,3))*sqrt(M_PI));
    coeffs[1][1] = ellE*((32*(10*pow(M,2) - 7*M*r0 + pow(r0,2))*sqrt((l + pow(l,2))/((M*M_PI + 2*l*M*M_PI)*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))))/(l*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(r0,1.5)) - (32*sqrt((l + pow(l,2))/((M*M_PI + 2*l*M*M_PI)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(16*pow(M,4) - 46*pow(M,3)*r0 + 35*pow(M,2)*pow(r0,2) - 10*M*pow(r0,3) + pow(r0,4)))/(1.*l*(-1 + 2*l)*(3 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*(3*M - r0)*pow(r0,2.5)) - (1680*sqrt((l*(1 + l)*(-2*M + r0))/(M + 2*l*M))*(578*pow(M,4) - 791*pow(M,3)*r0 + 420*pow(M,2)*pow(r0,2) - 93*M*pow(r0,3) + 6*pow(r0,4)))/(1.*l*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*sqrt(M_PI)*pow(r0*(-3*M + r0),2.5)) - (48*sqrt((l*(1 + l))/((M + 2*l*M)*(-2*M + r0)))*(1060*pow(M,5) - 2608*pow(M,4)*r0 + 2391*pow(M,3)*pow(r0,2) - 1022*pow(M,2)*pow(r0,3) + 201*M*pow(r0,4) - 14*pow(r0,5)))/(1.*l*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*sqrt(M_PI)*pow(r0*(-3*M + r0),2.5)) - (8*(-72*(-7595 - 270*l - 206*pow(l,2) + 128*pow(l,3) + 64*pow(l,4))*pow(M,5) - 6*(-11655 + 154502*l + 113278*pow(l,2) - 79568*pow(l,3) - 32584*pow(l,4) + 8640*pow(l,5) + 2880*pow(l,6))*pow(M,4)*r0 + 9*(-106155 + 155646*l + 123190*pow(l,2) - 65232*pow(l,3) - 33160*pow(l,4) + 64*pow(l,5) + 1216*pow(l,6) + 1024*pow(l,7) + 256*pow(l,8))*pow(M,3)*pow(r0,2) - 48*(-13335 + 13231*l + 12680*pow(l,2) - 1886*pow(l,3) - 2823*pow(l,4) - 2032*pow(l,5) - 304*pow(l,6) + 320*pow(l,7) + 80*pow(l,8))*pow(M,2)*pow(r0,3) + (-154035 + 92574*l + 135958*pow(l,2) + 65200*pow(l,3) - 19528*pow(l,4) - 57536*pow(l,5) - 10816*pow(l,6) + 7168*pow(l,7) + 1792*pow(l,8))*M*pow(r0,4) - 8*(-1575 + 102*l + 1304*pow(l,2) + 1988*pow(l,3) - 14*pow(l,4) - 1120*pow(l,5) - 224*pow(l,6) + 128*pow(l,7) + 32*pow(l,8))*pow(r0,5)))/(3.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(-2 + l + pow(l,2))*sqrt(M_PI)*pow(r0,4)*pow(-3*M + r0,2)*sqrt((l*(1 + 3*l + 2*pow(l,2))*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))/pow(r0,3)))) + ellK*((-32*(12*pow(M,2) - 7*M*r0 + pow(r0,2))*sqrt((l + pow(l,2))/((M*M_PI + 2*l*M*M_PI)*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))))/(1.*l*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(r0,1.5)) - (16*sqrt((l*(1 + l))/((M + 2*l*M)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(2*pow(M,3) + 31*pow(M,2)*r0 - 15*M*pow(r0,2) + 2*pow(r0,3)))/(1.*l*(-1 + 2*l)*(3 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*sqrt(M_PI)*pow(r0,2.5)) + (3360*sqrt((l*(1 + l)*(-2*M + r0))/(M + 2*l*M))*(-109*pow(M,3) + 123*pow(M,2)*r0 - 39*M*pow(r0,2) + 3*pow(r0,3)))/(1.*l*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*sqrt(M_PI)*pow(r0,2.5)*pow(-3*M + r0,1.5)) + (48*(442*pow(M,4) - 931*pow(M,3)*r0 + 640*pow(M,2)*pow(r0,2) - 169*M*pow(r0,3) + 14*pow(r0,4)))/(1.*l*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*sqrt(M_PI)*(3*M - r0)*pow(r0,2.5)*sqrt(((M + 2*l*M)*(-3*M + r0)*(-2*M + r0))/(1.*l*(1 + l)))) + (8*(72*(7945 - 542*l - 534*pow(l,2) + 16*pow(l,3) + 8*pow(l,4))*pow(M,5) - 6*(9975 + 149786*l + 106090*pow(l,2) - 84224*pow(l,3) - 34192*pow(l,4) + 9504*pow(l,5) + 3168*pow(l,6))*pow(M,4)*r0 + 3*(-274575 + 446798*l + 377038*pow(l,2) - 150176*pow(l,3) - 100192*pow(l,4) - 25824*pow(l,5) - 1440*pow(l,6) + 6144*pow(l,7) + 1536*pow(l,8))*pow(M,3)*pow(r0,2) - 12*(-50925 + 51106*l + 53550*pow(l,2) + 120*pow(l,3) - 11412*pow(l,4) - 12512*pow(l,5) - 2080*pow(l,6) + 1792*pow(l,7) + 448*pow(l,8))*pow(M,2)*pow(r0,3) + (-154665 + 91698*l + 143858*pow(l,2) + 79520*pow(l,3) - 20192*pow(l,4) - 66208*pow(l,5) - 12512*pow(l,6) + 8192*pow(l,7) + 2048*pow(l,8))*M*pow(r0,4) - 8*(-1575 + 102*l + 1304*pow(l,2) + 1988*pow(l,3) - 14*pow(l,4) - 1120*pow(l,5) - 224*pow(l,6) + 128*pow(l,7) + 32*pow(l,8))*pow(r0,5)))/(3.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(-2 + l + pow(l,2))*sqrt(M_PI)*(2*M - r0)*(3*M - r0)*pow(r0,4)*sqrt((l*(1 + 3*l + 2*pow(l,2))*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))/pow(r0,3))));
    coeffs[1][2] = (2*(ellE*((4*(M - r0))/((-2 - l + 2*pow(l,2) + pow(l,3))*pow(r0/(1 + 2*l),1.5)*sqrt((l*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))/(1 + l))) - (12*sqrt((l*(1 + l)*M)/(1.*(1 + 2*l)*pow(r0,3)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(30*pow(M,2) - 23*M*r0 + 2*pow(r0,2)))/(l*(-2 - l + 2*pow(l,2) + pow(l,3))*M*r0) + (6*sqrt((l*(1 + l)*M)/(1.*(1 + 2*l)*pow(r0,3)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(50*pow(M,3) - 101*pow(M,2)*r0 + 59*M*pow(r0,2) - 10*pow(r0,3)))/(1.*l*(-1 + 2*l)*(3 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*M*(3*M - r0)*r0) - (90*sqrt((l*(1 + l)*M)/(1.*(1 + 2*l)*pow(r0,3)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(230*pow(M,3) - 257*pow(M,2)*r0 + 85*M*pow(r0,2) - 8*pow(r0,3)))/(1.*l*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*M*(3*M - r0)*r0) + (6*(-4335 - 3166*l - 2230*pow(l,2) + 1872*pow(l,3) + 936*pow(l,4))*pow(M,4) + 3*(23115 + 9062*l - 3810*pow(l,2) - 23184*pow(l,3) - 5192*pow(l,4) + 7680*pow(l,5) + 2560*pow(l,6))*pow(M,3)*r0 - 2*(28845 + 4290*l - 12670*pow(l,2) - 30048*pow(l,3) - 5344*pow(l,4) + 11616*pow(l,5) + 3872*pow(l,6))*pow(M,2)*pow(r0,2) + 3*(5625 + 850*l - 3014*pow(l,2) - 6832*pow(l,3) - 1176*pow(l,4) + 2688*pow(l,5) + 896*pow(l,6))*M*pow(r0,3) - 16*(90 + 42*l - 49*pow(l,2) - 162*pow(l,3) - 31*pow(l,4) + 60*pow(l,5) + 20*pow(l,6))*pow(r0,4))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(-2 + l + pow(l,2))*(2*M - r0)*(3*M - r0)*pow(r0,4)*sqrt((l*(1 + 3*l + 2*pow(l,2))*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))/pow(r0,3)))) + ellK*((-4*(3*M - r0))/((-2 - l + 2*pow(l,2) + pow(l,3))*pow(r0/(1 + 2*l),1.5)*sqrt((l*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))/(1 + l))) + (180*sqrt((l*(1 + l)*M)/((1 + 2*l)*pow(r0,3)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(47*pow(M,2) - 33*M*r0 + 4*pow(r0,2)))/(1.*l*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*M*r0) - (6*sqrt((l*(1 + l)*M)/(1.*(1 + 2*l)*pow(r0,3)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(26*pow(M,3) - 83*pow(M,2)*r0 + 57*M*pow(r0,2) - 10*pow(r0,3)))/(1.*l*(-1 + 2*l)*(3 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*M*(2*M - r0)*r0) + (6*sqrt((l*(1 + l)*M)/(1.*(1 + 2*l)*pow(r0,3)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(162*pow(M,3) - 179*pow(M,2)*r0 + 55*M*pow(r0,2) - 4*pow(r0,3)))/(l*(-2 - l + 2*pow(l,2) + pow(l,3))*M*(2*M - r0)*r0) + (-12*(-1875 - 902*l - 614*pow(l,2) + 576*pow(l,3) + 288*pow(l,4))*pow(M,4) - 24*(2535 + 536*l - 484*pow(l,2) - 1832*pow(l,3) - 396*pow(l,4) + 624*pow(l,5) + 208*pow(l,6))*pow(M,3)*r0 + (51705 + 3738*l - 22390*pow(l,2) - 46176*pow(l,3) - 7888*pow(l,4) + 18240*pow(l,5) + 6080*pow(l,6))*pow(M,2)*pow(r0,2) + (-16065 - 1962*l + 8470*pow(l,2) + 18432*pow(l,3) + 3136*pow(l,4) - 7296*pow(l,5) - 2432*pow(l,6))*M*pow(r0,3) + 16*(90 + 42*l - 49*pow(l,2) - 162*pow(l,3) - 31*pow(l,4) + 60*pow(l,5) + 20*pow(l,6))*pow(r0,4))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(-2 + l + pow(l,2))*pow(r0,4)*pow(-2*M + r0,2)*sqrt((l*(1 + 3*l + 2*pow(l,2))*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))/pow(r0,3))))))/(3.*sqrt(M_PI));
    coeffs[1][3] = ellK*((-4*pow(1 + 2*l,2)*(-21*pow(M,2) + 16*M*r0 - 3*pow(r0,2))*sqrt((l + pow(l,2))/((M*M_PI + 2*l*M*M_PI)*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))))/(3.*l*(-2 - l + 2*pow(l,2) + pow(l,3))*(2*M - r0)*pow(r0,2.5)) + (10*(1074*pow(M,3) - 1367*pow(M,2)*r0 + 517*M*pow(r0,2) - 48*pow(r0,3)))/(1.*l*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*sqrt(M_PI)*pow(r0,3.5)*(-2*M + r0)*sqrt(((1 + 2*l)*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))/(1.*l*(1 + l)))) - (8*sqrt((l + pow(l,2))/((M*M_PI + 2*l*M*M_PI)*(-3*M + r0)))*(288*pow(M,4) - 486*pow(M,3)*r0 + 286*pow(M,2)*pow(r0,2) - 62*M*pow(r0,3) + 3*pow(r0,4)))/(3.*l*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(r0,3.5)*pow(-2*M + r0,2.5)) + (2*sqrt((l*(1 + l))/((M + 2*l*M)*(-3*M + r0)))*(996*pow(M,4) - 2476*pow(M,3)*r0 + 2091*pow(M,2)*pow(r0,2) - 737*M*pow(r0,3) + 90*pow(r0,4)))/(3.*l*(-3 + 4*l + 4*pow(l,2))*(-2 - l + 2*pow(l,2) + pow(l,3))*sqrt(M_PI)*pow(r0,3.5)*pow(-2*M + r0,2.5)) - (2*(-12*(4125 - 3622*l - 2854*pow(l,2) + 1536*pow(l,3) + 768*pow(l,4))*pow(M,4) + 8*(11595 - 8800*l - 8284*pow(l,2) + 1368*pow(l,3) + 1524*pow(l,4) + 1008*pow(l,5) + 336*pow(l,6))*pow(M,3)*r0 + (-61815 + 39642*l + 44906*pow(l,2) + 7136*pow(l,3) - 4912*pow(l,4) - 10176*pow(l,5) - 3392*pow(l,6))*pow(M,2)*pow(r0,2) + (16695 - 8026*l - 12026*pow(l,2) - 6592*pow(l,3) + 224*pow(l,4) + 4224*pow(l,5) + 1408*pow(l,6))*M*pow(r0,3) - 48*(30 - 6*l - 21*pow(l,2) - 26*pow(l,3) - 3*pow(l,4) + 12*pow(l,5) + 4*pow(l,6))*pow(r0,4)))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(-2 + l + pow(l,2))*sqrt(M_PI)*pow(r0,5)*pow(-2*M + r0,2)*sqrt((l*(1 + 3*l + 2*pow(l,2))*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))/pow(r0,3)))) + ellE*((-4*pow(1 + 2*l,2)*sqrt((l + pow(l,2))/((M*M_PI + 2*l*M*M_PI)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(16*pow(M,2) - 13*M*r0 + 3*pow(r0,2)))/(3.*l*(-2 - l + 2*pow(l,2) + pow(l,3))*(2*M - r0)*pow(r0,2.5)) + (2*(402*pow(M,3) - 551*pow(M,2)*r0 + 221*M*pow(r0,2) - 12*pow(r0,3)))/(3.*l*(-2 - l + 2*pow(l,2) + pow(l,3))*sqrt(M_PI)*(2*M - r0)*pow(r0,3.5)*sqrt(((M + 2*l*M)*(-3*M + r0)*(-2*M + r0))/(1.*l*(1 + l)))) + (10*sqrt((l + pow(l,2))/((M*M_PI + 2*l*M*M_PI)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(2652*pow(M,4) - 4454*pow(M,3)*r0 + 2652*pow(M,2)*pow(r0,2) - 634*M*pow(r0,3) + 48*pow(r0,4)))/(1.*l*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*(2*M - r0)*(3*M - r0)*pow(r0,3.5)) - (2*sqrt((l + pow(l,2))/(M*M_PI + 2*l*M*M_PI))*(1446*pow(M,4) - 3119*pow(M,3)*r0 + 2380*pow(M,2)*pow(r0,2) - 773*M*pow(r0,3) + 90*pow(r0,4)))/(3.*l*(-3 + 4*l + 4*pow(l,2))*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(r0,3.5)*pow(-3*M + r0,1.5)*pow(-2*M + r0,1.5)) + (2*(-6*(9645 - 7718*l - 6110*pow(l,2) + 3216*pow(l,3) + 1608*pow(l,4))*pow(M,4) + (105195 - 75034*l - 71650*pow(l,2) + 9840*pow(l,3) + 12600*pow(l,4) + 9216*pow(l,5) + 3072*pow(l,6))*pow(M,3)*r0 - 2*(33765 - 21262*l - 23886*pow(l,2) - 3488*pow(l,3) + 2656*pow(l,4) + 5280*pow(l,5) + 1760*pow(l,6))*pow(M,2)*pow(r0,2) + (17325 - 8422*l - 12350*pow(l,2) - 6448*pow(l,3) + 296*pow(l,4) + 4224*pow(l,5) + 1408*pow(l,6))*M*pow(r0,3) - 48*(30 - 6*l - 21*pow(l,2) - 26*pow(l,3) - 3*pow(l,4) + 12*pow(l,5) + 4*pow(l,6))*pow(r0,4)))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(-2 + l + pow(l,2))*sqrt(M_PI)*(2*M - r0)*(3*M - r0)*pow(r0,5)*sqrt((l*(1 + 3*l + 2*pow(l,2))*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))/pow(r0,3))));

    abscoeffs[1][0] = (-4*sqrt(M_PI)*sqrt((l*(1 + 3*l + 2*pow(l,2))*M)/(-3*M + r0)))/(l*(-2 - l + 2*pow(l,2) + pow(l,3)));
    abscoeffs[1][1] = (4*sqrt(M_PI)*((2*sqrt(-((l*(1 + 3*l + 2*pow(l,2))*M)/(3*M - r0))))/r0 - ((1 + 2*l)*sqrt((l*(1 + l)*M)/(1.*(1 + 2*l)*pow(r0,3)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(18*pow(M,3)*sqrt(r0*(-2*M + r0)) - 21*pow(M,2)*r0*sqrt(r0*(-2*M + r0)) + 8*M*sqrt(pow(r0,5)*(-2*M + r0)) - sqrt(pow(r0,7)*(-2*M + r0))))/((2*M - r0)*pow(-3*M + r0,2))))/(l*(-2 - l + 2*pow(l,2) + pow(l,3)));
    abscoeffs[1][2] = -((1 + 2*l)*sqrt(M_PI)*sqrt((l*(1 + l)*M)/(1.*(1 + 2*l)*pow(r0,4)*(-3*M + r0)))*(32*pow(M,2) - (14 + 9*l + 9*pow(l,2))*M*r0 + 4*l*(1 + l)*pow(r0,2)))/(2.*l*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(-2*M + r0,2));
    abscoeffs[1][3] = (sqrt(M + 2*l*M)*sqrt(M_PI)*(-32*pow(M,2) + (14 + 9*l + 9*pow(l,2))*M*r0 - 4*l*(1 + l)*pow(r0,2)))/(2.*(-2 + l + pow(l,2))*pow(r0,3)*sqrt(l*(1 + l)*(-3*M + r0))*pow(-2*M + r0,2)) - (sqrt(M_PI)*sqrt((l*(1 + l)*(1 + 2*l)*M)/pow(r0,3))*(-32*pow(M,2)*sqrt(r0*(-3*M + r0)) - 4*l*(1 + l)*sqrt(pow(r0,5)*(-3*M + r0)) + M*(14*r0*sqrt(r0*(-3*M + r0)) + 9*l*(1 + l)*sqrt(pow(r0,3)*(-3*M + r0)))))/(3.*l*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(r0,2)*(-3*M + r0)*pow(-2*M + r0,2));

    dcoeffs_dt[1][0] = 0.0;
    dcoeffs_dt[1][1] = 0.0;
    dcoeffs_dt[1][2] = 0.0;
    dcoeffs_dt[1][3] = 0.0;

    dabscoeffs_dt[1][0] = 0.0;
    dabscoeffs_dt[1][1] = 0.0;
    dabscoeffs_dt[1][2] = 0.0;
    dabscoeffs_dt[1][3] = 0.0;

    d2coeffs_dt2[1][0] = 0.0;
    d2coeffs_dt2[1][1] = 0.0;
    d2coeffs_dt2[1][2] = 0.0;
    d2coeffs_dt2[1][3] = 0.0;

    d2abscoeffs_dt2[1][0] = 0.0;
    d2abscoeffs_dt2[1][1] = 0.0;
    d2abscoeffs_dt2[1][2] = 0.0;
    d2abscoeffs_dt2[1][3] = 0.0;

    /* Coefficients of powers of r-r0: m'=2 */
    coeffs[2][0] = (32*m*(1 - (2*M)/r0)*(ellK*((2*(5*M - 2*r0)*(3*M - r0))/sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2)) + (10395*sqrt(-2*M + r0)*(455*pow(M,3) - 694*pow(M,2)*r0 + 303*M*pow(r0,2) - 40*pow(r0,3)))/(1.*(-7 + 2*l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(9 + 2*l)*r0*pow(-3*M + r0,1.5)) + (-52*pow(M,3) + 81*pow(M,2)*r0 - 33*M*pow(r0,2) + 4*pow(r0,3))/(1.*(-1 + 2*l)*(3 + 2*l)*r0*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2))) - (27*(2600*pow(M,4) - 4621*pow(M,3)*r0 + 2780*pow(M,2)*pow(r0,2) - 691*M*pow(r0,3) + 60*pow(r0,4)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*r0*pow(-3*M + r0,1.5)*sqrt(-2*M + r0)) + (4862*pow(M,4) - 13495*pow(M,3)*r0 + 10940*pow(M,2)*pow(r0,2) - 3495*M*pow(r0,3) + 388*pow(r0,4))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*r0*pow(-3*M + r0,1.5)*sqrt(-2*M + r0))) + ellE*((-4*(2*M - r0)*(3*M - r0))/sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2)) + (142*pow(M,4) - 259*pow(M,3)*r0 + 164*pow(M,2)*pow(r0,2) - 43*M*pow(r0,3) + 4*pow(r0,4))/(1.*(-1 + 2*l)*(3 + 2*l)*(3*M - r0)*r0*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2))) + (10395*sqrt(-2*M + r0)*(1148*pow(M,4) - 2181*pow(M,3)*r0 + 1444*pow(M,2)*pow(r0,2) - 403*M*pow(r0,3) + 40*pow(r0,4)))/(1.*(-7 + 2*l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(9 + 2*l)*r0*pow(-3*M + r0,2.5)) + (12332*pow(M,5) - 38388*pow(M,4)*r0 + 40587*pow(M,3)*pow(r0,2) - 19616*pow(M,2)*pow(r0,3) + 4465*M*pow(r0,4) - 388*pow(r0,5))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*r0*pow(-3*M + r0,2.5)*sqrt(-2*M + r0)) + (27*(-6722*pow(M,5) + 14175*pow(M,4)*r0 - 11493*pow(M,3)*pow(r0,2) + 4485*pow(M,2)*pow(r0,3) - 841*M*pow(r0,4) + 60*pow(r0,5)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*r0*pow(-3*M + r0,2.5)*sqrt(-2*M + r0)))))/(1.*(-1 + l)*(2 + l)*sqrt(M_PI)*(2*M - r0)*sqrt((l*(-2 - 5*l + 5*pow(l,3) + 2*pow(l,4))*M)/r0));
    coeffs[2][1] = (512*m*(-((ellE*sqrt((l*(1 + l)*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))/(-2 - 3*l + 3*pow(l,2) + 2*pow(l,3)))*(2*(-10458 - 16199*l - 15347*pow(l,2) + 1704*pow(l,3) + 852*pow(l,4))*pow(M,5) + (-1707615 + 671034*l + 602342*pow(l,2) - 135224*pow(l,3) - 62212*pow(l,4) + 6480*pow(l,5) + 2160*pow(l,6))*pow(M,4)*r0 - 16*(-164241 + 61113*l + 54695*pow(l,2) - 12629*pow(l,3) - 5797*pow(l,4) + 621*pow(l,5) + 207*pow(l,6))*pow(M,3)*pow(r0,2) + (-1430100 + 520469*l + 465645*pow(l,2) - 107872*pow(l,3) - 49496*pow(l,4) + 5328*pow(l,5) + 1776*pow(l,6))*pow(M,2)*pow(r0,3) + (331317 - 118371*l - 105983*pow(l,2) + 24376*pow(l,3) + 11188*pow(l,4) - 1200*pow(l,5) - 400*pow(l,6))*M*pow(r0,4) + 4*(-6930 + 2423*l + 2173*pow(l,2) - 492*pow(l,3) - 226*pow(l,4) + 24*pow(l,5) + 8*pow(l,6))*pow(r0,5)))/(1.*l*(1 + l)*M*(2*M - r0)*pow(3*M - r0,3))) + (ellK*(52*(-630 - 443*l - 419*pow(l,2) + 48*pow(l,3) + 24*pow(l,4))*pow(M,5) + 2*(-691047 + 271018*l + 243426*pow(l,2) - 54320*pow(l,3) - 25000*pow(l,4) + 2592*pow(l,5) + 864*pow(l,6))*pow(M,4)*r0 - 2*(-1131291 + 420047*l + 376059*pow(l,2) - 86560*pow(l,3) - 39740*pow(l,4) + 4248*pow(l,5) + 1416*pow(l,6))*pow(M,3)*pow(r0,2) + (-1302399 + 474009*l + 424105*pow(l,2) - 98192*pow(l,3) - 45056*pow(l,4) + 4848*pow(l,5) + 1616*pow(l,6))*pow(M,2)*pow(r0,3) + (317457 - 113525*l - 101637*pow(l,2) + 23392*pow(l,3) + 10736*pow(l,4) - 1152*pow(l,5) - 384*pow(l,6))*M*pow(r0,4) + 4*(-6930 + 2423*l + 2173*pow(l,2) - 492*pow(l,3) - 226*pow(l,4) + 24*pow(l,5) + 8*pow(l,6))*pow(r0,5)))/(sqrt(l*(-2 - 5*l + 5*pow(l,3) + 2*pow(l,4))*M)*pow(6*pow(M,2) - 5*M*r0 + pow(r0,2),1.5))))/(1.*(-7 + 2*l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(9 + 2*l)*sqrt(M_PI)*pow(r0,2.5));
    coeffs[2][2] = (64*m*((2*ellE*(18*(3150 - 479*l - 457*pow(l,2) + 44*pow(l,3) + 22*pow(l,4))*pow(M,4) - 3*(59675 - 6373*l - 8423*pow(l,2) - 3860*pow(l,3) - 1330*pow(l,4) + 720*pow(l,5) + 240*pow(l,6))*pow(M,3)*r0 + 3*(51345 - 5209*l - 7453*pow(l,2) - 4232*pow(l,3) - 1476*pow(l,4) + 768*pow(l,5) + 256*pow(l,6))*pow(M,2)*pow(r0,2) + 2*(-24360 + 2172*l + 3401*pow(l,2) + 2322*pow(l,3) + 821*pow(l,4) - 408*pow(l,5) - 136*pow(l,6))*M*pow(r0,3) + 8*(630 - 39*l - 77*pow(l,2) - 72*pow(l,3) - 26*pow(l,4) + 12*pow(l,5) + 4*pow(l,6))*pow(r0,4)))/(6*pow(M,2) - 5*M*r0 + pow(r0,2)) - (ellK*(18*(5355 - 877*l - 833*pow(l,2) + 88*pow(l,3) + 44*pow(l,4))*pow(M,4) + 3*(-101255 + 11497*l + 14693*pow(l,2) + 6008*pow(l,3) + 2044*pow(l,4) - 1152*pow(l,5) - 384*pow(l,6))*pow(M,3)*r0 + 3*(91385 - 9467*l - 13371*pow(l,2) - 7360*pow(l,3) - 2560*pow(l,4) + 1344*pow(l,5) + 448*pow(l,6))*pow(M,2)*pow(r0,2) - 4*(23100 - 2094*l - 3247*pow(l,2) - 2178*pow(l,3) - 769*pow(l,4) + 384*pow(l,5) + 128*pow(l,6))*M*pow(r0,3) + 16*(630 - 39*l - 77*pow(l,2) - 72*pow(l,3) - 26*pow(l,4) + 12*pow(l,5) + 4*pow(l,6))*pow(r0,4)))/pow(-2*M + r0,2)))/(3.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*sqrt(M_PI)*pow(r0,4)*sqrt((l*(-2 - 5*l + 5*pow(l,3) + 2*pow(l,4))*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))/pow(r0,3)));
    coeffs[2][3] = (64*m*(-(ellK*(18*(13545 - 4543*l - 4187*pow(l,2) + 712*pow(l,3) + 356*pow(l,4))*pow(M,4) + (-474775 + 174300*l + 157176*pow(l,2) - 33912*pow(l,3) - 16116*pow(l,4) + 1008*pow(l,5) + 336*pow(l,6))*pow(M,3)*r0 + (332185 - 134680*l - 118548*pow(l,2) + 31720*pow(l,3) + 14500*pow(l,4) - 1632*pow(l,5) - 544*pow(l,6))*pow(M,2)*pow(r0,2) + (-97720 + 42833*l + 37021*pow(l,2) - 11384*pow(l,3) - 5092*pow(l,4) + 720*pow(l,5) + 240*pow(l,6))*M*pow(r0,3) - 4*(-2520 + 1177*l + 1003*pow(l,2) - 340*pow(l,3) - 150*pow(l,4) + 24*pow(l,5) + 8*pow(l,6))*pow(r0,4))) + (ellE*(24*(24570 - 8243*l - 7597*pow(l,2) + 1292*pow(l,3) + 646*pow(l,4))*pow(M,5) + 2*(-700490 + 252557*l + 228793*pow(l,2) - 47144*pow(l,3) - 22612*pow(l,4) + 1152*pow(l,5) + 384*pow(l,6))*pow(M,4)*r0 - 2*(-645505 + 252219*l + 223991*pow(l,2) - 55632*pow(l,3) - 25756*pow(l,4) + 2472*pow(l,5) + 824*pow(l,6))*pow(M,3)*pow(r0,2) + (-574490 + 240781*l + 210277*pow(l,2) - 59872*pow(l,3) - 27096*pow(l,4) + 3408*pow(l,5) + 1136*pow(l,6))*pow(M,2)*pow(r0,3) + (122920 - 54603*l - 47051*pow(l,2) + 14784*pow(l,3) + 6592*pow(l,4) - 960*pow(l,5) - 320*pow(l,6))*M*pow(r0,4) + 4*(-2520 + 1177*l + 1003*pow(l,2) - 340*pow(l,3) - 150*pow(l,4) + 24*pow(l,5) + 8*pow(l,6))*pow(r0,5)))/(3*M - r0)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*sqrt(M_PI)*pow(r0,5)*pow(-2*M + r0,2)*sqrt((l*(-2 - 5*l + 5*pow(l,3) + 2*pow(l,4))*M*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))/pow(r0,3)));

    abscoeffs[2][0] = 0.0;
    abscoeffs[2][1] = 0.0;
    abscoeffs[2][2] = 0.0;
    abscoeffs[2][3] = (-2*m*M*sqrt(M_PI)*sqrt((l*(1 + l)*(1 + 2*l)*M*(-3*M + r0))/((-2 + l + pow(l,2))*pow(r0,4))))/(3.*l*(1 + l)*pow(-2*M + r0,3));

    dcoeffs_dt[2][0] = 0.0;
    dcoeffs_dt[2][1] = 0.0;
    dcoeffs_dt[2][2] = 0.0;
    dcoeffs_dt[2][3] = 0.0;

    dabscoeffs_dt[2][0] = 0.0;
    dabscoeffs_dt[2][1] = 0.0;
    dabscoeffs_dt[2][2] = 0.0;
    dabscoeffs_dt[2][3] = 0.0;

    d2coeffs_dt2[2][0] = 0.0;
    d2coeffs_dt2[2][1] = 0.0;
    d2coeffs_dt2[2][2] = 0.0;
    d2coeffs_dt2[2][3] = 0.0;

    d2abscoeffs_dt2[2][0] = 0.0;
    d2abscoeffs_dt2[2][1] = 0.0;
    d2abscoeffs_dt2[2][2] = 0.0;
    d2abscoeffs_dt2[2][3] = 0.0;

    /* Coefficients of powers of r-r0: m'=3 */
    coeffs[3][0] = (2*sqrt((l*(1 + l))/(12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5)))*(-(ellE*((8*pow(1 + 2*l,2)*(23*M - 8*r0)*r0*pow(-2*M + r0,2))/sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2)) + (108*pow(-2*M + r0,2)*(8571*pow(M,3) - 11181*pow(M,2)*r0 + 4864*M*pow(r0,2) - 688*pow(r0,3)))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(3*M - r0)*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2))) - (35640*pow(-2*M + r0,2)*(234*pow(M,3) - 2403*pow(M,2)*r0 + 1819*M*pow(r0,2) - 352*pow(r0,3)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(3*M - r0)*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2))) + (6486480*pow((2*M - r0)/(3*M - r0),1.5)*(237*pow(M,3) - 660*pow(M,2)*r0 + 431*M*pow(r0,2) - 80*pow(r0,3)))/(1.*(-7 + 2*l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(9 + 2*l)) + (12*pow(-2*M + r0,2)*(951*pow(M,3) - 837*pow(M,2)*r0 + 220*M*pow(r0,2) - 16*pow(r0,3)))/((3*M - r0)*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2))) - (2*pow(1 + 2*l,2)*r0*(-398*pow(M,3) + 477*pow(M,2)*r0 - 187*M*pow(r0,2) + 24*pow(r0,3)))/sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2)) - (24*pow(-2*M + r0,2)*(1191*pow(M,3) + 7851*pow(M,2)*r0 - 6368*M*pow(r0,2) + 1208*pow(r0,3)))/(1.*(-1 + 2*l)*(3 + 2*l)*(3*M - r0)*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2))) + (810810*sqrt(-2*M + r0)*(13470*pow(M,5) - 50913*pow(M,4)*r0 + 55976*pow(M,3)*pow(r0,2) - 26699*pow(M,2)*pow(r0,3) + 5838*M*pow(r0,4) - 480*pow(r0,5)))/(1.*(-7 + 2*l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(9 + 2*l)*pow(-3*M + r0,2.5)) + (12*(6996*pow(M,5) - 11732*pow(M,4)*r0 + 7457*pow(M,3)*pow(r0,2) - 2200*pow(M,2)*pow(r0,3) + 289*M*pow(r0,4) - 12*pow(r0,5)))/((3*M - r0)*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2))) - (3*(254988*pow(M,6) + 507996*pow(M,5)*r0 - 1532197*pow(M,4)*pow(r0,2) + 1297920*pow(M,3)*pow(r0,3) - 510803*pow(M,2)*pow(r0,4) + 97224*M*pow(r0,5) - 7248*pow(r0,6)))/((-3 + 4*l + 4*pow(l,2))*pow(-3*M + r0,2.5)*sqrt(-2*M + r0)) + (4455*(33204*pow(M,6) + 221952*pow(M,5)*r0 - 504253*pow(M,4)*pow(r0,2) + 403078*pow(M,3)*pow(r0,3) - 154619*pow(M,2)*pow(r0,4) + 28886*M*pow(r0,5) - 2112*pow(r0,6)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*pow(-3*M + r0,2.5)*sqrt(-2*M + r0)) + (216*(109482*pow(M,6) - 242877*pow(M,5)*r0 + 222282*pow(M,4)*pow(r0,2) - 107669*pow(M,3)*pow(r0,3) + 29223*pow(M,2)*pow(r0,4) - 4237*M*pow(r0,5) + 258*pow(r0,6)))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*pow(-3*M + r0,2.5)*sqrt(-2*M + r0)))) + (8*ellK*(108*(60639075 - 46956024*l - 33807101*pow(l,2) + 24885406*pow(l,3) + 8959139*pow(l,4) - 4047176*pow(l,5) - 1127224*pow(l,6) + 190144*pow(l,7) + 47536*pow(l,8))*pow(M,5)*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2)) + 224*pow(r0,4)*(289575*sqrt(pow(2*M - r0,3)*(3*M - r0)) + (209385 - 291888*l - 293988*pow(l,2) + 10972*pow(l,3) + 41528*pow(l,4) + 38028*pow(l,5) + 4164*pow(l,6) - 6912*pow(l,7) - 1248*pow(l,8) + 320*pow(l,9) + 64*pow(l,10))*r0*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2))) - 6*M*pow(r0,3)*(85270185*sqrt(pow(2*M - r0,3)*(3*M - r0)) + (68878215 - 95781792*l - 102256852*pow(l,2) - 6286808*pow(l,3) + 12737628*pow(l,4) + 16906672*pow(l,5) + 2128464*pow(l,6) - 2851712*pow(l,7) - 519968*pow(l,8) + 128640*pow(l,9) + 25728*pow(l,10))*r0*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2))) + 18*pow(M,4)*(26351325*sqrt(pow(2*M - r0,3)*(3*M - r0)) + (-392953680 + 274431000*l + 147170809*pow(l,2) - 231661862*pow(l,3) - 60319943*pow(l,4) + 62074760*pow(l,5) + 13248888*pow(l,6) - 6136768*pow(l,7) - 1230832*pow(l,8) + 202240*pow(l,9) + 40448*pow(l,10))*r0*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2))) - 3*pow(M,3)*r0*(506350845*sqrt(pow(2*M - r0,3)*(3*M - r0)) + (-459023355 - 13373088*l - 472699388*pow(l,2) - 785391320*pow(l,3) - 71891900*pow(l,4) + 350752928*pow(l,5) + 60998880*pow(l,6) - 45729280*pow(l,7) - 8680960*pow(l,8) + 1834240*pow(l,9) + 366848*pow(l,10))*r0*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2))) + pow(M,2)*pow(r0,2)*(1404728325*sqrt(pow(2*M - r0,3)*(3*M - r0)) + (902065815 - 1510072560*l - 1879786170*pow(l,2) - 556880516*pow(l,3) + 158473286*pow(l,4) + 470447280*pow(l,5) + 68927568*pow(l,6) - 71612544*pow(l,7) - 13252896*pow(l,8) + 3100160*pow(l,9) + 620032*pow(l,10))*r0*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2)))))/(1.*(-7 + 2*l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(9 + 2*l)*(3*M - r0)*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))))/(3.*l*(-2 - l + 2*pow(l,2) + pow(l,3))*sqrt(M_PI)*pow(M*r0,1.5));
    coeffs[3][1] = (sqrt((l*(1 + l))/(12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5)))*((ellK*((2*(-1 + 2*l)*pow(1 + 2*l,2)*(3 + 2*l)*(75*pow(M,2) - 49*M*r0 + 8*pow(r0,2)))/sqrt(pow(M,3)*r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2))) + (623700*(2763*pow(M,3) - 3359*pow(M,2)*r0 + 1322*M*pow(r0,2) - 168*pow(r0,3)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*pow(M,1.5)*r0*sqrt(r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) + (20*pow(1 + 2*l,2)*sqrt(M/(pow(r0,3)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(492*pow(M,3) - 545*pow(M,2)*r0 + 199*M*pow(r0,2) - 24*pow(r0,3)))/pow(M,2) + (40*pow(1 + 2*l,2)*(-492*pow(M,3) + 545*pow(M,2)*r0 - 199*M*pow(r0,2) + 24*pow(r0,3)))/(pow(M*r0,1.5)*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2))) - (30*(11538*pow(M,4) - 14441*pow(M,3)*r0 + 6219*pow(M,2)*pow(r0,2) - 1036*M*pow(r0,3) + 48*pow(r0,4)))/sqrt(pow(M,3)*pow(r0,5)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))) + (5*pow(1 + 2*l,2)*(9222*pow(M,4) - 14693*pow(M,3)*r0 + 8535*pow(M,2)*pow(r0,2) - 2134*M*pow(r0,3) + 192*pow(r0,4)))/(pow(M,1.5)*r0*(-2*M + r0)*sqrt(r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) + (32432400*sqrt(-2*M + r0)*(2775*pow(M,4) - 9315*pow(M,3)*r0 + 7597*pow(M,2)*pow(r0,2) - 2319*M*pow(r0,3) + 240*pow(r0,4)))/(1.*(-7 + 2*l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(9 + 2*l)*pow(r0,2.5)*pow(M*(-3*M + r0),1.5)) + (16216200*sqrt(M*(-2*M + r0))*(2775*pow(M,4) - 9315*pow(M,3)*r0 + 7597*pow(M,2)*pow(r0,2) - 2319*M*pow(r0,3) + 240*pow(r0,4)))/(1.*(-7 + 2*l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(9 + 2*l)*pow(M,2)*(3*M - r0)*pow(r0,2)*sqrt(r0*(-3*M + r0))) - (30*(23970*pow(M,4) - 46121*pow(M,3)*r0 + 31854*pow(M,2)*pow(r0,2) - 9415*M*pow(r0,3) + 1008*pow(r0,4)))/(pow(M,1.5)*r0*(-2*M + r0)*sqrt(r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) - (315*(50990*pow(M,4) - 81177*pow(M,3)*r0 + 47143*pow(M,2)*pow(r0,2) - 11770*M*pow(r0,3) + 1056*pow(r0,4)))/(1.*(-1 + 2*l)*(3 + 2*l)*pow(M,1.5)*r0*(-2*M + r0)*sqrt(r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) - (3780*(66330*pow(M,4) - 111857*pow(M,3)*r0 + 69333*pow(M,2)*pow(r0,2) - 18670*M*pow(r0,3) + 1836*pow(r0,4)))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*pow(M,1.5)*r0*(-2*M + r0)*sqrt(r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) + (540*(361182*pow(M,5) - 647913*pow(M,4)*r0 + 459990*pow(M,3)*pow(r0,2) - 162331*pow(M,2)*pow(r0,3) + 28736*M*pow(r0,4) - 2064*pow(r0,5)))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*pow(r0,2.5)*pow(M*(-3*M + r0),1.5)*sqrt(-2*M + r0)) + (178200*(6864*pow(M,5) + 48636*pow(M,4)*r0 - 83243*pow(M,3)*pow(r0,2) + 48000*pow(M,2)*pow(r0,3) - 11803*M*pow(r0,4) + 1056*pow(r0,5)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*pow(r0,2.5)*pow(M*(-3*M + r0),1.5)*sqrt(-2*M + r0)) - (89100*(6864*pow(M,5) + 48636*pow(M,4)*r0 - 83243*pow(M,3)*pow(r0,2) + 48000*pow(M,2)*pow(r0,3) - 11803*M*pow(r0,4) + 1056*pow(r0,5)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*pow(M,1.5)*pow(r0,2)*(-3*M + r0)*sqrt(r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) + (270*(-361182*pow(M,5) + 647913*pow(M,4)*r0 - 459990*pow(M,3)*pow(r0,2) + 162331*pow(M,2)*pow(r0,3) - 28736*M*pow(r0,4) + 2064*pow(r0,5)))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*pow(M,1.5)*pow(r0,2)*(-3*M + r0)*sqrt(r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) - (120*(52578*pow(M,5) + 127113*pow(M,4)*r0 - 261686*pow(M,3)*pow(r0,2) + 157201*pow(M,2)*pow(r0,3) - 39552*M*pow(r0,4) + 3624*pow(r0,5)))/((-3 + 4*l + 4*pow(l,2))*pow(r0,2.5)*pow(M*(-3*M + r0),1.5)*sqrt(-2*M + r0)) + (60*(52578*pow(M,5) + 127113*pow(M,4)*r0 - 261686*pow(M,3)*pow(r0,2) + 157201*pow(M,2)*pow(r0,3) - 39552*M*pow(r0,4) + 3624*pow(r0,5)))/((-3 + 4*l + 4*pow(l,2))*pow(M,1.5)*pow(r0,2)*(-3*M + r0)*sqrt(r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))))/(-2 - l + 2*pow(l,2) + pow(l,3)) + ellE*((-2*(-1 + 2*l)*pow(1 + 2*l,2)*(3 + 2*l)*(61*pow(M,2) - 45*M*r0 + 8*pow(r0,2)))/((-2 - l + 2*pow(l,2) + pow(l,3))*pow(M,1.5)*sqrt(r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) + (10*pow(1 + 2*l,2)*(1862*pow(M,3) - 2206*pow(M,2)*r0 + 827*M*pow(r0,2) - 96*pow(r0,3)))/((-2 - l + 2*pow(l,2) + pow(l,3))*pow(M,1.5)*r0*sqrt(r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) + (20*pow(1 + 2*l,2)*sqrt(M/(pow(r0,3)*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(398*pow(M,3) - 477*pow(M,2)*r0 + 187*M*pow(r0,2) - 24*pow(r0,3)))/((2 + l - 2*pow(l,2) - pow(l,3))*pow(M,2)) + (40*pow(1 + 2*l,2)*(398*pow(M,3) - 477*pow(M,2)*r0 + 187*M*pow(r0,2) - 24*pow(r0,3)))/((-2 - l + 2*pow(l,2) + pow(l,3))*pow(M*r0,1.5)*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2))) + (240*sqrt(M*r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))*(3498*pow(M,4) - 4117*pow(M,3)*r0 + 1670*pow(M,2)*pow(r0,2) - 265*M*pow(r0,3) + 12*pow(r0,4)))/((-2 - l + 2*pow(l,2) + pow(l,3))*pow(M,2)*pow(r0,3)*pow(-3*M + r0,2)) + (311850*(13398*pow(M,4) - 21985*pow(M,3)*r0 + 13265*pow(M,2)*pow(r0,2) - 3484*M*pow(r0,3) + 336*pow(r0,4)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(M,1.5)*r0*(-3*M + r0)*sqrt(r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) + (630*(30898*pow(M,4) - 46886*pow(M,3)*r0 + 25887*pow(M,2)*pow(r0,2) - 6149*M*pow(r0,3) + 528*pow(r0,4)))/(1.*(-1 + 2*l)*(3 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(M,1.5)*r0*(-3*M + r0)*sqrt(r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) + (15*(58086*pow(M,4) - 107447*pow(M,3)*r0 + 70729*pow(M,2)*pow(r0,2) - 19838*M*pow(r0,3) + 2016*pow(r0,4)))/((-2 - l + 2*pow(l,2) + pow(l,3))*pow(M,1.5)*r0*(-3*M + r0)*sqrt(r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) + (945*(321594*pow(M,4) - 518323*pow(M,3)*r0 + 305951*pow(M,2)*pow(r0,2) - 78352*M*pow(r0,3) + 7344*pow(r0,4)))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(M,1.5)*r0*(-3*M + r0)*sqrt(r0*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))) + (16216200*M*sqrt(-2*M + r0)*(13470*pow(M,5) - 50913*pow(M,4)*r0 + 55976*pow(M,3)*pow(r0,2) - 26699*pow(M,2)*pow(r0,3) + 5838*M*pow(r0,4) - 480*pow(r0,5)))/(1.*(-7 + 2*l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(9 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(-(M*(3*M - r0)*r0),2.5)) - (8108100*sqrt(M*(-2*M + r0))*(13470*pow(M,5) - 50913*pow(M,4)*r0 + 55976*pow(M,3)*pow(r0,2) - 26699*pow(M,2)*pow(r0,3) + 5838*M*pow(r0,4) - 480*pow(r0,5)))/((-7 + 2*l)*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(9 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(M,2)*pow(r0*(-3*M + r0),2.5)) - (4320*sqrt(-(M*(2*M - r0)*r0))*(54741*pow(M,5) - 94068*pow(M,4)*r0 + 64107*pow(M,3)*pow(r0,2) - 21781*pow(M,2)*pow(r0,3) + 3721*M*pow(r0,4) - 258*pow(r0,5)))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(M,2)*pow(r0,3)*pow(-3*M + r0,2.5)) - (120*(6996*pow(M,5) - 11732*pow(M,4)*r0 + 7457*pow(M,3)*pow(r0,2) - 2200*pow(M,2)*pow(r0,3) + 289*M*pow(r0,4) - 12*pow(r0,5)))/(1.*(2 + l - 2*pow(l,2) - pow(l,3))*pow(M,1.5)*pow(r0,2.5)*(-3*M + r0)*sqrt(6*pow(M,2) - 5*M*r0 + pow(r0,2))) - (89100*sqrt(-(M*(2*M - r0)*r0))*(16602*pow(M,5) + 119277*pow(M,4)*r0 - 192488*pow(M,3)*pow(r0,2) + 105295*pow(M,2)*pow(r0,3) - 24662*M*pow(r0,4) + 2112*pow(r0,5)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(M,2)*pow(r0,3)*pow(-3*M + r0,2.5)) + (60*sqrt(-(M*(2*M - r0)*r0))*(127494*pow(M,5) + 317745*pow(M,4)*r0 - 607226*pow(M,3)*pow(r0,2) + 345347*pow(M,2)*pow(r0,3) - 82728*M*pow(r0,4) + 7248*pow(r0,5)))/(1.*(-1 + 2*l)*(3 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(M,2)*pow(r0,3)*pow(-3*M + r0,2.5)) + (30*sqrt(-(M/(2*M - r0)))*(254988*pow(M,6) + 507996*pow(M,5)*r0 - 1532197*pow(M,4)*pow(r0,2) + 1297920*pow(M,3)*pow(r0,3) - 510803*pow(M,2)*pow(r0,4) + 97224*M*pow(r0,5) - 7248*pow(r0,6)))/((-3 + 4*l + 4*pow(l,2))*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(M,2)*pow(r0*(-3*M + r0),2.5)) - (44550*sqrt(-(M/(2*M - r0)))*(33204*pow(M,6) + 221952*pow(M,5)*r0 - 504253*pow(M,4)*pow(r0,2) + 403078*pow(M,3)*pow(r0,3) - 154619*pow(M,2)*pow(r0,4) + 28886*M*pow(r0,5) - 2112*pow(r0,6)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(M,2)*pow(r0*(-3*M + r0),2.5)) - (2160*sqrt(-(M/(2*M - r0)))*(109482*pow(M,6) - 242877*pow(M,5)*r0 + 222282*pow(M,4)*pow(r0,2) - 107669*pow(M,3)*pow(r0,3) + 29223*pow(M,2)*pow(r0,4) - 4237*M*pow(r0,5) + 258*pow(r0,6)))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(M,2)*pow(r0*(-3*M + r0),2.5)))))/(15.*l*sqrt(M_PI));
    coeffs[3][2] = (ellK*((831600*sqrt((l*(1 + l)*M)/((12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5))*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(2763*pow(M,3) - 3359*pow(M,2)*r0 + 1322*M*pow(r0,2) - 168*pow(r0,3)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*pow(r0,2.5)) + (2*(-1 + 2*l)*pow(1 + 2*l,2)*(3 + 2*l)*sqrt((l*(1 + l)*M)/((12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5))*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(705*pow(M,3) - 832*pow(M,2)*r0 + 319*M*pow(r0,2) - 40*pow(r0,3)))/(pow(r0,1.5)*(-2*M + r0)) + (831600*sqrt((l*(1 + l)*M)/((12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5))*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(-2763*pow(M,3) + 3359*pow(M,2)*r0 - 1322*M*pow(r0,2) + 168*pow(r0,3)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*pow(r0,2.5)) + (311850*sqrt((l*(1 + l)*M)/((12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5))*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(26502*pow(M,4) - 50373*pow(M,3)*r0 + 34363*pow(M,2)*pow(r0,2) - 9980*M*pow(r0,3) + 1040*pow(r0,4)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*pow(r0,2.5)*(-2*M + r0)) - (945*sqrt((l*(1 + l)*M)/((12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5))*(-3*M + r0)))*(1107588*pow(M,5) - 2659368*pow(M,4)*r0 + 2477453*pow(M,3)*pow(r0,2) - 1119449*pow(M,2)*pow(r0,3) + 245016*M*pow(r0,4) - 20720*pow(r0,5)))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*pow(r0*(-2*M + r0),2.5)) - (15*sqrt((l*(1 + l)*M)/((12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5))*(-3*M + r0)))*(120012*pow(M,5) - 447892*pow(M,4)*r0 + 537357*pow(M,3)*pow(r0,2) - 288871*pow(M,2)*pow(r0,3) + 72354*M*pow(r0,4) - 6880*pow(r0,5)))/pow(r0*(-2*M + r0),2.5) - (210*sqrt((l*(1 + l)*M)/((12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5))*(-3*M + r0)))*(270648*pow(M,5) - 599498*pow(M,4)*r0 + 514993*pow(M,3)*pow(r0,2) - 213114*pow(M,2)*pow(r0,3) + 42081*M*pow(r0,4) - 3120*pow(r0,5)))/((-3 + 4*l + 4*pow(l,2))*pow(r0*(-2*M + r0),2.5)) + (10*pow(1 + 2*l,2)*sqrt((l*(1 + l)*M)/((12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5))*(-3*M + r0)))*(8376*pow(M,5) - 21770*pow(M,4)*r0 + 20981*pow(M,3)*pow(r0,2) - 9478*pow(M,2)*pow(r0,3) + 2013*M*pow(r0,4) - 160*pow(r0,5)))/pow(r0*(-2*M + r0),2.5)))/(20.*l*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(M,2)*sqrt(M_PI)) + (ellE*((40*pow(1 + 2*l,2)*sqrt((l*(1 + l)*M)/((12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5))*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(1862*pow(M,3) - 2206*pow(M,2)*r0 + 827*M*pow(r0,2) - 96*pow(r0,3)))/((-2 - l + 2*pow(l,2) + pow(l,3))*pow(r0,2.5)) - (6*(-1 + 2*l)*pow(1 + 2*l,2)*(3 + 2*l)*sqrt((l*(1 + l)*M)/((12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5))*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(568*pow(M,3) - 725*pow(M,2)*r0 + 299*M*pow(r0,2) - 40*pow(r0,3)))/((-2 - l + 2*pow(l,2) + pow(l,3))*pow(r0,1.5)*(-2*M + r0)) + (40*pow(1 + 2*l,2)*sqrt((l*(1 + l)*M)/((12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5))*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(-1862*pow(M,3) + 2206*pow(M,2)*r0 - 827*M*pow(r0,2) + 96*pow(r0,3)))/((-2 - l + 2*pow(l,2) + pow(l,3))*pow(r0,2.5)) + (15*pow(1 + 2*l,2)*sqrt((l*(1 + l)*M)/((12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5))*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(6770*pow(M,4) - 14847*pow(M,3)*r0 + 10911*pow(M,2)*pow(r0,2) - 3226*M*pow(r0,3) + 320*pow(r0,4)))/((-2 - l + 2*pow(l,2) + pow(l,3))*pow(r0,2.5)*(-2*M + r0)) + (315*sqrt((l*(1 + l)*M)/((12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5))*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(656250*pow(M,5) - 1404377*pow(M,4)*r0 + 1160912*pow(M,3)*pow(r0,2) - 460899*pow(M,2)*pow(r0,3) + 87282*M*pow(r0,4) - 6240*pow(r0,5)))/((-3 + 4*l + 4*pow(l,2))*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(r0,2.5)*(-3*M + r0)*(-2*M + r0)) + (1871100*sqrt((l*(1 + l)*M)/((12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5))*(6*pow(M,2) - 5*M*r0 + pow(r0,2))))*(32130*pow(M,5) - 74721*pow(M,4)*r0 + 67448*pow(M,3)*pow(r0,2) - 29559*pow(M,2)*pow(r0,3) + 6290*M*pow(r0,4) - 520*pow(r0,5)))/(1.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(r0,2.5)*(-3*M + r0)*(-2*M + r0)) + (22680*sqrt((l*(1 + l)*M)/(12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5)))*(335700*pow(M,6) - 948654*pow(M,5)*r0 + 1091686*pow(M,4)*pow(r0,2) - 654980*pow(M,3)*pow(r0,3) + 216013*pow(M,2)*pow(r0,4) - 37102*M*pow(r0,5) + 2590*pow(r0,6)))/(1.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(2 + l - 2*pow(l,2) - pow(l,3))*pow(-3*M + r0,1.5)*pow(r0*(-2*M + r0),2.5)) + (90*sqrt((l*(1 + l)*M)/(12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5)))*(145500*pow(M,6) - 604844*pow(M,5)*r0 + 881501*pow(M,4)*pow(r0,2) - 624605*pow(M,3)*pow(r0,3) + 234233*pow(M,2)*pow(r0,4) - 44777*M*pow(r0,5) + 3440*pow(r0,6)))/((2 + l - 2*pow(l,2) - pow(l,3))*pow(-3*M + r0,1.5)*pow(r0*(-2*M + r0),2.5))))/(60.*l*pow(M,2)*sqrt(M_PI));
    coeffs[3][3] = (4*(-6 + l + pow(l,2))*(-((ellK*sqrt((l*(1 + l)*(6*pow(M,2) - 5*M*r0 + pow(r0,2)))/(12 + 16*l - 23*pow(l,2) - 12*pow(l,3) + 5*pow(l,4) + 2*pow(l,5)))*(-60*(-610365 - 231692*l - 209356*pow(l,2) + 44672*pow(l,3) + 22336*pow(l,4))*pow(M,5) + 40*(-985215 - 1149538*l - 1088402*pow(l,2) + 126784*pow(l,3) + 74672*pow(l,4) + 13536*pow(l,5) + 4512*pow(l,6))*pow(M,4)*r0 - (2134125 - 51869020*l - 50255644*pow(l,2) + 3529984*pow(l,3) + 2523072*pow(l,4) + 909696*pow(l,5) + 303232*pow(l,6))*pow(M,3)*pow(r0,2) + 5*(3415125 - 5335732*l - 5248820*pow(l,2) + 211456*pow(l,3) + 199808*pow(l,4) + 112896*pow(l,5) + 37632*pow(l,6))*pow(M,2)*pow(r0,3) - 12*(601125 - 535870*l - 533502*pow(l,2) + 8992*pow(l,3) + 15136*pow(l,4) + 12768*pow(l,5) + 4256*pow(l,6))*M*pow(r0,4) + 160*(5775 - 3680*l - 3704*pow(l,2) - 16*pow(l,3) + 72*pow(l,4) + 96*pow(l,5) + 32*pow(l,6))*pow(r0,5)))/(1.*l*(1 + l)*pow(M,1.5)*pow(2*M - r0,3)*(3*M - r0))) + (ellE*(-150*(-295995 - 112360*l - 101528*pow(l,2) + 21664*pow(l,3) + 10832*pow(l,4))*pow(M,5) + (-44461725 - 54494840*l - 51644104*pow(l,2) + 5919584*pow(l,3) + 3505072*pow(l,4) + 654336*pow(l,5) + 218112*pow(l,6))*pow(M,4)*r0 - 2*(2366175 - 29625420*l - 28731644*pow(l,2) + 1963104*pow(l,3) + 1420432*pow(l,4) + 526656*pow(l,5) + 175552*pow(l,6))*pow(M,3)*pow(r0,2) + (19585125 - 29194680*l - 28741352*pow(l,2) + 1114272*pow(l,3) + 1076176*pow(l,4) + 622848*pow(l,5) + 207616*pow(l,6))*pow(M,2)*pow(r0,3) - 4*(1918875 - 1681210*l - 1674586*pow(l,2) + 26656*pow(l,3) + 46848*pow(l,4) + 40224*pow(l,5) + 13408*pow(l,6))*M*pow(r0,4) + 160*(5775 - 3680*l - 3704*pow(l,2) - 16*pow(l,3) + 72*pow(l,4) + 96*pow(l,5) + 32*pow(l,6))*pow(r0,5)))/(sqrt(l*(12 + 28*l - 7*pow(l,2) - 35*pow(l,3) - 7*pow(l,4) + 7*pow(l,5) + 2*pow(l,6)))*pow(M*(6*pow(M,2) - 5*M*r0 + pow(r0,2)),1.5))))/(15.*(-5 + 2*l)*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*(7 + 2*l)*sqrt(M_PI)*pow(r0,3.5));

    abscoeffs[3][0] = 0.0;
    abscoeffs[3][1] = 0.0;
    abscoeffs[3][2] = -((-6 + l + pow(l,2))*M*sqrt(M_PI)*sqrt(-((l*(1 + l)*(1 + 2*l)*M)/((12 - 8*l - 7*pow(l,2) + 2*pow(l,3) + pow(l,4))*(3*M - r0)))))/(2.*l*(1 + l)*r0*pow(-2*M + r0,2));
    abscoeffs[3][3] = -((-6 + l + pow(l,2))*M*sqrt(M_PI)*sqrt(-((l*(1 + l)*(1 + 2*l)*M)/((12 - 8*l - 7*pow(l,2) + 2*pow(l,3) + pow(l,4))*(3*M - r0)*pow(r0,4)))))/(6.*l*(1 + l)*pow(-2*M + r0,2));

    dcoeffs_dt[3][0] = 0.0;
    dcoeffs_dt[3][1] = 0.0;
    dcoeffs_dt[3][2] = 0.0;
    dcoeffs_dt[3][3] = 0.0;

    dabscoeffs_dt[3][0] = 0.0;
    dabscoeffs_dt[3][1] = 0.0;
    dabscoeffs_dt[3][2] = 0.0;
    dabscoeffs_dt[3][3] = 0.0;

    d2coeffs_dt2[3][0] = 0.0;
    d2coeffs_dt2[3][1] = 0.0;
    d2coeffs_dt2[3][2] = 0.0;
    d2coeffs_dt2[3][3] = 0.0;

    d2abscoeffs_dt2[3][0] = 0.0;
    d2abscoeffs_dt2[3][1] = 0.0;
    d2abscoeffs_dt2[3][2] = 0.0;
    d2abscoeffs_dt2[3][3] = 0.0;
  }

if( parity == EvenParity) {
}
  /* Time derivatives of worldline quantities */
  rt = 0.0;
  rtt = 0.0;
  phit = sqrt(M*pow(r0,-3));
  phitt = 0.0;

  c = 0.0;
  ct = 0.0;
  ctt = 0.0;
}


void RWZ_EffectiveSource::operator()(const int n, const double r[],
  const double W[], const double dW_dr[], const double d2W_dr2[],
  double src_re[], double src_im[])
{
  const double md = m;
  const complex<double> I(0.0, 1.0);

  for (int i=0; i<n; i++) {
    /* Distance from the worldline */
    const double dr = r[i]-r0;

    /* The puncture field is given by:
     * \Psi = D^l_{mm'} e^{i m [- \phi_0(t) + c(t) (r - r_0(t))]}
     *        \sum_i a_i(t) [r-r0(t)]^i + b_i(t) |r-r_0(t)| [r-r_0(t)]^i
     *
     * First, compute the term without the phase factor
     * \Psi = D^l_{mm'} \sum_i a_i(t) [r-r0(t)]^i + b_i(t) |r-r_0(t)| [r-r_0(t)]^i
     */
    const double Psi_l0 =
      coeffs[0][0] + dr*(coeffs[0][1] + dr*(coeffs[0][2] + dr*coeffs[0][3]))
      + sgn(dr)*(abscoeffs[0][0] + dr*(abscoeffs[0][1] + dr*(abscoeffs[0][2] + dr*abscoeffs[0][3])));
    const double Psi_l1 =
      coeffs[1][0] + dr*(coeffs[1][1] + dr*(coeffs[1][2] + dr*coeffs[1][3]))
      + sgn(dr)*(abscoeffs[1][0] + dr*(abscoeffs[1][1] + dr*(abscoeffs[1][2] + dr*abscoeffs[1][3])));
    const double Psi_l2 =
      coeffs[2][0] + dr*(coeffs[2][1] + dr*(coeffs[2][2] + dr*coeffs[2][3]))
      + sgn(dr)*(abscoeffs[2][0] + dr*(abscoeffs[2][1] + dr*(abscoeffs[2][2] + dr*abscoeffs[2][3])));
    const double Psi_l3 = l < 3 ? 0 :
      coeffs[3][0] + dr*(coeffs[3][1] + dr*(coeffs[3][2] + dr*coeffs[3][3]))
      + sgn(dr)*(abscoeffs[3][0] + dr*(abscoeffs[3][1] + dr*(abscoeffs[3][2] + dr*abscoeffs[3][3])));

    /* Radial derivative */
    const double dPsi_dr_l0 =
      coeffs[0][1] + dr*(2.0*coeffs[0][2] + dr*3.0*coeffs[0][3])
      + sgn(dr)*(abscoeffs[0][1] + dr*(2.0*abscoeffs[0][2] + dr*3.0*abscoeffs[0][3]));
    const double dPsi_dr_l1 =
      coeffs[1][1] + dr*(2.0*coeffs[1][2] + dr*3.0*coeffs[1][3])
      + sgn(dr)*(abscoeffs[1][1] + dr*(2.0*abscoeffs[1][2] + dr*3.0*abscoeffs[1][3]));
    const double dPsi_dr_l2 =
      coeffs[2][1] + dr*(2.0*coeffs[2][2] + dr*3.0*coeffs[2][3])
      + sgn(dr)*(abscoeffs[2][1] + dr*(2.0*abscoeffs[2][2] + dr*3.0*abscoeffs[2][3]));
    const double dPsi_dr_l3 = l < 3 ? 0 :
      coeffs[3][1] + dr*(2.0*coeffs[3][2] + dr*3.0*coeffs[3][3] )
      + sgn(dr)*(abscoeffs[3][1] + dr*(2.0*abscoeffs[3][2] + dr*3.0*abscoeffs[3][3]));

    /* Second radial derivative */
    const double d2Psi_dr2_l0 =
      2.0*coeffs[0][2] + dr*6.0*coeffs[0][3]
      + sgn(dr)*(2.0*abscoeffs[0][2] + dr*6.0*abscoeffs[0][3]);
    const double d2Psi_dr2_l1 =
      2.0*coeffs[1][2] + dr*6.0*coeffs[1][3]
      + sgn(dr)*(2.0*abscoeffs[1][2] + dr*6.0*abscoeffs[1][3]);
    const double d2Psi_dr2_l2 =
      2.0*coeffs[2][2] + dr*6.0*coeffs[2][3]
      + sgn(dr)*(2.0*abscoeffs[2][2] + dr*6.0*abscoeffs[2][3]);
    const double d2Psi_dr2_l3 = l < 3 ? 0 :
      2.0*coeffs[3][2] + dr*6.0*coeffs[3][3]
      + sgn(dr)*(2.0*abscoeffs[3][2] + dr*6.0*abscoeffs[3][3]);

    /* The time derivative hits the coefficients and powers of dr = r - r_0(t).
     * First compute the derivative of the coefficients.
     */
    const double dPsi_dt_l0 =
      dcoeffs_dt[0][0] + dr*(dcoeffs_dt[0][1] + dr*(dcoeffs_dt[0][2] + dr*dcoeffs_dt[0][3]))
      + sgn(dr)*(dabscoeffs_dt[0][0] + dr*(dabscoeffs_dt[0][1] + dr*(dabscoeffs_dt[0][2] + dr*dabscoeffs_dt[0][3])));
    const double dPsi_dt_l1 =
      dcoeffs_dt[1][0] + dr*(dcoeffs_dt[1][1] + dr*(dcoeffs_dt[1][2] + dr*dcoeffs_dt[1][3]))
      + sgn(dr)*(dabscoeffs_dt[1][0] + dr*(dabscoeffs_dt[1][1] + dr*(dabscoeffs_dt[1][2] + dr*dabscoeffs_dt[1][3])));
    const double dPsi_dt_l2 =
      dcoeffs_dt[2][0] + dr*(dcoeffs_dt[2][1] + dr*(dcoeffs_dt[2][2] + dr*dcoeffs_dt[2][3]))
      + sgn(dr)*(dabscoeffs_dt[2][0] + dr*(dabscoeffs_dt[2][1] + dr*(dabscoeffs_dt[2][2] + dr*dabscoeffs_dt[2][3])));
    const double dPsi_dt_l3 = l < 3 ? 0 :
      dcoeffs_dt[3][0] + dr*(dcoeffs_dt[3][1] + dr*(dcoeffs_dt[3][2] + dr*dcoeffs_dt[3][3]))
      + sgn(dr)*(dabscoeffs_dt[3][0] + dr*(dabscoeffs_dt[3][1] + dr*(dabscoeffs_dt[3][2] + dr*dabscoeffs_dt[3][3])));

    /* And the second derivative of the coefficients */
    const double d2Psi_dt2_l0 =
      d2coeffs_dt2[0][0] + dr*(d2coeffs_dt2[0][1] + dr*(d2coeffs_dt2[0][2] + dr*d2coeffs_dt2[0][3]))
      + sgn(dr)*(d2abscoeffs_dt2[0][0] + dr*(d2abscoeffs_dt2[0][1] + dr*(d2abscoeffs_dt2[0][2] + dr*d2abscoeffs_dt2[0][3])));
    const double d2Psi_dt2_l1 =
      d2coeffs_dt2[1][0] + dr*(d2coeffs_dt2[1][1] + dr*(d2coeffs_dt2[1][2] + dr*d2coeffs_dt2[1][3]))
      + sgn(dr)*(d2abscoeffs_dt2[1][0] + dr*(d2abscoeffs_dt2[1][1] + dr*(d2abscoeffs_dt2[1][2] + dr*d2abscoeffs_dt2[1][3])));
    const double d2Psi_dt2_l2 =
      d2coeffs_dt2[2][0] + dr*(d2coeffs_dt2[2][1] + dr*(d2coeffs_dt2[2][2] + dr*d2coeffs_dt2[2][3]))
      + sgn(dr)*(d2abscoeffs_dt2[2][0] + dr*(d2abscoeffs_dt2[2][1] + dr*(d2abscoeffs_dt2[2][2] + dr*d2abscoeffs_dt2[2][3])));
    const double d2Psi_dt2_l3 = l < 3 ? 0 :
      d2coeffs_dt2[3][0] + dr*(d2coeffs_dt2[3][1] + dr*(d2coeffs_dt2[3][2] + dr*d2coeffs_dt2[3][3]))
      + sgn(dr)*(d2abscoeffs_dt2[3][0] + dr*(d2abscoeffs_dt2[3][1] + dr*(d2abscoeffs_dt2[3][2] + dr*d2abscoeffs_dt2[3][3])));

    /* The r-derivative of the time derivative of the coefficients */
    const double d2Psi_dtr_l0 =
      dcoeffs_dt[0][1] + dr*(2.0*dcoeffs_dt[0][2] + dr*3.0*dcoeffs_dt[0][3])
      + sgn(dr)*(dabscoeffs_dt[0][1] + dr*(2.0*dabscoeffs_dt[0][2] + dr*3.0*dabscoeffs_dt[0][3]));
    const double d2Psi_dtr_l1 =
      dcoeffs_dt[1][1] + dr*(2.0*dcoeffs_dt[1][2] + dr*3.0*dcoeffs_dt[1][3])
      + sgn(dr)*(dabscoeffs_dt[1][1] + dr*(2.0*dabscoeffs_dt[1][2] + dr*3.0*dabscoeffs_dt[1][3]));
    const double d2Psi_dtr_l2 =
      dcoeffs_dt[2][1] + dr*(2.0*dcoeffs_dt[2][2] + dr*3.0*dcoeffs_dt[2][3])
      + sgn(dr)*(dabscoeffs_dt[2][1] + dr*(2.0*dabscoeffs_dt[2][2] + dr*3.0*dabscoeffs_dt[2][3]));
    const double d2Psi_dtr_l3 = l < 3 ? 0 :
      dcoeffs_dt[3][1] + dr*(2.0*dcoeffs_dt[3][2] + dr*3.0*dcoeffs_dt[3][3])
      + sgn(dr)*(dabscoeffs_dt[3][1] + dr*(2.0*dabscoeffs_dt[3][2] + dr*3.0*dabscoeffs_dt[3][3]));

    complex<double> Psi, dPsi_dr, d2Psi_dr2, dPsi_dt, d2Psi_dt2, d2Psi_dtr;
    if (parity == EvenParity) {
      Psi = w0 * Psi_l0 + (w1p - w1m) * I * Psi_l1 + (w2p + w2m) * Psi_l2 + (w3p - w3m) * I * Psi_l3;
      dPsi_dr = w0 * dPsi_dr_l0 + (w1p - w1m) * I * dPsi_dr_l1 + (w2p + w2m) * dPsi_dr_l2 + (w3p - w3m) * I * dPsi_dr_l3;
      d2Psi_dr2 = w0 * d2Psi_dr2_l0 + (w1p - w1m) * I * d2Psi_dr2_l1 + (w2p + w2m) * d2Psi_dr2_l2 + (w3p - w3m) * I * d2Psi_dr2_l3;
      dPsi_dt = w0 * dPsi_dt_l0 + (w1p - w1m) * I * dPsi_dt_l1 + (w2p + w2m) * dPsi_dt_l2 + (w3p - w3m) * I * dPsi_dt_l3;
      d2Psi_dt2 = w0 * d2Psi_dt2_l0 + (w1p - w1m) * I * d2Psi_dt2_l1 + (w2p + w2m) * d2Psi_dt2_l2 + (w3p - w3m) * I * d2Psi_dt2_l3;
      d2Psi_dtr = w0 * d2Psi_dtr_l0 + (w1p - w1m) * I * d2Psi_dtr_l1 + (w2p + w2m) * d2Psi_dtr_l2 + (w3p - w3m) * I * d2Psi_dtr_l3;
    } else if (parity == OddParity) {
      Psi = (w1p + w1m) * I * Psi_l1 + (w2p - w2m) * Psi_l2 + (w3p + w3m) * I * Psi_l3;
      dPsi_dr = (w1p + w1m) * I * dPsi_dr_l1 + (w2p - w2m) * dPsi_dr_l2 + (w3p + w3m) * I * dPsi_dr_l3;
      d2Psi_dr2 = (w1p + w1m) * I * d2Psi_dr2_l1 + (w2p - w2m) * d2Psi_dr2_l2 + (w3p + w3m) * I * d2Psi_dr2_l3;
      dPsi_dt = (w1p + w1m) * I * dPsi_dt_l1 + (w2p - w2m) * dPsi_dt_l2 + (w3p + w3m) * I * dPsi_dt_l3;
      d2Psi_dt2 = (w1p + w1m) * I * d2Psi_dt2_l1 + (w2p - w2m) * d2Psi_dt2_l2 + (w3p + w3m) * I * d2Psi_dt2_l3;
      d2Psi_dtr = (w1p + w1m) * I * d2Psi_dtr_l1 + (w2p - w2m) * d2Psi_dtr_l2 + (w3p + w3m) * I * d2Psi_dtr_l3;
    }

    /* Add in derivative of time dependence of dr = r - r_0(t) */
    d2Psi_dt2 += - dPsi_dr*rtt + d2Psi_dr2*rt*rt - 2.0*d2Psi_dtr*rt;
    dPsi_dt += -dPsi_dr*rt;

    /* Now add in phase - note the order here is important */
    complex<double> phase = exp(complex<double>(0,md*(-phi0 + c*dr)));
    d2Psi_dt2 = phase*(
      - pow(md*(-ct*dr + phit + c*rt),2)*Psi
      - I*md*Psi*(2.0*ct*rt - ctt*dr + c*rtt + phitt)
      + 2.0*I*md*dPsi_dt*(ct*dr - phit - c*rt)
      + d2Psi_dt2);
    dPsi_dt = phase*(I*md*Psi*(ct*dr - phit - c*rt) + dPsi_dt);
    d2Psi_dr2 = phase*(d2Psi_dr2 + 2.0*I*c*md*dPsi_dr - c*c*md*md*Psi);
    dPsi_dr = phase*(dPsi_dr + I*c*md*Psi);
    Psi = phase*Psi;

    /* Add in factor of the window function: Psi -> W*Psi */
    d2Psi_dt2 = W[i] * d2Psi_dt2;
    dPsi_dt = W[i] * dPsi_dt;
    d2Psi_dr2 = W[i] * d2Psi_dr2 + 2.0*dW_dr[i]*dPsi_dr + d2W_dr2[i]*Psi;
    dPsi_dr = W[i] * dPsi_dr + dW_dr[i]*Psi;
    Psi = W[i] * Psi;

    /* Add in factor of time window function */
    const auto TPsi = T*Psi;
    const auto dTPsi_dr = T*dPsi_dr;
    const auto d2TPsi_dr2 = T*d2Psi_dr2;
    const auto d2TPsi_dt2 = T*d2Psi_dt2 + 2.0*dT_dt*dPsi_dt + d2T_dt2*Psi;

    /* Effective source = - D2[Psi] */
    const auto f = (r[i]-2.0*M)/r[i];
    double V;
    if(parity == EvenParity) {
      const double mu = (l-1)*(l+2);
      const double Lambda = mu + 6.0*M/r[i];
      /* Martel-Poisson Eq. (4.26) times a factor of f */
      V = f*pow(Lambda,-2)*(pow(mu,2)*((mu+2)*pow(r[i],-2) + 6.0*M*pow(r[i],-3)) + 36.0*pow(M,2)*pow(r[i],-4)*(mu+2.0*M/r[i]));
    } else if (parity == OddParity) {
      /* Martel-Poisson Eq. (4.26) */
      V = f*(l*(l+1)*pow(r[i],-2) - 6.0*M*pow(r[i],-3));
    } else {
      V = 1.e300;
    }

    const complex<double> src = - d2TPsi_dt2  + 2.0*M*f/pow(r[i],2)*dTPsi_dr
      + pow(f,2)*d2TPsi_dr2 - V*TPsi;

    src_re[i] = src.real();
    src_im[i] = src.imag();
  }
}


complex<double> RWZ_EffectiveSource::operator()(const double r)
{
  /* Window function and its derivatives */
  double W = 1.0, dW_dr = 0.0, d2W_dr2 = 0.0;

  double sre, sim;
  this->operator()(1, &r, &W, &dW_dr, &d2W_dr2, &sre, &sim);
  return complex<double>(sre,sim);
}


complex<double> RWZ_EffectiveSource::Psi(const double r, const int sign)
{
  const complex<double> I(0.0, 1.0);
  const double md = m;

  /* Distance from the worldline */
  const double dr = r-r0;

  /* The puncture field is given by:
   * \Psi = D^l_{mm'} e^{i m [- \phi_0(t) + c(t) (r - r_0(t))]}
   *        \sum_i a_i(t) [r-r0(t)]^i + b_i(t) |r-r_0(t)| [r-r_0(t)]^i
   *
   * First, compute the term without the phase factor
   * \Psi = D^l_{mm'} \sum_i a_i(t) [r-r0(t)]^i + b_i(t) |r-r_0(t)| [r-r_0(t)]^i
   */
  const double Psi_l0 =
    coeffs[0][0] + dr*(coeffs[0][1] + dr*(coeffs[0][2] + dr*coeffs[0][3]))
    + sign*(abscoeffs[0][0] + dr*(abscoeffs[0][1] + dr*(abscoeffs[0][2] + dr*abscoeffs[0][3])));
  const double Psi_l1 =
    coeffs[1][0] + dr*(coeffs[1][1] + dr*(coeffs[1][2] + dr*coeffs[1][3]))
    + sign*(abscoeffs[1][0] + dr*(abscoeffs[1][1] + dr*(abscoeffs[1][2] + dr*abscoeffs[1][3])));
  const double Psi_l2 =
    coeffs[2][0] + dr*(coeffs[2][1] + dr*(coeffs[2][2] + dr*coeffs[2][3]))
    + sign*(abscoeffs[2][0] + dr*(abscoeffs[2][1] + dr*(abscoeffs[2][2] + dr*abscoeffs[2][3])));
  const double Psi_l3 = l < 3 ? 0 :
    coeffs[3][0] + dr*(coeffs[3][1] + dr*(coeffs[3][2] + dr*coeffs[3][3]))
    + sign*(abscoeffs[3][0] + dr*(abscoeffs[3][1] + dr*(abscoeffs[3][2] + dr*abscoeffs[3][3])));

  complex<double> Psi;
  if (parity == EvenParity) {
    Psi = w0 * Psi_l0 + (w1p - w1m) * I * Psi_l1 + (w2p + w2m) * Psi_l2 + (w3p - w3m) * I * Psi_l3;
  } else if (parity == OddParity) {
    Psi = (w1p + w1m) * I * Psi_l1 + (w2p - w2m) * Psi_l2 + (w3p + w3m) * I * Psi_l3;
  }

  /* Now add in phase - note the order here is important */
  complex<double> phase = exp(complex<double>(0,md*(-phi0 + c*dr)));
  Psi = phase*Psi;

  /* Add in factor of time window function */
  const auto TPsi = T*Psi;

  return TPsi;
}


complex<double> RWZ_EffectiveSource::dPsi_dr(const double r, const int sign)
{
  const std::complex<double> I(0.0, 1.0);
  const double md = m;
  /* Distance from the worldline */
  const double dr = r-r0;

  /* The puncture field is given by:
   * \Psi = D^l_{mm'} e^{i m [- \phi_0(t) + c(t) (r - r_0(t))]}
   *        \sum_i a_i(t) [r-r0(t)]^i + b_i(t) |r-r_0(t)| [r-r_0(t)]^i
   *
   * First, compute the term without the phase factor
   * \Psi = D^l_{mm'} \sum_i a_i(t) [r-r0(t)]^i + b_i(t) |r-r_0(t)| [r-r_0(t)]^i
   */
  const double Psi_l0 =
    coeffs[0][0] + dr*(coeffs[0][1] + dr*(coeffs[0][2] + dr*coeffs[0][3]))
    + sign*(abscoeffs[0][0] + dr*(abscoeffs[0][1] + dr*(abscoeffs[0][2] + dr*abscoeffs[0][3])));
  const double Psi_l1 =
    coeffs[1][0] + dr*(coeffs[1][1] + dr*(coeffs[1][2] + dr*coeffs[1][3]))
    + sign*(abscoeffs[1][0] + dr*(abscoeffs[1][1] + dr*(abscoeffs[1][2] + dr*abscoeffs[1][3])));
  const double Psi_l2 =
    coeffs[2][0] + dr*(coeffs[2][1] + dr*(coeffs[2][2] + dr*coeffs[2][3]))
    + sign*(abscoeffs[2][0] + dr*(abscoeffs[2][1] + dr*(abscoeffs[2][2] + dr*abscoeffs[2][3])));
  const double Psi_l3 = l < 3 ? 0 :
    coeffs[3][0] + dr*(coeffs[3][1] + dr*(coeffs[3][2] + dr*coeffs[3][3]))
    + sign*(abscoeffs[3][0] + dr*(abscoeffs[3][1] + dr*(abscoeffs[3][2] + dr*abscoeffs[3][3])));

  /* Radial derivative */
  const double dPsi_dr_l0 =
    coeffs[0][1] + dr*(2.0*coeffs[0][2] + dr*3.0*coeffs[0][3])
    + sign*(abscoeffs[0][1] + dr*(2.0*abscoeffs[0][2] + dr*3.0*abscoeffs[0][3]));
  const double dPsi_dr_l1 =
    coeffs[1][1] + dr*(2.0*coeffs[1][2] + dr*3.0*coeffs[1][3])
    + sign*(abscoeffs[1][1] + dr*(2.0*abscoeffs[1][2] + dr*3.0*abscoeffs[1][3]));
  const double dPsi_dr_l2 =
    coeffs[2][1] + dr*(2.0*coeffs[2][2] + dr*3.0*coeffs[2][3])
    + sign*(abscoeffs[2][1] + dr*(2.0*abscoeffs[2][2] + dr*3.0*abscoeffs[2][3]));
  const double dPsi_dr_l3 = l < 3 ? 0 :
    coeffs[3][1] + dr*(2.0*coeffs[3][2] + dr*3.0*coeffs[3][3])
    + sign*(abscoeffs[3][1] + dr*(2.0*abscoeffs[3][2] + dr*3.0*abscoeffs[3][3]));

  complex<double> Psi, dPsi_dr;
  if (parity == EvenParity) {
    Psi = w0 * Psi_l0 + (w1p - w1m) * I * Psi_l1 + (w2p + w2m) * Psi_l2 + (w3p - w3m) * I * Psi_l3;
    dPsi_dr = w0 * dPsi_dr_l0 + (w1p - w1m) * I * dPsi_dr_l1 + (w2p + w2m) * dPsi_dr_l2 + (w3p - w3m) * I * dPsi_dr_l3;
  } else if (parity == OddParity) {
    Psi = (w1p + w1m) * I * Psi_l1 + (w2p - w2m) * Psi_l2 + (w3p + w3m) * I * Psi_l3;
    dPsi_dr = (w1p + w1m) * I * dPsi_dr_l1 + (w2p - w2m) * dPsi_dr_l2 + (w3p + w3m) * I * dPsi_dr_l3;
  }

  /* Now add in phase - note the order here is important */
  complex<double> phase = exp(complex<double>(0,md*(-phi0 + c*dr)));
  dPsi_dr = phase*(dPsi_dr + I*c*md*Psi);

  /* Account for time window */
  const auto dTPsi_dr = dPsi_dr*T;

  return dTPsi_dr;
}


complex<double> RWZ_EffectiveSource::dPsi_dt(const double r, const int sign)
{
  const std::complex<double> I(0.0, 1.0);
  const double md = m;

  /* Distance from the worldline */
  const double dr = r-r0;

  /* The puncture field is given by:
   * \Psi = D^l_{mm'} e^{i m [- \phi_0(t) + c(t) (r - r_0(t))]}
   *        \sum_i a_i(t) [r-r0(t)]^i + b_i(t) |r-r_0(t)| [r-r_0(t)]^i
   *
   * First, compute the term without the phase factor
   * \Psi = D^l_{mm'} \sum_i a_i(t) [r-r0(t)]^i + b_i(t) |r-r_0(t)| [r-r_0(t)]^i
   */
  const double Psi_l0 =
    coeffs[0][0] + dr*(coeffs[0][1] + dr*(coeffs[0][2] + dr*coeffs[0][3]))
    + sign*(abscoeffs[0][0] + dr*(abscoeffs[0][1] + dr*(abscoeffs[0][2] + dr*abscoeffs[0][3])));
  const double Psi_l1 =
    coeffs[1][0] + dr*(coeffs[1][1] + dr*(coeffs[1][2] + dr*coeffs[1][3]))
    + sign*(abscoeffs[1][0] + dr*(abscoeffs[1][1] + dr*(abscoeffs[1][2] + dr*abscoeffs[1][3])));
  const double Psi_l2 =
    coeffs[2][0] + dr*(coeffs[2][1] + dr*(coeffs[2][2] + dr*coeffs[2][3]))
    + sign*(abscoeffs[2][0] + dr*(abscoeffs[2][1] + dr*(abscoeffs[2][2] + dr*abscoeffs[2][3])));
  const double Psi_l3 = l < 3 ? 0 :
    coeffs[3][0] + dr*(coeffs[3][1] + dr*(coeffs[3][2] + dr*coeffs[3][3]))
    + sign*(abscoeffs[3][0] + dr*(abscoeffs[3][1] + dr*(abscoeffs[3][2] + dr*abscoeffs[3][3])));

  /* Radial derivative */
  const double dPsi_dr_l0 =
    coeffs[0][1] + dr*(2.0*coeffs[0][2] + dr*3.0*coeffs[0][3])
    + sign*(abscoeffs[0][1] + dr*(2.0*abscoeffs[0][2] + dr*3.0*abscoeffs[0][3]));
  const double dPsi_dr_l1 =
    coeffs[1][1] + dr*(2.0*coeffs[1][2] + dr*3.0*coeffs[1][3])
    + sign*(abscoeffs[1][1] + dr*(2.0*abscoeffs[1][2] + dr*3.0*abscoeffs[1][3]));
  const double dPsi_dr_l2 =
    coeffs[2][1] + dr*(2.0*coeffs[2][2] + dr*3.0*coeffs[2][3])
    + sign*(abscoeffs[2][1] + dr*(2.0*abscoeffs[2][2] + dr*3.0*abscoeffs[2][3]));
  const double dPsi_dr_l3 = l < 3 ? 0 :
    coeffs[3][1] + dr*(2.0*coeffs[3][2] + dr*3.0*coeffs[3][3])
    + sign*(abscoeffs[3][1] + dr*(2.0*abscoeffs[3][2] + dr*3.0*abscoeffs[3][3]));

  /* The time derivative hits the coefficients and powers of dr = r - r_0(t).
   * First compute the derivative of the coefficients.
   */
  const double dPsi_dt_l0 =
    dcoeffs_dt[0][0] + dr*(dcoeffs_dt[0][1] + dr*(dcoeffs_dt[0][2] + dr*dcoeffs_dt[0][3]))
    + sign*(dabscoeffs_dt[0][0] + dr*(dabscoeffs_dt[0][1] + dr*(dabscoeffs_dt[0][2] + dr*dabscoeffs_dt[0][3])));
  const double dPsi_dt_l1 =
    dcoeffs_dt[1][0] + dr*(dcoeffs_dt[1][1] + dr*(dcoeffs_dt[1][2] + dr*dcoeffs_dt[1][3]))
    + sign*(dabscoeffs_dt[1][0] + dr*(dabscoeffs_dt[1][1] + dr*(dabscoeffs_dt[1][2] + dr*dabscoeffs_dt[1][3])));
  const double dPsi_dt_l2 =
    dcoeffs_dt[2][0] + dr*(dcoeffs_dt[2][1] + dr*(dcoeffs_dt[2][2] + dr*dcoeffs_dt[2][3]))
    + sign*(dabscoeffs_dt[2][0] + dr*(dabscoeffs_dt[2][1] + dr*(dabscoeffs_dt[2][2] + dr*dabscoeffs_dt[2][3])));
  const double dPsi_dt_l3 = l < 3 ? 0 :
    dcoeffs_dt[3][0] + dr*(dcoeffs_dt[3][1] + dr*(dcoeffs_dt[3][2] + dr*dcoeffs_dt[3][3]))
    + sign*(dabscoeffs_dt[3][0] + dr*(dabscoeffs_dt[3][1] + dr*(dabscoeffs_dt[3][2] + dr*dabscoeffs_dt[3][3])));

  complex<double> Psi, dPsi_dr,dPsi_dt;
  if (parity == EvenParity) {
    Psi = w0 * Psi_l0 + (w1p - w1m) * I * Psi_l1 + (w2p + w2m) * Psi_l2 + (w3p - w3m) * I * Psi_l3;
    dPsi_dr = w0 * dPsi_dr_l0 + (w1p - w1m) * I * dPsi_dr_l1 + (w2p + w2m) * dPsi_dr_l2 + (w3p - w3m) * I * dPsi_dr_l3;
    dPsi_dt = w0 * dPsi_dt_l0 + (w1p - w1m) * I * dPsi_dt_l1 + (w2p + w2m) * dPsi_dt_l2 + (w3p - w3m) * I * dPsi_dt_l3;
  } else if (parity == OddParity) {
    Psi = (w1p + w1m) * I * Psi_l1 + (w2p - w2m) * Psi_l2 + (w3p + w3m) * I * Psi_l3;
    dPsi_dr = (w1p + w1m) * I * dPsi_dr_l1 + (w2p - w2m) * dPsi_dr_l2 + (w3p + w3m) * I * dPsi_dr_l3;
    dPsi_dt = (w1p + w1m) * I * dPsi_dt_l1 + (w2p - w2m) * dPsi_dt_l2 + (w3p + w3m) * I * dPsi_dt_l3;
  }

  /* Add in derivative of time dependence of dr = r - r_0(t) */
  dPsi_dt += -dPsi_dr*rt;

  /* Now add in phase - note the order here is important */
  complex<double> phase = exp(complex<double>(0,md*(-phi0 + c*dr)));
  dPsi_dt = phase*(I*md*Psi*(ct*dr - phit - c*rt) + dPsi_dt);

  /* Add factor of time window */
  const auto dTPsi_dt = T*dPsi_dt + dT_dt*Psi;

  return dTPsi_dt;
}

int RWZ_EffectiveSource::get_l()
{
  return l;

}
