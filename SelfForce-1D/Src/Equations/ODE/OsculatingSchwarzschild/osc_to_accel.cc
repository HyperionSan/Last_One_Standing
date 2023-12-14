#include <math.h>
#include <stdio.h>

extern "C" {

double compute_w(double alpha, double beta)
{
  double w = atan2(alpha, beta);
  return w;
}

double compute_e(double alpha, double beta)
{
  double e = sqrt(alpha*alpha + beta*beta);
  return e;
}

double compute_v(double chi, double w)
{
  double v = chi - w;
  return v;
}

double compute_dwdt(double alpha, double beta, double dalphadt, double dbetadt)
{
  double dwdt = (beta*dalphadt - alpha*dbetadt)/(alpha*alpha + beta*beta);
  return dwdt;
}

double compute_dedt(double alpha, double beta, double dalphadt, double dbetadt)
{
  double dedt = (alpha*dalphadt + beta*dbetadt)/sqrt(alpha*alpha + beta*beta);
  return dedt;
}

double compute_d2wdt2(double alpha, double beta, double dalphadt, double dbetadt, double d2alphadt2, double d2betadt2)
{
  double d2wdt2 = (-(pow(alpha,3)*d2betadt2) + pow(beta,2)*(beta*d2alphadt2 - 2*dalphadt*dbetadt) + pow(alpha,2)*(beta*d2alphadt2 + 2*dalphadt*dbetadt) - alpha*beta*(beta*d2betadt2 + 2*pow(dalphadt,2) - 2*pow(dbetadt,2)))/pow(pow(alpha,2) + pow(beta,2),2);
  return d2wdt2;
}

double compute_d2edt2(double alpha, double beta, double dalphadt, double dbetadt, double d2alphadt2, double d2betadt2)
{
  double d2edt2 = (pow(beta,3)*d2betadt2 + pow(beta,2)*(alpha*d2alphadt2 + pow(dalphadt,2)) + alpha*beta*(alpha*d2betadt2 - 2*dalphadt*dbetadt) + pow(alpha,2)*(alpha*d2alphadt2 + pow(dbetadt,2)))/pow(pow(alpha,2) + pow(beta,2),1.5);
  return d2edt2;
}

double compute_d3wdt3(double alpha, double beta, double dalphadt, double dbetadt, double d2alphadt2, double d2betadt2, double d3alphadt3, double d3betadt3)
{
  double d3wdt3 = (-(pow(alpha,5)*d3betadt3) + pow(alpha,4)*(beta*d3alphadt3 + 3*d2betadt2*dalphadt + 3*d2alphadt2*dbetadt) + 2*pow(alpha,2)*beta*(pow(beta,2)*d3alphadt3 + 3*pow(dalphadt,3) - 9*dalphadt*pow(dbetadt,2)) - 2*pow(alpha,3)*(pow(beta,2)*d3betadt3 + 3*beta*d2alphadt2*dalphadt - 3*beta*d2betadt2*dbetadt + 3*pow(dalphadt,2)*dbetadt - pow(dbetadt,3)) - alpha*pow(beta,2)*(pow(beta,2)*d3betadt3 + 6*beta*d2alphadt2*dalphadt - 6*beta*d2betadt2*dbetadt - 18*pow(dalphadt,2)*dbetadt + 6*pow(dbetadt,3)) + pow(beta,3)*(-2*pow(dalphadt,3) + beta*(beta*d3alphadt3 - 3*d2alphadt2*dbetadt) + dalphadt*(-3*beta*d2betadt2 + 6*pow(dbetadt,2))))/ pow(pow(alpha,2) + pow(beta,2),3);
  return d3wdt3;
}

double compute_d3edt3(double alpha, double beta, double dalphadt, double dbetadt, double d2alphadt2, double d2betadt2, double d3alphadt3, double d3betadt3)
{
  double d3edt3 = (3*pow(alpha*dalphadt + beta*dbetadt,3) + pow(pow(alpha,2) + pow(beta,2),2)* (alpha*d3alphadt3 + beta*d3betadt3 + 3*d2alphadt2*dalphadt + 3*d2betadt2*dbetadt) - 3*(pow(alpha,2) + pow(beta,2))*(alpha*dalphadt + beta*dbetadt)* (alpha*d2alphadt2 + beta*d2betadt2 + pow(dalphadt,2) + pow(dbetadt,2)))/pow(pow(alpha,2) + pow(beta,2),2.5);
  return d3edt3;
}
          
double compute_dardt(double d2edt2, double d2wdt2, double d2chidt2, double dpdt, double dedt, double dwdt, double dchidt, double e, double p, double v, double M)
{
//  printf("d2edt2 = %20.15e, d2wdt2 = %20.15e, d2chidt2 = %20.15e, dpdt = %20.15e, dedt = %20.15e, dwdt = %20.15e, dchidt = %20.15e, e = %20.15e, p = %20.15e, v = %20.15e, mass = %20.15e\n", d2edt2, d2wdt2, d2chidt2, dpdt, dedt, dwdt, dchidt, e, p, v, M);
  double dardt = (dchidt*(3 + pow(e,2) - p)*p*(1 + e*cos(v))*(-2*dwdt*pow(e,2)*(4*dedt*e + dpdt*(-5 + p)) - 2*dedt*dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) - d2wdt2*pow(e,2)*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 2*dpdt*dwdt*e*(6 + 2*pow(e,2) - p)*cos(v) + 2*dwdt*e*(-dpdt + 4*dedt*e)*(-3 + p)*cos(v) + 2*dedt*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p)*cos(v) + 2*d2wdt2*e*(6 + 2*pow(e,2) - p)*(-3 + p)*cos(v) + dpdt*dwdt*pow(e,2)*(6 + 2*pow(e,2) - p)*pow(cos(v),2) + dwdt*pow(e,2)*(-dpdt + 4*dedt*e)*(-6 + p)*pow(cos(v),2) + 2*dedt*dwdt*e*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(cos(v),2) + d2wdt2*pow(e,2)*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(cos(v),2) - 2*dwdt*pow(e,3)*(-dpdt + 4*dedt*e)*pow(cos(v),3) - 6*dedt*dwdt*pow(e,2)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) - 2*d2wdt2*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) - dedt*(-dchidt + dwdt)*cos(v)*(2*(18 - 9*p + pow(p,2)) + e*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - 2*pow(e,2)*(-6 + p)*pow(cos(v),2)) + 12*dedt*(-dchidt + dwdt)*e*p*cos(2*v) + 2*dwdt*(-dchidt + dwdt)*e*(6 + 2*pow(e,2) - p)*(-3 + p)*sin(v) - 6*dwdt*(-dchidt + dwdt)*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),2)*sin(v) + d2edt2*(2*(18 - 9*p + pow(p,2)) + e*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - 2*pow(e,2)*(-6 + p)*pow(cos(v),2))*sin(v) - 6*dedt*dpdt*e*sin(2*v) + dwdt*(-dchidt + dwdt)*pow(e,2)*(6 + 2*pow(e,2) - p)*(-6 + p)*sin(2*v) - 6*pow(dedt,2)*p*sin(2*v) - 6*d2edt2*e*p*sin(2*v) + dedt*sin(v)*(2*dpdt*(-9 + 2*p) + (2*dpdt*e*p + dedt*(36 - 12*pow(e,2) + pow(p,2)))*cos(v) - 2*e*(dpdt*e + 2*dedt*(-6 + p))*pow(cos(v),2) - (-dchidt + dwdt)*e*(-36 + 4*pow(e,2) - pow(p,2))*sin(v) - 2*(-dchidt + dwdt)*pow(e,2)*(-6 + p)*sin(2*v))) - 2*dchidt*dpdt*(3 + pow(e,2) - p)*(1 + e*cos(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))) - 2*dchidt*(-dpdt + 2*dedt*e)*p*(1 + e*cos(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))) - d2chidt2*(3 + pow(e,2) - p)*p*(1 + e*cos(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))) + dchidt*(3 + pow(e,2) - p)*p*(dedt*cos(v) + (-dchidt + dwdt)*e*sin(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))))/(2.*pow(dchidt,2)*M*pow(3 + pow(e,2) - p,3)*pow(p,3));
//  printf("dardt = %20.15e\n", dardt);
//  fflush(stdout);
  return dardt;
}

double compute_daphidt(double d2edt2, double d2wdt2, double d2chidt2, double dpdt, double dedt, double dwdt, double dchidt, double e, double p, double v, double M)
{
  double daphidt = (pow(1 + e*cos(v),2)*(-5*dchidt*dpdt*(3 + pow(e,2) - p)*(-6 + p - 2*e*cos(v))*(1 + e*cos(v))*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) - 4*dchidt*(-dpdt + 2*dedt*e)*p*(-6 + p - 2*e*cos(v))*(1 + e*cos(v))*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) - 2*d2chidt2*(3 + pow(e,2) - p)*p*(-6 + p - 2*e*cos(v))*(1 + e*cos(v))*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) + dchidt*(3 + pow(e,2) - p)*p*(1 + e*cos(v))*(dpdt - 2*dedt*cos(v) - 2*(-dchidt + dwdt)*e*sin(v))*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) + 6*dchidt*(3 + pow(e,2) - p)*p*(-6 + p - 2*e*cos(v))*(dedt*cos(v) + (-dchidt + dwdt)*e*sin(v))*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) + 2*dchidt*(3 + pow(e,2) - p)*p*(-6 + p - 2*e*cos(v))*(1 + e*cos(v))*(dwdt*(-dchidt + dwdt)*e*(6 + 2*pow(e,2) - p)*cos(v) + d2edt2*(2*e + (-6 + p)*cos(v)) + dwdt*e*(dpdt - 4*dedt*e)*sin(v) + dedt*dwdt*(-6 - 2*pow(e,2) + p)*sin(v) + d2wdt2*e*(-6 - 2*pow(e,2) + p)*sin(v) + dedt*(2*dedt + dpdt*cos(v) + (-dchidt + dwdt)*(-6 + p)*sin(v)))))/(4.*pow(dchidt,2)*pow(M,2)*pow(3 + pow(e,2) - p,3)*pow(p,3.5)*sqrt(-6 + p - 2*e*cos(v)));
  return daphidt;
}

double compute_d2ardt2(double d3edt3, double d3wdt3, double d3chidt3, double d2pdt2, double d2edt2, double d2wdt2, double d2chidt2, double dpdt, double dedt, double dwdt, double dchidt, double e, double p, double v, double M)
{
  double d2ardt2 = (pow(dchidt,2)*pow(3 + pow(e,2) - p,2)*pow(p,2)*(1 + e*cos(v))*(-2*dwdt*pow(e,2)*(4*pow(dedt,2) + pow(dpdt,2) + 4*d2edt2*e + d2pdt2*(-5 + p)) - 8*dedt*dwdt*e*(4*dedt*e + dpdt*(-5 + p)) - 4*d2wdt2*pow(e,2)*(4*dedt*e + dpdt*(-5 + p)) - 2*pow(dedt,2)*dwdt*(12 + 4*pow(e,2) - 10*p + pow(p,2)) - 4*d2wdt2*dedt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) - 2*d2edt2*dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) - d3wdt3*pow(e,2)*(12 + 4*pow(e,2) - 10*p + pow(p,2)) - 4*dpdt*dwdt*e*(dpdt - 4*dedt*e)*cos(v) + 4*dedt*dpdt*dwdt*(6 + 2*pow(e,2) - p)*cos(v) + 4*d2wdt2*dpdt*e*(6 + 2*pow(e,2) - p)*cos(v) + 2*d2pdt2*dwdt*e*(6 + 2*pow(e,2) - p)*cos(v) + 4*dedt*dwdt*(-dpdt + 4*dedt*e)*(-3 + p)*cos(v) + 4*d2wdt2*e*(-dpdt + 4*dedt*e)*(-3 + p)*cos(v) + 2*dwdt*e*(-d2pdt2 + 4*(pow(dedt,2) + d2edt2*e))*(-3 + p)*cos(v) + 4*d2wdt2*dedt*(6 + 2*pow(e,2) - p)*(-3 + p)*cos(v) + 2*d2edt2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p)*cos(v) + 2*d3wdt3*e*(6 + 2*pow(e,2) - p)*(-3 + p)*cos(v) - 2*dwdt*pow(-dchidt + dwdt,2)*e*(6 + 2*pow(e,2) - p)*(-3 + p)*cos(v) - 2*dpdt*dwdt*pow(e,2)*(dpdt - 4*dedt*e)*pow(cos(v),2) + 4*dedt*dpdt*dwdt*e*(6 + 2*pow(e,2) - p)*pow(cos(v),2) + 2*d2wdt2*dpdt*pow(e,2)*(6 + 2*pow(e,2) - p)*pow(cos(v),2) + d2pdt2*dwdt*pow(e,2)*(6 + 2*pow(e,2) - p)*pow(cos(v),2) + 4*dedt*dwdt*e*(-dpdt + 4*dedt*e)*(-6 + p)*pow(cos(v),2) + 2*d2wdt2*pow(e,2)*(-dpdt + 4*dedt*e)*(-6 + p)*pow(cos(v),2) + dwdt*pow(e,2)*(-d2pdt2 + 4*(pow(dedt,2) + d2edt2*e))*(-6 + p)*pow(cos(v),2) + 2*pow(dedt,2)*dwdt*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(cos(v),2) + 4*d2wdt2*dedt*e*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(cos(v),2) + 2*d2edt2*dwdt*e*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(cos(v),2) + d3wdt3*pow(e,2)*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(cos(v),2) - 2*dwdt*pow(-dchidt + dwdt,2)*pow(e,2)*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(cos(v),2) - 12*dedt*dwdt*pow(e,2)*(-dpdt + 4*dedt*e)*pow(cos(v),3) - 4*d2wdt2*pow(e,3)*(-dpdt + 4*dedt*e)*pow(cos(v),3) - 2*dwdt*pow(e,3)*(-d2pdt2 + 4*(pow(dedt,2) + d2edt2*e))*pow(cos(v),3) - 12*pow(dedt,2)*dwdt*e*(6 + 2*pow(e,2) - p)*pow(cos(v),3) - 12*d2wdt2*dedt*pow(e,2)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) - 6*d2edt2*dwdt*pow(e,2)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) - 2*d3wdt3*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 6*dwdt*pow(-dchidt + dwdt,2)*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) - (-d2chidt2 + d2wdt2)*dedt*cos(v)*(2*(18 - 9*p + pow(p,2)) + e*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - 2*pow(e,2)*(-6 + p)*pow(cos(v),2)) - 2*d2edt2*(-dchidt + dwdt)*cos(v)*(2*(18 - 9*p + pow(p,2)) + e*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - 2*pow(e,2)*(-6 + p)*pow(cos(v),2)) + 24*dedt*dpdt*(-dchidt + dwdt)*e*cos(2*v) + 24*pow(dedt,2)*(-dchidt + dwdt)*p*cos(2*v) + 12*(-d2chidt2 + d2wdt2)*dedt*e*p*cos(2*v) + 24*d2edt2*(-dchidt + dwdt)*e*p*cos(2*v) + 4*dpdt*dwdt*(-dchidt + dwdt)*e*(6 + 2*pow(e,2) - p)*sin(v) + 4*dwdt*(-dchidt + dwdt)*e*(-dpdt + 4*dedt*e)*(-3 + p)*sin(v) + 4*dedt*dwdt*(-dchidt + dwdt)*(6 + 2*pow(e,2) - p)*(-3 + p)*sin(v) + 2*(-d2chidt2 + d2wdt2)*dwdt*e*(6 + 2*pow(e,2) - p)*(-3 + p)*sin(v) + 4*d2wdt2*(-dchidt + dwdt)*e*(6 + 2*pow(e,2) - p)*(-3 + p)*sin(v) - 12*dwdt*(-dchidt + dwdt)*pow(e,3)*(-dpdt + 4*dedt*e)*pow(cos(v),2)*sin(v) - 36*dedt*dwdt*(-dchidt + dwdt)*pow(e,2)*(6 + 2*pow(e,2) - p)*pow(cos(v),2)*sin(v) - 6*(-d2chidt2 + d2wdt2)*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),2)*sin(v) - 12*d2wdt2*(-dchidt + dwdt)*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),2)*sin(v) + d3edt3*(2*(18 - 9*p + pow(p,2)) + e*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - 2*pow(e,2)*(-6 + p)*pow(cos(v),2))*sin(v) - dedt*pow(-dchidt + dwdt,2)*(2*(18 - 9*p + pow(p,2)) + e*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - 2*pow(e,2)*(-6 + p)*pow(cos(v),2))*sin(v) + 2*dwdt*pow(-dchidt + dwdt,2)*pow(e,2)*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(sin(v),2) - 12*pow(dedt,2)*dpdt*sin(2*v) - 6*d2pdt2*dedt*e*sin(2*v) - 12*d2edt2*dpdt*e*sin(2*v) + 2*dpdt*dwdt*(-dchidt + dwdt)*pow(e,2)*(6 + 2*pow(e,2) - p)*sin(2*v) + 2*dwdt*(-dchidt + dwdt)*pow(e,2)*(-dpdt + 4*dedt*e)*(-6 + p)*sin(2*v) + 4*dedt*dwdt*(-dchidt + dwdt)*e*(6 + 2*pow(e,2) - p)*(-6 + p)*sin(2*v) + (-d2chidt2 + d2wdt2)*dwdt*pow(e,2)*(6 + 2*pow(e,2) - p)*(-6 + p)*sin(2*v) + 2*d2wdt2*(-dchidt + dwdt)*pow(e,2)*(6 + 2*pow(e,2) - p)*(-6 + p)*sin(2*v) - 18*d2edt2*dedt*p*sin(2*v) - 6*d3edt3*e*p*sin(2*v) + 24*dedt*pow(-dchidt + dwdt,2)*e*p*sin(2*v) - 6*dwdt*pow(-dchidt + dwdt,2)*pow(e,3)*(6 + 2*pow(e,2) - p)*sin(v)*sin(2*v) + dedt*sin(v)*(4*pow(dpdt,2) + 2*d2pdt2*(-9 + 2*p) + 2*e*(-4*pow(dedt,2) + pow(dpdt,2) - 4*d2edt2*e + d2pdt2*p)*cos(v) + 2*dedt*(-8*dedt*e + 2*dpdt*p)*cos(v) + d2edt2*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - pow(-dchidt + dwdt,2)*e*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - 8*dedt*dpdt*e*pow(cos(v),2) - 2*d2pdt2*pow(e,2)*pow(cos(v),2) - 4*pow(dedt,2)*(-6 + p)*pow(cos(v),2) - 4*d2edt2*e*(-6 + p)*pow(cos(v),2) + 4*pow(-dchidt + dwdt,2)*pow(e,2)*(-6 + p)*pow(cos(v),2) + 2*(-dchidt + dwdt)*e*(-8*dedt*e + 2*dpdt*p)*sin(v) + 2*dedt*(-dchidt + dwdt)*(36 - 4*pow(e,2) + pow(p,2))*sin(v) + (-d2chidt2 + d2wdt2)*e*(36 - 4*pow(e,2) + pow(p,2))*sin(v) - 4*pow(-dchidt + dwdt,2)*pow(e,2)*(-6 + p)*pow(sin(v),2) - 4*dpdt*(-dchidt + dwdt)*pow(e,2)*sin(2*v) - 8*dedt*(-dchidt + dwdt)*e*(-6 + p)*sin(2*v) - 2*(-d2chidt2 + d2wdt2)*pow(e,2)*(-6 + p)*sin(2*v)) - 2*dedt*(-dchidt + dwdt)*cos(v)*(2*dpdt*(-9 + 2*p) + (2*dpdt*e*p + dedt*(36 - 12*pow(e,2) + pow(p,2)))*cos(v) - 2*e*(dpdt*e + 2*dedt*(-6 + p))*pow(cos(v),2) - (-dchidt + dwdt)*e*(-36 + 4*pow(e,2) - pow(p,2))*sin(v) - 2*(-dchidt + dwdt)*pow(e,2)*(-6 + p)*sin(2*v)) + 2*d2edt2*sin(v)*(2*dpdt*(-9 + 2*p) + (2*dpdt*e*p + dedt*(36 - 12*pow(e,2) + pow(p,2)))*cos(v) - 2*e*(dpdt*e + 2*dedt*(-6 + p))*pow(cos(v),2) - (-dchidt + dwdt)*e*(-36 + 4*pow(e,2) - pow(p,2))*sin(v) - 2*(-dchidt + dwdt)*pow(e,2)*(-6 + p)*sin(2*v))) - 4*pow(dchidt,2)*dpdt*pow(3 + pow(e,2) - p,2)*p*(1 + e*cos(v))*(-2*dwdt*pow(e,2)*(4*dedt*e + dpdt*(-5 + p)) - 2*dedt*dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) - d2wdt2*pow(e,2)*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 2*dpdt*dwdt*e*(6 + 2*pow(e,2) - p)*cos(v) + 2*dwdt*e*(-dpdt + 4*dedt*e)*(-3 + p)*cos(v) + 2*dedt*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p)*cos(v) + 2*d2wdt2*e*(6 + 2*pow(e,2) - p)*(-3 + p)*cos(v) + dpdt*dwdt*pow(e,2)*(6 + 2*pow(e,2) - p)*pow(cos(v),2) + dwdt*pow(e,2)*(-dpdt + 4*dedt*e)*(-6 + p)*pow(cos(v),2) + 2*dedt*dwdt*e*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(cos(v),2) + d2wdt2*pow(e,2)*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(cos(v),2) - 2*dwdt*pow(e,3)*(-dpdt + 4*dedt*e)*pow(cos(v),3) - 6*dedt*dwdt*pow(e,2)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) - 2*d2wdt2*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) - dedt*(-dchidt + dwdt)*cos(v)*(2*(18 - 9*p + pow(p,2)) + e*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - 2*pow(e,2)*(-6 + p)*pow(cos(v),2)) + 12*dedt*(-dchidt + dwdt)*e*p*cos(2*v) + 2*dwdt*(-dchidt + dwdt)*e*(6 + 2*pow(e,2) - p)*(-3 + p)*sin(v) - 6*dwdt*(-dchidt + dwdt)*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),2)*sin(v) + d2edt2*(2*(18 - 9*p + pow(p,2)) + e*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - 2*pow(e,2)*(-6 + p)*pow(cos(v),2))*sin(v) - 6*dedt*dpdt*e*sin(2*v) + dwdt*(-dchidt + dwdt)*pow(e,2)*(6 + 2*pow(e,2) - p)*(-6 + p)*sin(2*v) - 6*pow(dedt,2)*p*sin(2*v) - 6*d2edt2*e*p*sin(2*v) + dedt*sin(v)*(2*dpdt*(-9 + 2*p) + (2*dpdt*e*p + dedt*(36 - 12*pow(e,2) + pow(p,2)))*cos(v) - 2*e*(dpdt*e + 2*dedt*(-6 + p))*pow(cos(v),2) - (-dchidt + dwdt)*e*(-36 + 4*pow(e,2) - pow(p,2))*sin(v) - 2*(-dchidt + dwdt)*pow(e,2)*(-6 + p)*sin(2*v))) - 4*pow(dchidt,2)*(-dpdt + 2*dedt*e)*(3 + pow(e,2) - p)*pow(p,2)*(1 + e*cos(v))*(-2*dwdt*pow(e,2)*(4*dedt*e + dpdt*(-5 + p)) - 2*dedt*dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) - d2wdt2*pow(e,2)*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 2*dpdt*dwdt*e*(6 + 2*pow(e,2) - p)*cos(v) + 2*dwdt*e*(-dpdt + 4*dedt*e)*(-3 + p)*cos(v) + 2*dedt*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p)*cos(v) + 2*d2wdt2*e*(6 + 2*pow(e,2) - p)*(-3 + p)*cos(v) + dpdt*dwdt*pow(e,2)*(6 + 2*pow(e,2) - p)*pow(cos(v),2) + dwdt*pow(e,2)*(-dpdt + 4*dedt*e)*(-6 + p)*pow(cos(v),2) + 2*dedt*dwdt*e*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(cos(v),2) + d2wdt2*pow(e,2)*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(cos(v),2) - 2*dwdt*pow(e,3)*(-dpdt + 4*dedt*e)*pow(cos(v),3) - 6*dedt*dwdt*pow(e,2)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) - 2*d2wdt2*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) - dedt*(-dchidt + dwdt)*cos(v)*(2*(18 - 9*p + pow(p,2)) + e*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - 2*pow(e,2)*(-6 + p)*pow(cos(v),2)) + 12*dedt*(-dchidt + dwdt)*e*p*cos(2*v) + 2*dwdt*(-dchidt + dwdt)*e*(6 + 2*pow(e,2) - p)*(-3 + p)*sin(v) - 6*dwdt*(-dchidt + dwdt)*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),2)*sin(v) + d2edt2*(2*(18 - 9*p + pow(p,2)) + e*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - 2*pow(e,2)*(-6 + p)*pow(cos(v),2))*sin(v) - 6*dedt*dpdt*e*sin(2*v) + dwdt*(-dchidt + dwdt)*pow(e,2)*(6 + 2*pow(e,2) - p)*(-6 + p)*sin(2*v) - 6*pow(dedt,2)*p*sin(2*v) - 6*d2edt2*e*p*sin(2*v) + dedt*sin(v)*(2*dpdt*(-9 + 2*p) + (2*dpdt*e*p + dedt*(36 - 12*pow(e,2) + pow(p,2)))*cos(v) - 2*e*(dpdt*e + 2*dedt*(-6 + p))*pow(cos(v),2) - (-dchidt + dwdt)*e*(-36 + 4*pow(e,2) - pow(p,2))*sin(v) - 2*(-dchidt + dwdt)*pow(e,2)*(-6 + p)*sin(2*v))) - 2*d2chidt2*dchidt*pow(3 + pow(e,2) - p,2)*pow(p,2)*(1 + e*cos(v))*(-2*dwdt*pow(e,2)*(4*dedt*e + dpdt*(-5 + p)) - 2*dedt*dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) - d2wdt2*pow(e,2)*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 2*dpdt*dwdt*e*(6 + 2*pow(e,2) - p)*cos(v) + 2*dwdt*e*(-dpdt + 4*dedt*e)*(-3 + p)*cos(v) + 2*dedt*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p)*cos(v) + 2*d2wdt2*e*(6 + 2*pow(e,2) - p)*(-3 + p)*cos(v) + dpdt*dwdt*pow(e,2)*(6 + 2*pow(e,2) - p)*pow(cos(v),2) + dwdt*pow(e,2)*(-dpdt + 4*dedt*e)*(-6 + p)*pow(cos(v),2) + 2*dedt*dwdt*e*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(cos(v),2) + d2wdt2*pow(e,2)*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(cos(v),2) - 2*dwdt*pow(e,3)*(-dpdt + 4*dedt*e)*pow(cos(v),3) - 6*dedt*dwdt*pow(e,2)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) - 2*d2wdt2*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) - dedt*(-dchidt + dwdt)*cos(v)*(2*(18 - 9*p + pow(p,2)) + e*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - 2*pow(e,2)*(-6 + p)*pow(cos(v),2)) + 12*dedt*(-dchidt + dwdt)*e*p*cos(2*v) + 2*dwdt*(-dchidt + dwdt)*e*(6 + 2*pow(e,2) - p)*(-3 + p)*sin(v) - 6*dwdt*(-dchidt + dwdt)*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),2)*sin(v) + d2edt2*(2*(18 - 9*p + pow(p,2)) + e*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - 2*pow(e,2)*(-6 + p)*pow(cos(v),2))*sin(v) - 6*dedt*dpdt*e*sin(2*v) + dwdt*(-dchidt + dwdt)*pow(e,2)*(6 + 2*pow(e,2) - p)*(-6 + p)*sin(2*v) - 6*pow(dedt,2)*p*sin(2*v) - 6*d2edt2*e*p*sin(2*v) + dedt*sin(v)*(2*dpdt*(-9 + 2*p) + (2*dpdt*e*p + dedt*(36 - 12*pow(e,2) + pow(p,2)))*cos(v) - 2*e*(dpdt*e + 2*dedt*(-6 + p))*pow(cos(v),2) - (-dchidt + dwdt)*e*(-36 + 4*pow(e,2) - pow(p,2))*sin(v) - 2*(-dchidt + dwdt)*pow(e,2)*(-6 + p)*sin(2*v))) + 2*pow(dchidt,2)*pow(3 + pow(e,2) - p,2)*pow(p,2)*(dedt*cos(v) + (-dchidt + dwdt)*e*sin(v))*(-2*dwdt*pow(e,2)*(4*dedt*e + dpdt*(-5 + p)) - 2*dedt*dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) - d2wdt2*pow(e,2)*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 2*dpdt*dwdt*e*(6 + 2*pow(e,2) - p)*cos(v) + 2*dwdt*e*(-dpdt + 4*dedt*e)*(-3 + p)*cos(v) + 2*dedt*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p)*cos(v) + 2*d2wdt2*e*(6 + 2*pow(e,2) - p)*(-3 + p)*cos(v) + dpdt*dwdt*pow(e,2)*(6 + 2*pow(e,2) - p)*pow(cos(v),2) + dwdt*pow(e,2)*(-dpdt + 4*dedt*e)*(-6 + p)*pow(cos(v),2) + 2*dedt*dwdt*e*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(cos(v),2) + d2wdt2*pow(e,2)*(6 + 2*pow(e,2) - p)*(-6 + p)*pow(cos(v),2) - 2*dwdt*pow(e,3)*(-dpdt + 4*dedt*e)*pow(cos(v),3) - 6*dedt*dwdt*pow(e,2)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) - 2*d2wdt2*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) - dedt*(-dchidt + dwdt)*cos(v)*(2*(18 - 9*p + pow(p,2)) + e*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - 2*pow(e,2)*(-6 + p)*pow(cos(v),2)) + 12*dedt*(-dchidt + dwdt)*e*p*cos(2*v) + 2*dwdt*(-dchidt + dwdt)*e*(6 + 2*pow(e,2) - p)*(-3 + p)*sin(v) - 6*dwdt*(-dchidt + dwdt)*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),2)*sin(v) + d2edt2*(2*(18 - 9*p + pow(p,2)) + e*(36 - 4*pow(e,2) + pow(p,2))*cos(v) - 2*pow(e,2)*(-6 + p)*pow(cos(v),2))*sin(v) - 6*dedt*dpdt*e*sin(2*v) + dwdt*(-dchidt + dwdt)*pow(e,2)*(6 + 2*pow(e,2) - p)*(-6 + p)*sin(2*v) - 6*pow(dedt,2)*p*sin(2*v) - 6*d2edt2*e*p*sin(2*v) + dedt*sin(v)*(2*dpdt*(-9 + 2*p) + (2*dpdt*e*p + dedt*(36 - 12*pow(e,2) + pow(p,2)))*cos(v) - 2*e*(dpdt*e + 2*dedt*(-6 + p))*pow(cos(v),2) - (-dchidt + dwdt)*e*(-36 + 4*pow(e,2) - pow(p,2))*sin(v) - 2*(-dchidt + dwdt)*pow(e,2)*(-6 + p)*sin(2*v))) + 6*pow(dchidt,2)*pow(dpdt,2)*pow(3 + pow(e,2) - p,2)*(1 + e*cos(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))) + 8*pow(dchidt,2)*dpdt*(-dpdt + 2*dedt*e)*(3 + pow(e,2) - p)*p*(1 + e*cos(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))) - 2*d2pdt2*pow(dchidt,2)*pow(3 + pow(e,2) - p,2)*p*(1 + e*cos(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))) + 4*d2chidt2*dchidt*dpdt*pow(3 + pow(e,2) - p,2)*p*(1 + e*cos(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))) + 6*pow(dchidt,2)*pow(dpdt - 2*dedt*e,2)*pow(p,2)*(1 + e*cos(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))) + 4*d2chidt2*dchidt*(-dpdt + 2*dedt*e)*(3 + pow(e,2) - p)*pow(p,2)*(1 + e*cos(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))) - 2*pow(dchidt,2)*(-d2pdt2 + 2*(pow(dedt,2) + d2edt2*e))*(3 + pow(e,2) - p)*pow(p,2)*(1 + e*cos(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))) + 2*pow(d2chidt2,2)*pow(3 + pow(e,2) - p,2)*pow(p,2)*(1 + e*cos(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))) - d3chidt3*dchidt*pow(3 + pow(e,2) - p,2)*pow(p,2)*(1 + e*cos(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))) - 4*pow(dchidt,2)*dpdt*pow(3 + pow(e,2) - p,2)*p*(dedt*cos(v) + (-dchidt + dwdt)*e*sin(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))) - 4*pow(dchidt,2)*(-dpdt + 2*dedt*e)*(3 + pow(e,2) - p)*pow(p,2)*(dedt*cos(v) + (-dchidt + dwdt)*e*sin(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))) - 2*d2chidt2*dchidt*pow(3 + pow(e,2) - p,2)*pow(p,2)*(dedt*cos(v) + (-dchidt + dwdt)*e*sin(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))) + pow(dchidt,2)*pow(3 + pow(e,2) - p,2)*pow(p,2)*((d2edt2 - pow(-dchidt + dwdt,2)*e)*cos(v) + (2*dedt*(-dchidt + dwdt) + (-d2chidt2 + d2wdt2)*e)*sin(v))*(-2*dwdt*pow(e,3)*(6 + 2*pow(e,2) - p)*pow(cos(v),3) + 2*dedt*(18 - 9*p + pow(p,2))*sin(v) + pow(e,2)*(-6 + p)*pow(cos(v),2)*(dwdt*(6 + 2*pow(e,2) - p) - 2*dedt*sin(v)) + e*cos(v)*(2*dwdt*(6 + 2*pow(e,2) - p)*(-3 + p) + dedt*(36 - 4*pow(e,2) + pow(p,2))*sin(v)) - e*(dwdt*e*(12 + 4*pow(e,2) - 10*p + pow(p,2)) + 6*dedt*p*sin(2*v))))/(2.*pow(dchidt,3)*M*pow(3 + pow(e,2) - p,4)*pow(p,4));
  return d2ardt2;
}

double compute_d2aphidt2(double d3edt3, double d3wdt3, double d3chidt3, double d2pdt2, double d2edt2, double d2wdt2, double d2chidt2, double dpdt, double dedt, double dwdt, double dchidt, double e, double p, double v, double M)
{
  double d2aphidt2 = ((1 + e*cos(v))*(4*pow(dchidt,2)*pow(3 + pow(e,2) - p,2)*pow(p,2)*pow(-6 + p - 2*e*cos(v),2)*pow(1 + e*cos(v),2)*(6*d2edt2*dedt + 2*d3edt3*e + d2pdt2*dedt*cos(v) + 2*d2edt2*dpdt*cos(v) - 2*dwdt*(-dchidt + dwdt)*e*(dpdt - 4*dedt*e)*cos(v) + 2*dedt*dwdt*(-dchidt + dwdt)*(6 + 2*pow(e,2) - p)*cos(v) + (-d2chidt2 + d2wdt2)*dwdt*e*(6 + 2*pow(e,2) - p)*cos(v) + 2*d2wdt2*(-dchidt + dwdt)*e*(6 + 2*pow(e,2) - p)*cos(v) + d3edt3*(-6 + p)*cos(v) - dedt*pow(-dchidt + dwdt,2)*(-6 + p)*cos(v) + 2*dedt*dpdt*(-dchidt + dwdt)*sin(v) + 2*dedt*dwdt*(dpdt - 4*dedt*e)*sin(v) + 2*d2wdt2*e*(dpdt - 4*dedt*e)*sin(v) + dwdt*e*(d2pdt2 - 4*(pow(dedt,2) + d2edt2*e))*sin(v) + dwdt*pow(-dchidt + dwdt,2)*e*(6 + 2*pow(e,2) - p)*sin(v) + (-d2chidt2 + d2wdt2)*dedt*(-6 + p)*sin(v) + 2*d2edt2*(-dchidt + dwdt)*(-6 + p)*sin(v) + 2*d2wdt2*dedt*(-6 - 2*pow(e,2) + p)*sin(v) + d2edt2*dwdt*(-6 - 2*pow(e,2) + p)*sin(v) + d3wdt3*e*(-6 - 2*pow(e,2) + p)*sin(v)) + 35*pow(dchidt,2)*pow(dpdt,2)*pow(3 + pow(e,2) - p,2)*pow(-6 + p - 2*e*cos(v),2)*pow(1 + e*cos(v),2)*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) + 40*pow(dchidt,2)*dpdt*(-dpdt + 2*dedt*e)*(3 + pow(e,2) - p)*p*pow(-6 + p - 2*e*cos(v),2)*pow(1 + e*cos(v),2)*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) - 10*d2pdt2*pow(dchidt,2)*pow(3 + pow(e,2) - p,2)*p*pow(-6 + p - 2*e*cos(v),2)*pow(1 + e*cos(v),2)*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) + 20*d2chidt2*dchidt*dpdt*pow(3 + pow(e,2) - p,2)*p*pow(-6 + p - 2*e*cos(v),2)*pow(1 + e*cos(v),2)*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) + 24*pow(dchidt,2)*pow(dpdt - 2*dedt*e,2)*pow(p,2)*pow(-6 + p - 2*e*cos(v),2)*pow(1 + e*cos(v),2)*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) + 16*d2chidt2*dchidt*(-dpdt + 2*dedt*e)*(3 + pow(e,2) - p)*pow(p,2)*pow(-6 + p - 2*e*cos(v),2)*pow(1 + e*cos(v),2)*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) - 8*pow(dchidt,2)*(-d2pdt2 + 2*(pow(dedt,2) + d2edt2*e))*(3 + pow(e,2) - p)*pow(p,2)*pow(-6 + p - 2*e*cos(v),2)*pow(1 + e*cos(v),2)*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) + 8*pow(d2chidt2,2)*pow(3 + pow(e,2) - p,2)*pow(p,2)*pow(-6 + p - 2*e*cos(v),2)*pow(1 + e*cos(v),2)*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) - 4*d3chidt3*dchidt*pow(3 + pow(e,2) - p,2)*pow(p,2)*pow(-6 + p - 2*e*cos(v),2)*pow(1 + e*cos(v),2)*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) - 10*pow(dchidt,2)*dpdt*pow(3 + pow(e,2) - p,2)*p*(-6 + p - 2*e*cos(v))*pow(1 + e*cos(v),2)*(dpdt - 2*dedt*cos(v) - 2*(-dchidt + dwdt)*e*sin(v))*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) - 8*pow(dchidt,2)*(-dpdt + 2*dedt*e)*(3 + pow(e,2) - p)*pow(p,2)*(-6 + p - 2*e*cos(v))*pow(1 + e*cos(v),2)*(dpdt - 2*dedt*cos(v) - 2*(-dchidt + dwdt)*e*sin(v))*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) - 4*d2chidt2*dchidt*pow(3 + pow(e,2) - p,2)*pow(p,2)*(-6 + p - 2*e*cos(v))*pow(1 + e*cos(v),2)*(dpdt - 2*dedt*cos(v) - 2*(-dchidt + dwdt)*e*sin(v))*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) - pow(dchidt,2)*pow(3 + pow(e,2) - p,2)*pow(p,2)*pow(1 + e*cos(v),2)*pow(dpdt - 2*dedt*cos(v) - 2*(-dchidt + dwdt)*e*sin(v),2)*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) - 60*pow(dchidt,2)*dpdt*pow(3 + pow(e,2) - p,2)*p*pow(-6 + p - 2*e*cos(v),2)*(1 + e*cos(v))*(dedt*cos(v) + (-dchidt + dwdt)*e*sin(v))*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) - 48*pow(dchidt,2)*(-dpdt + 2*dedt*e)*(3 + pow(e,2) - p)*pow(p,2)*pow(-6 + p - 2*e*cos(v),2)*(1 + e*cos(v))*(dedt*cos(v) + (-dchidt + dwdt)*e*sin(v))*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) - 24*d2chidt2*dchidt*pow(3 + pow(e,2) - p,2)*pow(p,2)*pow(-6 + p - 2*e*cos(v),2)*(1 + e*cos(v))*(dedt*cos(v) + (-dchidt + dwdt)*e*sin(v))*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) + 12*pow(dchidt,2)*pow(3 + pow(e,2) - p,2)*pow(p,2)*(-6 + p - 2*e*cos(v))*(1 + e*cos(v))*(dpdt - 2*dedt*cos(v) - 2*(-dchidt + dwdt)*e*sin(v))*(dedt*cos(v) + (-dchidt + dwdt)*e*sin(v))*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) + 24*pow(dchidt,2)*pow(3 + pow(e,2) - p,2)*pow(p,2)*pow(-6 + p - 2*e*cos(v),2)*pow(dedt*cos(v) + (-dchidt + dwdt)*e*sin(v),2)*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) + 2*pow(dchidt,2)*pow(3 + pow(e,2) - p,2)*pow(p,2)*(-6 + p - 2*e*cos(v))*pow(1 + e*cos(v),2)*(d2pdt2 - 2*(d2edt2 - pow(-dchidt + dwdt,2)*e)*cos(v) + (4*dchidt*dedt - 4*dedt*dwdt + 2*d2chidt2*e - 2*d2wdt2*e)*sin(v))*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) + 12*pow(dchidt,2)*pow(3 + pow(e,2) - p,2)*pow(p,2)*pow(-6 + p - 2*e*cos(v),2)*(1 + e*cos(v))*((d2edt2 - pow(-dchidt + dwdt,2)*e)*cos(v) + (2*dedt*(-dchidt + dwdt) + (-d2chidt2 + d2wdt2)*e)*sin(v))*(2*dedt*e + dedt*(-6 + p)*cos(v) + dwdt*e*(-6 - 2*pow(e,2) + p)*sin(v)) - 20*pow(dchidt,2)*dpdt*pow(3 + pow(e,2) - p,2)*p*pow(-6 + p - 2*e*cos(v),2)*pow(1 + e*cos(v),2)*(dwdt*(-dchidt + dwdt)*e*(6 + 2*pow(e,2) - p)*cos(v) + d2edt2*(2*e + (-6 + p)*cos(v)) + dwdt*e*(dpdt - 4*dedt*e)*sin(v) + dedt*dwdt*(-6 - 2*pow(e,2) + p)*sin(v) + d2wdt2*e*(-6 - 2*pow(e,2) + p)*sin(v) + dedt*(2*dedt + dpdt*cos(v) + (-dchidt + dwdt)*(-6 + p)*sin(v))) - 16*pow(dchidt,2)*(-dpdt + 2*dedt*e)*(3 + pow(e,2) - p)*pow(p,2)*pow(-6 + p - 2*e*cos(v),2)*pow(1 + e*cos(v),2)*(dwdt*(-dchidt + dwdt)*e*(6 + 2*pow(e,2) - p)*cos(v) + d2edt2*(2*e + (-6 + p)*cos(v)) + dwdt*e*(dpdt - 4*dedt*e)*sin(v) + dedt*dwdt*(-6 - 2*pow(e,2) + p)*sin(v) + d2wdt2*e*(-6 - 2*pow(e,2) + p)*sin(v) + dedt*(2*dedt + dpdt*cos(v) + (-dchidt + dwdt)*(-6 + p)*sin(v))) - 8*d2chidt2*dchidt*pow(3 + pow(e,2) - p,2)*pow(p,2)*pow(-6 + p - 2*e*cos(v),2)*pow(1 + e*cos(v),2)*(dwdt*(-dchidt + dwdt)*e*(6 + 2*pow(e,2) - p)*cos(v) + d2edt2*(2*e + (-6 + p)*cos(v)) + dwdt*e*(dpdt - 4*dedt*e)*sin(v) + dedt*dwdt*(-6 - 2*pow(e,2) + p)*sin(v) + d2wdt2*e*(-6 - 2*pow(e,2) + p)*sin(v) + dedt*(2*dedt + dpdt*cos(v) + (-dchidt + dwdt)*(-6 + p)*sin(v))) + 4*pow(dchidt,2)*pow(3 + pow(e,2) - p,2)*pow(p,2)*(-6 + p - 2*e*cos(v))*pow(1 + e*cos(v),2)*(dpdt - 2*dedt*cos(v) - 2*(-dchidt + dwdt)*e*sin(v))*(dwdt*(-dchidt + dwdt)*e*(6 + 2*pow(e,2) - p)*cos(v) + d2edt2*(2*e + (-6 + p)*cos(v)) + dwdt*e*(dpdt - 4*dedt*e)*sin(v) + dedt*dwdt*(-6 - 2*pow(e,2) + p)*sin(v) + d2wdt2*e*(-6 - 2*pow(e,2) + p)*sin(v) + dedt*(2*dedt + dpdt*cos(v) + (-dchidt + dwdt)*(-6 + p)*sin(v))) + 24*pow(dchidt,2)*pow(3 + pow(e,2) - p,2)*pow(p,2)*pow(-6 + p - 2*e*cos(v),2)*(1 + e*cos(v))*(dedt*cos(v) + (-dchidt + dwdt)*e*sin(v))*(dwdt*(-dchidt + dwdt)*e*(6 + 2*pow(e,2) - p)*cos(v) + d2edt2*(2*e + (-6 + p)*cos(v)) + dwdt*e*(dpdt - 4*dedt*e)*sin(v) + dedt*dwdt*(-6 - 2*pow(e,2) + p)*sin(v) + d2wdt2*e*(-6 - 2*pow(e,2) + p)*sin(v) + dedt*(2*dedt + dpdt*cos(v) + (-dchidt + dwdt)*(-6 + p)*sin(v)))))/(8.*pow(dchidt,3)*pow(M,2)*pow(3 + pow(e,2) - p,4)*pow(p,4.5)*pow(-6 + p - 2*e*cos(v),1.5));
  return d2aphidt2;
}

}
