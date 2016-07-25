#include <Rcpp.h>

using namespace Rcpp;

/*
  Author: Dan Taranu
 
  Semi-analytic integration for a Gaussian profile.
 
  I call this semi-analytic because there is an analytical solution for the integral
  of a 2D Gaussian over one dimension of a rectangular area, which is basically just 
  the product of an exponential and error function (not too surprising). Unfortunately 
  Wolfram Alpha wasn't able to integrate it over the second dimension. Perhaps a more
  clever person could find a useful normal integral as from here:
 
  http://www.tandfonline.com/doi/pdf/10.1080/03610918008812164

  xmod = xmid * INVREX + ymid * INVREY;
  ymod = (xmid * INVREY - ymid * INVREX)*INVAXRAT;
  return exp(-BN*(xmod*xmod + ymod * ymod -1));
  -> integral of
  exp(-BN*((x*INVREX + y*INVREY)^2 + (x*INVREY - y*INVREX)^2*INVAXRAT^2 -1));
  constants: BN = b, INVREX=c, INVREY=d, INVAXRAT^2 = a
  
  exp(-b*((x*c + y*d)^2 + a*(x*d - y*c)^2) - 1)
  
  Integrate[Exp[-(b ((x c + y d)^2 + a (x d - y c)^2)) - 1], x]
  = (Sqrt[Pi] Erf[(Sqrt[b] (c^2 x + a d^2 x - (-1 + a) c d y))/Sqrt[c^2 + a d^2]])/(2 Sqrt[b] Sqrt[c^2 + a d^2] E^((c^2 + a d^2 + a b c^4 y^2 + 2 a b c^2 d^2 y^2 + a b d^4 y^2)/(c^2 + a d^2)))
  = Erf((sqrt(b)*(c^2*x + a*d^2*x - (-1 + a)*c*d*y))/sqrt(c^2 + a*d^2))*exp((c^2 + a*d^2 + a*b*c^4*y^2 + 2*a*b*c^2*d^2*y^2 + a*b*d^4*y^2)/(c^2 + a*d^2)) dy * const
  = erf((sqb*(c2pad2tx - am1tcd*y))*isqc2pad2)*
    exp(-(c2pad2 + abc4*y2 + 2*abc2d2*y2 + abd4*y^2)/(c2pad2)) dy * const
  = erf((sqb*(c2pad2tx - am1tcd*y))*isqc2pad2)*
    exp(-1 - y^2*(abc4 + 2*abc2d2 + abd4)/c2pad2) dy * const
  
  where:
    
  const = Sqrt[Pi]/(2 Sqrt[b] Sqrt[c^2 + a d^2])
        = sqrt(pi)/(2*sqb)*isqc2pad2
        
  And:
  
  Integrate[Exp[-(b ((x c + y d)^2 + a (x d - y c)^2)) - 1], y]
  = (Sqrt[Pi] Erf[(Sqrt[b] (-((-1 + a) c d x) + a c^2 y + d^2 y))/Sqrt[a c^2 + d^2]])/(2 Sqrt[b] Sqrt[a c^2 + d^2] E^((d^2 + a (b c^4 x^2 + b d^4 x^2 + c^2 (1 + 2 b d^2 x^2)))/(a c^2 + d^2)))
  = Erf[(Sqrt[b] (-((-1 + a) c d x) + a c^2 y + d^2 y))/Sqrt[a c^2 + d^2]]*E^(-(d^2 + a (b c^4 x^2 + b d^4 x^2 + c^2 (1 + 2 b d^2 x^2)))/(a c^2 + d^2))*const
  = erf((sqb*(d2pac2ty - am1tcd*x))*isqd2pac2)*
    exp(-1 - x^2*(abc4 + 2*abc2d2 + abd4)/d2pac2) dx * const
  
  where:
  
  const = Sqrt[Pi]/(2 Sqrt[b] Sqrt[d^2 + a c^2])
        = sqrt(pi)/(2*sqb)*isqd2pac2
        
  (In other words, exactly the same but replacing c2pad2 with d2pac2 and isqc2pad2 with isqd2pac2)
 
 TODO: Integrate this into libprofit
*/

class Profit2DGaussianIntegrator{
  private:
    double a, b, c, d;
    double sqb, x2fac, y2fac, am1tcd, c2pad2, isqc2pad2, d2pac2, isqd2pac2;
    short i;
    double acc[2];
    double xintfac, yintfac;
    std::vector<double> rval;
    
    template <bool isx> inline double value(double x, double axfac1, double axfac2);
    
    template <bool x>
    inline double integral(double x1, double x2, double dx2, 
      double axfac1, double axfac2, double leftval=0,
      double rightval = 0);
    
  public:
    const std::vector<double> & integral(double x1, double x2, double y1, 
      double y2, double iacc=1e-5, double bottomval = 0, double topval = 0,
      double leftval = 0, double rightval = 0);
    
    // constants: BN = b, INVREX=c, INVREY=d, INVRAT = a
    Profit2DGaussianIntegrator(const double BN, const double INVREX, const double INVREY,
      const double INVAXRAT) : a(INVAXRAT*INVAXRAT),b(BN),c(INVREX),d(INVREY)
    {
      rval.resize(5);
      sqb = sqrt(b);
      double ab = a*b;
      double c2 = c*c;
      double d2 = d*d;
      double cd = c*d;
      am1tcd = (-1 + a)*cd;
      c2pad2 = c2+a*d2;
      d2pac2 = d2+a*c2;
      isqc2pad2 = 1.0/sqrt(c2pad2);
      isqd2pac2 = 1.0/sqrt(d2pac2);
      y2fac = (ab*pow(c,4) + 2*ab*c2*d2 + ab*pow(d,4));
      x2fac = y2fac/d2pac2;
      y2fac /= c2pad2;
      yintfac = sqrt(M_PI)*isqc2pad2/(2.*sqb);
      xintfac = sqrt(M_PI)*isqd2pac2/(2.*sqb);
    }
};

template <> inline double Profit2DGaussianIntegrator::value<false>(
  double x, double axfac1, double axfac2)
{
  double zi = am1tcd*x;
  return (erf(sqb*(axfac2 - zi)*isqc2pad2) - erf(sqb*(axfac1 - zi)*isqc2pad2))*
      exp(-1.0 - y2fac*x*x);
}

template<> inline double Profit2DGaussianIntegrator::value<true>(
  double x, double axfac1, double axfac2)
{
  double zi = am1tcd*x;
  return (erf(sqb*(axfac2 - zi)*isqd2pac2) - erf(sqb*(axfac1 - zi)*isqd2pac2))*
      exp(-1.0 - x2fac*x*x);
}

template <bool isx> inline double Profit2DGaussianIntegrator::integral(
  double x1, double x2, double dx2, double axfac1, double axfac2, 
  double leftval, double rightval)
{
  double xmid = x1+dx2;
  if(leftval == 0)
  {
    leftval = value<isx>(x1, axfac1, axfac2);
  }
  double cenval = value<isx>(xmid, axfac1, axfac2);
  if(rightval == 0)
  {
    rightval = value<isx>(x2, axfac1, axfac2);
  }
  double z1 = 0.5*(leftval+cenval)*dx2;
  if(z1 > acc[isx]) z1 = integral<isx>(x1,xmid,dx2/2.,axfac1,axfac2,
    leftval,cenval);
  double z2 = 0.5*(rightval+cenval)*dx2;
  if(z2 > acc[isx]) z2 = integral<isx>(xmid,x2,dx2/2.,axfac1,axfac2,
    cenval,rightval);
  //std::cout << y1 << "-" << y2 << "," << ymid << ":" << (z1+z2) << " " << dy2 << " " << leftval << " " << cenval << " " << rightval << std::endl;
  return(z1+z2);
}

const std::vector<double> & Profit2DGaussianIntegrator::integral(
  double x1, double x2, double y1, double y2, double iacc, 
  double bottomval, double topval, double leftval, double rightval)
{
  acc[true] = iacc*xintfac;
  acc[false] = iacc*yintfac;
  double d2pac2tx1 = x1*d2pac2;
  double d2pac2tx2 = x2*d2pac2;
  
  double c2pad2ty1 = y1*c2pad2;
  double c2pad2ty2 = y2*c2pad2;
  
  if(bottomval == 0)
  {
    bottomval = value<false>(y1, d2pac2tx1, d2pac2tx2);
    rval[1] = bottomval;
  }
  if(topval == 0)
  {
    topval = value<false>(y2, d2pac2tx1, d2pac2tx2);
    rval[2] = topval;
  }
  rval[0] = integral<false>(y1,y2,(y2-y1)/2.,d2pac2tx1,d2pac2tx2,bottomval,topval)*yintfac;
  if(leftval == 0)
  {
    leftval = value<true>(x1, c2pad2ty1, c2pad2ty2);
    rval[3] = leftval;
  }
  if(rightval == 0)
  {
    rightval = value<true>(x2, c2pad2ty1, c2pad2ty2);
    rval[4] = rightval;
  }
  rval[0] = (integral<true>(x1,x2,(x2-x1)/2.,c2pad2ty1,c2pad2ty2,leftval,rightval)*xintfac + rval[0])/2.;
  return rval;
}

NumericMatrix profitMakeGaussian(
    const double XCEN, const double YCEN, const double MAG, const double RE, 
    const double ANG, const double AXRAT, const double BOX, const double MAGZERO, 
    const NumericVector & XLIM, const NumericVector & YLIM, const IntegerVector & DIM, 
    const double ACC) {
  //const double BN=R::qgamma(0.5, 2 * NSER,1,1,0);  
  const double BN = 0.69314718055994528623;
  const double RBOX=PI*(BOX+2.)/(4.*R::beta(1./(BOX+2.),1+1./(BOX+2.)));
  const double LUMTOT = pow(10,(-0.4*(MAG-MAGZERO)));
  // TODO: Figure out why this empirical factor of 2*exp(1) is required to normalize properly
  const double Ie=2.*exp(1.)*LUMTOT/(RE*RE*AXRAT*PI*((exp(BN))/BN)*R::gammafn(1)/RBOX);
  const double INVRE = 1.0/RE;

  NumericMatrix mat(DIM(0), DIM(1));
  double x,y,xhi,yhi,angmod;
  double xbin=(XLIM(1)-XLIM(0))/DIM(0);
  double ybin=(YLIM(1)-YLIM(0))/DIM(1);

  angmod = std::fmod(ANG+90.,360.);
  if(angmod > 180.) angmod -= 180.;
  const double PX = cos(angmod*M_PI/180.);
  const double INVREY = sin(angmod*M_PI/180.)*INVRE*pow(-1,angmod < PI);
  const double INVREX = PX*INVRE;
  const double INVAXRAT = 1.0/AXRAT;
  Profit2DGaussianIntegrator gauss2(BN, INVREX, INVREY, INVAXRAT);
  // Want each pixels' integral to stop recursing only when less than 1e-3
  // But keep in mind we don't include the ie term, so the total over the
  // image won't be lumtot but lumtot/ie
  const double ACC2 = ACC*(LUMTOT/Ie);
  
  double bottomval=0;
  std::vector<double> leftvals(DIM(1));
  int i=0,j=0;
  x=XLIM(0)-XCEN; xhi = x+xbin;
  for(i = 0; i < DIM(0); i++) {
    y=YLIM(0)-YCEN; yhi = y+ybin;
    for(j = 0; j < DIM(1); j++) {
      const std::vector<double> & rval = gauss2.integral(x,xhi,y,yhi,ACC2,bottomval,0,leftvals[j]);
      bottomval = rval[2];
      leftvals[j] = rval[4];
      mat(i,j)=rval[0]*Ie;
      y = yhi;
      yhi += ybin;
      bottomval = 0;
    }
    x = xhi;
    xhi += xbin;
  }
  return(mat);
}
