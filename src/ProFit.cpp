#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

/*
  Some fast(er than pure R) functions for integrating a Sersic function
  and brute force convolution.
*/

/*
 The functions below MUST call this function first to get the value of
 n_ser to pass in, because in the non-boxy case, we want to pass in
 1/(2*nser), whereas in the boxiness it must be 1/nser.
 
 I've called the variable NSERFAC instead of INVNSER to emphasize this..
 */
inline double nserfac(double nser, double box)
{
  return 1.0/(nser*(1.0 + (box == 0)));
}

enum nsertype {general=0, exponent=1, two=2, three=3, four=4, gauss=5};

inline nsertype getnsertype(double nser)
{
  if(nser == -0.5) return gauss;
  if(nser == 1) return exponent;
  if(nser == 2) return two;
  if(nser == 3) return three;
  if(nser == 4) return four;
  return general;
}

inline double boxyre(const double & XMOD, const double & YMOD, const double & BOX)
{
  return pow(pow(std::abs(XMOD),2.+BOX)+pow(std::abs(YMOD),2.+BOX),1./(2.+BOX));
}

// The general case. Specific optimizations follow
template<bool hasbox, nsertype t>
inline double profitEvalSersic(const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  return exp(-BN*(pow(boxyre(XMOD,YMOD,BOX),NSERFAC)-1.));
}

template<> inline double profitEvalSersic<true, gauss>(
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  return exp(-BN*(boxyre(XMOD,YMOD,BOX)-1.));
}

template<> inline double profitEvalSersic<true, exponent>(
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  return exp(-BN*(sqrt(boxyre(XMOD,YMOD,BOX))-1.));
}

template<> inline double profitEvalSersic<true, two>(
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  return exp(-BN*(sqrt(sqrt(boxyre(XMOD,YMOD,BOX)))-1.));
}

template<> inline double profitEvalSersic<true, three>(
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  return exp(-BN*(cbrt(sqrt(boxyre(XMOD,YMOD,BOX)))-1.));
}

template<> inline double profitEvalSersic<true, four>(
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  return exp(-BN*(sqrt(sqrt(sqrt(boxyre(XMOD,YMOD,BOX))))-1.));
}

template<> inline double profitEvalSersic<false, general> (
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  double rdivresq = XMOD*XMOD + YMOD*YMOD;
  return exp(-BN*(pow(rdivresq,NSERFAC)-1.));
}

template<> inline double profitEvalSersic<false, gauss> (
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  double rdivresq = XMOD*XMOD + YMOD*YMOD;
  return exp(-BN*(rdivresq-1.));
}

template<> inline double profitEvalSersic<false, exponent> (
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  double rdivresq = sqrt(XMOD*XMOD + YMOD*YMOD);
  return exp(-BN*(rdivresq-1.));
}

template<> inline double profitEvalSersic<false, two> (
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  double rdivresq = sqrt(XMOD*XMOD + YMOD*YMOD);
  return exp(-BN*(sqrt(rdivresq)-1.));
}

template<> inline double profitEvalSersic<false, three> (
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  double rdivresq = sqrt(XMOD*XMOD + YMOD*YMOD);
  return exp(-BN*(cbrt(rdivresq)-1.));
}

template<> inline double profitEvalSersic<false, four> (
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  double rdivresq = sqrt(XMOD*XMOD + YMOD*YMOD);
  return exp(-BN*(sqrt(sqrt(rdivresq))-1.));
}

/*

Integral for an elliptical Gaussian via Wolfram Alpha:

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
/*
  Integrates the flux in a single pixel from a Sersic function, determining whether to recurse from
  an extra evaluation of the gradient along the minor axis (where it is maximized).
*/
template<bool hasbox, nsertype t>
double profitSumPixMinorAxisGrad(double XCEN, double YCEN, const NumericVector & XLIM, const NumericVector & YLIM,
                    const double INVREX, const double INVREY, const double INVAXRAT, const double INVNSER,
                    const double BOX, const double BN, const int UPSCALE=8L, const int RECURLEVEL=0L,
                    const int MAXDEPTH=2, const double ACC=0.1) {
  const bool RECURSE = (RECURLEVEL < MAXDEPTH) && (UPSCALE > 1);

  double x,y,xmid,ymid,xmod,ymod,testvaly;
  double sumpixel=0, addval=0;
  NumericVector xlim2(2),ylim2(2);
  
  // Reverse the order of integration to maintain symmetry
  const bool XREV = (XLIM(0) + (XLIM(1) - XLIM(0)/2.0)) > XCEN;
  const short XSIGN = XREV ? -1 : 1;
  const double XBIN = XSIGN*(XLIM(1)-XLIM(0))/double(UPSCALE);
  const double DX = XBIN/2.0;
  const double XINIT = XREV ? XLIM(1) : XLIM(0);
  const int XI0 = XREV ? 1 : 0;
  const int XI1 = 1 - XI0;
  
  const bool YREV = (YLIM(0) + (YLIM(1) - YLIM(0)/2.0)) > YCEN;
  const double YSIGN = YREV ? -1 : 1;
  const double YBIN = YSIGN * (YLIM(1)-YLIM(0))/double(UPSCALE);
  const double DY = YBIN/2.0;
  const double YINIT = YREV ? YLIM(1) : YLIM(0);
  const int YI0 = YREV ? 1 : 0;
  const int YI1 = 1 - YI0;
  //Rcout << XINIT << " " << XBIN << " " << YINIT << " " <<YBIN << std::endl;
  int i,j;
  x = XINIT;
  
  //int uprescount=0;
  
  for(i = 0; i < UPSCALE; ++i) {
    xmid = x + DX - XCEN;
    xlim2(XI0) = x;
    xlim2(XI1) = x + XBIN;
    y = YINIT;
    
    for(j = 0; j < UPSCALE; ++j) {
      ymid = y + DY - YCEN;
      // Project (xmid, ymid) onto the projection vector of the major axis,
      // pre-normalized to 1/Re to avoid extraneous division
      xmod = xmid * INVREX + ymid * INVREY;
      ymod = (xmid * INVREY - ymid * INVREX)*INVAXRAT;
      addval = profitEvalSersic<hasbox,t>(xmod, ymod, BN, BOX, INVNSER);
      //std::cout << "(" << xmid << "," << ymid << ") -> (" << xmod << "," << ymod << "): " << 
      //  addval << " for (" << BN << "," << INVNSER << "," << BOX << "), RECURSE=" << RECURSE << std::endl;
      if(RECURSE){
          testvaly = profitEvalSersic<hasbox,t>(xmod, std::abs(ymod)+std::abs(YBIN*INVREY*INVAXRAT), BN, BOX, INVNSER);
          if(std::abs(testvaly/addval - 1.0) > ACC){
            //Rcpp::Rcout << testvaly/addval<< " " << ACC << std::endl;
            ylim2(YI0) = y;
            ylim2(YI1) = y + YBIN;
            //uprescount++;
            addval=profitSumPixMinorAxisGrad<hasbox,t>(XCEN,YCEN,xlim2,ylim2,INVREX,INVREY,INVAXRAT,
              INVNSER, BOX, BN, UPSCALE, RECURLEVEL+1,MAXDEPTH,ACC);
          }
      }
      sumpixel+=addval;
      y += YBIN;
    }
    x += XBIN;
  }
  return(sumpixel/double(UPSCALE*UPSCALE));
}

/*
  Integrates the flux in a single pixel from a Sersic function, determining whether to recurse from
  the gradient relative to the nearest neighbour along the grid (saving function evaluations, 
  but potentially giving a less accurate answer than along the minor axis).
*/
template<bool hasbox, nsertype t>
double profitSumPixGridGrad(double XCEN, double YCEN, const NumericVector & XLIM, const NumericVector & YLIM,
                    const double INVREX, const double INVREY, const double INVAXRAT, const double NSERFAC,
                    const double BOX, const double BN, const int UPSCALE, const int RECURLEVEL=0, const int MAXDEPTH=3, 
                    const double ACC=0.1) {
  const bool RECURSE = (RECURLEVEL < MAXDEPTH) && (UPSCALE > 1);

  double x=0,y=0,xmid=0,ymid=0,xmod=0,ymod=0;
  double sumpixel=0, addval=0, oldaddval=1, olderaddval = 0;
  NumericVector xlim2(2),ylim2(2);
  
  // Reverse the order of integration to maintain symmetry
  const bool XREV = (XLIM(0) + (XLIM(1) - XLIM(0)/2.0)) > XCEN;
  const short XSIGN = XREV ? -1 : 1;
  const double XBIN = XSIGN*(XLIM(1)-XLIM(0))/double(UPSCALE);
  const double DX = XBIN/2.0;
  const double XINIT = XREV ? XLIM(1) : XLIM(0);
  const int XI0 = XREV ? 1 : 0;
  const int XI1 = 1 - XI0;
  
  const bool YREV = (YLIM(0) + (YLIM(1) - YLIM(0)/2.0)) > YCEN;
  const double YSIGN = YREV ? -1 : 1;
  const double YBIN = YSIGN * (YLIM(1)-YLIM(0))/double(UPSCALE);
  const double DY = YBIN/2.0;
  const double YINIT = YREV ? YLIM(1) : YLIM(0);
  const int YI0 = YREV ? 1 : 0;
  const int YI1 = 1 - YI0;

  // TODO: Separate this into another inline method
  if(UPSCALE > 1)
  {
    xmid = x + DX - XCEN;
    ymid = y + DY - YCEN - YBIN;
    xmod = xmid * INVREX + ymid * INVREY;
    ymod = (xmid * INVREY - ymid * INVREX)*INVAXRAT;
    oldaddval = profitEvalSersic<hasbox,t>(xmod, ymod, BN, BOX, NSERFAC);
  }
  
  int i,j;
  x = XINIT;
  
  for(i = 0; i < UPSCALE; ++i) {
    xmid = x + DX - XCEN;
    xlim2(XI0) = x;
    xlim2(XI1) = x + XBIN;
    y = YINIT;
    
    for(j = 0; j < UPSCALE; ++j) {
      ymid = y + DY - YCEN;
      // Project (xmid, ymid) onto the projection vector of the major axis,
      // pre-normalized to 1/Re to avoid extraneous division
      xmod = xmid * INVREX + ymid * INVREY;
      ymod = (xmid * INVREY - ymid * INVREX)*INVAXRAT;
      addval = profitEvalSersic<hasbox,t>(xmod, ymod, BN, BOX, NSERFAC);
      //std::cout << "(" << xmid << "," << ymid << ") -> (" << xmod << "," << ymod << "): " << 
      //  addval << " for (" << BN << "," << NSERFAC << "," << BOX << "), RECURSE=" << RECURSE << std::endl;
      if(RECURSE && (std::abs(addval/oldaddval - 1.0) > ACC)) {
        ylim2(YI0) = y;
        ylim2(YI1) = y + YBIN;
        addval=profitSumPixMinorAxisGrad<hasbox,t>(XCEN,YCEN,xlim2,ylim2,INVREX,INVREY,INVAXRAT, NSERFAC, BOX, BN,
          UPSCALE, RECURLEVEL+1,MAXDEPTH,ACC);
      }
      sumpixel+=addval;
      oldaddval=addval;
      if(j == 0) olderaddval = oldaddval;
      y += YBIN;
    }
    // Reset oldaddval to the value of the first pixel in the row
    // when jumping to the next row
    oldaddval = olderaddval;
    x += XBIN;
  }
  return(sumpixel/double(UPSCALE*UPSCALE));
}

/*
  This function computes the integrated flux of a Sersic profile on a grid, determining whether to recurse.
  TODO: Consider leaving the recursion logic entirely to profitSumPix.
*/
template<bool hasbox, nsertype t>
NumericMatrix profitMakeBoxySersic(const IntegerMatrix CALCREGION,
    const double XCEN=0, const double YCEN=0, const double MAG=15, 
    const double RE=1, const double NSER=1, const double ANG=0,
    const double AXRAT=1, const double BOX=0, 
    const double MAGZERO=0, const bool ROUGH=false,
    const NumericVector & XLIM = NumericVector::create(-100,100),
    const NumericVector & YLIM = NumericVector::create(-100,100),
    const IntegerVector & DIM = IntegerVector::create(200,200),
    const int UPSCALE=8L, const int MAXDEPTH=2L, const double RESWITCH=1,
    const double ACC=0.1, const bool DOCALCREGION=false, const double REMAX=10) {
  // Precompute things we only need to do once.
  const double BN=R::qgamma(0.5, 2 * NSER,1,1,0);  
  const double RBOX=PI*(BOX+2.)/(4.*R::beta(1./(BOX+2.),1+1./(BOX+2.)));
  const double LUMTOT = pow(10,(-0.4*(MAG-MAGZERO)));
  const double Ie=LUMTOT/(RE*RE*AXRAT*2.*PI*NSER*((exp(BN))/pow(BN,2*NSER))*R::gammafn(2*NSER)/RBOX);
  const double INVRE = 1.0/RE;
  // Do not change this! Read the function definitions for justification
  const double NSERFAC = nserfac(NSER,BOX);
  NumericMatrix mat(DIM(0), DIM(1));
  double x,y,xmid,ymid,xmod,ymod,rdivre,angmod;
  double xbin=(XLIM(1)-XLIM(0))/DIM(0);
  double ybin=(YLIM(1)-YLIM(0))/DIM(1);
  NumericVector xlim2(2),ylim2(2);
  //Get things into GALFIT's angle format (to make comparison easier)
  angmod = std::fmod(ANG+90.,360.);
  if(angmod > 180.) angmod -= 180.;
  const double PX = cos(angmod*M_PI/180.);
  const double INVREY = sin(angmod*M_PI/180.)*INVRE*pow(-1,angmod < PI);
  const double INVREX = PX*INVRE;
  const double INVAXRAT = 1.0/AXRAT;
  const double IEPIX = xbin*ybin*Ie;
  // End of precompute block
  
  // std::cout << RESWITCH << " " << UPSCALE << " " << MAXDEPTH << " " << ACC << std::endl;
  int i=0,j=0;
  x=XLIM(0);
  for(i = 0; i < DIM(0); i++) {
    xmid = x+xbin/2. - XCEN;
    y=YLIM(0);
    for(j = 0; j < DIM(1); j++) {
      mat(i,j)=0;
      if(DOCALCREGION==FALSE || CALCREGION(i,j)==1){
        ymid = y+ybin/2.-YCEN;
        xmod = xmid * INVREX + ymid * INVREY;
        ymod = (xmid * INVREY - ymid * INVREX)*INVAXRAT;
        rdivre = sqrt(xmod*xmod + ymod*ymod);
        if(rdivre<REMAX){
          if(rdivre>RESWITCH || ROUGH){
            mat(i,j)=profitEvalSersic<hasbox,t>(xmod, ymod, BN, BOX, NSERFAC);
          }
          else{
            xlim2(0)=x;
            xlim2(1)=x+xbin;
            ylim2(0)=y;
            ylim2(1)=y+ybin;
            if((std::abs(xmid)<1.0) && (std::abs(ymid)<1.0)) {
              mat(i,j)=profitSumPixMinorAxisGrad<hasbox,t>(XCEN,YCEN,xlim2,ylim2,INVREX,INVREY,INVAXRAT,
              NSERFAC, BOX,BN,8,0,10,ACC);
             }else{
               //Rcout << UPSCALE << std::endl;
            mat(i,j)=profitSumPixMinorAxisGrad<hasbox,t>(XCEN,YCEN,xlim2,ylim2,INVREX,INVREY,INVAXRAT,
              NSERFAC, BOX,BN,UPSCALE,0,MAXDEPTH,ACC);
            }
          }
        }
        mat(i,j)*=IEPIX;
      }
      y=y+ybin;
    }
    x=x+xbin;
  }
  return(mat);
}

template<> NumericMatrix profitMakeBoxySersic<false,gauss>(
    const IntegerMatrix CALCREGION, const double XCEN, const double YCEN, 
    const double MAG, const double RE, const double NSER, const double ANG,
    const double AXRAT, const double BOX, const double MAGZERO, const bool ROUGH,
    const NumericVector & XLIM, const NumericVector & YLIM, const IntegerVector & DIM, 
    const int UPSCALE, const int MAXDEPTH, const double RESWITCH, const double ACC, 
    const bool DOCALCREGION, const double REMAX) {
  //const double BN=R::qgamma(0.5, 2 * NSER,1,1,0);  
  const double BN = 0.69314718055994528623;
  const double RBOX=PI*(BOX+2.)/(4.*R::beta(1./(BOX+2.),1+1./(BOX+2.)));
  const double LUMTOT = pow(10,(-0.4*(MAG-MAGZERO)));
  // TODO: Figure out why this empirical factor of 2*exp(1) is required to normalize properly
  const double Ie=2.*exp(1.)*LUMTOT/(RE*RE*AXRAT*2.*PI*NSER*((exp(BN))/pow(BN,2*NSER))*R::gammafn(2*NSER)/RBOX);
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
      if(DOCALCREGION==FALSE || CALCREGION(i,j)==1){
        const std::vector<double> & rval = gauss2.integral(x,xhi,y,yhi,ACC2,bottomval,0,leftvals[j]);
        bottomval = rval[2];
        leftvals[j] = rval[4];
        mat(i,j)=rval[0]*Ie;
      } else {
        mat(i,j)=0;
        bottomval = 0;
        leftvals[j] = 0;
      }
      y = yhi;
      yhi += ybin;
      bottomval = 0;
    }
    x = xhi;
    xhi += xbin;
  }
  return(mat);
}

/*
  This is the exported, publicly available function for integrating a Sersic profile.
  It calls the appropriate templated method. Unfortunately, exporting templated
  functions is difficult, and anyway this ensures that the correct template
  is always called.
*/
// [[Rcpp::export]]
NumericMatrix profitMakeSersic(const IntegerMatrix & CALCREGION,
    const double XCEN=0, const double YCEN=0, const double MAG=15, 
    const double RE=1, const double NSER=1, const double ANG=0,
    const double AXRAT=1, const double BOX=0, 
    const double MAGZERO=0, const bool ROUGH=false,
    const NumericVector & XLIM = NumericVector::create(-100,100),
    const NumericVector & YLIM = NumericVector::create(-100,100),
    const IntegerVector & DIM = IntegerVector::create(200,200),
    const int UPSCALE=8L, const int MAXDEPTH=2L, const double RESWITCH=2,
    const double ACC=0.1, const bool DOCALCREGION=false, const double REMAX=10)
{
  //Rcout << UPSCALE << std::endl;
  if(BOX == 0) 
  {
    if(NSER == -0.5) return profitMakeBoxySersic<false,gauss>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
      AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC*1e-4, DOCALCREGION, REMAX);
    else if(NSER == 1) return profitMakeBoxySersic<false,exponent>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
      AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
    else if(NSER == 2) return profitMakeBoxySersic<false,two>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
      AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
    else if(NSER == 3) return profitMakeBoxySersic<false,three>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
      AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
    else if(NSER == 4) return profitMakeBoxySersic<false,four>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
      AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
    return profitMakeBoxySersic<false,general>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
      AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
  }
  if(NSER == -0.5) return profitMakeBoxySersic<true,gauss>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
    AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
  else if(NSER == 1) return profitMakeBoxySersic<true,exponent>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
    AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
  else if(NSER == 2) return profitMakeBoxySersic<true,two>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
    AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
  else if(NSER == 3) return profitMakeBoxySersic<true,three>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
    AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
  else if(NSER == 4) return profitMakeBoxySersic<true,four>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
    AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
  return profitMakeBoxySersic<true,general>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
    AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
}

inline void profitBruteConvRowCol(const double & IMGVAL, const NumericMatrix & PSF,
  NumericMatrix & output, const int & X_K, const int & Y_K, const int & ROW, const int & COL)
{
  double * output_row_col_j; 
  int i; int j;
  for (j = 0; j < Y_K; j++) {
    output_row_col_j = & output(ROW,COL+j) ;
    const double * psf_j = &PSF(0,j) ;
    for (i = 0; i < X_K; i++) {
      output_row_col_j[i] += IMGVAL*psf_j[i] ;
    }
  }
}

/*
  Brute force convolution of an image with a PSF.
  For the moment, it seems like convolving onto a padded image and then copying
  the result into a smaller (unpadded) result is faster than any kind of bounds
  checking (definitely faster than if statements and faster even than setting 
  appropriate limits for each loop).
*/
// [[Rcpp::export]]
NumericMatrix profitBruteConv(const NumericMatrix & IMG, const NumericMatrix & PSF,
  const IntegerMatrix & CALCREGION, const bool DOCALCREGION=false){
  const int X_S = IMG.nrow(), X_K = PSF.nrow();
  const int Y_S = IMG.ncol(), Y_K = PSF.ncol();
  const int EXTRAX= X_K % 2;
  const int EXTRAY= Y_K % 2;
    
  NumericMatrix output(X_S + X_K - EXTRAX, Y_S + Y_K - EXTRAY);
  //double* output_row_col_j ;
  //double image_row_col = 0.0 ;
  
  int col=0,row=0;
  if(DOCALCREGION)
  {
    for (col = 0; col < Y_S; col++) {
      const int * CALC_C = &CALCREGION(0,col);
      const double * IMG_C = &IMG(0,col);
      for (row = 0; row < X_S; row++) {
        if(CALC_C[row]) profitBruteConvRowCol(IMG_C[row], PSF, output, X_K, Y_K, row, col);
      }
    }
  } else {
    for (col = 0; col < Y_S; col++) {
      const double * IMG_C = &IMG(0,col);
      for (row = 0; row < X_S; row++) {
        //std::cout << col << " " << row << std::endl;
        profitBruteConvRowCol(IMG_C[row], PSF, output, X_K, Y_K, row, col);
      }
    } 
  }
  NumericMatrix cutout(X_S, Y_S);
  for (int row = 0; row < X_S; row++) {
    for (int col = 0; col < Y_S; col++) {
      cutout(row,col)=output(row+(X_K-EXTRAX)/2,col+(Y_K-EXTRAY)/2);
    }
  }
  return cutout;
}

// Downsamples an image. Probably not optimally implemented, but does it matter?
// [[Rcpp::export]]
NumericMatrix profitDownsample(const NumericMatrix & IMG, const int DOWNSAMPLEFAC){
  if(DOWNSAMPLEFAC == 1) return IMG;
  if(!DOWNSAMPLEFAC > 1) return NumericMatrix(0,0);
  const int X_S = IMG.nrow(), Y_S = IMG.ncol();
  const int X_D = X_S / DOWNSAMPLEFAC, Y_D = Y_S / DOWNSAMPLEFAC;
  NumericMatrix d(X_D,Y_D);
  int row,col,rowd, cold;
  for(col = 0; col < Y_S; ++col)
  {
    cold = col / DOWNSAMPLEFAC;
    const double * image_row = &IMG(0,col);
    double * ds_row = &d(0,cold);
    for(row = 0; row < X_S; ++row)
    {
      rowd = row / DOWNSAMPLEFAC;
      ds_row[rowd] += image_row[row];
    }
  }
  return d;
}

// Upsamples an image without interpolation. Probably not optimally implemented, but does it matter?
// [[Rcpp::export]]
NumericMatrix profitUpsample(const NumericMatrix & IMG, const int UPSAMPLEFAC){
  if(UPSAMPLEFAC == 1) return IMG;
  if(!UPSAMPLEFAC > 1) return NumericMatrix(0,0);
  const int X_S = IMG.nrow(), Y_S = IMG.ncol();
  const int X_U = X_S * UPSAMPLEFAC, Y_U = Y_S * UPSAMPLEFAC;
  NumericMatrix u(X_U,Y_U);
  int row,col,rowu, colu;
  for(colu = 0; colu < Y_U; ++colu)
  {
    col = colu / UPSAMPLEFAC;
    const double * image_row = &IMG(0,col);
    double * us_row = &u(0,colu);
    for(rowu = 0; rowu < X_U; ++rowu)
    {
      row = rowu / UPSAMPLEFAC;
      us_row[rowu] += image_row[row];
    }
  }
  return u;
}