/**
 * King profile implementation
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Aaron Robotham
 *
 * This file is part of libprofit.
 *
 * libprofit is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libprofit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libprofit.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>

#include <R.h>
#include <R_ext/Applic.h>

#include "profit/king.h"
#include "profit/utils.h"

using namespace std;

namespace profit
{

/**
 * The evaluation of the king profile at king coordinates (x,y).
 *
 * The king profile has this form:
 *
 *    (1+(r/rb)^(-a))^(b/a)*
 *        exp(-bn*(((r^a+rb^a)/re^a))^(1/(nser*a)))
 *
 * with::
 *
 *           r = (x^{2+B} + y^{2+B})^{1/(2+B)}
 *           B = box parameter
 */
static
double _king_for_xy_r(RadialProfile *sp,
                        double x, double y,
                        double r, bool reuse_r) {

	KingProfile *kp = static_cast<KingProfile *>(sp);
	double r_factor;
	if( reuse_r && kp->box == 0 ) {
		r_factor = r;
	}
	else if( kp->box == 0 ) {
		r_factor = sqrt(x*x + y*y);
	}
	else {
		double box = kp->box + 2.;
		r_factor = pow( pow(abs(x), box) + pow(abs(y), box), 1./box);
	}
  if( r_factor < kp->rt ) {
		double temp=1/pow(1+pow(kp->rt/kp->rc,2),1/kp->a);
    return pow(1-temp,-kp->a)*
          pow(1/pow(1+pow(r_factor/kp->rc,2),1/kp->a)-temp,kp->a);
	}
	else {
		return 0;
	}
}

eval_function_t KingProfile::get_evaluation_function() {
	return &_king_for_xy_r;
}

void king_int(double *x, int n, void *ex) {
  KingProfile *kp = (KingProfile *)ex;
  double temp=1/pow(1+pow(kp->rt/kp->rc,2),1/kp->a);
  for(auto i=0; i!=n; i++) {
    double r = x[i];
    if( r < kp->rt ) {
      x[i] = r*pow(1-temp,-kp->a)*
          pow(1/pow(1+pow(r/kp->rc,2),1/kp->a)-temp,kp->a);
    }
    else {
      x[i] = 0;
    }
  }
}


double KingProfile::get_lumtot(double r_box) {
  /* 
   * Need correct numerical integral here. In R we do:
   * 
   * if(box!=0){Rbox=pi*(box+2)/(4*beta(1/(box+2),1+1/(box+2)))}else{Rbox=1}
   * lumtot = 2*pi*axrat*integrate(.profitKingR, 0, Inf, re=re, rb=rb, nser=nser, a=a, b=b, bn=bn)$value/Rbox
   * magtot = -2.5 * log10(lumtot)
   * return(1/(10^(0.4 * (mag - magtot))))
   * 
   * Can somehow replace R's stats::integrate with GSL version.
   * 
   */
  
  int neval, ier, last;
  int limit = 100;
  int lenw = 4 * limit;
  int *iwork = new int[limit];
  double *work = new double[lenw];
  double result, abserr;
  double a= 0, b= this->rt, epsabs = 1e-4, epsrel = 1e-4;
  ::Rdqags(&king_int, this, &a, &b,
            &epsabs, &epsrel, &result, &abserr, &neval, &ier,
            &limit, &lenw, &last,
            iwork, work);

  delete [] iwork;
  delete [] work;
	return 2*M_PI * axrat * result/r_box;
}

double KingProfile::adjust_rscale_switch() {
	return 0.5;
}

double KingProfile::adjust_rscale_max() {
	return 1;
}

double KingProfile::get_rscale() {
	return this->rt;
}

double KingProfile::adjust_acc() {
	return this->acc;
}

KingProfile::KingProfile() :
	RadialProfile(),
	rc(1), rt(3), a(2)
{
	// this profile defaults to a different accuracy
	
}

} /* namespace profit */
