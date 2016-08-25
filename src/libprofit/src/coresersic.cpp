/**
 * CoreSersic profile implementation
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

#include "profit/coresersic.h"
#include "profit/utils.h"

using namespace std;

namespace profit
{

/**
 * The evaluation of the coresersic profile at coresersic coordinates (x,y).
 *
 * The coresersic profile has this form:
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
double _coresersic_for_xy_r(RadialProfile *sp,
                        double x, double y,
                        double r, bool reuse_r) {

	CoreSersicProfile *csp = static_cast<CoreSersicProfile *>(sp);
	double r_factor;
	if( reuse_r && csp->box == 0 ) {
		r_factor = r;
	}
	else if( csp->box == 0 ) {
		r_factor = sqrt(x*x + y*y);
	}
	else {
		double box = csp->box + 2.;
		r_factor = pow( pow(abs(x), box) + pow(abs(y), box), 1./box);
	}

	return pow(1+pow(r_factor/csp->rb,-csp->a),csp->b/csp->a)*
         exp(-csp->_bn*pow(((pow(r_factor,csp->a)+pow(csp->rb,csp->a))/pow(csp->re,csp->a)),1/(csp->nser*csp->a)));

}

eval_function_t CoreSersicProfile::get_evaluation_function() {
	return &_coresersic_for_xy_r;
}

void coresersic_int(double *x, int n, void *ex) {
  CoreSersicProfile *csp = (CoreSersicProfile *)ex;
  for(auto i=0; i!=n; i++) {
    double r = x[i];
    x[i] = r*pow(1+pow(r/csp->rb,-csp->a),csp->b/csp->a)*
           exp(-csp->_bn*pow(((pow(r,csp->a)+pow(csp->rb,csp->a))/pow(csp->re,csp->a)),1/(csp->nser*csp->a)));
  }
}


double CoreSersicProfile::get_lumtot(double r_box) {
  /* 
   * Need correct numerical integral here. In R we do:
   * 
   * if(box!=0){Rbox=pi*(box+2)/(4*beta(1/(box+2),1+1/(box+2)))}else{Rbox=1}
   * lumtot = 2*pi*axrat*integrate(.profitCoreSersicR, 0, Inf, re=re, rb=rb, nser=nser, a=a, b=b, bn=bn)$value/Rbox
   * magtot = -2.5 * log10(lumtot)
   * return(1/(10^(0.4 * (mag - magtot))))
   * 
   * Can somehow replace R's stats::integrate with GSL version.
   * 
   */
  
  int neval, ier, last, inf = 1;
  int limit = 100;
  int lenw = 4 * limit;
  int *iwork = new int[limit];
  double *work = new double[lenw];
  double result, abserr;
  double a= 0, epsabs = 1e-4, epsrel = 1e-4;
  ::Rdqagi(&coresersic_int, this, &a, &inf,
            &epsabs, &epsrel, &result, &abserr, &neval, &ier,
            &limit, &lenw, &last,
            iwork, work);
  
  delete [] iwork;
  delete [] work;
	return 2*M_PI * axrat * result/r_box;
}

void CoreSersicProfile::initial_calculations() {

	/*
	 * bn needs to be calculated before calling the super method
	 * because it's used to calculate the total luminosity
	 */
	this->_bn = qgamma(0.5, 2*this->nser);

	/* Common calculations first */
	RadialProfile::initial_calculations();

}

double CoreSersicProfile::adjust_rscale_switch() {
	return this->re*4;
}

double CoreSersicProfile::adjust_rscale_max() {
	return this->re*100;
}

double CoreSersicProfile::get_rscale() {
	return this->re;
}

double CoreSersicProfile::adjust_acc() {
	return this->acc;
}

CoreSersicProfile::CoreSersicProfile() :
	RadialProfile(),
	re(1), rb(1), nser(4), a(1), b(1)
{
	// this profile defaults to a different accuracy
	
}

} /* namespace profit */
