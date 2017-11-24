/**
 * Utility routines for libprofit
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Rodrigo Tobar
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

#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "profit/common.h"
#include "profit/config.h"
#include "profit/utils.h"

/*
 * We use either GSL or R to provide the low-level
 * beta, gamma and pgamma and qgamma functions needed by some profiles.
 * If neither is given the compilation should fail
 */
#if defined(PROFIT_USES_GSL)
	#include <gsl/gsl_errno.h>
	#include <gsl/gsl_cdf.h>
	#include <gsl/gsl_sf_gamma.h>
	#include <gsl/gsl_integration.h>
#elif defined(PROFIT_USES_R)
	#include <Rmath.h>
	#include <R_ext/Applic.h>

	/*
	 * R and Rmath play a tricky game where functions with "fully qualified
	 * names" (e.g., Rf_qgamma) are present only when compiling code with
	 * undef(MATHLIB_STANDALONE) || undef(R_NO_REMAP_RMATH).
	 * We currently use the same names in some of our functions. On the other
	 * hand we compile this code both using the stand-alone Rmath library and
	 * as part of an R extension. This means that:
	 *
	 * * When compiling as standalone we can't access the Rf_* names and
	 * * When compiling as an R extension our function names are replaced during
	 *   pre-processing with Rf_* names, and our namespace ends up containing
	 *   functions like profit::Rf_qgamma.
	 *
	 * The code below solves these two problems by avoiding name clashes and
	 * letting us use the Rf_* names in our code seamessly. In the future we may
	 * want to change our function names and be safer
	 */
	#if defined(MATHLIB_STANDALONE)
	#define Rf_qgamma   qgamma
	#define Rf_pgamma   pgamma
	#define Rf_gammafn  gammafn
	#define Rf_beta     beta
	#else
	#undef qgamma
	#undef pgamma
	#undef gammafn
	#undef beta
	#endif

#else
	#error("No high-level library (GSL or R) provided")
#endif


namespace profit {

bool almost_equals(double x, double y, double e) {
	return std::abs(x - y) < std::abs(e);
}

/*
 * GSL-based functions
 */
#if defined(PROFIT_USES_GSL)
double qgamma(double p, double shape) {
	return gsl_cdf_gamma_Pinv(p, shape, 1);
}

double pgamma(double q, double shape) {
	return gsl_cdf_gamma_P(q, shape, 1);
}

double gammafn(double x) {

	gsl_sf_result result;
	int status = gsl_sf_gamma_e(x, &result);
	if( status ) {
		if( status == GSL_EUNDRFLW ) {
			return 0.;
		}
		else if( status == GSL_EOVRFLW && x > 0 ) {
			return std::numeric_limits<double>::infinity();
		}
		return std::numeric_limits<double>::quiet_NaN();
	}

	return result.val;
}

double beta(double a, double b) {

	if( a < 0. || b < 0. ) {
		return std::numeric_limits<double>::quiet_NaN();
	}
	if( a == 0. || b == 0. ) {
		return std::numeric_limits<double>::infinity();
	}

	gsl_sf_result result;
	int status = gsl_sf_beta_e(a, b, &result);
	if( status ) {
		if( status == GSL_EUNDRFLW ) {
			return 0.;
		}
		return std::numeric_limits<double>::quiet_NaN();
	}

	return result.val;
}

static
double __gsl_integrate_qag(integration_func_t f, void *params,
                           double a, double b, bool to_infinity) {

	size_t limit = 100;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc (limit);
	gsl_function F;
	F.function = f;
	F.params = params;
	double epsabs = 1e-4, epsrel = 1e-4;
	double result, abserr;

	if( to_infinity ) {
		gsl_integration_qagiu(&F, a, epsabs, epsrel, limit, w, &result, &abserr);
	}
	else {
		gsl_integration_qags(&F, a, b, epsabs, epsrel, limit, w, &result, &abserr);
	}
	gsl_integration_workspace_free (w);

	return result;
}

double integrate_qagi(integration_func_t f, double a, void *params) {
	return __gsl_integrate_qag(f, params, a, 0, true);
}

double integrate_qags(integration_func_t f, double a, double b, void *params) {
	return __gsl_integrate_qag(f, params, a, b, false);
}

/*
 * R-based functions
 */
#elif defined(PROFIT_USES_R)

double qgamma(double p, double shape) {
	return ::Rf_qgamma(p, shape, 1, 1, 0);
}

double pgamma(double q, double shape) {
	return ::Rf_pgamma(q, shape, 1, 1, 0);
}

double gammafn(double x) {
	return ::Rf_gammafn(x);
}

double beta(double a, double b) {
	return ::Rf_beta(a, b);
}

struct __r_integrator_args {
	integration_func_t f;
	void *params;
};

static
void __r_integrator(double *x, int n, void *ex) {
	struct __r_integrator_args *int_args = (struct __r_integrator_args *)ex;
	for(auto i=0; i<n; i++) {
		x[i] = int_args->f(x[i], int_args->params);
	}
}

static
double __r_integrate_qag(integration_func_t f, void *params,
                       double a, double b, bool to_infinity) {

	int neval, ier, last;
	int limit = 100;
	int lenw = 4 * limit;
	std::vector<int> iwork(limit);
	std::vector<double> work(lenw);
	double result, abserr;
	double epsabs = 1e-4, epsrel = 1e-4;
	struct __r_integrator_args int_args = {f, params};

	if( to_infinity ) {
		int inf = 1;
		::Rdqagi(&__r_integrator, &int_args, &a, &inf,
		         &epsabs, &epsrel, &result, &abserr, &neval, &ier,
		         &limit, &lenw, &last,
		         iwork.data(), work.data());
	}
	else {
		::Rdqags(&__r_integrator, &int_args, &a, &b,
		         &epsabs, &epsrel, &result, &abserr, &neval, &ier,
		         &limit, &lenw, &last,
		         iwork.data(), work.data());
	}

	return result;
}

double integrate_qagi(integration_func_t f, double a, void *params) {
	return __r_integrate_qag(f, params, a, 0, true);
}

double integrate_qags(integration_func_t f, double a, double b, void *params) {
	return __r_integrate_qag(f, params, a, b, false);
}

#endif

} /* namespace profit */
