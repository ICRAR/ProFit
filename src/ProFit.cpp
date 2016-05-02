#include <profit.h>
#include <sersic.h>
#include <sky.h>

#include <math.h>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

extern "C" {

double _qgamma_wrapper(double a, double b, double c) {
	return Rf_qgamma(a, b, c, 1, 0);
}

SEXP _get_list_element(SEXP list, const char *name) {
	SEXP names = getAttrib(list, R_NamesSymbol);
	for(unsigned int i=0; i!=length(list); i++) {
		if( !strcmp(name, CHAR(STRING_ELT(names, i))) ) {
			return VECTOR_ELT(list, i);
		}
	}
	return R_NilValue;
}

void _read_real(SEXP list, const char *name, unsigned int idx, double *target) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		*target = REAL(element)[idx];
	}
}

profit_profile **_read_sky_profiles(SEXP sky_list, unsigned int *count) {

	unsigned int i;

	SEXP bg = _get_list_element(sky_list, "bg");
	if( bg == R_NilValue ) {
		*count = 0;
		return NULL;
	}

	/* OK, we know now how many there are... */
	*count = length(bg);
	profit_profile **all = (profit_profile **)malloc(sizeof(profit_profile*) * *count);

	/* Create that many profiles and then start reading the info if available */
	for(i=0; i!=*count; i++) {
		profit_profile *p = profit_get_profile("sky");
		all[i] = p;
		profit_sky_profile *sp = (profit_sky_profile *)p;

		_read_real(sky_list, "bg",  i, &(sp->bg));
	}

	return all;
}

profit_profile **_read_sersic_profiles(SEXP sersic_list, unsigned int *count) {

	unsigned int i;

	SEXP xcen = _get_list_element(sersic_list, "xcen");
	if( xcen == R_NilValue ) {
		*count = 0;
		return NULL;
	}

	/* OK, we know now how many there are... */
	*count = length(xcen);
	profit_profile **all = (profit_profile **)malloc(sizeof(profit_profile*) * *count);

	/* Create that many profiles and then start reading the info if available */
	for(i=0; i!=*count; i++) {
		profit_profile *p = profit_get_profile("sersic");
		all[i] = p;
		profit_sersic_profile *sp = (profit_sersic_profile *)p;

		_read_real(sersic_list, "xcen",  i, &(sp->xcen));
		_read_real(sersic_list, "ycen",  i, &(sp->ycen));
		_read_real(sersic_list, "mag",   i, &(sp->mag));
		_read_real(sersic_list, "re",    i, &(sp->re));
		_read_real(sersic_list, "nser",  i, &(sp->nser));
		_read_real(sersic_list, "ang",   i, &(sp->ang));
		_read_real(sersic_list, "axrat", i, &(sp->axrat));
		_read_real(sersic_list, "box",   i, &(sp->box));
		sp->_qgamma = &_qgamma_wrapper;
		sp->_gammafn = &Rf_gammafn;
		sp->_beta = &Rf_beta;
	}

	return all;
}

SEXP R_profit_make_model(SEXP model_list, SEXP magzero, SEXP dim) {

	ssize_t size;
	double *dim_l;
	unsigned int n_sersic = 0, n_sky = 0, n_profiles = 0;
	unsigned int i, p;
	profit_profile **sersic_profiles = NULL;
	profit_profile **sky_profiles = NULL;
	profit_profile **all_profiles = NULL;

	/* Inspect which profiles we have specified and creates them */
	SEXP sersic_list = _get_list_element(model_list, "sersic");
	if( sersic_list != R_NilValue ) {
		sersic_profiles = _read_sersic_profiles(sersic_list, &n_sersic);
		n_profiles += n_sersic;
	}
	SEXP sky_list = _get_list_element(model_list, "sky");
	if( sky_list != R_NilValue ) {
		sky_profiles = _read_sky_profiles(sky_list, &n_sky);
		n_profiles += n_sky;
	}

	if( !n_profiles ) {
		error("No profiles found in incoming model");
		return R_NilValue;
	}

	/* Combine all profiles into a single list */
	all_profiles = (profit_profile **)malloc(sizeof(profit_profile *) * n_profiles);
	profit_profile **dst = all_profiles;
	memcpy(dst,             sersic_profiles, sizeof(profit_profile *) * n_sersic);
	memcpy(dst += n_sersic, sky_profiles,    sizeof(profit_profile *) * n_sky);
	free(sersic_profiles);
	free(sky_profiles);

	/* Create the model */
	profit_model *m = (profit_model *)malloc(sizeof(profit_model));
	m->n_profiles = n_profiles;
	m->profiles = all_profiles;
	m->magzero = asReal(magzero);
	dim_l = REAL(dim);
	m->width  = m->res_x = (unsigned int)dim_l[0];
	m->height = m->res_y = (unsigned int)dim_l[1];
	size = m->width * m->height;

	/* Go, go, go! */
	if( profit_make_model(m) ) {
		Rprintf("Error while calculating model :(");
		return R_NilValue;
	}

	/* Copy the image, clean up, and good bye */
	SEXP image = PROTECT(allocVector(REALSXP, size));
	memcpy(REAL(image), m->image, sizeof(double) * size);

	for(i=0; i!=m->n_profiles; i++) {
		free(m->profiles[i]);
	}
	free(m->profiles);
	free(m->image);

	UNPROTECT(1);
	return image;
}

}

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix profitBruteConv(NumericMatrix image, NumericMatrix psf){
    int x_s = image.nrow(), x_k = psf.nrow();
    int y_s = image.ncol(), y_k = psf.ncol();
    int extrax= x_k % 2;
    int extray= y_k % 2;

    NumericMatrix output(x_s + x_k - extrax, y_s + y_k - extray);
    double* output_row_col_j ;
    double image_row_col = 0.0 ;
    double* psf_j ;

    for (int row = 0; row < x_s; row++) {
      for (int col = 0; col < y_s; col++) {
        image_row_col = image(row,col) ;
        for (int j = 0; j < y_k; j++) {
          output_row_col_j = & output(row,col+j) ;
          psf_j = &psf(0,j) ;
          for (int i = 0; i < x_k; i++) {
           output_row_col_j[i] += image_row_col*psf_j[i] ;
          }
        }
      }
    }
    NumericMatrix cutout(x_s, y_s);
    for (int row = 0; row < x_s; row++) {
      for (int col = 0; col < y_s; col++) {
      cutout(row,col)=output(row+(x_k-extrax)/2,col+(y_k-extray)/2);
      }
    }
   return cutout;
}
