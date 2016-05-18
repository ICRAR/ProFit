#include <math.h>

/* Use the cannonical Rf_* names */
#define R_NO_REMAP
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <profit.h>
#include <sersic.h>
#include <sky.h>
#include <psf.h>


extern "C" {

static
SEXP _get_list_element(SEXP list, const char *name) {
	SEXP names = Rf_getAttrib(list, R_NamesSymbol);
	for(unsigned int i=0; i!=Rf_length(list); i++) {
		if( !strcmp(name, CHAR(STRING_ELT(names, i))) ) {
			return VECTOR_ELT(list, i);
		}
	}
	return R_NilValue;
}

static
void _read_real(SEXP list, const char *name, unsigned int idx, double *target) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		*target = REAL(element)[idx];
	}
}

static
void list_to_sersic(SEXP sersic_list, profit_profile *p, unsigned int idx) {
	profit_sersic_profile *sp = (profit_sersic_profile *)p;
	_read_real(sersic_list, "xcen",  idx, &(sp->xcen));
	_read_real(sersic_list, "ycen",  idx, &(sp->ycen));
	_read_real(sersic_list, "mag",   idx, &(sp->mag));
	_read_real(sersic_list, "re",    idx, &(sp->re));
	_read_real(sersic_list, "nser",  idx, &(sp->nser));
	_read_real(sersic_list, "ang",   idx, &(sp->ang));
	_read_real(sersic_list, "axrat", idx, &(sp->axrat));
	_read_real(sersic_list, "box",   idx, &(sp->box));
}

static
void list_to_sky(SEXP sky_list, profit_profile *p, unsigned int idx) {
	profit_sky_profile *sp = (profit_sky_profile *)p;
	_read_real(sky_list, "bg",  idx, &(sp->bg));
}

static
void list_to_psf(SEXP psf_list, profit_profile *p, unsigned int idx) {
	profit_psf_profile *psf = (profit_psf_profile *)p;
	_read_real(psf_list, "xcen",  idx, &(psf->xcen));
	_read_real(psf_list, "ycen",  idx, &(psf->ycen));
	_read_real(psf_list, "mag",   idx, &(psf->mag));
}


static
void _read_profiles(profit_model *model, SEXP model_list,
                    const char *profile_name, const char* example_property,
                    void (*list_to_profile)(SEXP, profit_profile *, unsigned int)) {

	/* Is the profile specified in the model? */
	SEXP profile_list = _get_list_element(model_list, profile_name);
	if( profile_list == R_NilValue ) {
		return;
	}

	/* This property should exist and should tell us the number of profiles there are */
	SEXP property = _get_list_element(profile_list, example_property);
	if( property == R_NilValue ) {
		return;
	}

	/* OK, we know now how many there are... */
	unsigned int i, count = Rf_length(property);

	/* Create that many profiles and then start reading the info if available */
	for(i=0; i!=count; i++) {
		profit_profile *p = profit_create_profile(profile_name);
		profit_add_profile(model, p);
		list_to_profile(profile_list, p, i);
	}
}

static
void _read_sersic_profiles(profit_model *model, SEXP model_list) {
	_read_profiles(model, model_list, "sersic", "xcen", &list_to_sersic);
}

static
void _read_sky_profiles(profit_model *model, SEXP model_list) {
	_read_profiles(model, model_list, "sky", "bg", &list_to_sky);
}

static
void _read_psf_profiles(profit_model *model, SEXP model_list) {
	_read_profiles(model, model_list, "psf", "xcen", &list_to_psf);
}

SEXP R_profit_make_model(SEXP model_list, SEXP magzero, SEXP dim) {

	ssize_t size;
	double *dim_l = REAL(dim);
	unsigned int i, p;
	char *error;

	/* Read model parameters */
	profit_model *m = profit_create_model();
	m->magzero = Rf_asReal(magzero);
	m->width  = m->res_x = (unsigned int)dim_l[0];
	m->height = m->res_y = (unsigned int)dim_l[1];

	/* Read profiles and parameters and append them to the model */
	_read_sersic_profiles(m, model_list);
	_read_sky_profiles(m, model_list);
	_read_psf_profiles(m, model_list);
	if( !m->n_profiles ) {
		Rf_error("No profiles found in incoming model");
		profit_cleanup(m);
		return R_NilValue;
	}

	/* Go, go, go! */
	profit_eval_model(m);
	error = profit_get_error(m);
	if( error ) {
		Rprintf("Error while calculating model: %s", error);
		profit_cleanup(m);
		return R_NilValue;
	}

	/* Copy the image, clean up, and good bye */
	size = m->width * m->height;
	SEXP image = PROTECT(Rf_allocVector(REALSXP, size));
	memcpy(REAL(image), m->image, sizeof(double) * size);

	profit_cleanup(m);

	UNPROTECT(1);
	return image;
}

}
