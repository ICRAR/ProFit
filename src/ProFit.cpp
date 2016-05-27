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
void _read_bool(SEXP list, const char *name, unsigned int idx, bool *target) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		if( TYPEOF(element) == LGLSXP ) {
			*target = (bool)LOGICAL(element)[idx];
		}
		else if( TYPEOF(element) == INTSXP ) {
			*target = (bool)INTEGER(element)[idx];
		}
		else {
			Rf_error("Parameter %s[%u] should be of boolean or integer type", name, idx);
		}
	}
}

static
void _read_unsigned_int(SEXP list, const char *name, unsigned int idx, unsigned int *target) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		if( TYPEOF(element) == INTSXP ) {
			*target = (unsigned int)INTEGER(element)[idx];
		}
		else if( TYPEOF(element) == LGLSXP ) {
			*target = (unsigned int)LOGICAL(element)[idx];
		}
		else if( TYPEOF(element) == REALSXP ) {
			*target = (unsigned int)REAL(element)[idx];
		}
		else {
			Rf_error("Parameter %s[%u] should be of numeric type", name, idx);
		}
	}
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

	_read_bool(sersic_list, "rough",   idx, &(sp->rough));
	_read_real(sersic_list, "acc",   idx, &(sp->acc));
	_read_real(sersic_list, "re_switch",   idx, &(sp->re_switch));
	_read_unsigned_int(sersic_list, "resolution",   idx, &(sp->resolution));
	_read_unsigned_int(sersic_list, "max_recursions",   idx, &(sp->max_recursions));
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
void _read_profiles(profit_model *model, SEXP profiles_list,
                    const char *profile_name, const char* example_property,
                    void (*list_to_profile)(SEXP, profit_profile *, unsigned int)) {

	/* Is the profile specified in the model? */
	SEXP profile_list = _get_list_element(profiles_list, profile_name);
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

		_read_bool(profile_list, "convolve", i, &(p->convolve));
	}
}

static
void _read_sersic_profiles(profit_model *model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "sersic", "xcen", &list_to_sersic);
}

static
void _read_sky_profiles(profit_model *model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "sky", "bg", &list_to_sky);
}

static
void _read_psf_profiles(profit_model *model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "psf", "xcen", &list_to_psf);
}

SEXP R_profit_make_model(SEXP model_list) {

	ssize_t size;
	unsigned int i, p;
	char *error;
	double psf_width = 0, psf_height = 0;
	double *psf = NULL;

	SEXP width = _get_list_element(model_list, "width");
	if( width == R_NilValue ) {
		Rf_error("No width provided in the model\n");
		return R_NilValue;
	}
	SEXP height = _get_list_element(model_list, "height");
	if( height == R_NilValue ) {
		Rf_error("No height provided in the model\n");
		return R_NilValue;
	}
	SEXP magzero = _get_list_element(model_list, "magzero");
	if( magzero == R_NilValue ) {
		Rf_error("No magzero provided in the model\n");
		return R_NilValue;
	}
	SEXP r_psf = _get_list_element(model_list, "psf");
	if( r_psf != R_NilValue ) {
		psf_width = Rf_nrows(r_psf);
		psf_height = Rf_ncols(r_psf);
		psf = (double *)malloc(sizeof(double) * psf_width * psf_height);
		memcpy(psf, REAL(r_psf), sizeof(double) * psf_width * psf_height);
	}

	/* Read model parameters */
	profit_model *m = profit_create_model();
	m->width  = m->res_x = (unsigned int)Rf_asReal(width);
	m->height = m->res_y = (unsigned int)Rf_asReal(height);
	m->magzero = Rf_asReal(magzero);
	m->psf = psf;
	m->psf_width = psf_width;
	m->psf_height = psf_height;

	/* Read profiles and parameters and append them to the model */
	SEXP profiles = _get_list_element(model_list, "profiles");
	if( profiles == R_NilValue ) {
		Rf_error("No profiles provided in the model\n");
		return R_NilValue;
	}
	_read_sersic_profiles(m, profiles);
	_read_sky_profiles(m, profiles);
	_read_psf_profiles(m, profiles);
	if( !m->n_profiles ) {
		Rf_error("No valid profiles found in profiles list\n");
		profit_cleanup(m);
		return R_NilValue;
	}

	/* Go, go, go! */
	profit_eval_model(m);
	error = profit_get_error(m);
	if( error ) {
		Rprintf("Error while calculating model: %s\n", error);
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
