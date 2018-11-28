#include <cmath>
#include <memory>
#include <sstream>
#include <vector>

#include <profit/profit.h>

/* Use the cannonical Rf_* names */
#define R_NO_REMAP
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

using namespace profit;
using namespace std;

static
Image _read_image(SEXP r_image) {

	unsigned int im_width = Rf_nrows(r_image);
	unsigned int im_height = Rf_ncols(r_image);
	unsigned int size = im_width * im_height;

	double *image_real_ptr;
	int *image_int_ptr;
	vector<double> image;
	switch (TYPEOF(r_image)) {

		case REALSXP:
			image_real_ptr = REAL(r_image);
			image = vector<double>(image_real_ptr, image_real_ptr + size);
			break;

		case INTSXP:
			image_int_ptr = INTEGER(r_image);
			image = vector<double>(size);
			std::copy(image_int_ptr, image_int_ptr + size, image.begin());
			break;

		case LGLSXP:
			image_int_ptr = LOGICAL(r_image);
			image = vector<double>(size);
			std::copy(image_int_ptr, image_int_ptr + size, image.begin());
			break;


		default:
			Rf_error("Image not in one of the supported formats (logical, integer, double)");
			image = {};
			break;
	}

	return Image(std::move(image), im_width, im_height);
}

static
Mask _read_mask(SEXP r_mask) {
	int *r_raw_mask = LOGICAL(r_mask);
	unsigned int m_width = Rf_nrows(r_mask);
	unsigned int m_height = Rf_ncols(r_mask);
	vector<bool> mask(r_raw_mask, r_raw_mask + (m_width * m_height));
	return Mask(std::move(mask), m_width, m_height);
}

static
SEXP _get_list_element(SEXP list, const char *name) {
	SEXP names = Rf_getAttrib(list, R_NamesSymbol);
	for(int i=0; i < Rf_length(list); i++) {
		if( !strcmp(name, CHAR(STRING_ELT(names, i))) ) {
			return VECTOR_ELT(list, i);
		}
	}
	return R_NilValue;
}

static
void _read_bool(shared_ptr<Profile> p, SEXP list, const char *name, unsigned int idx) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		if( TYPEOF(element) == LGLSXP ) {
			p->parameter(name, (bool)LOGICAL(element)[idx]);
		}
		else if( TYPEOF(element) == INTSXP ) {
			p->parameter(name, (bool)INTEGER(element)[idx]);
		}
		else if( TYPEOF(element) == REALSXP ) {
			p->parameter(name, (bool)REAL(element)[idx]);
		}
		else {
			Rf_error("Parameter %s[%u] should be of boolean or integer type", name, idx);
		}
	}
}

static
void _read_uint(shared_ptr<Profile> p, SEXP list, const char *name, unsigned int idx) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		if( TYPEOF(element) == INTSXP ) {
			p->parameter(name, (unsigned int)INTEGER(element)[idx]);
		}
		else if( TYPEOF(element) == LGLSXP ) {
			p->parameter(name, (unsigned int)LOGICAL(element)[idx]);
		}
		else if( TYPEOF(element) == REALSXP ) {
			p->parameter(name, (unsigned int)REAL(element)[idx]);
		}
		else {
			Rf_error("Parameter %s[%u] should be of numeric type", name, idx);
		}
	}
}

static
void _read_real(shared_ptr<Profile> p, SEXP list, const char *name, unsigned int idx) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		p->parameter(name, REAL(element)[idx]);
	}
}

static
void list_to_radial(SEXP radial_list, shared_ptr<Profile> p, unsigned int idx) {
	// p->parameter("adjust", true);
	_read_real(p, radial_list, "xcen",  idx);
	_read_real(p, radial_list, "ycen",  idx);
	_read_real(p, radial_list, "mag",   idx);
	_read_real(p, radial_list, "ang",   idx);
	_read_real(p, radial_list, "axrat", idx);
	_read_real(p, radial_list, "box",   idx);

	_read_bool(p, radial_list, "rough",          idx);
	_read_bool(p, radial_list, "adjust",         idx);
	_read_real(p, radial_list, "acc",            idx);
	_read_real(p, radial_list, "rscale_switch",  idx);
	_read_uint(p, radial_list, "resolution",     idx);
	_read_uint(p, radial_list, "max_recursions", idx);

	_read_real(p, radial_list, "rscale_max", idx);
}

static
void list_to_sersic(SEXP sersic_list, shared_ptr<Profile> p, unsigned int idx) {
	list_to_radial(sersic_list, p, idx);
	_read_real(p, sersic_list, "re",           idx);
	_read_real(p, sersic_list, "nser",         idx);
	_read_bool(p, sersic_list, "rescale_flux", idx);
}

static
void list_to_moffat(SEXP moffat_list, shared_ptr<Profile> p, unsigned int idx) {
	list_to_radial(moffat_list, p, idx);
	_read_real(p, moffat_list, "fwhm", idx);
	_read_real(p, moffat_list, "con",  idx);
}

static
void list_to_ferrer(SEXP ferrer_list, shared_ptr<Profile> p, unsigned int idx) {
	list_to_radial(ferrer_list, p, idx);
	_read_real(p, ferrer_list, "rout",  idx);
	_read_real(p, ferrer_list, "a",     idx);
	_read_real(p, ferrer_list, "b",     idx);
}

static
void list_to_coresersic(SEXP coresersic_list, shared_ptr<Profile> p, unsigned int idx) {
	list_to_radial(coresersic_list, p, idx);
	_read_real(p, coresersic_list, "re",   idx);
	_read_real(p, coresersic_list, "rb",   idx);
	_read_real(p, coresersic_list, "nser", idx);
	_read_real(p, coresersic_list, "a",    idx);
	_read_real(p, coresersic_list, "b",    idx);
}

static
void list_to_king(SEXP king_list, shared_ptr<Profile> p, unsigned int idx) {
	list_to_radial(king_list, p, idx);
	_read_real(p, king_list, "rc", idx);
	_read_real(p, king_list, "rt", idx);
	_read_real(p, king_list, "a",  idx);
}

static
void list_to_brokenexponential(SEXP brokenexponential_list, shared_ptr<Profile> p, unsigned int idx) {
	list_to_radial(brokenexponential_list, p, idx);
	_read_real(p, brokenexponential_list, "h1", idx);
	_read_real(p, brokenexponential_list, "h2", idx);
	_read_real(p, brokenexponential_list, "rb", idx);
	_read_real(p, brokenexponential_list, "a", idx);
}

static
void list_to_sky(SEXP sky_list, shared_ptr<Profile> p, unsigned int idx) {
	_read_real(p, sky_list, "bg", idx);
}

static
void list_to_psf(SEXP psf_list, shared_ptr<Profile> p, unsigned int idx) {
	_read_real(p, psf_list, "xcen",  idx);
	_read_real(p, psf_list, "ycen",  idx);
	_read_real(p, psf_list, "mag",   idx);
}


static
void _read_profiles(Model &model, SEXP profiles_list,
                    const char *profile_name, const char* example_property,
                    void (*list_to_profile)(SEXP, shared_ptr<Profile>, unsigned int)) {

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
		try {
			shared_ptr<Profile> p = model.add_profile(profile_name);
			_read_bool(p, profile_list, "convolve", i);
			list_to_profile(profile_list, p, i);
		} catch (invalid_parameter &e) {
			ostringstream os;
			os << "Error while creating profile '" << profile_name << "': " << e.what();
			Rf_error("%s\n", os.str().c_str());
		}
	}
}

static
void _read_sersic_profiles(Model &model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "sersic", "xcen", &list_to_sersic);
}

static
void _read_moffat_profiles(Model &model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "moffat", "xcen", &list_to_moffat);
}

static
void _read_ferrer_profiles(Model &model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "ferrer", "xcen", &list_to_ferrer);
  _read_profiles(model, profiles_list, "ferrers", "xcen", &list_to_ferrer);
}

static
void _read_coresersic_profiles(Model &model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "coresersic", "xcen", &list_to_coresersic);
}

static
void _read_king_profiles(Model &model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "king", "xcen", &list_to_king);
}

static
void _read_brokenexponential_profiles(Model &model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "brokenexp", "xcen", &list_to_brokenexponential);
}

static
void _read_sky_profiles(Model &model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "sky", "bg", &list_to_sky);
}

static
void _read_psf_profiles(Model &model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "psf", "xcen", &list_to_psf);
}


/*
 * OpenCL-related functionality follows
 * ----------------------------------------------------------------------------
 */
static
SEXP _R_profit_openclenv_info() {

	map<int, OpenCL_plat_info> clinfo;
	try {
		clinfo = get_opencl_info();
	} catch (const exception &e) {
		ostringstream os;
		os << "Error while querying OpenCL environment: " << e.what();
		Rf_error(os.str().c_str());
		return R_NilValue;
	}

	unsigned int protections = 0;
	SEXP r_platinfo_names = PROTECT(Rf_allocVector(STRSXP, 3));
	SEXP r_devinfo_names = PROTECT(Rf_allocVector(STRSXP, 2));
	SET_STRING_ELT(r_platinfo_names, 0, Rf_mkChar("name"));
	SET_STRING_ELT(r_platinfo_names, 1, Rf_mkChar("opencl_version"));
	SET_STRING_ELT(r_platinfo_names, 2, Rf_mkChar("devices"));
	SET_STRING_ELT(r_devinfo_names, 0, Rf_mkChar("name"));
	SET_STRING_ELT(r_devinfo_names, 1, Rf_mkChar("supports_double"));
	protections += 2;

	SEXP r_clinfo = PROTECT(Rf_allocVector(VECSXP, clinfo.size()));
	protections += 1;
	unsigned int plat = 0;
	for(auto platform_info: clinfo) {

		auto plat_info = std::get<1>(platform_info);

		unsigned int dev = 0;
		SEXP r_devsinfo = PROTECT(Rf_allocVector(VECSXP, plat_info.dev_info.size()));
		protections += 1;
		for(auto device_info: plat_info.dev_info) {

			auto dev_info = std::get<1>(device_info);
			SEXP r_double_support = PROTECT(Rf_ScalarLogical(dev_info.double_support ? TRUE : FALSE));
			SEXP r_dev_name = PROTECT(Rf_mkString(dev_info.name.c_str()));
			SEXP r_devinfo = PROTECT(Rf_allocVector(VECSXP, 2));
			Rf_setAttrib(r_devinfo, R_NamesSymbol, r_devinfo_names);
			SET_VECTOR_ELT(r_devinfo, 0, r_dev_name);
			SET_VECTOR_ELT(r_devinfo, 1, r_double_support);
			SET_VECTOR_ELT(r_devsinfo, dev++, r_devinfo);
			protections += 3;
		}

		SEXP r_plat_clver = PROTECT(Rf_ScalarReal(plat_info.supported_opencl_version/100.));
		SEXP r_plat_name = PROTECT(Rf_mkString(plat_info.name.c_str()));
		SEXP r_platinfo = PROTECT(Rf_allocVector(VECSXP, 3));
		Rf_setAttrib(r_platinfo, R_NamesSymbol, r_platinfo_names);
		SET_VECTOR_ELT(r_platinfo, 0, r_plat_name);
		SET_VECTOR_ELT(r_platinfo, 1, r_plat_clver);
		SET_VECTOR_ELT(r_platinfo, 2, r_devsinfo);
		protections += 3;

		SET_VECTOR_ELT(r_clinfo, plat++, r_platinfo);
	}

	UNPROTECT(protections);
	return r_clinfo;
}

struct openclenv_wrapper {
	OpenCLEnvPtr env;
};

static
void _R_profit_openclenv_finalizer(SEXP ptr) {

	if(!R_ExternalPtrAddr(ptr)) {
		return;
	}

	openclenv_wrapper *wrapper = reinterpret_cast<openclenv_wrapper *>(R_ExternalPtrAddr(ptr));
	wrapper->env.reset();
	delete wrapper;
	R_ClearExternalPtr(ptr); /* not really needed */

}

static
OpenCLEnvPtr unwrap_openclenv(SEXP openclenv) {

	if( TYPEOF(openclenv) != EXTPTRSXP ) {
		Rf_error("Given openclenv not of proper type\n");
		return nullptr;
	}

	openclenv_wrapper *wrapper = reinterpret_cast<openclenv_wrapper *>(R_ExternalPtrAddr(openclenv));
	if( !wrapper ) {
		Rf_error("No OpenCL environment found in openclenv\n");
		return nullptr;
	}

	return wrapper->env;
}

static
SEXP _R_profit_openclenv(SEXP plat_idx, SEXP dev_idx, SEXP use_dbl) {

	unsigned int platform_idx = INTEGER(plat_idx)[0];
	unsigned int device_idx = INTEGER(dev_idx)[0];
	bool use_double = static_cast<bool>(INTEGER(use_dbl)[0]);

	OpenCLEnvPtr env;
	try {
		env = get_opencl_environment(platform_idx, device_idx, use_double, false);
	} catch (const opencl_error &e) {
		ostringstream os;
		os << "Error while creating OpenCL environment for plat/dev/double " << platform_idx << "/" << device_idx << "/" << use_double << ": " << e.what();
		Rf_error(os.str().c_str());
		return R_NilValue;
	} catch (const invalid_parameter &e) {
		ostringstream os;
		os << "Error while creating OpenCL environment, invalid parameter: " << e.what();
		Rf_error(os.str().c_str());
		return R_NilValue;
	}

	openclenv_wrapper *wrapper = new openclenv_wrapper();
	wrapper->env = env;
	SEXP r_openclenv = R_MakeExternalPtr(wrapper, Rf_install("OpenCL_env"), R_NilValue);
	PROTECT(r_openclenv);
	R_RegisterCFinalizerEx(r_openclenv, _R_profit_openclenv_finalizer, TRUE);
	UNPROTECT(1);
	return r_openclenv;
}

/*
 * OpenMP-related functionality follows
 * ----------------------------------------------------------------------------
 */
static
SEXP _R_profit_has_openmp() {
	return Rf_ScalarLogical(profit::has_openmp());
}

/*
 * FFTW-related functionality follows
 * ----------------------------------------------------------------------------
 */
static
SEXP _R_profit_has_fftw() {
	return Rf_ScalarLogical(profit::has_fftw());
}

/*
 * Convolver exposure support follows
 * ----------------------------------------------------------------------------
 */
struct convolver_wrapper {
	ConvolverPtr convolver;
};

static
void _R_profit_convolver_finalizer(SEXP ptr) {

	if(!R_ExternalPtrAddr(ptr)) {
		return;
	}

	convolver_wrapper *wrapper = reinterpret_cast<convolver_wrapper *>(R_ExternalPtrAddr(ptr));
	wrapper->convolver.reset();
	delete wrapper;
	R_ClearExternalPtr(ptr); /* not really needed */

}

static
ConvolverPtr unwrap_convolver(SEXP convolver)
{
	if( TYPEOF(convolver) != EXTPTRSXP ) {
		Rf_error("Given convolver not of proper type\n");
		return nullptr;
	}

	convolver_wrapper *wrapper = reinterpret_cast<convolver_wrapper *>(R_ExternalPtrAddr(convolver));
	if (!wrapper) {
		Rf_error("No Convolver found in convolver object");
		return nullptr;
	}

	return wrapper->convolver;
}

static
SEXP _R_profit_convolvers()
{
	vector<string> convolvers {"brute"};
	if (profit::has_opencl()) {
		convolvers.push_back("opencl");
	}
	if (profit::has_fftw()) {
		convolvers.push_back("fftw");
	}

	const size_t n_convolvers = convolvers.size();
	SEXP convolvers_r = PROTECT(Rf_allocVector(STRSXP, n_convolvers));
	for(size_t i = 0; i != n_convolvers; i++) {
		SET_STRING_ELT(convolvers_r, i, Rf_mkChar(convolvers[i].c_str()));
	}

	UNPROTECT(1);
	return convolvers_r;
}

static
SEXP _R_profit_make_convolver(SEXP type, SEXP image_dimensions, SEXP psf,
                              SEXP reuse_psf_fft, SEXP fft_effort, SEXP omp_threads,
                              SEXP openclenv)
{

	ConvolverCreationPreferences conv_prefs;
	conv_prefs.src_dims = {(unsigned int)INTEGER(image_dimensions)[0], (unsigned int)INTEGER(image_dimensions)[1]};
	conv_prefs.krn_dims = _read_image(psf).getDimensions();

	if( omp_threads != R_NilValue ) {
		conv_prefs.omp_threads = (unsigned int)Rf_asInteger(omp_threads);
	}
	if (reuse_psf_fft != R_NilValue ) {
		conv_prefs.reuse_krn_fft = (bool)Rf_asLogical(reuse_psf_fft);
	}
	if (fft_effort != R_NilValue) {
		conv_prefs.effort = effort_t((unsigned int)Rf_asInteger(fft_effort));
	}
	if (openclenv != R_NilValue ) {
		if((conv_prefs.opencl_env = unwrap_openclenv(openclenv)) == nullptr) {
			return R_NilValue;
		}
	}

	convolver_wrapper *wrapper = new convolver_wrapper();
	std::string error;
	try {
		wrapper->convolver = create_convolver(CHAR(STRING_ELT(type, 0)), conv_prefs);
	} catch (std::exception &e) {
		Rf_error(e.what());
		return R_NilValue;
	}

	SEXP r_convolver = R_MakeExternalPtr(wrapper, Rf_install("Convolver"), R_NilValue);
	PROTECT(r_convolver);
	R_RegisterCFinalizerEx(r_convolver, _R_profit_convolver_finalizer, TRUE);
	UNPROTECT(1);
	return r_convolver;
}

/*
 * Public exported functions follow now
 * ----------------------------------------------------------------------------
 */
static
SEXP _R_profit_make_model(SEXP model_list) {

	unsigned int img_w, img_h;
	double scale_x = 1, scale_y = 1;
	string error;

	SEXP dimensions = _get_list_element(model_list, "dimensions");
	img_w = (unsigned int)INTEGER(dimensions)[0];
	img_h = (unsigned int)INTEGER(dimensions)[1];

	SEXP magzero = _get_list_element(model_list, "magzero");
	if( magzero == R_NilValue ) {
		Rf_error("No magzero provided in the model\n");
		return R_NilValue;
	}

	SEXP r_scale_x = _get_list_element(model_list, "scale_x");
	if( r_scale_x != R_NilValue ) {
		scale_x = Rf_asReal(r_scale_x);
	}
	SEXP r_scale_y = _get_list_element(model_list, "scale_y");
	if( r_scale_y != R_NilValue ) {
		scale_y = Rf_asReal(r_scale_y);
	}

	/* Read model parameters */
	Model m {img_w, img_h};
	m.set_image_pixel_scale({scale_x, scale_y});
	m.set_magzero(Rf_asReal(magzero));
	m.set_finesampling(static_cast<unsigned int>(Rf_asInteger(_get_list_element(model_list, "finesampling"))));
	m.set_return_finesampled(static_cast<bool>(Rf_asLogical(_get_list_element(model_list, "returnfine"))));
	m.set_crop(static_cast<bool>(Rf_asLogical(_get_list_element(model_list, "crop"))));

	SEXP r_psf = _get_list_element(model_list, "psf");
	if( r_psf != R_NilValue ) {
		m.set_psf(_read_image(r_psf));
		m.set_psf_pixel_scale({scale_x, scale_y});
	}

	SEXP r_calcregion = _get_list_element(model_list, "calcregion");
	if( r_calcregion != R_NilValue ) {
		auto mask = _read_mask(r_calcregion);
		if (mask.getDimensions() != Dimensions{img_w, img_h}) {
			Rf_error("Calc region has different dimensions than the image");
			return R_NilValue;
		}
		m.set_mask(std::move(mask));
	}

	/* A convolver, if any */
	SEXP r_convolver = _get_list_element(model_list, "convolver");
	if (r_convolver != R_NilValue) {
		auto convolver = unwrap_convolver(r_convolver);
		if (!convolver) {
			return R_NilValue;
		}
		m.set_convolver(convolver);
	}

	/* OpenCL environment, if any */
	SEXP openclenv = _get_list_element(model_list, "openclenv");
	if( openclenv != R_NilValue ) {

		if( TYPEOF(openclenv) != EXTPTRSXP ) {
			Rf_error("Given openclenv not of proper type\n");
			return R_NilValue;
		}

		openclenv_wrapper *wrapper = reinterpret_cast<openclenv_wrapper *>(R_ExternalPtrAddr(openclenv));
		if( !wrapper ) {
			Rf_error("No OpenCL environment found in openclenv\n");
			return R_NilValue;
		}

		m.set_opencl_env(wrapper->env);
	}

	/* Number of OpenMP threads, if any */
	SEXP omp_threads = _get_list_element(model_list, "omp_threads");
	if( omp_threads != R_NilValue ) {
		m.set_omp_threads((unsigned int)Rf_asInteger(omp_threads));
	}

	/* Read profiles and parameters and append them to the model */
	SEXP profiles = _get_list_element(model_list, "profiles");
	if( profiles == R_NilValue ) {
		Rf_error("No profiles provided in the model\n");
		return R_NilValue;
	}
	_read_sersic_profiles(m, profiles);
	_read_moffat_profiles(m, profiles);
	_read_ferrer_profiles(m, profiles);
	_read_coresersic_profiles(m, profiles);
	_read_king_profiles(m, profiles);
	_read_brokenexponential_profiles(m, profiles);
	_read_sky_profiles(m, profiles);
	_read_psf_profiles(m, profiles);
	if( !m.has_profiles() ) {
		Rf_error("No valid profiles found in profiles list\n");
		return R_NilValue;
	}

	/* Go, go, go! */
	Image image;
	Point offset;
	try {
		image = m.evaluate(offset);
	} catch (const std::exception &e) {
		stringstream ss;
		ss << "Error while calculating model: " << e.what() << endl;
		Rf_error(ss.str().c_str());
		return R_NilValue;
	}

	/* Copy the image into a matrix, pack together with the offset, and bye */
	SEXP r_image = PROTECT(Rf_allocMatrix(REALSXP, image.getWidth(), image.getHeight()));
	SEXP r_offset = PROTECT(Rf_list2(Rf_ScalarInteger(offset.x + 1), Rf_ScalarInteger(offset.y + 1)));
	memcpy(REAL(r_image), image.data(), sizeof(double) * image.size());
	SEXP ret = PROTECT(Rf_list2(r_image, r_offset));
	UNPROTECT(3);
	return ret;
}

static
SEXP _R_profit_convolve(SEXP r_convolver, SEXP r_image, SEXP r_psf, SEXP r_mask) {

	ConvolverPtr convolver = unwrap_convolver(r_convolver);
	if (!convolver) {
		return R_NilValue;
	}

	auto src = _read_image(r_image);
	auto psf = _read_image(r_psf);
	Mask mask;
	if (r_mask != R_NilValue) {
		mask = _read_mask(r_mask);
	}

	auto image = convolver->convolve(src, psf, mask);
	SEXP ret_image = PROTECT(Rf_allocVector(REALSXP, image.size()));
	memcpy(REAL(ret_image), image.data(), sizeof(double) * image.size());

	UNPROTECT(1);
	return ret_image;

}

static
SEXP _R_profit_downsample(SEXP img, SEXP factor)
{
	auto image = _read_image(img);
	if (!image) {
		return R_NilValue;
	}

	try {
		image = image.downsample((unsigned int)Rf_asInteger(factor), Image::DownsamplingMode::SUM);
	} catch (const std::exception &e) {
		stringstream ss;
		ss << "Error while downsampling image: " << e.what() << endl;
		Rf_error(ss.str().c_str());
		return R_NilValue;
	}

	/* Copy the image, clean up, and good bye */
	SEXP r_image = PROTECT(Rf_allocMatrix(REALSXP, image.getWidth(), image.getHeight()));
	memcpy(REAL(r_image), image.data(), sizeof(double) * image.size());

	UNPROTECT(1);
	return r_image;
}

static
SEXP _R_profit_upsample(SEXP img, SEXP factor)
{
	auto image = _read_image(img);
	if (!image) {
		return R_NilValue;
	}

	try {
		image = image.upsample((unsigned int)Rf_asInteger(factor), Image::UpsamplingMode::COPY);
	} catch (const std::exception &e) {
		stringstream ss;
		ss << "Error while downsampling image: " << e.what() << endl;
		Rf_error(ss.str().c_str());
		return R_NilValue;
	}

	/* Copy the image, clean up, and good bye */
	SEXP r_image = PROTECT(Rf_allocMatrix(REALSXP, image.getWidth(), image.getHeight()));
	memcpy(REAL(r_image), image.data(), sizeof(double) * image.size());

	UNPROTECT(1);
	return r_image;
}

// TODO: to be moved down to libprofit at some point
namespace detail {

int setenv(const char *name, const char *value)
{
#ifdef _WIN32
	return ::_putenv_s(name, value);
#else
	return ::setenv(name, value, 1);
#endif
}

} // namespace detail

extern "C" {
	SEXP R_profit_make_model(SEXP model_list) {
		return _R_profit_make_model(model_list);
	}

	SEXP R_profit_convolvers() {
		return _R_profit_convolvers();
	}

	SEXP R_profit_make_convolver(SEXP type, SEXP image_dimensions, SEXP psf,
	                             SEXP reuse_psf_fft, SEXP fft_effort, SEXP omp_threads,
	                             SEXP openclenv) {
		return _R_profit_make_convolver(type, image_dimensions, psf, reuse_psf_fft,
		                                fft_effort, omp_threads, openclenv);
	}

	SEXP R_profit_convolve(SEXP convolver, SEXP r_image, SEXP r_psf, SEXP mask) {
		return _R_profit_convolve(convolver, r_image, r_psf, mask);
	}

	SEXP R_profit_has_openmp() {
		return _R_profit_has_openmp();
	}

	SEXP R_profit_has_fftw() {
		return _R_profit_has_fftw();
	}

	SEXP R_profit_openclenv_info() {
		return _R_profit_openclenv_info();
	}

	SEXP R_profit_openclenv(SEXP plat_idx, SEXP dev_idx, SEXP use_dbl) {
		return _R_profit_openclenv(plat_idx, dev_idx, use_dbl);
	}

	SEXP R_profit_downsample(SEXP img, SEXP factor) {
		return _R_profit_downsample(img, factor);
	}

	SEXP R_profit_upsample(SEXP img, SEXP factor) {
		return _R_profit_upsample(img, factor);
	}

	SEXP R_profit_clear_cache() {
		profit::clear_cache();
		return R_NilValue;
	}

	/*
	 * Registering the methods above at module loading time
	 * This should speed symbol lookup, and anyway it's considered a better
	 * practice.
	 */
	static const R_CallMethodDef callMethods[]  = {

		/* Defined in this module */
		{"R_profit_make_model",     (DL_FUNC) &R_profit_make_model,     1},
		{"R_profit_convolvers",     (DL_FUNC) &R_profit_convolvers,     0},
		{"R_profit_make_convolver", (DL_FUNC) &R_profit_make_convolver, 7},
		{"R_profit_convolve",       (DL_FUNC) &R_profit_convolve,       4},
		{"R_profit_has_openmp",     (DL_FUNC) &R_profit_has_openmp,     0},
		{"R_profit_has_fftw",       (DL_FUNC) &R_profit_has_fftw,       0},
		{"R_profit_openclenv_info", (DL_FUNC) &R_profit_openclenv_info, 0},
		{"R_profit_openclenv",      (DL_FUNC) &R_profit_openclenv,      3},
		{"R_profit_downsample",     (DL_FUNC) &R_profit_downsample,     2},
		{"R_profit_upsample",       (DL_FUNC) &R_profit_upsample,       2},
		{"R_profit_clear_cache",    (DL_FUNC) &R_profit_clear_cache,    0},

		/* Sentinel */
		{NULL, NULL, 0}
	};

	void R_init_ProFit(DllInfo *dll) {

		/*
		 * Changing PROFIT_HOME's environment variable to point to R's session
		 * temporary directory. It was originally suggested that this could be
		 * done in the .onLoad() function, but that would be too late, as this
		 * shared library would have been loaded already, which in turn means
		 * that the profit::init function (the one actually creating the
		 * PROFIT_HOME directory if non-existing) would have executed already.
		 *
		 * This is to avoid violating CRAN's policy that packages should not
		 * write anywhere on the filesystem apart from the R session's tmp dir
		 * Obviously this completely misses the point of having a cache in the
		 * first place, so in the future we'll have to revisit this and figure
		 * out how to cache things in R.
		 *
		 * For the life of me, I also couldn't find an entry point in R's C API
		 * that would return the R session's temporary directory. This is
		 * required by both R_tmpnam and R_tmpnam2, which surprisingly do not
		 * calculate this themselves. The two solutions I could find were: a)
		 * using the R_TempDir symbol declared in Rembedded.h (but is not meant
		 * to be used, it issues a check warning), or b) call the tempdir R
		 * function, which returns exactly what we need. Here I implemented the
		 * latter.
		 */
		SEXP tempdir_func, r_session_tmpdir;
		PROTECT(tempdir_func = Rf_lang1(Rf_install("tempdir")));
		PROTECT(r_session_tmpdir = Rf_eval(tempdir_func, R_GlobalEnv));
		auto tmpdir = R_tmpnam("profit", CHAR(STRING_ELT(r_session_tmpdir, 0)));
		detail::setenv("PROFIT_HOME", tmpdir);
		free(tmpdir);
		UNPROTECT(2);

		auto success = profit::init();
		if (!success) {
			std::ostringstream os;
			os << "Error when initializing libprofit: " << profit::init_diagnose();
			Rf_error(os.str().c_str());
		}

		/* Using registered symbols only from now on */
		R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
		R_useDynamicSymbols(dll, FALSE);
	}

	void R_unload_ProFit(DllInfo *info) {
		profit::finish();
		auto diagnose = profit::finish_diagnose();
		if (!diagnose.empty()) {
			std::ostringstream os;
			os << "Warning when finishing libprofit: " << diagnose;
			Rf_warning(os.str().c_str());
		}
	}

}
