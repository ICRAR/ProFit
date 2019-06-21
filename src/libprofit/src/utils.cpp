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
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <cerrno>
#include <functional>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef _WIN32
# include <windows.h>
#else
# include <dirent.h>
# include <sys/stat.h>
# include <sys/time.h>
# include <sys/types.h>
# include <unistd.h>
#endif // _WIN32

#include "profit/common.h"
#include "profit/config.h"
#include "profit/exceptions.h"
#include "profit/utils.h"

/*
 * We use either GSL or R to provide the low-level
 * beta, gamma and pgamma and qgamma functions needed by some profiles.
 * If neither is given the compilation should fail
 *
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
#if defined(PROFIT_USES_GSL)
#  include <gsl/gsl_errno.h>
#  include <gsl/gsl_cdf.h>
#  include <gsl/gsl_sf_gamma.h>
#  include <gsl/gsl_integration.h>
#elif defined(PROFIT_USES_R)
#  include <Rmath.h>
#  include <R_ext/Applic.h>
#  if defined(MATHLIB_STANDALONE)
#    define Rf_qgamma   qgamma
#    define Rf_pgamma   pgamma
#    define Rf_gammafn  gammafn
#    define Rf_beta     beta
#  else
#    undef qgamma
#    undef pgamma
#    undef gammafn
#    undef beta
#  endif
#else
#  error("No high-level library (GSL or R) provided")
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
	double epsabs = 1e-4;
	double epsrel = 1e-4;
	double result;
	double abserr;

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
#endif // defined(PROFIT_USES_GSL) / defined(PROFIT_USES_R)

#ifdef _WIN32
static inline
bool path_exists(const std::string &path, bool dir_expected)
{
	auto attrs = GetFileAttributes(path.c_str());
	if (attrs == INVALID_FILE_ATTRIBUTES) {

		auto last_error = GetLastError();
		if (last_error == ERROR_FILE_NOT_FOUND || last_error == ERROR_PATH_NOT_FOUND ||
		    last_error == ERROR_INVALID_NAME || last_error == ERROR_INVALID_DRIVE ||
		    last_error == ERROR_BAD_PATHNAME) {
			return false;
		}

		std::ostringstream os;
		os << "Unexpected error found when inspecting " << path << ": " << last_error;
		throw std::runtime_error(os.str());
	}

	bool is_dir = attrs & FILE_ATTRIBUTE_DIRECTORY;
	return is_dir == dir_expected;
}
#else
static inline
bool inode_exists(const std::string &fname, mode_t expected_type, const char *type_name) {

	struct ::stat st;
	int result = ::stat(fname.c_str(), &st);

	if (result == -1) {

		// it doesn't exist
		if (errno == ENOENT) {
			return false;
		}

		// Another kind of unexpected error
		std::ostringstream os;
		os << "Unexpected error found when inspecting " << fname << ": ";
		os << strerror(errno);
		throw fs_error(os.str());
	}


	// it exists, but is it a the expected type
	bool is_expected = (st.st_mode & S_IFMT) == expected_type;
	if (!is_expected) {
		std::ostringstream os;
		os << fname << " exists but is not a " << type_name << ". Please remove it and try again";
		throw fs_error(os.str());
	}

	return true;
}
#endif // _WIN32

bool dir_exists(const std::string &fname)
{
#ifdef _WIN32
	return path_exists(fname, true);
#else
	return inode_exists(fname, S_IFDIR, "directory");
#endif // _WIN32
}

bool file_exists(const std::string &fname)
{
#ifdef _WIN32
	return path_exists(fname, false);
#else
	return inode_exists(fname, S_IFREG, "regular file");
#endif // _WIN32
}

static inline
void create_dir(const std::string &fname) {
#ifdef _WIN32
	CreateDirectory(fname.c_str(), NULL);
#else
	// mkdir with 755 permissions
	::mkdir(fname.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
#endif // _WIN32
}

std::string create_dirs(const std::string &at, const std::vector<std::string> &parts)
{
	std::string the_dir = at;
	for(auto &part: parts) {
#ifdef _WIN32
		the_dir += "\\" + part;
#else
		the_dir += "/" + part;
#endif
		if (!dir_exists(the_dir)) {
			create_dir(the_dir);
		}
	}
	return the_dir;
}

static
fs_error _removal_error(const char *path)
{
	std::ostringstream os;
	os << "Unexpected error found when removing " << path << ": ";
#ifdef _WIN32
	auto err = GetLastError();
	LPTSTR errormsg = 0;
	FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM, nullptr, err, 0, (LPTSTR)&errormsg, 0, NULL);
	os << errormsg;
	LocalFree(errormsg);
#else
	os << errno << " (" << strerror(errno) << ")";
#endif
	return fs_error(os.str());
}

static
void _recursive_remove(const char *path)
{
#ifdef _WIN32
	std::string pattern(path);
	pattern.append("\\*");

	WIN32_FIND_DATA data;
	HANDLE find_handle = FindFirstFile(pattern.c_str(), &data);
	if (find_handle == INVALID_HANDLE_VALUE) {
		throw _removal_error(path);
	}

	do {

		if (!strcmp(".", data.cFileName) || !strcmp("..", data.cFileName)) {
			continue;
		}

		if (data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
			std::ostringstream full_path;
			full_path << path << "\\" << data.cFileName;
			_recursive_remove(full_path.str().c_str());
		}
		else {
			std::ostringstream full_path;
			full_path << path << "\\" << data.cFileName;
			if (!DeleteFile(full_path.str().c_str())) {
				throw _removal_error(full_path.str().c_str());
			}
		}
	} while (FindNextFile(find_handle, &data) != 0);

	FindClose(find_handle);

	if (!RemoveDirectory(path)) {
		throw _removal_error(path);
	}
#else

	struct ::stat st;
	int result = ::stat(path, &st);
	if (result == -1) {
		throw _removal_error(path);
	}

	auto mode = st.st_mode & S_IFMT;
	if (mode == S_IFDIR) {

		// We open the directory, recursively remove its contents
		// (making sure we skip . and ..), close it, and finally remove it
		DIR *dir;
		if ((dir = ::opendir(path)) == nullptr) {
			throw _removal_error(path);
		}

		struct dirent *ent;
		while ((ent = ::readdir (dir)) != nullptr) {

			if (!strcmp(".", ent->d_name) || !strcmp("..", ent->d_name)) {
				continue;
			}

			std::ostringstream full_path;
			full_path << path << "/" << ent->d_name;
			_recursive_remove(full_path.str().c_str());
		}

		if (::closedir(dir) == -1) {
			throw _removal_error(path);
		}

		if (::rmdir(path) == -1) {
			throw _removal_error(path);
		}

		return;
	}

	// No recursion needed, just remove
	auto ret = ::unlink(path);
	if (ret == -1) {
		throw _removal_error(path);
	}
#endif // _WIN32
}

void recursive_remove(const std::string &path)
{
	_recursive_remove(path.c_str());
}


std::string get_profit_home()
{
	auto *profit_home = std::getenv("PROFIT_HOME");
	if (profit_home) {
		if (!dir_exists(profit_home)) {
			create_dir(profit_home);
		}
		return profit_home;
	}

#ifdef _WIN32
	constexpr const char *home_var = "APPDATA";
	constexpr const char *profit_basedir = "profit";
#else
	constexpr const char *home_var = "HOME";
	constexpr const char *profit_basedir = ".profit";
#endif // _WIN32

	auto *user_home = std::getenv(home_var);
	if (!user_home) {
		throw exception("User doesn't have a home :(");
	}

	return create_dirs(user_home, {std::string(profit_basedir)});
}

void setenv(const std::string &name, const std::string &value)
{
#ifdef _WIN32
	::_putenv_s(name.c_str(), value.c_str());
#else
	if (!value.empty()) {
		::setenv(name.c_str(), value.c_str(), 1);
	}
	else {
		::unsetenv(name.c_str());
	}
#endif // _WIN32
}

// Adapted from: http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
std::vector<std::string> split(const std::string &s, const std::string &delims)
{
	auto lastPos = s.find_first_not_of(delims, 0);
	auto pos = s.find_first_of(delims, lastPos);

	std::vector<std::string> tokens;
	while (std::string::npos != pos || std::string::npos != lastPos) {
		tokens.push_back(s.substr(lastPos, pos - lastPos));
		lastPos = s.find_first_not_of(delims, pos);
		pos = s.find_first_of(delims, lastPos);
	}
	return tokens;
}

// trim from start
static inline std::string &ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(),
	        std::not1(std::ptr_fun<int, int>(std::isspace))));
	return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(),
	        std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	return s;
}

// trim from both ends
std::string &trim(std::string &s) {
	return ltrim(rtrim(s));
}

std::string trim(const std::string &s) {
	std::string s_copy(s);
	ltrim(rtrim(s_copy));
	return s_copy;
}

} /* namespace profit */
