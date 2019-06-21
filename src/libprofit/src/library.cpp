/**
 * Implementation of library-related routines for libprofit
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2017
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

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <sstream>

#include "profit/config.h"
#include "profit/library.h"
#include "profit/utils.h"
#include "profit/fft_impl.h"

namespace profit {

static std::string _version = std::to_string(PROFIT_VERSION_MAJOR) + "." +
                              std::to_string(PROFIT_VERSION_MINOR) + "." +
                              std::to_string(PROFIT_VERSION_PATCH) +
                              (std::string("") != PROFIT_VERSION_SUFFIX ? std::string("-") + PROFIT_VERSION_SUFFIX : "");

std::string version()
{
	return _version;
}

unsigned short version_major()
{
	return PROFIT_VERSION_MAJOR;
}

unsigned short version_minor()
{
	return PROFIT_VERSION_MINOR;
}

unsigned short version_patch()
{
	return PROFIT_VERSION_PATCH;
}

std::string version_suffix()
{
	return PROFIT_VERSION_SUFFIX;
}

#ifdef PROFIT_FFTW
static inline
std::string get_fftw_wisdom_filename()
{
	auto fftw_cache_dir = create_dirs(get_profit_home(), {std::string("fftw_cache")});
#ifdef PROFIT_FFTW_OPENMP
	auto fftw_wisdom_fname = fftw_cache_dir + "/threaded-wisdom";
#else
	auto fftw_wisdom_fname = fftw_cache_dir + "/unthreaded-wisdom";
#endif
	return fftw_wisdom_fname + "_" + fftw_version;
}
#endif // PROFIT_FFTW

static std::string _init_diagnose;
static std::string _finish_diagnose;

std::string init_diagnose()
{
	return _init_diagnose;
}

std::string finish_diagnose()
{
	return _finish_diagnose;
}

bool init()
{
	// Initialize FFTW library, including its OpenMP support
	// It is important to configure the OpenMP support before reading the wisdom,
	// otherwise the plans will fail to import
#ifdef PROFIT_FFTW
	std::lock_guard<std::mutex> guard(fftw_mutex);
#ifdef PROFIT_FFTW_OPENMP
	int res = fftw_init_threads();
	if (!res) {
		std::ostringstream os;
		os << "Error while initializing FFTW threads support, errno = " << res;
		_init_diagnose = os.str();
		return false;
	}
#endif // PROFIT_FFTW_OPENMP

	auto fftw_wisdom_filename = get_fftw_wisdom_filename();
	if (file_exists(fftw_wisdom_filename)) {
		auto *fftw_wisdom_file = fopen(fftw_wisdom_filename.c_str(), "r");
		if (!fftw_wisdom_file) {
			std::ostringstream os;
			os << "Opening fftw wisdom from " << fftw_wisdom_filename << " failed: " << strerror(errno);
			_init_diagnose = os.str();
		}
		else {
			auto import_status = fftw_import_wisdom_from_file(fftw_wisdom_file);
			if (import_status == 0) {
				std::ostringstream os;
				os << "Importing fftw wisdom from " << fftw_wisdom_filename << " failed: " << import_status;
				_init_diagnose = os.str();
			}
			fclose(fftw_wisdom_file);
		}
	}

#if !(defined(__WIN32__) || defined(WIN32) || defined(_WINDOWS))
	if (file_exists("/etc/fftw/wisdom") && fftw_import_system_wisdom() == 0) {
		std::ostringstream os;
		os << _init_diagnose << '\n';
		os << "Importing fftw system wisdom failed (returned 0)";
		_init_diagnose = os.str();
	}
#endif // !WIN32
#endif // PROFIT_FFTW

	return true;
}

void finish()
{
#ifdef PROFIT_FFTW
	std::lock_guard<std::mutex> guard(fftw_mutex);
	auto fftw_wisdom_filename = get_fftw_wisdom_filename();
	auto *fftw_wisdom_file = fopen(fftw_wisdom_filename.c_str(), "w");
	if (!fftw_wisdom_file) {
		std::ostringstream os;
		os << "Error when exporting fftw wisdom from " << fftw_wisdom_file << ": " << strerror(errno);
		_finish_diagnose = os.str();
	}
	else {
		fftw_export_wisdom_to_file(fftw_wisdom_file);
		fclose(fftw_wisdom_file);
	}
#ifdef PROFIT_FFTW_OPENMP
	fftw_cleanup_threads();
#endif /* PROFIT_FFTW_OPENMP */
	fftw_cleanup();
#endif // PROFIT_FFTW
}

bool has_openmp()
{
#ifdef PROFIT_OPENMP
	return true;
#else
	return false;
#endif // PROFIT_OPENMP
}

bool has_fftw()
{
#ifdef PROFIT_FFTW
	return true;
#else
	return false;
#endif // PROFIT_FFTW
}

bool has_fftw_with_openmp()
{
#ifdef PROFIT_FFTW_OPENMP
	return true;
#else
	return false;
#endif
}

bool has_opencl()
{
#ifdef PROFIT_OPENCL
	return true;
#else
	return false;
#endif // PROFIT_OPENCL
}

unsigned short opencl_version_major()
{
#ifdef PROFIT_OPENCL
	return PROFIT_OPENCL_MAJOR;
#else
	return 0;
#endif // PROFIT_OPENCL
}

unsigned short opencl_version_minor()
{
#ifdef PROFIT_OPENCL
	return PROFIT_OPENCL_MINOR;
#else
	return 0;
#endif // PROFIT_OPENCL
}

bool has_simd_instruction_set(simd_instruction_set instruction_set)
{
	if (instruction_set == simd_instruction_set::AUTO || instruction_set == simd_instruction_set::NONE) {
		return true;
	}

	// AVX implies SSE2
#ifdef PROFIT_HAS_SSE2
	if (instruction_set == simd_instruction_set::SSE2) {
		return true;
	}
#endif // PROFIT_HAS_SSE2
#ifdef PROFIT_HAS_AVX
	if (instruction_set == simd_instruction_set::AVX) {
		return true;
	}
#endif // PROFIT_HAS_AVX

	return false;
}

bool has_avx()
{
#ifdef PROFIT_HAS_AVX
	return true;
#else
	return false;
#endif // PROFIT_HAS_AVX
}

void clear_cache()
{
	auto profit_home = get_profit_home();

#ifdef PROFIT_FFTW
	std::lock_guard<std::mutex> guard(fftw_mutex);
	fftw_forget_wisdom();
	auto fftw_cache = profit_home + "/fftw_cache";
	if (dir_exists(fftw_cache)) {
		recursive_remove(fftw_cache);
	}
#endif

#ifdef PROFIT_OPENCL
	auto opencl_cache = profit_home + "/opencl_cache";
	if (dir_exists(opencl_cache)) {
		recursive_remove(opencl_cache);
	}
#endif
}

} // namespace profit