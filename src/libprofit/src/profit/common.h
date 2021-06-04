/**
 * Header file with common definitions for libprofit
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

#ifndef PROFIT_COMMON_H
#define PROFIT_COMMON_H

#include <cmath>
#include <chrono>
#include <iosfwd>

/* M_PI is not part of C/C++, but usually there */
#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

/* The override keyword is not supported until gcc 4.7 */
#if defined(__GNUG__) && (__GNUC__ <= 4) && (__GNUC_MINOR__ < 7)
# define override
#endif

// Proper function exporting/importing under Windows
#define PROFIT_API

/* Sometimes we don't use all arguments */
#define UNUSED(x) do { (void)x; } while(0)

namespace profit {

	/// A type to hold nanosecond values
	typedef std::chrono::nanoseconds::rep nsecs_t;

	/// Trait naming different types
	template <typename T>
	struct type_info {};

	template <>
	struct type_info<bool> {
		constexpr static const char *name = "bool";
	};

	template <>
	struct type_info<unsigned int> {
		constexpr static const char *name = "unsigned int";
	};

	template <>
	struct type_info<float> {
		constexpr static const char *name = "float";
	};

	template <>
	struct type_info<double> {
		constexpr static const char *name = "double";
	};

	/// Trait describing specific float and double floating types
	template <typename FT>
	struct float_traits {
		const static bool is_float = false;
		const static bool is_double = false;
	};

	template <>
	struct float_traits<float> : type_info<float> {
		const static bool is_float = true;
		const static bool is_double = false;
	};

	template <>
	struct float_traits<double> : type_info<double> {
		const static bool is_float = false;
		const static bool is_double = true;
	};

	/**
	 * SIMD instruction sets choosers can choose from
	 */
	enum simd_instruction_set {
		/// Automatically choose the best available SIMD instruction set
		AUTO = 0,
		/// No SIMD instruction set
		NONE,
		/// The SSE2 instruction set
		SSE2,
		/// The AVX instruction set
		AVX
	};

	template <typename T, typename CharT>
	std::basic_ostream<T, CharT> &operator<<(std::basic_ostream<T, CharT> &os, simd_instruction_set instruction_set) {
		if (instruction_set == AUTO) {
			os << "AUTO";
		}
		else if (instruction_set == NONE) {
			os << "NONE";
		}
		else if (instruction_set == SSE2) {
			os << "SSE2";
		}
		else if (instruction_set == AVX) {
			os << "AVX";
		}
		else {
			os << "unknown";
		}
		return os;
	}

}

#endif /* PROFIT_COMMON_H */
