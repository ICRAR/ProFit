/**
 * Dot product function implementations for libprofit
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Aaron Robotham, Dan Taranu, Rodrigo Tobar
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

#ifndef PROFIT_DOT_PRODUCT_H_
#define PROFIT_DOT_PRODUCT_H_

#include "profit/config.h"
#include "profit/common.h"

/* Convenience macro */
#if defined(PROFIT_HAS_SSE2) || defined(PROFIT_HAS_AVX)
# define PROFIT_HAS_INTRINSICS
#else
# undef PROFIT_HAS_INTRINSICS
#endif

#ifdef PROFIT_HAS_SSE2
#include <emmintrin.h>
#endif // PROFIT_HAS_SSE2

#ifdef PROFIT_HAS_AVX
#include <immintrin.h>
#endif // PROFIT_HAS_AVX

namespace profit {

/*
 * =============================================================================
 * Software-based dot product implementation follows
 * =============================================================================
 */

//
// _dot_sw<N_Elements> performs the dot product of src and krn
// for the constant given number of elements.
// We first define it for 8, 4, 2 and 1 and build the cases for 7, 6, 5 and 3
// based on the former
//
template <unsigned int N_Elements>
static inline
double _dot_sw(const double *src, const double *krn);

template <>
double _dot_sw<1>(const double *src, const double *krn)
{
	return src[0] * krn[0];
}

template <>
double _dot_sw<2>(const double *src, const double *krn)
{
	double tmp1 = src[0] * krn[0];
	double tmp2 = src[1] * krn[1];
	return tmp1 + tmp2;
}

template <>
double _dot_sw<3>(const double *src, const double *krn)
{
	auto t1 = _dot_sw<2>(src, krn);
	auto t2 = _dot_sw<1>(src + 2, krn + 2);
	return t1 + t2;
}

template <>
double _dot_sw<4>(const double *src, const double *krn)
{
	double tmp1 = src[0] * krn[0];
	double tmp2 = src[1] * krn[1];
	double tmp3 = src[2] * krn[2];
	double tmp4 = src[3] * krn[3];
	return (tmp1 + tmp3) + (tmp2 + tmp4);
}

//
// _dot_remainder_sw<Batch_Size> performs the dot product of src and krn
// for a variable given number of elements n, which is known to be less than
// Batch_Size.
// We define it for 4 and for 8 (using the definition for 4), which are the
// two definitions we need.
template <unsigned int Batch_Size>
static inline
double _dot_remainder_sw(const double *src, const double *krn, std::size_t n);

template <>
double _dot_remainder_sw<2>(const double *src, const double *krn, std::size_t n)
{
	double buf = 0;
	if (n == 1) {
		buf += _dot_sw<1>(src, krn);
	}
	return buf;
}

template <>
double _dot_remainder_sw<4>(const double *src, const double *krn, std::size_t n)
{
	double buf = 0;
	if (n == 3) {
		buf += _dot_sw<3>(src, krn);
	}
	else if (n == 2) {
		buf += _dot_sw<2>(src, krn);
	}
	buf += _dot_remainder_sw<2>(src, krn, n);
	return buf;
}

//
// _dot_remainder_sw<Batch_Size> performs the dot product of src and krn
// for a variable given number of elements n, trying to perform most of the
// dot products using a batch size Batch_Size. The elements at the end of the
// arrays that don't make a full batch are processed via _dot_remainder_sw
//
template <unsigned int Batch_Size>
static inline
double _dot_sw(const double *src, const double *krn, std::size_t n)
{
	double buf = 0;
	for (size_t k = 0; k < n / Batch_Size; k++) {
		buf += _dot_sw<Batch_Size>(src + Batch_Size * k, krn + Batch_Size * k);
	}
	auto rem = n % Batch_Size;
	auto src_rem = src + n - rem;
	auto krn_rem = krn + n - rem;
	buf += _dot_remainder_sw<Batch_Size>(src_rem, krn_rem, rem);
	return buf;
}

//
// The sw-based function that is finally exposed to users
//
static inline
double dot_sw(const double *src, const double *krn, std::size_t n)
{
	return _dot_sw<4>(src, krn, n);
}


/*
 * =============================================================================
 * Intrinsics-based dot product implementation follows
 * =============================================================================
 */
#ifdef PROFIT_HAS_INTRINSICS
template <int>
struct intrinsic_traits;

template <simd_instruction_set Intrinsic>
static inline
typename intrinsic_traits<Intrinsic>::accum_type _zero_intrinsic();

template <simd_instruction_set Intrinsic>
static inline
double _extract_intrinsic(typename intrinsic_traits<Intrinsic>::accum_type final);

template <simd_instruction_set Intrinsic, unsigned short N>
static inline
typename intrinsic_traits<Intrinsic>::accum_type
_dot_intrinsic(const double *src, const double *krn, typename intrinsic_traits<Intrinsic>::accum_type accum);


#ifdef PROFIT_HAS_SSE2
template <>
struct intrinsic_traits<SSE2> {
	typedef __m128d accum_type;
};

template <>
__m128d _zero_intrinsic<SSE2>()
{
	return _mm_setzero_pd();
}

template <>
double _extract_intrinsic<SSE2>(__m128d final)
{
	auto x = _mm_add_pd(final, _mm_shuffle_pd(final, final, 0x01));
	return _mm_cvtsd_f64(x);
}


template <>
__m128d _dot_intrinsic<SSE2, 1>(const double *src, const double *krn, __m128d accum)
{
	auto _src = _mm_load_sd(src);
	auto _krn = _mm_load_sd(krn);
	auto dp = _mm_mul_pd(_src, _krn);
	return _mm_add_pd(dp, accum);
}

template <>
__m128d _dot_intrinsic<SSE2, 2>(const double *src, const double *krn, __m128d accum)
{
	auto _src = _mm_loadu_pd(src);
	auto _krn = _mm_loadu_pd(krn);
	auto dp = _mm_mul_pd(_src, _krn);
	return _mm_add_pd(dp, accum);
}

template <>
__m128d _dot_intrinsic<SSE2, 4>(const double *src, const double *krn, __m128d accum)
{
	auto src_1 = _mm_loadu_pd(src);
	auto src_2 = _mm_loadu_pd(src + 2);

	auto krn_1 = _mm_loadu_pd(krn);
	auto krn_2 = _mm_loadu_pd(krn + 2);

	auto dp_1 = _mm_mul_pd(src_1, krn_1);
	auto dp_2 = _mm_mul_pd(src_2, krn_2);

	return _mm_add_pd(_mm_add_pd(dp_1, dp_2), accum);
}
#endif // PROFIT_HAS_SSE2

#ifdef PROFIT_HAS_AVX
template <>
struct intrinsic_traits<AVX> {
	typedef __m256d accum_type;
};

template <>
__m256d _zero_intrinsic<AVX>()
{
	return _mm256_setzero_pd();
}

template <>
double _extract_intrinsic<AVX>(__m256d final)
{
	auto x = _mm256_hadd_pd(final, final);
	auto y = _mm256_extractf128_pd(x, 1);
	y = _mm_add_sd(y, _mm256_castpd256_pd128(x));
	return _mm_cvtsd_f64(y);
}

template <>
__m256d _dot_intrinsic<AVX, 1>(const double *src, const double *krn, __m256d accum)
{
	auto dot128 = _dot_intrinsic<SSE2, 1>(src, krn, _zero_intrinsic<SSE2>());
	auto t = _mm256_insertf128_pd(_zero_intrinsic<AVX>(), dot128, 0);
	return _mm256_add_pd(t, accum);
}

template <>
__m256d _dot_intrinsic<AVX, 2>(const double *src, const double *krn, __m256d accum)
{
	auto dot128 = _dot_intrinsic<SSE2, 2>(src, krn, _zero_intrinsic<SSE2>());
	auto t = _mm256_insertf128_pd(_zero_intrinsic<AVX>(), dot128, 0);
	return _mm256_add_pd(t, accum);
}

template <>
__m256d _dot_intrinsic<AVX, 4>(const double *src, const double *krn, __m256d accum)
{
	auto _src = _mm256_loadu_pd(src);
	auto _krn = _mm256_loadu_pd(krn);
	auto mul = _mm256_mul_pd(_src, _krn);
	return _mm256_add_pd(mul, accum);
}

template <>
__m256d _dot_intrinsic<AVX, 8>(const double *src, const double *krn, __m256d accum)
{
	auto src_1 = _mm256_loadu_pd(src);
	auto src_2 = _mm256_loadu_pd(src + 4);

	auto krn_1 = _mm256_loadu_pd(krn);
	auto krn_2 = _mm256_loadu_pd(krn + 4);

	auto mul_1 = _mm256_mul_pd(src_1, krn_1);
	auto mul_2 = _mm256_mul_pd(src_2, krn_2);

	return _mm256_add_pd(_mm256_add_pd(mul_1, mul_2), accum);
}
#endif // PROFIT_HAS_AVX

template <simd_instruction_set Intrinsic, unsigned int Batch_Size>
class _dot_remainder_instrinsic_calculator;

// We implement this using classes because C++ function templates cannot be
// partially specialized while classes can
template <simd_instruction_set Intrinsic, unsigned int Batch_Size>
static inline
typename intrinsic_traits<Intrinsic>::accum_type
_dot_remainder_intrinsic(const double *src, const double *krn, const std::size_t n, typename intrinsic_traits<Intrinsic>::accum_type buf)
{
	return _dot_remainder_instrinsic_calculator<Intrinsic, Batch_Size>::_(src, krn, n, buf);
}

template <simd_instruction_set Intrinsic>
class _dot_remainder_instrinsic_calculator<Intrinsic, 2> {
public:
	typedef typename intrinsic_traits<Intrinsic>::accum_type accum_type;
	static accum_type _(const double *src, const double *krn, const std::size_t n, typename intrinsic_traits<Intrinsic>::accum_type buf)
	{
		if (n == 1) {
			buf = _dot_intrinsic<Intrinsic, 1>(src, krn, buf);
		}
		return buf;
	}
};

template <simd_instruction_set Intrinsic>
class _dot_remainder_instrinsic_calculator<Intrinsic, 4> {
public:
	typedef typename intrinsic_traits<Intrinsic>::accum_type accum_type;
	static accum_type _(const double *src, const double *krn, const std::size_t n, typename intrinsic_traits<Intrinsic>::accum_type buf)
	{
		if (n == 3) {
			buf = _dot_intrinsic<Intrinsic, 2>(src, krn, buf);
			buf = _dot_intrinsic<Intrinsic, 1>(src + 2, krn + 2, buf);
		}
		else if (n == 2) {
			buf = _dot_intrinsic<Intrinsic, 2>(src, krn, buf);
		}
		return _dot_remainder_intrinsic<Intrinsic, 2>(src, krn, n, buf);
	}
};

template <simd_instruction_set Intrinsic>
class _dot_remainder_instrinsic_calculator<Intrinsic, 8> {
public:
	typedef typename intrinsic_traits<Intrinsic>::accum_type accum_type;
	static accum_type _(const double *src, const double *krn, const std::size_t n, typename intrinsic_traits<Intrinsic>::accum_type buf)
	{
		if (n == 7) {
			buf = _dot_intrinsic<Intrinsic, 4>(src, krn, buf);
			buf = _dot_intrinsic<Intrinsic, 2>(src + 4, krn + 4, buf);
			buf = _dot_intrinsic<Intrinsic, 1>(src + 6, krn + 6, buf);
		}
		else if (n == 6) {
			buf = _dot_intrinsic<Intrinsic, 4>(src, krn, buf);
			buf = _dot_intrinsic<Intrinsic, 2>(src + 4, krn + 4, buf);
		}
		else if (n == 5) {
			buf = _dot_intrinsic<Intrinsic, 4>(src, krn, buf);
			buf = _dot_intrinsic<Intrinsic, 1>(src + 4, krn + 4, buf);
		}
		else if (n == 4) {
			buf = _dot_intrinsic<Intrinsic, 4>(src, krn, buf);
		}
		return _dot_remainder_intrinsic<Intrinsic, 4>(src, krn, n, buf);
	}
};

template <simd_instruction_set Intrinsic, unsigned short Batch_Size>
static inline
double _dot_intrinsic(const double * src, const double * krn, std::size_t n)
{
	typedef typename intrinsic_traits<Intrinsic>::accum_type accum_type;
	accum_type accum_buf = _zero_intrinsic<Intrinsic>();
	for (size_t k = 0; k < n / Batch_Size; k++) {
		accum_buf = _dot_intrinsic<Intrinsic, Batch_Size>(src + Batch_Size * k, krn + Batch_Size * k, accum_buf);
	}
	auto rem = n % Batch_Size;
	auto src_rem = src + n - rem;
	auto krn_rem = krn + n - rem;
	accum_buf = _dot_remainder_intrinsic<Intrinsic, Batch_Size>(src_rem, krn_rem, rem, accum_buf);
	return _extract_intrinsic<Intrinsic>(accum_buf);
}
#endif // PROFIT_HAS_INTRINSICS

template <simd_instruction_set SIMD>
static inline
double dot_product(const double * src, const double * krn, std::size_t n)
{
	return dot_sw(src, krn, n);
}

#ifdef PROFIT_HAS_SSE2
template <>
double dot_product<SSE2>(const double *src, const double *krn, std::size_t n)
{
	return _dot_intrinsic<SSE2, 4>(src, krn, n);
}
#endif // PROFIT_HAS_SSE2

#ifdef PROFIT_HAS_AVX
template <>
double dot_product<AVX>(const double * src, const double * krn, std::size_t n)
{
	return _dot_intrinsic<AVX, 8>(src, krn, n);
}
#endif // PROFIT_HAS_AVX

template <>
double dot_product<AUTO>(const double *src, const double *krn, std::size_t n)
{
	// Choose the best if requested
#if defined(PROFIT_HAS_AVX)
	return _dot_intrinsic<AVX, 8>(src, krn, n);
#elif defined(PROFIT_HAS_SSE2)
	return _dot_intrinsic<SSE2, 4>(src, krn, n);
#else
	return dot_sw(src, krn, n);
#endif
}

}  // namespace profit

#endif /* PROFIT_DOT_PRODUCT_H_ */