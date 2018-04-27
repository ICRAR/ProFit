/**
 * Image class implementation
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

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <utility>

#include "profit/image.h"

namespace profit {

Mask::Mask(unsigned int width, unsigned int height) :
	surface({width, height})
{
}

Mask::Mask(Dimensions dimensions) :
	surface(dimensions)
{
}

Mask::Mask(const std::vector<bool>& data, unsigned int width, unsigned int height) :
	surface(data, {width, height})
{
}

Mask::Mask(const std::vector<bool>& data, Dimensions dimensions) :
	surface(data, dimensions)
{
}

Mask::Mask(std::vector<bool>&& data, unsigned int width, unsigned int height) :
	surface(std::move(data), {width, height})
{
}

Mask::Mask(std::vector<bool>&& data, Dimensions dimensions) :
	surface(std::move(data), dimensions)
{
}

Mask::Mask(const Mask& other) :
	surface(other)
{
}

Mask::Mask(Mask&& other) :
	surface(std::move(other))
{
}

Image::Image(unsigned int width, unsigned int height) :
	surface({width, height})
{
}

Image::Image(Dimensions dimensions) :
	surface(dimensions)
{
}

Image::Image(const std::vector<double>& data, unsigned int width, unsigned int height) :
	surface(data, {width, height})
{
}

Image::Image(const std::vector<double>& data, Dimensions dimensions) :
	surface(data, dimensions)
{
}

Image::Image(std::vector<double> &&data, unsigned int width, unsigned int height) :
	surface(std::move(data), {width, height})
{
}

Image::Image(std::vector<double>&& data, Dimensions dimensions) :
	surface(std::move(data), dimensions)
{
}

Image::Image(const Image& other) :
	surface(other)
{
}

Image::Image(Image&& other) :
	surface(std::move(other))
{
}

double Image::getTotal() const {
	return std::accumulate(begin(), end(), 0.);
}

void Image::normalize()
{
	double sum = getTotal();
	if( sum > 0 ) {
		*this /= sum;
	}
}

Image Image::normalize() const
{
	Image normalized(*this);
	normalized.normalize();
	return normalized;
}

Image Image::upsample(unsigned int factor, UpsamplingMode mode) const
{
	if (factor == 0) {
		throw std::invalid_argument("upsampling factor is 0");
	}
	if (factor == 1) {
		return Image(*this);
	}

	auto up_dims = getDimensions() * factor;
	Image upsampled(up_dims);

	// Each pixel from this image is repeated `factor x factor` times
	double divide_factor = mode == SCALE ? (factor * factor) : 1;
	for (unsigned int row_u = 0; row_u != up_dims.y; row_u++) {
		auto row = row_u / factor;
		for (unsigned int col_u = 0; col_u != up_dims.x; col_u++) {
			auto col = col_u / factor;
			upsampled[col_u + row_u * up_dims.x] = this->operator[]({col, row}) / divide_factor;
		}
	}
	return upsampled;
}

Image Image::downsample(unsigned int factor, DownsamplingMode mode) const
{
	using std::ceil;

	if (factor == 0) {
		throw std::invalid_argument("downsampling factor is 0");
	}
	if (factor == 1) {
		return Image(*this);
	}

	auto dims = getDimensions();
	auto down_dims = Dimensions{static_cast<unsigned int>(ceil(dims.x / float(factor))),
	                            static_cast<unsigned int>(ceil(dims.y / float(factor)))};
	Image downsampled(down_dims);

	// The following implementations might not be ideal, but our images are not
	// that big usually.

	if (mode == SAMPLE) {
		// We loop over the dimensions of the downsampled image and copy the
		// value of the first pixel of this image that corresponds
		for (unsigned int row_d = 0; row_d != down_dims.y; row_d++) {
			auto row = row_d * factor;
			for (unsigned int col_d = 0; col_d != down_dims.x; col_d++) {
				auto col = col_d * factor;
				downsampled[col_d + row_d * down_dims.x] = this->operator[]({col, row});
			}
		}
	}
	else if (mode == SUM) {
		// Pixels of this image are summed into the corresponding target pixels
		for (unsigned int row = 0; row != dims.y; row++) {
			auto row_d = row / factor;
			for (unsigned int col = 0; col != dims.x; col++) {
				auto col_d = col / factor;
				downsampled[col_d + row_d * down_dims.x] += this->operator[]({col, row});
			}
		}
	}
	else { // mode == AVERAGE
		// Pixels of this image are averaged into the corresponding target pixels
		for (unsigned int row_d = 0; row_d != down_dims.y; row_d++) {

			auto row_0 = row_d * factor;
			auto row_last = std::min(row_0 + factor, dims.y);

			for (unsigned int col_d = 0; col_d != down_dims.x; col_d++) {

				auto col_0 = col_d * factor;
				auto col_last = std::min(col_0 + factor, dims.x);

				// Accumulate the total flux from all the pixels that correspond
				// to this downsampled pixel and divide by the count
				double total = 0;
				unsigned int count = 0;
				for (unsigned int row = row_0; row != row_last; row++) {
					for (unsigned int col = col_0; col != col_last; col++) {
						total += this->operator[]({col, row});
						count++;
					}
				}
				downsampled[col_d + row_d * down_dims.x] = total / count;

			}
		}
	}

	return downsampled;
}


Image &Image::operator&=(const Mask &mask)
{
	// Don't apply empty masks
	if (mask.empty()) {
		return *this;
	}

	std::transform(begin(), end(), mask.begin(), begin(),
		[](const double i, const bool m) {
			return m ? i : 0.;
	});
	return *this;
}

const Image Image::operator&(const Mask &mask) const
{
	Image masked(*this);
	masked &= mask;
	return masked;
}

Image &Image::operator+=(const Image& rhs)
{
	std::transform(begin(), end(), rhs.begin(), begin(), std::plus<double>());
	return *this;
}

Image Image::operator+(const Image& rhs) const
{
	Image sum(*this);
	sum += rhs;
	return sum;

}

Image &Image::operator/=(double denominator)
{
	using std::placeholders::_1;
	std::transform(begin(), end(), begin(),
	               std::bind(std::divides<double>(), _1, denominator));
	return *this;
}

Image &Image::operator*=(double multiplier)
{
	using std::placeholders::_1;
	std::transform(begin(), end(), begin(),
	               std::bind(std::multiplies<double>(), _1, multiplier));
	return *this;
}

Image Image::operator/(double denominator) const
{
	Image sum(*this);
	sum /= denominator;
	return sum;
}

Image Image::operator/(int denominator) const
{
	return operator/(static_cast<double>(denominator));
}

Image Image::operator/(unsigned int denominator) const
{
	return operator/(static_cast<double>(denominator));
}

}