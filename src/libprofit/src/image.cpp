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
#include <functional>
#include <numeric>
#include <utility>

#include "profit/image.h"

namespace profit {

Image::Image(unsigned int width, unsigned int height) :
	_2ddata(width, height)
{
}

Image::Image(const std::vector<double>& data, unsigned int width,
		unsigned int height) :
	_2ddata(data, width, height)
{
}

Image::Image(std::vector<double>&& data, unsigned int width,
		unsigned int height) :
	_2ddata(std::move(data), width, height)
{
}

Image::Image(const Image& other) :
	_2ddata(other)
{
}

Image::Image(Image&& other) :
	_2ddata(std::move(other))
{
}

double Image::getTotal() const {
	const auto &data = getData();
	return std::accumulate(data.begin(), data.end(), 0.);
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

Image Image::extend(unsigned int new_width, unsigned int new_height,
                    unsigned int start_x, unsigned int start_y) const
{
	auto width = getWidth();
	auto height = getHeight();

	if (new_width < width) {
		throw std::invalid_argument("new_width should be >= width");
	}
	if (new_height < height) {
		throw std::invalid_argument("new_height should be >= height");
	}
	if (start_x + width > new_width) {
		throw std::invalid_argument("start_x + new_width should be <= width");
	}
	if (start_y + height > new_height) {
		throw std::invalid_argument("start_y + new_height <= image.height");
	}

	Image extended(new_width, new_height);
	const auto &data = getData();
	auto &extended_data = extended.getData();
	for(unsigned int j = 0; j < height; j++) {
		for(unsigned int i = 0; i < width; i++) {
			extended_data[(i+start_x) + (j+start_y)*new_width] = data[i + j*width];
		}
	}
	return extended;
}

Image Image::crop(unsigned int new_width, unsigned int new_height,
                  unsigned int start_x, unsigned int start_y) const
{
	auto width = getWidth();
	auto height = getHeight();

	if (new_width > width) {
		throw std::invalid_argument("new_width should be <= width");
	}
	if (new_height > height) {
		throw std::invalid_argument("new_height should be <= height");
	}
	if (start_x + new_width > width) {
		throw std::invalid_argument("start_x + new_width should be <= image.width");
	}
	if (start_y + new_height > height) {
		throw std::invalid_argument("start_y + new_height should be <= image.height");
	}

	Image crop(new_width, new_height);
	const auto &data = getData();
	auto &crop_data = crop.getData();
	for(unsigned int j = 0; j < new_height; j++) {
		for(unsigned int i = 0; i < new_width; i++) {
			crop_data[i + j * new_width] = data[(i + start_x) + (j + start_y) * width];
		}
	};
	return crop;
}

Image &Image::operator&=(const Mask &mask)
{
	// Don't apply empty masks
	if (mask.empty()) {
		return *this;
	}

	auto &data = getData();
	const auto &mask_data = mask.getData();
	std::transform(data.begin(), data.end(), mask_data.begin(), data.begin(),
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

Image &Image::operator=(Image &&rhs)
{
	_2ddata<double>::operator=(std::move(rhs));
	return *this;
}

Image &Image::operator+=(const Image& rhs)
{
	auto &data = getData();
	const auto &other_data = rhs.getData();
	std::transform(data.begin(), data.end(), other_data.begin(), data.begin(), std::plus<double>());
	return *this;
}

const Image Image::operator+(const Image& rhs) const
{
	Image sum(*this);
	sum += rhs;
	return sum;

}

Image &Image::operator/=(double denominator)
{
	using std::placeholders::_1;
	auto &data = getData();
	std::transform(data.begin(), data.end(), data.begin(),
	               std::bind(std::divides<double>(), _1, denominator));
	return *this;
}

const Image Image::operator/(double denominator) const
{
	Image sum(*this);
	sum /= denominator;
	return sum;
}

}