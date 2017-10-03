/**
 * Image and related classes definition
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


#ifndef PROFIT_IMAGE_H
#define PROFIT_IMAGE_H

#include <stdexcept>
#include <vector>

namespace profit {

/**
 * Base class for 2D-organized data
 */
template <typename T>
class _2ddata {

public:

	_2ddata(unsigned int width = 0, unsigned int height = 0) :
		data(width * height),
		width(width),
		height(height)
	{
		check_size();
	}

	_2ddata(const std::vector<T> &data, unsigned int width, unsigned int height) :
		data(data.begin(), data.end()),
		width(width),
		height(height)
	{
		check_size();
	}

	_2ddata(std::vector<T> &&data, unsigned int width, unsigned int height) :
		data(std::move(data)),
		width(width),
		height(height)
	{
		if (width * height != this->data.size()) {
			data = std::move(this->data);
			throw std::invalid_argument("data.size() != weight * height");
		}
	}

	/**
	 * Copy constructor
	 * @param other A different image
	 */
	_2ddata(const _2ddata &other) :
		data(other.data),
		width(other.width),
		height(other.height)
	{
		check_size();
	}

	/**
	 * Move constructor
	 * @param other A different image
	 */
	_2ddata(_2ddata &&other) :
		data(std::move(other.data)),
		width(other.width),
		height(other.height)
	{
		other.width = 0;
		other.height = 0;
	}

	bool empty() const {
		return width == 0 || height == 0;
	}

	const std::vector<T>& getData() const {
		return data;
	}

	std::vector<T>& getData() {
		return data;
	}

	unsigned int getHeight() const {
		return height;
	}

	unsigned int getSize() const {
		return width * height;
	}

	unsigned int getWidth() const {
		return width;
	}

	bool operator==(const _2ddata &other) const {
		return width == other.width &&
		       height == other.height &&
		       data == other.data;
	}

	_2ddata &operator=(_2ddata &&rhs) {
		data = std::move(rhs.data);
		width = rhs.width;
		height = rhs.height;
		rhs.width = 0;
		rhs.height = 0;
		return *this;
	}

protected:
	void setWidth(unsigned int width) {
		this->width = width;
	}

	void setHeight(unsigned int height) {
		this->height = height;
	}

	void setData(std::vector<double> &&data) {
		this->data = data;
	}

private:
	std::vector<T> data;
	unsigned int width;
	unsigned int height;

	void check_size()
	{
		if (width * height != data.size()) {
			throw std::invalid_argument("data.size() != weight * height");
		}
	}
};

/**
 * A mask is 2D area of bools
 */
typedef _2ddata<bool> Mask;

/**
 * An image is a 2D area of doubles.
 */
class Image : public _2ddata<double> {

public:

	// Constructors that look that those from _2ddata
	//
	// C++11 supports inheriting ctors like this:
	//
	// using _2ddata::_2ddata;
	//
	// Sadly, gcc supports this only since 4.8. Since we still want to keep
	// supporting gcc >= 4.6.3 for a while, we still need to keep this verbose
	// code that simply forwards each of the arguments to the parent ctor
	Image(unsigned int width = 0, unsigned int height = 0);
	Image(const std::vector<double> &data, unsigned int width, unsigned int height);
	Image(std::vector<double> &&data, unsigned int width, unsigned int height);
	Image(const Image &other);
	Image(Image &&other);

	/**
	 * Returns the sum of the image pixel's values
	 *
	 * @return The sum of the image pixel's values
	 */
	double getTotal() const;

	/**
	 * Normalized this image; i.e., rescales its values so the sum of all its
	 * pixels' values is 1. If all pixels are 0 the image is not changed.
	 */
	void normalize();

	/**
	 * Returns a normalized version of this image; i.e., one where the sum of
	 * all pixels' values is 1. If all pixels are 0 the returned image is
	 * identical to this image.
	 */
	Image normalize() const;

	/**
	 * Creates a new image that is an extension of this image. The new
	 * dimensions must be greater or equal to the current dimensions.
	 * Additionally, the contents of this image are placed on the new image
	 * at ``(start_x, start_y)`` relative to the new image's dimensions; the
	 * rest is initialized with zeroes.
	 *
	 * @param new_width The width of the new image, it should be ``>= width``
	 * @param new_height The height of the new image, it should be ``>= height``
	 * @param start_x The X start of the new image relative to the new image
	 * @param start_y The Y start of the new image relative to the new image
	 * @return The new extended image
	 */

	Image extend(unsigned int new_width, unsigned int new_height,
	             unsigned int start_x, unsigned int start_y) const;

	/**
	 * Creates a new image that is a crop of this image. The cropped image
	 * starts at ``(start_x, start_y)`` (relative to this image) and its
	 * dimensions are ``new_width x new_height``.
	 *
	 * @param new_width The width of the new image. It should be ``<= width``
	 * @param new_height The height of the new image. It should be ``<= height``
	 * @param start_x The X start of the new image relative to this image.`
	 * @param start_y The Y start of the new image relative to this image.
	 * @return The new cropped image
	 */
	Image crop(unsigned int new_width, unsigned int new_height,
	           unsigned int start_x, unsigned int start_y) const;


	// Support for operator overloading

	/// Move assignment of an Image
	Image &operator=(Image &&rhs);

	/// Addition assignment of another Image
	Image &operator+=(const Image &rhs);

	/// Addition of another image
	const Image operator+(const Image &rhs) const;

	/// Division assignment against a double denominator
	Image &operator/=(double denominator);

	/// Division against a double denominator
	const Image operator/(double denominator) const;

	Image &operator|=(const Mask &mask);

	/// Bitwise AND assignment with a Mask (applies the mask to the image).
	Image &operator&=(const Mask &mask);

	/// Bitwise AND with a Mask (applies the mask to the image).
	const Image operator&(const Mask &mask) const;

};

}  // namespace profit

#endif // PROFIT_IMAGE_H