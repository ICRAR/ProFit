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

#include <algorithm>
#include <ostream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "profit/common.h"
#include "profit/coordinates.h"
#include "profit/point.h"

namespace profit {

/// @typedef A discrete 2-dimensional dimension definition
using Dimensions = discrete_2d_coordinate;

/// A pair of points
typedef std::pair<Point, Point> PointPair;

/// A box defined by the lowest (inclusive) and highest (exclusive) 2D points
/// that contain it.
class Box : public PointPair {
public:

	/// Creates an empty box
	Box() : PointPair({0, 0}, {0, 0})
	{ }

	/**
	 * Creates a box starting at @p lb (inclusive) and ending at @p ub
	 * (exclusive).
	 * @param lb The lower boundaries of the box
	 * @param ub The upper boundaries (exclusive in both dimensions) of the box
	 */
	Box(Point lb, Point ub) : PointPair(lb, ub)
	{
		if (!(first <=second)) {
			throw std::invalid_argument("box's lower boundary must be <= than high boundary");
		}
	}

	/// Whether this object represents an empty box
	bool empty() const {
		return first == second;
	}

	Box operator*(unsigned int n)
	{
		return {first * n, second * n};
	}
};

template <typename CharT>
std::basic_ostream<CharT> &operator<<(std::basic_ostream<CharT> &os, const Box &box)
{
	os << '[' << box.first << ", " << box.second << ']';
	return os;
}

///
/// Non-templated code common to 2D surface classes, holding only the dimensions
/// of the surface but not its data.
///
/// Note that when surface_base is moved, its dimensions are explicitly set to
/// (0, 0) so the moved-from image acts like an empty image. We still want
class surface_base {

public:

	/// Default constructor, creates an empty surface_base
	surface_base() = default;

	/// Copy constructor
	surface_base(const surface_base &other) = default;

	/// Move constructor, sets the moved-from surface_base as empty
	surface_base(surface_base &&other) :
		dimensions(std::move(other.dimensions))
	{
		other.dimensions = {0, 0};
	}

	/// Copy assignment operator
	surface_base &operator=(const surface_base &other) = default;

	/// Move assignment operator, sets the moved-from surface_base as empty
	surface_base &operator=(surface_base &&other)
	{
		dimensions = std::move(other.dimensions);
		other.dimensions = {0, 0};
		return *this;
	}

	/// Constructs a surface_base with the given dimensions
	explicit surface_base(Dimensions dimensions) :
		dimensions(dimensions)
	{
		// no-op
	}

	unsigned int getHeight() const {
		return dimensions.y;
	}

	unsigned int size() const {
		return dimensions.x * dimensions.y;
	}

	unsigned int getWidth() const {
		return dimensions.x;
	}

	Dimensions getDimensions() const {
		return dimensions;
	}

	bool empty() const {
		return dimensions.x == 0 && dimensions.y == 0;
	}

	/// Surfaces are true if they have a dimension
	explicit operator bool() const {
		return dimensions.x > 0 && dimensions.y > 0;
	}

	/// Comparison operator
	bool operator==(const surface_base &other) const {
		return dimensions == other.dimensions;
	}

protected:

	void _extension_is_possible(const Dimensions &new_size, const Point &start) const {
		if (new_size.x < dimensions.x) {
			throw std::invalid_argument("new_width should be >= width");
		}
		if (new_size.y < dimensions.y) {
			throw std::invalid_argument("new_height should be >= height");
		}
		if (start.x + dimensions.x > new_size.x) {
			throw std::invalid_argument("start_x + new_width should be <= width");
		}
		if (start.y + dimensions.y > new_size.y) {
			throw std::invalid_argument("start_y + new_height <= image.height");
		}
	}

	void _crop_is_possible(const Dimensions &new_size, const Point &start) const
	{
		if (new_size.x > dimensions.x) {
			throw std::invalid_argument("new_width should be <= width");
		}
		if (new_size.y > dimensions.y) {
			throw std::invalid_argument("new_height should be <= height");
		}
		if (start.x + new_size.x > dimensions.x) {
			throw std::invalid_argument("start_x + new_width should be <= image.width");
		}
		if (start.y + new_size.y > dimensions.y) {
			throw std::invalid_argument("start_y + new_height should be <= image.height");
		}
	}

private:
	Dimensions dimensions;

};

namespace detail {

template <typename T>
struct is_bool : public std::true_type {};

template <>
struct is_bool<bool> : public std::false_type {};

template <typename T>
struct iterator_type
{
	using iterator = T *;
	using const_iterator = const T *;
};

template <>
struct iterator_type<bool>
{
	using iterator = std::vector<bool>::iterator;
	using const_iterator = std::vector<bool>::const_iterator;
};
}

/**
 * Base class for 2D-organized data
 */
template <typename T, typename D>
class surface : public surface_base {

public:

	typedef typename std::vector<T>::value_type value_type;
	typedef typename std::vector<T>::reference reference;
	typedef typename std::vector<T>::const_reference const_reference;
	typedef typename std::vector<T>::size_type size_type;
	using iterator = typename detail::iterator_type<T>::iterator;
	using const_iterator = typename detail::iterator_type<T>::const_iterator;

	surface() = default;

	explicit surface(Dimensions dimensions) :
		surface_base(dimensions),
		_managed_data(dimensions.x * dimensions.y)
	{
		// no-op
	}

	surface(const std::vector<T> &data, Dimensions dimensions) :
		surface_base(dimensions),
		_managed_data(data.begin(), data.end())
	{
		check_size();
	}

	surface(std::vector<T> &&data, Dimensions dimensions) :
		surface_base(dimensions),
		_managed_data(std::move(data))
	{
		if (dimensions.x * dimensions.y != this->_managed_data.size()) {
			data = std::move(this->_managed_data);
			throw std::invalid_argument("data.size() != weight * height");
		}
	}

	template <typename Void=std::enable_if<!detail::is_bool<T>::value>>
	surface(T *data, Dimensions dimensions) :
		surface_base(dimensions),
		_external_data(data)
	{
	}

	/**
	 * Assigns zero to all elements of this surface.
	 */
	void zero() {
		auto *_data = data();
		std::fill(_data, _data + size(), static_cast<T>(0));
	}

	/**
	 * Creates a new surface that is an extension of this object. The new
	 * dimensions must be greater or equal to the current dimensions.
	 * The current contents of this surface are placed at @p start, relative to
	 * the new surface's dimension.
	 *
	 * @param dimensions The dimensions of the new extended surface.
	 * @param start The starting point of the original surface relative to the new one
	 * @return The new extended surface
	 */
	D extend(Dimensions dimensions, Point start = Point()) const
	{
		_extension_is_possible(dimensions, start);
		D extended(dimensions);
		extend(extended, start);
		return extended;
	}

	/**
	 * Extends this object into the given surface. The new surface's
	 * dimensions must be greater or equal to the current dimensions.
	 * The current contents of this surface are placed at @p start, relative to
	 * the new surface's dimension.
	 *
	 * @param extended The new surface to hold the extended version of this image.
	 * Its dimensions mandate how much the current image should extend.
	 * @param start The starting point of the original surface relative to the new one
	 */
	void extend(D &extended, Point start = Point()) const
	{
		auto dimensions = extended.getDimensions();
		_extension_is_possible(dimensions, start);
		auto _data = data();
		for(unsigned int j = 0; j < getHeight(); j++) {
			for(unsigned int i = 0; i < getWidth(); i++) {
				extended[(i+start.x) + (j+start.y)*dimensions.x] = _data[i + j*getWidth()];
			}
		}
	}

	/**
	 * Creates a new image that is a crop of this image. The cropped image
	 * starts at @p start (relative to this image) and has new dimensions
	 * @p dimensions.
	 *
	 * @param dimensions The dimensions of the cropped image. They should be less or
	 *        equal than the dimensions of this image.
	 * @param start The start of the new image relative to this image.
	 * @return The new cropped image
	 */
	D crop(Dimensions dimensions, Point start = Point()) const
	{
		_crop_is_possible(dimensions, start);
		D crop(dimensions);
		auto _data = data();
		for(unsigned int j = 0; j < dimensions.y; j++) {
			for(unsigned int i = 0; i < dimensions.x; i++) {
				crop[i + j * dimensions.x] = _data[(i + start.x) + (j + start.y) * getWidth()];
			}
		}
		return crop;
	}

	/**
	 * Returns a copy of this surface with its underlying values in the reversed
	 * order, such that the top-right corner is now that bottom-left corner and
	 * vice-versa.
	 *
	 * @return A new object with reversed values
	 */
	D reverse() const
	{
		D reversed(static_cast<const D &>(*this));
		std::reverse(reversed.begin(), reversed.end());
		return reversed;
	}

	/**
	 * Returns a "value-interesting" bounding box for this surface; that is,
	 * the subset of this surface inside which all values are different from
	 * zero.
	 *
	 * @return The minimum bounding box within which all non-zero values of this
	 * surface are contained.
	 */
	Box bounding_box() const
	{
		Point lb = getDimensions();
		Point ub;
		auto _data = data();
		bool only_zeros = true;
		for(unsigned int j = 0; j < getHeight(); j++) {
			for(unsigned int i = 0; i < getWidth(); i++) {
				if (_data[i + j * getWidth()] == 0) {
					continue;
				}
				if (i < lb.x) {
					lb.x = i;
				}
				if (j < lb.y) {
					lb.y = j;
				}
				if ((i + 1) > ub.x) {
					ub.x = i + 1;
				}
				if ((j + 1) > ub.y) {
					ub.y = j + 1;
				}
				only_zeros = false;
			}
		}
		if (only_zeros) {
			return {};
		}
		return {lb, ub};
	}

	/// Comparison operator
	bool operator==(const surface &other) const {
		return surface_base::operator==(other) &&
		       _managed_data == other._managed_data &&
		       _external_data == other._external_data;
	}

	/// subscript operator
	reference operator[](const size_type idx)
	{
		return data()[idx];
	}

	/// subscript operator, const
	const_reference operator[](const size_type idx) const
	{
		return data()[idx];
	}

	/// [] operator that works with a Point
	reference operator[](const Point &p)
	{
		return data()[p.x + p.y * getWidth()];
	}

	const_reference operator[](const Point &p) const
	{
		return data()[p.x + p.y * getWidth()];
	}

	/// iterator to beginning of data
	iterator begin() { return data(); }
	const_iterator begin() const { return data(); }
	const_iterator cbegin() const { return data(); }

	/// iterator to end of data
	iterator end() { return data() + size(); }
	const_iterator end() const { return data() + size(); }
	const_iterator cend() const { return data() + size(); }

	/// type casting to std::vector<T>
	explicit operator std::vector<T>() const {
		return std::vector<T>(data(), data() + size());
	}

	iterator data()
	{
		return static_cast<D &>(*this).data_impl();
	}

	const_iterator data() const
	{
		return static_cast<const_iterator>(
			const_cast<surface *>(this)->data());
	}

protected:
	std::vector<T> &_get_managed_data() {
		return _managed_data;
	}

	const std::vector<T> &_get_managed_data() const {
		return _managed_data;
	}

	T *_get_external_data() {
		return _external_data;
	}

	const T *_get_external_data() const {
		return _external_data;
	}

	void _set_external_data(T *external_data)
	{
		_external_data = external_data;
	}

private:
	std::vector<T> _managed_data;
	T *_external_data = nullptr;

	void check_size()
	{
		if (getWidth() * getHeight() != _managed_data.size()) {
			throw std::invalid_argument("data.size() != weight * height");
		}
	}

};

namespace detail {

	template <typename ValueT>
	struct printing_return {
		typedef ValueT type;
	};

	template <>
	struct printing_return<bool> {
		typedef char type;
	};

	template <typename ValueT>
	inline
	typename printing_return<ValueT>::type surface_value_for_printing(ValueT value)
	{
		return value;
	}

	template <>
	inline
	char surface_value_for_printing<bool>(bool value)
	{
		if (value) {
			return 'T';
		}
		return 'F';
	}
} // namespace detail

template <typename CharT, typename ValueT, typename Derived>
std::basic_ostream<CharT> &operator<<(std::basic_ostream<CharT> &os, const surface<ValueT, Derived> &s)
{
	os << '[';
	auto width = s.getWidth();
	auto height = s.getHeight();
	for (unsigned int j = 0; j != height; j++) {
		os << '[';
		for (unsigned int i = 0; i != width; i++) {
			os << detail::surface_value_for_printing(s[Point{i, j}]);
			if (i < width - 1) {
				os << ", ";
			}
		}
		os << ']';
		if (j < height - 1) {
			os << ", ";
		}
	}
	os << ']';
	return os;
}

/**
 * A mask is surface of bools
 */
class PROFIT_API Mask : public surface<bool, Mask> {

public:

	// Constructors that look like those from _surface
	Mask() = default;
	Mask(unsigned int width, unsigned int height);
	Mask(bool value, unsigned int width, unsigned int height);
	Mask(bool value, Dimensions dimensions);
	explicit Mask(Dimensions dimensions);
	Mask(const std::vector<bool> &data, unsigned int width, unsigned int height);
	Mask(const std::vector<bool> &data, Dimensions dimensions);
	Mask(std::vector<bool> &&data, unsigned int width, unsigned int height);
	Mask(std::vector<bool> &&data, Dimensions dimensions);

	/**
	 * Returns a new Mask where the area covered by the new mask (i.e., where
	 * the new mask's value is @p true) is an "expanded" version of this mask.
	 * This is similar in nature to a convolution, but simpler as it is a
	 * simpler boolean operation that requires no additions or further scaling.
	 *
	 * @param pad the amount of cells to expand each input pixel on each dimension.
	 * @param threads threads to use to perform computation. Only valid if
	 * compiled with OpenMP support
	 */
	Mask expand_by(Dimensions pad, int threads=1) const;

	/**
	 * Upsamples this mask by the given factor.
	 *
	 * The resulting mask's dimensions will be the original mask's times the
	 * upsampling factor. The original mask's values are copied on the
	 * corresponding <pre>factor * factor</pre> cells of the upsampled
	 * mask.
	 *
	 * @param factor The upsampling factor. Must be greater than 0. If equals to
	 * 1, the upsampled mask is equals to the original mask.
	 * @return The upsampled mask
	 */
	Mask upsample(unsigned int factor) const;

	iterator data_impl()
	{
		return _get_managed_data().begin();
	}
};

/**
 * An image is a surface of doubles.
 */
class PROFIT_API Image : public surface<double, Image> {

public:

	// Constructors that look like those from _surface
	Image() = default;
	Image(const Image &other);
	Image(Image &&other) = default;
	Image(unsigned int width, unsigned int height);
	Image(double value, Dimensions dimensions);
	Image(double *values, Dimensions dimensions);
	Image(double value, unsigned int width, unsigned int height);
	explicit Image(Dimensions dimensions);
	Image(const std::vector<double> &data, unsigned int width, unsigned int height);
	Image(const std::vector<double> &data, Dimensions dimensions);
	Image(std::vector<double> &&data, unsigned int width, unsigned int height);
	Image(std::vector<double> &&data, Dimensions dimensions);

	Image &operator=(const Image &other);
	Image &operator=(Image &&other) = default;

	/** Available image upsampling modes */
	enum UpsamplingMode {

		/**
		 * Scales the value of the original pixel by `factor * factor` before
		 * copying it into each corresponding upsampled pixel. This has the
		 * effect of preserving the total flux of the original image
		 */
		SCALE = 0,

		/**
		 * Copies the value of the original pixel unmodified into each
		 * corresponding upsampled pixel
		 */
		COPY

	};

	/** Available image downsampling modes */
	enum DownsamplingMode {

		/**
		 * Pixel values on the resulting image are the average of the
		 * corresponding pixels on the original image.
		 */
		AVERAGE = 0,

		/**
		 * Pixel values on the resulting image are the sum of the corresponding
		 * pixels on the original image.
		 */
		SUM,

		/**
		 * Pixel values on the resulting image are samples from the original
		 * image. Samples are taken from the lowest placed pixel, in both
		 * dimensions, of the corresponding pixels of the original image.
		 */
		SAMPLE

	};

	/**
	 * Returns the sum of the image pixel's values (or "total flux").
	 *
	 * @return The sum of the image pixel's values
	 */
	double total() const;

	/**
	 * Upsamples this image by the given factor.
	 *
	 * The resulting image's dimensions will be the original image's times the
	 * upsampling factor. The particular upsampling method is determined by
	 * @p mode
	 *
	 * @param factor The upsampling factor. Must be greater than 0. If equals to
	 * 1, the upsampled image is equals to the original image.
	 * @param mode The upsampling mode to use
	 * @return An upsampled image, without interpolation.
	 */
	Image upsample(unsigned int factor, UpsamplingMode mode = SCALE) const;

	/**
	 * Upsamples this image by the given factor. The resulting image is written
	 * to @p target.
	 *
	 * The resulting image's dimensions will be the original image's times the
	 * upsampling factor. If @p target doesn't have these dimensions an error is
	 * raised. The particular upsampling method is determined by @p mode.
	 *
	 * @param target The target image where the result of the upsampling will be
	 * written to.
	 * @param factor The upsampling factor. Must be greater than 0. If equals to
	 * 1, the upsampled image is equals to the original image.
	 * @param mode The upsampling mode to use
	 * @return An upsampled image, without interpolation.
	 */
	void upsample(Image &target, unsigned int factor, UpsamplingMode mode = SCALE) const;

	/**
	 * Downsamples this image by the given factor.
	 *
	 * The resulting image's dimensions will be the ceiling of this image's
	 * divided by the downsampling factor. The particular downsampling method is
	 * determined by @p mode
	 *
	 * @param factor The downsampling factor. Must be greater than 0. If equals
	 * to 1, the upsampled image is equals to the original image.
	 * @param mode The downsampling mode to use
	 * @return A downsampled image, without interpolation.
	 */
	Image downsample(unsigned int factor, DownsamplingMode = SUM) const;

	/**
	 * Downsamples this image by the given factor. The resulting image is written
	 * to @p target.
	 *
	 * The resulting image's dimensions will be the ceiling of this image's
	 * divided by the downsampling factor. If @p target doesn't have these
	 * dimensions an error is raised. The particular downsampling method is
	 * determined by @p mode.
	 *
	 * @param target The target image where the result of the downsampling will
	 * be written to.
	 * @param factor The downsampling factor. Must be greater than 0. If equals
	 * to 1, the upsampled image is equals to the original image.
	 * @param mode The downsampling mode to use
	 * @return A downsampled image, without interpolation.
	 */
	void downsample(Image &downsampled, unsigned int factor, DownsamplingMode mode) const;

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
	 * Exposes the underlying data pointer in case it becomes necessary to
	 * access it directly.
	 *
	 * @return The underlying data pointer
	 */
	value_type *data_impl() {
		auto external_data = _get_external_data();
		if (external_data) {
			return external_data;
		}
		return _get_managed_data().data();
	}

	/**
	 * Exposes the underlying data pointer in case it becomes necessary to
	 * access it directly
	 *
	 * @return The underlying data pointer
	 */
	const value_type *data_impl() const {
		return const_cast<const value_type *>(
		    const_cast<Image *>(this)->data_impl()
		);
	}

	/// Addition assignment of another Image
	Image &operator+=(const Image &rhs);

	/// Addition of another image
	Image operator+(const Image &rhs) const;

	/// Division assignment against a double denominator
	Image &operator/=(double denominator);

	/// Multiplication assignment against a double multiplier
	Image &operator*=(double denominator);

	/// Division against a double denominator
	Image operator/(double denominator) const;
	Image operator/(int denominator) const;
	Image operator/(unsigned int denominator) const;

	Image &operator|=(const Mask &mask);

	/// Bitwise AND assignment with a Mask (applies the mask to the image).
	Image &operator&=(const Mask &mask);

	/// Bitwise AND with a Mask (applies the mask to the image).
	const Image operator&(const Mask &mask) const;

};

}  // namespace profit

#endif // PROFIT_IMAGE_H